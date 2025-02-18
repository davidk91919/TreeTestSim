library(data.table)

#' Generate a Complete k-ary Tree as a data.table
#'
#' @description
#' Constructs a complete k-ary tree with levels 0 to \code{max_level}.
#' The tree is stored in a data.table with one row per node and columns for:
#' \code{node} (the node's index),
#' \code{level} (its level in the tree, with 0 = root),
#' \code{parent} (the index of its parent; NA for the root),
#' and \code{nonnull} (a logical indicating whether the node is non-null).
#'
#' At the leaves (nodes at \code{level == max_level}), the \code{nonnull} flag
#' is sampled as \code{TRUE} with probability \code{t}. For internal nodes,
#' \code{nonnull} is set to \code{TRUE} if any of its children are non-null.
#'
#' @param max_level Integer. The maximum level of the tree (with the root at level 0).
#' @param k Integer. The branching factor (each internal node has exactly \code{k} children).
#' @param t Numeric in [0,1]. The probability that a leaf is non-null (has a treatment effect).
#'
#' @return A data.table with columns \code{node}, \code{level}, \code{parent}, and \code{nonnull}.
#'
#' @examples
#' dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.2)
#' print(dt)
#'
#' @export
generate_tree_DT <- function(max_level, k, t) {
  # Total number of nodes in a complete k-ary tree:
  total_nodes <- sum(k^(0:max_level))

  # Initialize the data.table with node IDs.
  treeDT <- data.table(node = 1:total_nodes)

  # Compute levels by assigning nodes in breadth-first order.
  levels <- integer(total_nodes)
  cum <- cumsum(k^(0:max_level))
  start <- 1
  for (l in 0:max_level) {
    end <- cum[l + 1]
    levels[start:end] <- l
    start <- end + 1
  }
  treeDT[, level := levels]

  # Compute parent indices.
  # In a complete k-ary tree stored in breadth-first order:
  # For node i (i > 1), parent = floor((i - 2) / k) + 1.
  treeDT[, parent := ifelse(node == 1, NA_integer_, floor((node - 2) / k) + 1)]

  # Allocate the "nonnull" flag.
  # For leaves (level == max_level), sample from Bernoulli(t).
  ## treeDT[level == max_level, nonnull := runif(.N) < t]
  # Complete Random assignment at leaves
  treeDT[, nonnull := NA]
  treeDT[level == max_level, nonnull := {
    n_leaves <- .N
    m_leaves <- floor(t * n_leaves)
    sample(rep(c(FALSE, TRUE), c(n_leaves - m_leaves, m_leaves)))
  }]

  # Propagate nonnull flags upward (for internal nodes).
  if (max_level > 0) {
    for (l in (max_level - 1):0) {
      # For nodes at level l, set nonnull = TRUE if any of its children are TRUE.
      treeDT[level == l, nonnull := {
        sapply(node, function(i) {
          # Compute children indices for node i:
          children <- ((i - 1) * k + 2):((i - 1) * k + k + 1)
          children <- children[children <= total_nodes]
          any(treeDT[children, nonnull], na.rm = TRUE)
        })
      }]
    }
  } else {
    # Only one node.
    treeDT[node == 1, nonnull := runif(1) < t]
  }

  return(treeDT)
}

#' Simulate the Testing Procedure on a k-ary Tree (data.table version)
#'
#' @description
#' Simulates the hierarchical testing procedure on a complete k-ary tree (given as a data.table).
#' The procedure works level-by-level. At the root the test is performed using parent p = 0.
#' For each active node (one that is "tested" because its parent's p-value was \(\le \alpha\)),
#' a local p-value is simulated:
#'
#' \itemize{
#'   \item For a non-null node:
#'     \[ p_{\text{local}} = \text{parent}_p + (1-\text{parent}_p) \times \mathrm{rbeta}(1,\, \beta_{\text{base}} \times \sqrt{effN/N_{\text{total}}},\, 1) \]
#'   \item For a null node:
#'     \[ p_{\text{local}} \sim \mathrm{Uniform}(\text{parent}_p,1). \]
#' }
#'
#' Then, each child’s final p-value is defined as
#' \[ p_{\text{child}} = \max(\text{parent}_p,\, p_{\text{local}}), \]
#' and for each parent the children are combined using a Simes test.
#' If the Simes p-value exceeds \(\alpha\), that branch is gated (the children’s p-values remain NA).
#'
#' @param treeDT A data.table as produced by \code{generate_tree_DT()}.
#' @param alpha Numeric. The significance level.
#' @param k Integer. The branching factor.
#' @param effN Numeric. The effective sample size at the root (assumed to decrease by a factor of k per level).
#' @param N_total Numeric. The total sample size at the root. For now effN should be same as N_total
#' @param beta_base Numeric. The base parameter for the beta distribution.
#' @param  adj_effN Logical. Whether the power of the test should decrease to simulate splitting the sample size.
#' @param  local_adj_p_fn Function. The name of a function that operates like \code{local_simes} or \code{local_hommel_all_ps} to adjust the p-values at a given node.
#' @param  global_adj Character. The name of a function that will adjust the p-values on the leaves for comparison with the bottom-up approaches
#' @param  return_details Logical. Whether the full data set that is simulated should be returned.
#'
#' @return A list with components:
#' \describe{
#'   \item{treeDT}{The data.table augmented with a column \code{p_val} for the simulated p-value.}
#'   \item{false_error}{A logical flag equal to TRUE if any truly null node (nonnull == FALSE)
#'          that was tested yielded \(\texttt{p_val} \le \alpha\).}
#' }
#'
#' @examples
#' dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.2)
#' res <- simulate_test_DT(dt, alpha = 0.05, k = 3, effN = 1000, N_total = 1000, beta_base = 0.1)
#' print(res$treeDT)
#'
#' @export
simulate_test_DT <- function(treeDT, alpha, k, effN, N_total, beta_base,
                             adj_effN = TRUE, local_adj_p_fn = local_simes, global_adj = "hommel", return_details = TRUE) {
  # Work on a copy so that the original tree is preserved.
  tree_sim <- copy(treeDT)
  tree_sim[, p_val := NA_real_]
  setkey(tree_sim, "node")

  # --- Level 0: Simulate for the root.
  tree_sim[node == 1, p_val := fifelse(
    nonnull,
    ## 0 + (1 - 0) * rbeta(1, beta_base * sqrt(effN / N_total), 1),
    rbeta(1, beta_base, 1),
    runif(1, min = 0, max = 1)
  )]

  max_level <- max(tree_sim$level)

  ######################################## The bottom up approach #####################
  ## For comparison with the bottom up approaches:  test in all of the leaves
  ## according to whether they are null or nonnull using the effective sample
  ## size (which here is capped at a power of about .22 but is otherwise
  ## proportional to N/number of leaves

  effN_leaves <- N_total / (k^max_level)
  beta_eff_leaves_raw <- beta_base * sqrt(N_total / effN_leaves)
  cap_leaves <- .5
  eff_beta_leaves <- min(cap_leaves, beta_eff_leaves_raw)
  tree_sim[
    level == max_level,
    p_sim := fifelse(
      nonnull,
      rbeta(.N, eff_beta_leaves, 1),
      runif(.N, min = 0, max = 1)
    )
  ]

  # Record a false error if any truly null node that was tested (i.e. has a non-NA p-value)
  # produces a p-value <= alpha.
  if (global_adj == "hommel") {
    global_adj_fn <- local_hommel_all_ps
  }
  tree_sim[level == max_level, bottom_up_p_adj := global_adj_fn(p_sim)]
  bottom_up_false_error <- any(tree_sim[level == max_level & nonnull == FALSE, bottom_up_p_adj] <= alpha)

  ## Record number of true discoveries --- non-null nodes tested and with p <= alpha
  bottom_up_true_discoveries <- sum(tree_sim[level == max_level & nonnull == TRUE, bottom_up_p_adj] <= alpha)

  ## Record number of nodes tested
  num_leaves <- k^max_level
  #######################################################

  # --- Process level-by-level.
  for (l in 1:max_level) {
    # Identify active parents at level l-1 (those that were tested and had p_val <= alpha).
    # For example, If we cannot reject at root (l-1)=0 but p_val > alpha,  the loop stops here
    active_parents <- tree_sim[level == (l - 1) & !is.na(p_val) & p_val <= alpha, node]
    if (length(active_parents) == 0) {
      break
    } # no branch continues.

    # Get children of active parents (at level l).
    child_rows <- tree_sim[parent %in% active_parents & level == l]
    if (nrow(child_rows) == 0) {
      next
    }

    if (adj_effN) {
      # Assume that the effective sample size is split equally. But put a lower
      # limit on power, say, no more than 1/100 the power of the root node
      ## mean(rbeta(10000,.5,1)<=.05)
      # Raw effective beta based on the loss in sample size:
      ## Basically we make equal splits
      effN <- effN / k^l
      ## And we know that power is a function of the sqrt of N.
      beta_eff_raw <- beta_base * sqrt(N_total / effN)
      # Compute the cap; note that log(0.05^(beta_base)) is beta_base*log(0.05)
      # cap <- (beta_base * log(alpha) - log(100)) / log(alpha)
      # mean(rbeta(10000,.9,1)<=.05)  is about .07 which is more powerful than
      # uniform but still low power mean(rbeta(10000,.5,1)<=.05)  is about .22
      # which is more powerful than uniform but still low power

      # TODO: we should actually be reducing power pbeta() etc. and then
      # stopping because in real life there are only so many blocks.

      cap <- .5
      effective_beta <- min(beta_eff_raw, cap)
    } else {
      effective_beta <- beta_base
    }
    # For each child, retrieve its parent's p-value.
    parents <- tree_sim[, .(parent = node, parent_p = p_val)]
    setkey(parents, "parent")
    setkey(child_rows, "parent")
    child_rows <- parents[child_rows]
    ## child_rows <- merge(child_rows, parents, by = "parent", all.x = TRUE, sort = FALSE)

    # draw local p-values. the beta p-values are rescaled to go from parent_p to 1.

    ## Don't replace the p_sim already drawn for the leaves but we will enforce
    ## monotonicity on those p-values below

    child_rows[is.na(p_sim), p_sim := fifelse(
      nonnull,
      parent_p + (1 - parent_p) * rbeta(.N, effective_beta, 1),
      runif(.N, min = parent_p, max = 1)
    )]

    ## add this column for debugging and testing for now.
    ## child_rows[, effective_beta := effective_beta]

    # Enforce global monotonicity.
    child_rows[, p_val := pmax(parent_p, p_sim)]

    # Compute the locally adjusted p-value for each group of children sharing the same parent.
    ## child_rows[, local_adj_p := {
    ##   m <- .N
    ##   sorted <- sort(p_val)
    ##   min(sorted * m / seq_len(m))
    ## }, by = parent]

    child_rows[, local_adj_p := local_adj_p_fn(p_val), by = parent]

    # Gate the branch: if local_adj_p p-value > alpha, do not continue testing the children.
    child_rows[local_adj_p > alpha, p_val := NA_real_]

    ## Reduce chance of false rejections later by further inflating/adjusting the local p-values
    child_rows[local_adj_p <= alpha, p_val := pmax(local_adj_p, p_val)]

    # Update the main table.
    ## tree_sim[node %in% child_rows$node, p_val := child_rows$p_val]
    tree_sim[child_rows$node, p_val := child_rows$p_val]
  }

  # Record a false error if any truly null node that was tested (i.e. has a non-NA p-value)
  # produces a p-value <= alpha.
  false_error <- any(tree_sim[nonnull == FALSE & !is.na(p_val), p_val] <= alpha)

  ## Record number of true discoveries --- non-null nodes tested and with p <= alpha
  ## Notice that a discovery on a node means that at least one of its leaves is nonnull
  true_discoveries <- sum(tree_sim[nonnull == TRUE & !is.na(p_val), p_val] <= alpha)

  ## Record number of nodes tested
  num_nodes_tested <- sum(!is.na(tree_sim$p_val))

  ## Collect results about the tests into a data.table
  sim_res <- data.table(
    false_error = false_error, true_discoveries = true_discoveries, num_nodes_tested = num_nodes_tested,
    bottom_up_false_error = bottom_up_false_error,
    bottom_up_true_discoveries = bottom_up_true_discoveries, num_leaves = num_leaves
  )

  if (return_details) {
    return(list(
      treeDT = tree_sim, sim_res = sim_res
    ))
  } else {
    return(sim_res)
  }
}

#' Simulate the Family-Wise Error Rate (FWER) using a Data.table-based Tree
#'
#' @description
#' Runs the hierarchical testing procedure on complete k-ary trees (generated via \code{generate_tree_DT()})
#' \code{n_sim} times and returns the estimated FWER. For each simulated tree the testing procedure is applied
#' using \code{simulate_test_DT()}.
#'
#' @param n_sim Integer. The number of simulation replicates.
#' @param t Numeric in [0,1]. The probability that a leaf is non-null.
#' @param k Integer. The branching factor.
#' @param max_level Integer. The maximum level (depth) of the tree.
#' @param alpha Numeric. The significance level.
#' @param N_total Numeric. The total sample size at the root.
#' @param beta_base Numeric. The base parameter for the Beta distribution at the root.
#' @param return_details Logical. For now this must be FALSE.
#'
#' @return A numeric value in [0,1] giving the estimated FWER.
#'
#' @examples
#' fwer <- simulate_FWER_DT(
#'   n_sim = 1000, t = 0.2, k = 3, max_level = 3,
#'   alpha = 0.05, N_total = 1000, beta_base = 0.1
#' )
#' print(fwer)
#'
#' @export
simulate_many_runs_DT <- function(n_sim, t, k, max_level, alpha, N_total, beta_base = 0.1,
                                  adj_effN = TRUE, local_adj_p_fn = local_simes, return_details = FALSE,
                                  global_adj = "hommel") {
  treeDT <- generate_tree_DT(max_level, k, t)
  res <- replicate(
    n_sim,
    simulate_test_DT(treeDT, alpha, k,
      effN = N_total, N_total = N_total, beta_base = beta_base,
      adj_effN = adj_effN, local_adj_p_fn = local_adj_p_fn, global_adj = global_adj, return_details = return_details
    ),
    simplify = FALSE
  )

  res_dt <- rbindlist(res)
  mean_res <- unlist(res_dt[, lapply(.SD, mean)])

  return(mean_res)
}
