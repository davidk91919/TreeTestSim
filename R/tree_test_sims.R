
#' @title Compute Level and Offset Information for a k-ary Tree
#'
#' @description Given a branching factor \code{k} and a depth \code{l}, this function
#'   computes the number of nodes at each level and the starting offsets for a
#'   level-by-level indexing scheme in a single vector.
#'
#' @param k Integer; the branching factor (i.e., each internal node has \code{k} children).
#' @param l Integer; the depth of the tree (0-based). If \code{l=0}, there is only a root.
#'
#' @return A list with the elements:
#' \describe{
#'   \item{\code{level_sizes}}{An integer vector of length \code{l+1}, where \code{level_sizes[d+1] = k^d}.}
#'   \item{\code{level_offsets}}{An integer vector of length \code{l+1}, where
#'       \code{level_offsets[d+1]} is the total number of nodes up to (but not including) level \code{d}.}
#'   \item{\code{n_tot}}{The total number of nodes in the tree, \eqn{\sum_{d=0}^{l} k^d}.}
#' }
#'
#' @examples
#' # For k=2, l=3, we have a binary tree of depth 3
#' info <- get_level_info(k = 2, l = 3)
#' info$level_sizes # c(1, 2, 4, 8)
#' info$level_offsets # c(0, 1, 3, 7)
#' info$n_tot # 15 total nodes
#' 
#' @export
get_level_info <- function(k, l) {
  level_sizes <- k^(0:l)
  level_offsets <- cumsum(c(0, level_sizes[1:l]))
  n_tot <- sum(level_sizes)
  list(
    level_sizes = level_sizes,
    level_offsets = level_offsets,
    n_tot = n_tot
  )
}

#' @title Get Children of a Node in a k-ary Tree
#'
#' @description Given a node's level \code{d} and local index \code{j} in a k-ary tree,
#'   return the global indices of its \code{k} children (if \code{d < l}).
#'
#' @param d Integer; the current level of the node (0-based).
#' @param j Integer; the node's local index within level \code{d} (1-based).
#' @param k Integer; the branching factor of the tree.
#' @param level_offsets Integer vector of length \code{l+1}, as returned by \code{\link{get_level_info}}.
#'
#' @details We assume a level-by-level storage scheme:
#'   \itemize{
#'     \item Level \eqn{d} has \eqn{k^d} nodes.
#'     \item The children of node \eqn{(d, j)} are those with local indices
#'           \eqn{k(j-1)+1, \ldots, k(j-1)+k} in level \eqn{d+1}.
#'     \item The function returns an empty integer vector if \eqn{d=l} (a leaf).
#'   }
#'
#' @return An integer vector of global indices of the children.
#'
#' @seealso \code{\link{get_level_info}}
#'
#' @examples
#' info <- get_level_info(k = 2, l = 3)
#' # Root is at level d=0, j=1 => children are global indices 2 and 3
#' children_indices(d = 0, j = 1, k = 2, level_offsets = info$level_offsets)
#' # Next look at level 1 where we have 2 nodes
#' # d=1, j=1 => global index node 2 => children are 4 and 5
#' children_indices(d = 1, j = 1, k = 2, level_offsets = info$level_offsets)
#' # d=1, j=2 ==> the global indx node 3 => children are nodes 6 and 7
#' children_indices(d = 1, j = 2, k = 2, level_offsets = info$level_offsets)

#' @export
children_indices <- function(d, j, k, level_offsets) {
  l <- length(level_offsets) - 1
  if (d == l) {
    return(integer(0))
  } # no children at the leaf level

  start_child <- k * (j - 1) + 1
  end_child <- k * (j - 1) + k
  child_locals <- seq.int(start_child, end_child)

  idx_global <- function(dd, localJ) level_offsets[dd + 1] + localJ

  child_globals <- idx_global(d + 1, child_locals)
  child_globals
}

#' @title Assign Null/Alternative Status in a k-ary Tree
#'
#' @description Randomly marks each leaf node as alternative (non-null) with probability
#'   \code{t}, then propagates that status upward so that all ancestors of a non-null leaf
#'   become non-null as well.
#'
#' @param k Integer; the branching factor.
#' @param l Integer; the depth of the tree (0-based).
#' @param t Numeric in \eqn{[0, 1]}; the probability that a leaf is non-null.
#' @param level_offsets Integer vector, part of the output from \code{\link{get_level_info}}.
#' @param level_sizes Integer vector, part of the output from \code{\link{get_level_info}}.
#'
#' @return A logical vector of length \code{sum(level_sizes)}, indicating for each node
#'   whether it is alternative (\code{TRUE}) or null (\code{FALSE}).
#'
#' @examples
#' info <- get_level_info(k = 2, l = 2)
#' alt_vec <- assign_alt(
#'   k = 2, l = 2, t = 0.5,
#'   level_offsets = info$level_offsets, level_sizes = info$level_sizes
#' )
#' alt_vec
#'
#' @export
assign_alt <- function(k, l, t, level_offsets, level_sizes) {
  n_tot <- sum(level_sizes)
  alt <- rep(FALSE, n_tot)

  # Identify leaves
  leaf_start <- level_offsets[l + 1] + 1
  leaf_end <- level_offsets[l + 1] + level_sizes[l + 1]
  leaf_inds <- seq(leaf_start, leaf_end)

  # Random assignment at leaves
  alt[leaf_inds] <- (runif(length(leaf_inds)) < t)

  # Propagate alt up
  for (dd in seq(l, 1, by = -1)) {
    offset_d <- level_offsets[dd + 1]
    n_nodes_d <- level_sizes[dd + 1]
    for (j in seq_len(n_nodes_d)) {
      child_global <- offset_d + j
      if (alt[child_global]) {
        # Mark parent alt
        parent_local <- ceiling(j / k)
        parent_global <- level_offsets[dd] + parent_local
        alt[parent_global] <- TRUE
      }
    }
  }

  alt
}

#' @title Draw a Single Node's P-value under Monotonicity
#'
#' @description Generates a final p-value for a node given its parent's final p-value,
#'   respecting the monotonic rule \eqn{p_\text{child} \ge p_\text{parent}}.
#'
#' @param isAlt Logical; \code{TRUE} if the node is non-null, \code{FALSE} if null.
#' @param parent_p Numeric in \eqn{[0,1]}; the parent's final p-value.
#' @param beta_params Numeric vector of length 2, e.g. \code{c(a, b)}, specifying
#'   the Beta(\code{a},\code{b}) distribution for alternative nodes. The actual
#'   p-value is then scaled to \eqn{[parent\_p, 1]}.
#'
#' @details
#' If \code{isAlt=FALSE}, we draw \eqn{\mathrm{Uniform}(parent\_p, 1)}.
#' If \code{isAlt=TRUE}, we draw \eqn{\mathrm{Beta}(a,b)} and linearly
#'   rescale the result into \eqn{[parent\_p, 1]}.
#'
#' @return A numeric p-value in \eqn{[parent\_p, 1]}.
#'
#' @examples
#' # Null node p-value if parent_p=0.3
#' draw_node_p_value(FALSE, parent_p = 0.3, beta_params = c(0.1, 1))
#'
#' # Alt node p-value if parent_p=0.3
#' draw_node_p_value(TRUE, parent_p = 0.3, beta_params = c(0.1, 1))
#'
#' @export
draw_node_p_value <- function(isAlt, parent_p, beta_params) {
  if (isAlt) {
    rawBeta <- rbeta(1, beta_params[1], beta_params[2])
    parent_p + (1 - parent_p) * rawBeta
  } else {
    runif(1, min = parent_p, max = 1)
  }
}

#' @title Compute local Simes p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the Simes p-value
#'   \eqn{\min_{i=1\ldots k} \{ (k/i) * p_{(i)} \}}, where \eqn{p_{(1)} \le \ldots \le p_{(k)}}.
#'
#' @param pvals_children Numeric vector of child p-values.
#'
#' @details The Simes p-value is a valid test for the intersection hypothesis
#'   \eqn{H_{\cap} : \text{all child hypotheses are null}}, under certain independence
#'   or positive dependence conditions among the p-values.
#'
#' @return A single numeric value: the Simes combination p-value.
#'
#' @examples
#' local_simes(c(0.01, 0.04, 0.10, 0.20))
#'
#' @export
local_simes <- function(pvals_children) {
  k <- length(pvals_children)
  sort_p <- sort(pvals_children)
  i_seq <- seq_len(k)
  simes_vals <- (k / i_seq) * sort_p
  min(simes_vals)
}
## NOTE maybe TODO: we could use hommel(pvals_children,simes=TRUE) from the hommel package here too.

#' @title Simulate a Single Top-Down Run of the Hierarchical Simes Procedure
#'
#' @description Performs one realization of the top-down testing with strict gating
#'   (\eqn{p > \alpha} stops further testing) plus a local Simes gate for each parent.
#'
#' @param k Integer; branching factor.
#' @param l Integer; depth of the tree (0-based).
#' @param alpha Numeric in \eqn{(0,1)}; the significance level threshold.
#' @param alt Logical vector of length \eqn{\sum_{d=0}^{l} k^d}, indicating which nodes
#'   are non-null (\code{TRUE}) vs. null (\code{FALSE}).
#' @param level_offsets Integer vector, part of \code{\link{get_level_info}} output.
#' @param level_sizes Integer vector, part of \code{\link{get_level_info}} output.
#' @param beta_params Numeric vector of length 2 for the Beta distribution
#'   parameters of alternative nodes (see \code{\link{draw_node_p_value}}).
#' @param local_adj_fn A function that takes the child p-values of a node and returns a single value (like the min or max p-value) that can be used to decide on a local gating procedure. For now the default is the \code{local_simes} function.
#' @details
#' \enumerate{
#'   \item The root node (level 0) is tested unconditionally:
#'         if \code{alt[1] = FALSE}, its p-value is \eqn{\mathrm{Unif}(0,1)};
#'         if \code{alt[1] = TRUE}, it is a scaled Beta(\code{beta_params}).
#'   \item If a node's final p-value \eqn{> \alpha}, we do not generate child p-values
#'         (gating).
#'   \item If \eqn{\le \alpha}, we generate the child's p-values.  Then we apply
#'         a \code{\link{local_simes}} test on those \eqn{k} child p-values.
#'         If Simes \eqn{> \alpha}, we forcibly set those child p-values above \eqn{\alpha}
#'         to block deeper testing.
#'   \item We track any null node with p-value \eqn{\le \alpha}; that is counted as
#'         a "false rejection."
#' }
#'
#' @return A logical scalar: \code{TRUE} if at least one null node was
#'   \eqn{\le \alpha} (false rejection) in this single run, \code{FALSE} otherwise.
#'
#' @seealso \code{\link{simulate_hier_simes_local_modular}}
#'
#' @examples
#' info <- get_level_info(k = 2, l = 2)
#' alt_vec <- assign_alt(
#'   k = 2, l = 2, t = 0.5,
#'   level_offsets = info$level_offsets, level_sizes = info$level_sizes
#' )
#' simulate_single_run(
#'   k = 2, l = 2, alpha = 0.05, alt_vec,
#'   level_offsets = info$level_offsets, level_sizes = info$level_sizes, beta_params = c(0.1, 1)
#' )
#'
#' @export
simulate_single_run <- function(k, l, alpha, alt, level_offsets, level_sizes, beta_params, local_adjust_fn = local_simes) {
  n_tot <- sum(level_sizes)
  pvals <- rep(NA_real_, n_tot)
  tested <- rep(FALSE, n_tot)

  false_reject <- FALSE

  # Root
  root_idx <- 1
  pvals[root_idx] <- draw_node_p_value(alt[root_idx], 0, beta_params)
  tested[root_idx] <- TRUE
  if (!alt[root_idx] && pvals[root_idx] <= alpha) {
    false_reject <- TRUE
  }

  # Traverse level by level
  for (d in seq_len(l)) {
    start_d <- level_offsets[d + 1] + 1
    end_d <- level_offsets[d + 1] + level_sizes[d + 1]
    nodeRange <- seq(start_d, end_d)

    # Parents that are tested and p <= alpha
    tested_parents <- nodeRange[tested[nodeRange] & (pvals[nodeRange] <= alpha)]
    if (!length(tested_parents)) next

    for (parent_idx in tested_parents) {
      # local index j for the parent at level d
      j <- parent_idx - level_offsets[d + 1]

      # identify children
      child_inds <- children_indices(d, j, k, level_offsets)
      if (!length(child_inds)) next

      # draw child p-values
      parent_p <- pvals[parent_idx]
      child_pvals <- numeric(k)
      for (iC in seq_len(k)) {
        c_idx <- child_inds[iC]
        child_pvals[iC] <- draw_node_p_value(alt[c_idx], parent_p, beta_params)
      }
      pvals[child_inds] <- child_pvals
      tested[child_inds] <- TRUE

      # check false rejections among children
      null_children <- (!alt[child_inds])
      if (any(null_children & (child_pvals <= alpha))) {
        false_reject <- TRUE
      }

      # local Simes test
      # simes_p <- local_simes(child_pvals)
      local_p <- local_adjust_fn(child_pvals)
      if (local_p > alpha) {
        # block deeper testing from these children
        pvals[child_inds] <- alpha + 1e-8
      }
    }
  }

  false_reject
}

#' @title Monte Carlo Simulation of a Hierarchical Simes-Gated Testing Procedure
#'
#' @description Repeats \code{n_sim} independent realizations of the following steps:
#'   \enumerate{
#'     \item Randomly assign half the leaves as alt with probability \code{t},
#'           and propagate alt status upward.
#'     \item Perform a top-down testing procedure with:
#'       \itemize{
#'         \item strict gating: if a node's p-value \eqn{> \alpha}, do not test children
#'         \item local Simes: among each parent's \eqn{k} children, if the Simes p-value
#'               \eqn{> \alpha}, block deeper descendants
#'         \item monotonic p-values: null => \eqn{\mathrm{Uniform}(\text{parent\_p},1)},
#'               alt => scaled \eqn{\mathrm{Beta}(a,b)} to \eqn{[\text{parent\_p},1]}
#'       \item track any null node \eqn{\le \alpha} as a false rejection
#'       }
#'     \item Estimate the familywise error rate (FWER) as
#'           the fraction of runs in which at least one null node was \eqn{\le \alpha}.
#'   }
#'
#' @param n_sim Number of simulation replicates.
#' @param k Integer; branching factor.
#' @param l Integer; depth of the tree.
#' @param t Numeric in \eqn{[0,1]}; probability that a leaf is non-null.
#' @param alpha Numeric in \eqn{(0,1)}; significance level threshold.
#' @param beta_params Numeric vector of length 2, \code{c(a,b)}, specifying Beta(a,b)
#'   for alt nodes. The final alt p-value is scaled to be in \eqn{[\text{parent\_p},1]}.
#'
#' @return A numeric estimate of the familywise error rate:
#'   the proportion of runs that had at least one null node p-value \eqn{\le \alpha}.
#'
#' @seealso \code{\link{simulate_single_run}}, \code{\link{assign_alt}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' fwer_est <- simulate_hier_simes_local_modular(
#'   n_sim = 2000,
#'   k = 2,
#'   l = 3,
#'   t = 0.5,
#'   alpha = 0.05,
#'   beta_params = c(0.1, 1)
#' )
#' fwer_est
#' }
#'
#' @export
simulate_hier_simes_local_modular <- function(n_sim = 10000,
                                              k = 2,
                                              l = 3,
                                              t = 0.5,
                                              alpha = 0.05,
                                              beta_params = c(0.1, 1),
                                              adjfn = local_simes) {
  info <- get_level_info(k, l)
  off <- info$level_offsets
  sizes <- info$level_sizes

  false_count <- 0

  for (rep in seq_len(n_sim)) {
    # 1) alt assignment
    alt_vec <- assign_alt(k, l, t, off, sizes)

    # 2) single run
    had_false_rej <- simulate_single_run(k, l, alpha, alt_vec, off, sizes, beta_params)

    if (had_false_rej) {
      false_count <- false_count + 1
    }
  }

  # fraction that had false rejection => FWER estimate
  fwer <- false_count / n_sim
  fwer
}
