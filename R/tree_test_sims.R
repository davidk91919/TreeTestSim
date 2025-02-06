#' @title Compute Level and Offset Information for a k-ary Tree
#'
#' @description Given a branching factor \code{k} and a depth \code{l}, this function
#'   computes the number of nodes at each level and the starting offsets for a
#'   level-by-level indexing scheme in a single vector.
#'
#' @param k Integer; the branching factor (i.e., each internal node has \code{k} children). Must be greater than 0.
#' @param l Integer; the depth of the tree (0-based). If \code{l=0}, there is only a root. If l=1, this means a root and one level below the root.
#'
#' @return A list with the elements:
#' \describe{
#'   \item{\code{k}{k}}
#'   \item{\code{l}{l}}
#'   \item{\code{level_sizes}}{An integer vector of length \code{l+1}, where \code{level_sizes[d+1] = k^d}.}
#'   \item{\code{level_offsets}}{An integer vector of length \code{l+1}, where
#'       \code{level_offsets[d+1]} is the total number of nodes up to (but not including) level \code{d}. Used for creating indices.}
#'   \item{\code{n_tot}}{The total number of nodes in the tree, \eqn{\sum_{d=0}^{l} k^d}.}
#' }
#'
#' @examples
#' # For k=2, l=3, we have a binary tree of depth 3
#' info <- get_level_info(k = 2, l = 3)
#' info$level_sizes # c(1, 2, 4, 8)
#' info$level_offsets # c(0, 1, 3, 7)
#' # Index of root is 0 in a 0-based system (we go between 0-based for mathematical clarity and 1-based for R
#' # Here we show both:
#' # for 1-based indexing the next level l=1 is level_offsets[2]+c(1,2)=c(2,3)
#' # for 0-based indexing the next level l=1 is level_offsets[2]+c(0,1)=c(1,2)
#' # For the next level, l=2 (or the third level in the tree if you are counting the root as a level)
#' # for the first node in level 1 the children in level 2 using: 1-based indexing level_offsets[3] + c(1,2)=c(4,5)
#' # for the second node in level 1 the children in level 2  1-based indexing level_offsets[3] + c(k*(2-1)+1,k*(2-1)+2) = 3+c(3,4)=c(6,7)
#' # for the first node in level 1 the children in level 2 0-based indexing level_offsets[3] + seq(k*j,k*j+(k-1))=3+c(0,1)=c(3,4)
#' # for the second node in level 1 the children in level 2 0-based indexing level_offsets[3] + seq(k*j,k*j + (k-1))=3+c(2,3)=c(5,6)
#' info$n_tot # 15 total nodes
#'
#' @export
get_level_info <- function(k, l) {
  stopifnot(k > 0)
  level_sizes <- k^(0:l)
  if (l == 0) {
    level_offsets <- 0
  } else {
    level_offsets <- cumsum(c(0, level_sizes[1:l]))
  }
  n_tot <- sum(level_sizes)
  list(
    k = k,
    l = l,
    level_sizes = level_sizes,
    level_offsets = level_offsets,
    n_tot = n_tot
  )
}

#' @title Get Children of a Node in a k-ary Tree
#'
#' @description Given a node's level \code{d} and local index \code{j} in a
#' k-ary tree, return the global indices of its \code{k} children (if \code{d <
#' l}) using a 1-based indexing system given the way R indexes. When d=0, j must be 1 since we only have one root node.
#'
#' @param d Integer; the current level or depth of the node (0-based).
#' @param j Integer; the node's local index within level \code{d} (1-based) -- a function of k (number of nodes per level)
#' @param tree_info A list such as the one returned by \code{\link{get_level_info}}
#'
#' @details We assume a level-by-level storage scheme:
#' \itemize{
#' \item Level \eqn{d} has \eqn{k^d} nodes.
#' \item The children of local node j in level d, \eqn{(d, j)}, are those with local indices \eqn{k(j-1)+1, \ldots, k(j-1)+k} in level
#' \eqn{d+1}.
#' \item The function returns an empty integer vector if \eqn{d=l} (a leaf). }
#'
#' @return An integer vector of global indices of the children.
#'
#' @seealso \code{\link{get_level_info}}
#'
#' @examples
#' info <- get_level_info(k = 2, l = 2)
#' # Root is at level d=0, j=1 => children are global indices 2 and 3
#' children_indices(d = 0, j = 1, tree_info = info)
#' # Next look at level 1 where we have 2 nodes # d=1, j=1 => global index node 2 => children are 4 and 5
#' children_indices(d = 1, j = 1, tree_info = info)
#' # d=1, j=2 ==> the global indx node 3 => children are nodes 6 and 7
#' children_indices(d = 1, j = 2, tree_info = info)

#' @export
children_indices <- function(d, j, tree_info) {
  level_offsets <- tree_info$level_offsets
  k <- tree_info$k

  l <- length(level_offsets) - 1
  if (d == l) {
    return(integer(0))
  } # no children at the leaf level
  if (d == 0) {
    stopifnot("There should only be one node at the root (d=0 requires j=1)" = j == 1)
  }
  stopifnot("There are no nodes labeled as zero" = j != 0)
  # stopifnot("The local node number must be between 1 and k"=(j>0 & j<=k))
  start_child <- k * (j - 1) + 1
  end_child <- k * (j - 1) + k
  child_locals <- seq.int(start_child, end_child)

  idx_global <- function(dd, localJ) {
    level_offsets[dd + 1] + localJ
  }

  child_globals <- idx_global(d + 1, child_locals)
  child_globals
}


#' @title Return a vector of all ancestors of a node
#'
#' @description Given a 1-based global node id and information about the number of nodes at a level, return the global indices of all ancestors of that node.
#'
#' @param node_id Integer; The global node id
#' @param k Integer; The number of nodes at each level of the k-ary tree
#'
#' @return a vector with the global ids of the ancestors of the node (always includes node 1 or the root node)
#'
#' @examples
#' ## This next returns c(3,1)
#' get_ancestors(node_id = 8, k = 3)
#' ## This next returns c(4,2,1)
#' get_ancestors(node_id = 8, k = 2)
#'
#' @export
get_ancestors <- function(node_id, k) {
  if (node_id < 2) {
    return(integer(0))
  } # the root has no ancestors
  # Set an upper bound on the number of steps needed.
  max_exp <- ceiling(log(node_id - 1, base = k)) + 1
  # Compute ancestors using vectorized arithmetic:
  # For exponent i, the ancestor is floor((node_id - 2)/k^i) + 1.
  anc <- floor((node_id - 2) / k^(1:max_exp)) + 1
  # The sequence will eventually hit the root (1) and then remain 1.
  # We take ancestors up to the first occurrence of 1.
  last <- which(anc == 1)[1]
  return(anc[1:last])
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
assign_alt <- function(t, tree_info) {
  k <- tree_info$k
  l <- tree_info$l
  level_offsets <- tree_info$level_offsets
  level_sizes <- tree_info$level_sizes
  n_tot <- tree_info$n_tot

  ## Initialize the vector indicating whether a node has a non-null efffect as FALSE
  alt <- rep(FALSE, n_tot)

  # Identify leaves
  ## This might be slightly slower than
  ## seq.int(level_offsets[l+1],level_offsets[l+1]+k^l)
  leaf_start <- level_offsets[l + 1] + 1
  leaf_end <- level_offsets[l + 1] + level_sizes[l + 1]
  leaf_inds <- seq.int(leaf_start, leaf_end)
  ## Stop this test once we are satified
  leaf_inds2 <- seq.int(level_offsets[l + 1] + 1, level_offsets[l + 1] + k^l)
  stopifnot(all.equal(leaf_inds, leaf_inds2))

  # Complete Random assignment at leaves
  ## alt[leaf_inds] <- (runif(length(leaf_inds)) < t)
  n_leaves <- length(leaf_inds)
  m_leaves <- floor(t * n_leaves)
  alt[leaf_inds] <- sample(rep(c(FALSE, TRUE), c(n_leaves - m_leaves, m_leaves)))

  ## Assign all ancestors of each TRUE leaf node to be TRUE
  for (i in leaf_inds[alt[leaf_inds]]) {
    # TODO: this could be made more efficient maybe since some leafs will share exactly the same ancestors
    # but I'm not sure how much more overhead would go into checking equality of the vectors or the math of get_ancestors
    leaf_ancestors <- get_ancestors(node_id = i, k = k)
    alt[leaf_ancestors] <- TRUE
  }

  alt
}

#' @title Draw a Single Node's P-value under Monotonicity
#'
#' @description Generates a final p-value for a node given its parent's final p-value,
#'   respecting the monotonic rule \eqn{p_\text{child} \ge p_\text{parent}}.
#'
#' @param is_alt Logical; \code{TRUE} if the node is non-null, \code{FALSE} if null.
#' @param parent_p Numeric in \eqn{[0,1]}; the parent's final p-value.
#' @param beta_params Numeric vector of length 2, e.g. \code{c(a, b)}, specifying
#'   the Beta(\code{a},\code{b}) distribution for alternative nodes. The actual
#'   p-value is then scaled to \eqn{[parent\_p, 1]}.
#'
#' @details
#' If \code{is_alt=FALSE}, we draw \eqn{\mathrm{Uniform}(parent\_p, 1)}.
#' If \code{is_alt=TRUE}, we draw \eqn{\mathrm{Beta}(a,b)} and linearly
#'   rescale the result into \eqn{[parent\_p, 1]}.
#'
#' @return A numeric p-value in \eqn{[parent\_p, 1]}.
#'
#' @examples
#' # Null node p-value if parent_p=0.05 all should be greater than or equal to .05
#' draw_node_p_value(FALSE, parent_p = 0.05, beta_params = c(0.1, 1))
#'
#' # Alt node p-value if parent_p=0.05 all should be greater than or equal to .05
#' draw_node_p_value(TRUE, parent_p = 0.05, beta_params = c(0.1, 1))
#'
#' @export
draw_node_p_value <- function(is_alt, parent_p, beta_params) {
  if (is_alt) {
    raw_beta_p <- rbeta(1, beta_params[1], beta_params[2])
    parent_p + (1 - parent_p) * raw_beta_p
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
#' TODO maybe this should just apply the local Simes gate rather than do the whole thing?
#'
#'  @param tree_info; a list object like that created by \code{\link{get_level_info}}.
#' @param alpha Numeric in \eqn{(0,1)}; the significance level threshold.
#' @param alt Logical vector of length \eqn{\sum_{d=0}^{l} k^d}, indicating which nodes
#'   are non-null (\code{TRUE}) vs. null (\code{FALSE}).
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
#'   \eqn{\le \alpha} (false rejection) in this single run, \code{FALSE} otherwise. TODO maybe track FDR later.
#'
#' @seealso \code{\link{simulate_hier_simes_local_modular}}
#'
#' @examples
#' info <- get_level_info(k = 2, l = 2)
#' alt_vec <- assign_alt(t = 0.5, tree_info = info)
#' simulate_single_run(tree_info = info, alpha = 0.05, alt = alt_vec, beta_params = c(0.1, 1))
#'
#' @export
simulate_single_run <- function(tree_info, alpha, alt, beta_params, local_adjust_fn = local_simes) {
  k <- tree_info$k
  l <- tree_info$l
  level_offsets <- tree_info$level_offsets
  level_sizes <- tree_info$level_sizes
  n_tot <- tree_info$n_tot

  ## Initialize some variables
  node_idx <- seq(1, n_tot)
  pvals <- rep(NA_real_, n_tot)
  tested <- rep(FALSE, n_tot)
  false_reject <- FALSE

  # Root
  root_idx <- 1
  pvals[root_idx] <- draw_node_p_value(is_alt = alt[root_idx], parent_p = 0, beta_params = beta_params)
  tested[root_idx] <- TRUE
  ## If the root is FALSE and p <= alpha, then we have a false_rejection
  if (!alt[root_idx] && pvals[root_idx] <= alpha) {
    false_reject <- TRUE
    ## We do not have to continue the simulation since FWER is **any** false rejection.
    return(false_reject)
  }

  ## If didn't reject the root but it was a non-null effect, this is a power
  ## problem but not a false rejection. Regardless, we cannot continue testing.
  if (!alt[root_idx] && pvals[root_idx] > alpha) {
    false_reject <- FALSE
    return(false_reject)
  }

  # Level 1
  idx_l1 <- children_indices(d = 0, j = 1, tree_info = tree_info)
  # Draw p-values for them all
  for (i in idx_l1) {
    pvals[i] <- draw_node_p_value(is_alt = alt[i], parent_p = pvals[root_idx], beta_params = beta_params)
  }
  local_adj_l1 <- local_adjust_fn(pvals[idx_l1])
  ## if local_adj_l1 <= alpha and at least one of the nodes is non-null go to next level otherwise return false_reject <- FALSE
  if (local_adj_l1 > alpha) {
    false_reject <- FALSE
    return(false_reject)
  }

  ## If all of the nodes are FALSE (i.e. null) but we reject using the local
  ## adjustment, then this is a false rejection (at least one false rejection)
  ## If we reject using the local adjustment and at least one of the nodes is
  ## non-null, then this is not a false rejection --- if that node p <= alpha,
  ## we will continue down the tree.

  if (local_adj_l1 <= alpha && all(!alt[idx_l1])) {
    false_reject <- TRUE
    return(false_reject)
  }

  ## If we have a global reject across the three nodes and node level rejection
  ## of a node that should not have been rejected -- alt=FALSE, then this is a
  ## false rejection. No need to continue the simulation to look for others.

  if (local_adj_l1 <= alpha & any(pvals[idx_l1] <= alpha & !alt[idx_l1])) {
    false_reject <- TRUE
    return(false_reject)
  }

  level_idxs <- c(0, rep(seq_len(l), k^seq_len(l)))

  ## Loop over the global indices of the nodes at level d>2
  node_gate <- list()
  for (d in seq(2, l)) {
    parent_level <- d - 1
    parent_idxs_raw <- node_idx[level_idxs == parent_level]
    ## You can only go on with testing if the node p <= alpha
    parent_idxs_non_gated <- parent_idxs_raw[pvals[parent_idxs_raw] <= alpha]
    ## This next is mostly because I'm avoiding names for memory reasons
    ## Not elegant to replace values like this
    node_gate[[d - 1]] <- parent_idxs_non_gated
    ## Loop over the nodes that have pvals <= alpha
    for (i in parent_idxs_non_gated) {
      local_idx <- i - tree_info$level_offsets[d]
      child_idxs <- children_indices(d = parent_level, j = local_idx, tree_info = tree_info)
      for (j in child_idxs) {
        pvals[j] <- draw_node_p_value(is_alt = alt[j], parent_p = pvals[i], beta_params = beta_params)
      }
      local_adj <- local_adjust_fn(pvals[child_idxs])

      ## Assign **all** pvals for these children to be a bit higher than alpha
      ## if we can't pass the local adjustment
      if (local_adj > alpha) {
        pvals[child_idxs] <- alpha + 1e-8
      }
    }
  }
  ## Continue testing for nodes with node_gate==0

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
