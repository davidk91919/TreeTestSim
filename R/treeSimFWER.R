#' Generate a k-ary Tree with Treatment Effects
#'
#' @description
#' Recursively generate a k-ary tree with levels 0 to \code{max_level}. At the leaves,
#' each node is assigned a treatment effect (\code{nonnull} flag) with probability \code{t}.
#' For internal nodes, the \code{nonnull} flag is set to \code{TRUE} if any descendant is non-null.
#'
#' @param level An integer specifying the current level in the tree (starting at 0).
#' @param max_level An integer specifying the maximum depth of the tree.
#' @param k An integer specifying the number of children per node.
#' @param t A numeric value in [0,1] giving the probability that a leaf is non-null (has a treatment effect).
#'
#' @return A list representing the tree with components \code{level}, \code{nonnull}, and \code{children}.
#'
#' @examples
#' tree <- generate_tree(0, max_level = 3, k = 3, t = 0.2)
#' print(tree)
#'
#' @export
generate_tree <- function(level, max_level, k, t) {
  if (level == max_level) {
    nonnull <- runif(1) < t
    return(list(level = level, nonnull = nonnull, children = NULL))
  } else {
    children <- vector("list", k)
    for (i in 1:k) {
      children[[i]] <- generate_tree(level + 1, max_level, k, t)
    }
    # Propagate the nonnull flag upward:
    nonnull <- any(sapply(children, function(child) child$nonnull))
    return(list(level = level, nonnull = nonnull, children = children))
  }
}

#' Simulate the Testing Procedure on a Tree Node
#'
#' @description
#' Recursively simulate the testing procedure along a tree branch. At each node, a p-value is
#' generated using a uniform distribution if the node is null or a rescaled beta distribution if non-null.
#' The function applies a global gate by computing a Simes test over the children.
#'
#' @param node A list representing a node (or subtree) created by \code{generate_tree()}.
#' @param parent_p A numeric value giving the lower bound for the p-value inherited from the parent.
#' @param alpha A numeric value for the significance level.
#' @param k An integer for the number of children per node.
#' @param effN A numeric value for the effective sample size available at the current node.
#' @param N_total A numeric value for the total sample size available at the root.
#' @param beta_base A numeric value for the base parameter of the beta distribution (used at the root).
#'
#' @return A list with components:
#' \describe{
#'   \item{p_val}{The simulated p-value at the current node.}
#'   \item{false_error}{A logical value: \code{TRUE} if a type I error occurred somewhere in the branch, otherwise \code{FALSE}.}
#' }
#'
#' @examples
#' tree <- generate_tree(0, max_level = 3, k = 3, t = 0.2)
#' res <- simulate_test(tree, parent_p = 0, alpha = 0.05, k = 3, effN = 1000, N_total = 1000, beta_base = 0.1)
#' print(res)
#'
#' @export
simulate_test <- function(node, parent_p, alpha, k, effN, N_total, beta_base) {
  if (node$nonnull) {
    effective_beta <- beta_base * sqrt(effN / N_total)
    p_val <- parent_p + (1 - parent_p) * rbeta(1, effective_beta, 1)
  } else {
    p_val <- runif(1, min = parent_p, max = 1)
  }

  false_error_here <- (!node$nonnull && p_val <= alpha)

  # Global gate: if the p-value is > alpha, stop testing further.
  if (p_val > alpha) {
    return(list(false_error = false_error_here, p_val = p_val))
  }

  # If this is a leaf, return the result.
  if (is.null(node$children)) {
    return(list(false_error = false_error_here, p_val = p_val))
  }

  n_children <- length(node$children)
  children_p <- numeric(n_children)
  children_effN <- numeric(n_children)
  child_errors <- logical(n_children)

  # Assume that the effective sample size is split equally.
  new_effN <- effN / k

  for (i in seq_len(n_children)) {
    res <- simulate_test(
      node = node$children[[i]],
      parent_p = p_val, alpha = alpha,
      k = k, effN = new_effN, N_total = N_total, beta_base = beta_base
    )
    children_p[i] <- res$p_val
    child_errors[i] <- res$false_error
  }

  # Compute the Simes p-value from the children.
  ## TODO: here is where we could use hommel(children_p)@adjusted if we want to adjust all of the p-values.
  sorted_p <- sort(children_p)
  simes_vals <- sorted_p * n_children / seq_along(sorted_p)
  simes_p <- min(simes_vals)
  adjusted_p <- max(simes_p, children_p)

  # Gate the branch: if the Simes p-value is > alpha, do not continue testing down the branch.
  if (simes_p > alpha) {
    return(list(false_error = false_error_here, p_val = adjusted_p))
  } else {
    return(list(false_error = (false_error_here || any(child_errors)), p_val = adjusted_p))
  }
}

#' Simulate the Family-Wise Error Rate (FWER)
#'
#' @description
#' Run the hierarchical testing procedure \code{n_sim} times with a given configuration and
#' estimate the family-wise error rate (FWER), defined as the proportion of runs in which a type I error occurred.
#'
#' @param n_sim An integer giving the number of simulation runs.
#' @param t A numeric value in [0,1] for the probability that a leaf is non-null (has a treatment effect).
#' @param k An integer specifying the branching factor (number of children per node).
#' @param max_level An integer for the maximum depth of the tree (levels 0 to \code{max_level}).
#' @param alpha A numeric value for the significance level.
#' @param N_total A numeric value for the total sample size available at the root.
#' @param beta_base A numeric value for the base parameter of the beta distribution at the root (default is 0.1).
#'
#' @return A numeric value between 0 and 1 representing the estimated FWER.
#'
#' @examples
#' fwer <- simulate_FWER(n_sim = 10000, t = 0.2, k = 3, max_level = 3, alpha = 0.05, N_total = 1000)
#' print(fwer)
#'
#' @export
simulate_FWER <- function(n_sim, t, k, max_level, alpha, N_total, beta_base = 0.1) {
  false_count <- numeric(length = n_sim)
  tree <- generate_tree(level = 0, max_level = max_level, k = k, t = t)
  for (sim in 1:n_sim) {
    res <- simulate_test(tree,
      parent_p = 0, alpha = alpha, k = k,
      effN = N_total, N_total = N_total, beta_base = beta_base
    )
    false_count[sim] <- res$false_error
  }
  return(mean(false_count))
}
