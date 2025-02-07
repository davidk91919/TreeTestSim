#' @title Compute local Simes p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the Simes p-value
#'   \eqn{\min_{i=1\ldots k} \{ (k/i) * p_{(i)} \}}, where \eqn{p_{(1)} \le \ldots \le p_{(k)}}.
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used in this function)
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
local_simes <- function(pvals_children, alpha = .05) {
  k <- length(pvals_children)
  sort_p <- sort(pvals_children)
  i_seq <- seq_len(k)
  simes_vals <- (k / i_seq) * sort_p
  min(simes_vals)
}
## NOTE maybe TODO: we could use hommel(pvals_children,simes=TRUE) from the hommel package here too.

#' @title Unadjusted local step
#'
#' @description Given \eqn{k} child p-values, return the highest p-value below alpha. This is the same as respecting the global gating criteria of not continuing to test along any path where p > alpha, but doesn't do a local, parent-level adjustment. This is mostly useful for simulation and study purposes although with k=2 and low l  (for example, less than 5 or so) this will also tend to control the FWER.
#'
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha
#'
#' @details The idea is to not adjust so that we can compare the consequences of not adjusting to adjusting
#'
#' @return A single numeric value: the highest p<=alpha among the children of the parent node or the minimum of the p-values of the children if they are all greater than alpha
#'
#' @examples
#' local_min_p(c(0.01, 0.04, 0.10, 0.20)) ## should return .04
#' local_min_p(c(0.10, 0.20)) ## should return .04
#'
#' @export
local_min_p <- function(pvals_children, alpha = .05) {
  p_le_alpha <- pvals_children[pvals_children <= alpha]
  if (length(p_le_alpha) == 0) {
    return(min(pvals_children))
  } else {
    return(max(p_le_alpha))
  }
}
