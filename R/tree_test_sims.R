# hier_simes_modular.R

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
#'   \item{\code{levelSizes}}{An integer vector of length \code{l+1}, where \code{levelSizes[d+1] = k^d}.}
#'   \item{\code{levelOffsets}}{An integer vector of length \code{l+1}, where
#'       \code{levelOffsets[d+1]} is the total number of nodes up to (but not including) level \code{d}.}
#'   \item{\code{nTot}}{The total number of nodes in the tree, \eqn{\sum_{d=0}^{l} k^d}.}
#' }
#'
#' @examples
#' # For k=2, l=3, we have a binary tree of depth 3
#' info <- get_level_info(k = 2, l = 3)
#' info$levelSizes # c(1, 2, 4, 8)
#' info$levelOffsets # c(0, 1, 3, 7)
#' info$nTot # 15 total nodes
#'
#' @export
get_level_info <- function(k, l) {
  levelSizes <- k^(0:l)
  levelOffsets <- cumsum(c(0, levelSizes[1:l]))
  nTot <- sum(levelSizes)
  list(
    levelSizes   = levelSizes,
    levelOffsets = levelOffsets,
    nTot         = nTot
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
#' @param levelOffsets Integer vector of length \code{l+1}, as returned by \code{\link{get_level_info}}.
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
#' info <- get_level_info(k = 2, l = 2)
#' # Root is at level d=0, j=1 => children are global indices 2 and 3
#' children_indices(d = 0, j = 1, k = 2, levelOffsets = info$levelOffsets)
#'
#' @export
children_indices <- function(d, j, k, levelOffsets) {
  l <- length(levelOffsets) - 1
  if (d == l) {
    return(integer(0))
  } # no children at the leaf level

  startChild <- k * (j - 1) + 1
  endChild <- k * (j - 1) + k
  childLocals <- seq.int(startChild, endChild)

  idx_global <- function(dd, localJ) levelOffsets[dd + 1] + localJ

  childGlobals <- idx_global(d + 1, childLocals)
  childGlobals
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
#' @param levelOffsets Integer vector, part of the output from \code{\link{get_level_info}}.
#' @param levelSizes Integer vector, part of the output from \code{\link{get_level_info}}.
#'
#' @return A logical vector of length \code{sum(levelSizes)}, indicating for each node
#'   whether it is alternative (\code{TRUE}) or null (\code{FALSE}).
#'
#' @examples
#' info <- get_level_info(k = 2, l = 2)
#' altVec <- assign_alt(
#'   k = 2, l = 2, t = 0.5,
#'   levelOffsets = info$levelOffsets, levelSizes = info$levelSizes
#' )
#' altVec
#'
#' @export
assign_alt <- function(k, l, t, levelOffsets, levelSizes) {
  nTot <- sum(levelSizes)
  alt <- rep(FALSE, nTot)

  # Identify leaves
  leafStart <- levelOffsets[l + 1] + 1
  leafEnd <- levelOffsets[l + 1] + levelSizes[l + 1]
  leafInds <- seq(leafStart, leafEnd)

  # Random assignment at leaves
  alt[leafInds] <- (runif(length(leafInds)) < t)

  # Propagate alt up
  for (dd in seq(l, 1, by = -1)) {
    offset_d <- levelOffsets[dd + 1]
    nNodes_d <- levelSizes[dd + 1]
    for (j in seq_len(nNodes_d)) {
      childGlobal <- offset_d + j
      if (alt[childGlobal]) {
        # Mark parent alt
        parentLocal <- ceiling(j / k)
        parentGlobal <- levelOffsets[dd] + parentLocal
        alt[parentGlobal] <- TRUE
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
#' @param betaParams Numeric vector of length 2, e.g. \code{c(a, b)}, specifying
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
#' draw_node_p_value(FALSE, parent_p = 0.3, betaParams = c(0.1, 1))
#'
#' # Alt node p-value if parent_p=0.3
#' draw_node_p_value(TRUE, parent_p = 0.3, betaParams = c(0.1, 1))
#'
#' @export
draw_node_p_value <- function(isAlt, parent_p, betaParams) {
  if (isAlt) {
    rawBeta <- rbeta(1, betaParams[1], betaParams[2])
    parent_p + (1 - parent_p) * rawBeta
  } else {
    runif(1, min = parent_p, max = 1)
  }
}

#' @title Compute Local Simes p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the Simes p-value
#'   \eqn{\min_{i=1\ldots k} \{ (k/i) * p_{(i)} \}}, where \eqn{p_{(1)} \le \ldots \le p_{(k)}}.
#'
#' @param pvalsChildren Numeric vector of child p-values.
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
local_simes <- function(pvalsChildren) {
  k <- length(pvalsChildren)
  sortP <- sort(pvalsChildren)
  iSeq <- seq_len(k)
  simesVals <- (k / iSeq) * sortP
  min(simesVals)
}

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
#' @param levelOffsets Integer vector, part of \code{\link{get_level_info}} output.
#' @param levelSizes Integer vector, part of \code{\link{get_level_info}} output.
#' @param betaParams Numeric vector of length 2 for the Beta distribution
#'   parameters of alternative nodes (see \code{\link{draw_node_p_value}}).
#'
#' @details
#' \enumerate{
#'   \item The root node (level 0) is tested unconditionally:
#'         if \code{alt[1] = FALSE}, its p-value is \eqn{\mathrm{Unif}(0,1)};
#'         if \code{alt[1] = TRUE}, it is a scaled Beta(\code{betaParams}).
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
#' altVec <- assign_alt(
#'   k = 2, l = 2, t = 0.5,
#'   levelOffsets = info$levelOffsets, levelSizes = info$levelSizes
#' )
#' simulate_single_run(
#'   k = 2, l = 2, alpha = 0.05, altVec,
#'   levelOffsets = info$levelOffsets, levelSizes = info$levelSizes, betaParams = c(0.1, 1)
#' )
#'
#' @export
simulate_single_run <- function(k, l, alpha, alt, levelOffsets, levelSizes, betaParams) {
  nTot <- sum(levelSizes)
  pvals <- rep(NA_real_, nTot)
  tested <- rep(FALSE, nTot)

  falseReject <- FALSE

  # Root
  rootIdx <- 1
  pvals[rootIdx] <- draw_node_p_value(alt[rootIdx], 0, betaParams)
  tested[rootIdx] <- TRUE
  if (!alt[rootIdx] && pvals[rootIdx] <= alpha) {
    falseReject <- TRUE
  }

  # Traverse level by level
  for (d in seq_len(l)) {
    start_d <- levelOffsets[d + 1] + 1
    end_d <- levelOffsets[d + 1] + levelSizes[d + 1]
    nodeRange <- seq(start_d, end_d)

    # Parents that are tested and p <= alpha
    testedParents <- nodeRange[tested[nodeRange] & (pvals[nodeRange] <= alpha)]
    if (!length(testedParents)) next

    for (parentIdx in testedParents) {
      # local index j for the parent at level d
      j <- parentIdx - levelOffsets[d + 1]

      # identify children
      childInds <- children_indices(d, j, k, levelOffsets)
      if (!length(childInds)) next

      # draw child p-values
      parent_p <- pvals[parentIdx]
      childPvals <- numeric(k)
      for (iC in seq_len(k)) {
        cIdx <- childInds[iC]
        childPvals[iC] <- draw_node_p_value(alt[cIdx], parent_p, betaParams)
      }
      pvals[childInds] <- childPvals
      tested[childInds] <- TRUE

      # check false rejections among children
      nullChildren <- (!alt[childInds])
      if (any(nullChildren & (childPvals <= alpha))) {
        falseReject <- TRUE
      }

      # local Simes test
      simesP <- local_simes(childPvals)
      if (simesP > alpha) {
        # block deeper testing from these children
        pvals[childInds] <- alpha + 1e-8
      }
    }
  }

  falseReject
}

#' @title Monte Carlo Simulation of a Hierarchical Simes-Gated Testing Procedure
#'
#' @description Repeats \code{nSim} independent realizations of the following steps:
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
#' @param nSim Number of simulation replicates.
#' @param k Integer; branching factor.
#' @param l Integer; depth of the tree.
#' @param t Numeric in \eqn{[0,1]}; probability that a leaf is non-null.
#' @param alpha Numeric in \eqn{(0,1)}; significance level threshold.
#' @param betaParams Numeric vector of length 2, \code{c(a,b)}, specifying Beta(a,b)
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
#'   nSim = 2000,
#'   k = 2,
#'   l = 3,
#'   t = 0.5,
#'   alpha = 0.05,
#'   betaParams = c(0.1, 1)
#' )
#' fwer_est
#' }
#'
#' @export
simulate_hier_simes_local_modular <- function(nSim = 10000,
                                              k = 2,
                                              l = 3,
                                              t = 0.5,
                                              alpha = 0.05,
                                              betaParams = c(0.1, 1)) {
  info <- get_level_info(k, l)
  off <- info$levelOffsets
  sizes <- info$levelSizes

  falseCount <- 0

  for (rep in seq_len(nSim)) {
    # 1) alt assignment
    altVec <- assign_alt(k, l, t, off, sizes)

    # 2) single run
    hadFalseRej <- simulate_single_run(k, l, alpha, altVec, off, sizes, betaParams)

    if (hadFalseRej) {
      falseCount <- falseCount + 1
    }
  }

  # fraction that had false rejection => FWER estimate
  fwer <- falseCount / nSim
  fwer
}
