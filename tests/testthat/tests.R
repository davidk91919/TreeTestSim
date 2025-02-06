library(testthat)

devtools::load_all()
# source("R/tree_test_sims.R") # Source the file with the code

info_k3_l2 <- get_level_info(k = 3, l = 2)
## k=3, l=2
## 1-based
## 1
##   2       3       4
## 5 6 7  8 9 10 11 12 13

## 0-based
## 0
##   1       2       3
## 4 5 6  7 8 9  10 11 12

test_that("get_level_info() returns correct sizes and offsets", {
  k <- 2
  l <- 3
  info <- get_level_info(k, l)
  expect_equal(info$k, k)
  expect_equal(info$l, l)
  expect_equal(info$level_sizes, c(1, 2, 4, 8)) # 2^(0:3)
  expect_equal(info$level_offsets, c(0, 1, 3, 7))
  expect_equal(info$n_tot, 15)
  ## We know that this is the formula for the total nodes in a k-ary tree
  total_nodes <- (k^(l + 1) - 1) / (k - 1)
  expect_equal(info$n_tot, total_nodes)
  ## k > 0
  expect_error(get_level_info(k = 0, l = 2))
  ## When l=0 we get just level_sizes=1 and level_offsets=0. ignores k since we have only a root node
  info_l0 <- get_level_info(k = 3, l = 0)
  expect_equal(info_l0$level_sizes, 1)
  expect_equal(info_l0$level_offsets, 0)
})

test_that("children_indices() works", {
  # levels: d=0 => 1 node, d=1 => 3 nodes, d=2 => 9 nodes, d=3 (leaves) => 3^3 = 27 leaves
  # total nodes = 1+3+9+27 = 40 = (3^(3+1) - 1)/(3-1)
  info <- get_level_info(k = 3, l = 3)
  # The root is level 0 => global idx=1 => local j=1 at d=0
  # children of root => level 1 => 3 nodes => global indices 2,3,4
  cidx_l1 <- children_indices(d = 0, j = 1, tree_info = info)
  expect_equal(cidx_l1, c(2, 3, 4))

  expect_error(children_indices(d = 0, j = 2, tree_info = info))
  expect_error(children_indices(d = 0, j = 0, tree_info = info))
  expect_error(children_indices(d = 1, j = 0, tree_info = info))

  # now for level=1, local j=1 => global idx=2 => children => global indices 5,6,7
  cidx_l2_node1 <- children_indices(d = 1, j = 1, tree_info = info)
  expect_equal(cidx_l2_node1, c(5, 6, 7))
  # now for level=1, local j=2 => global idx=3 => children => global indices 8,9,10
  cidx_l2_node2 <- children_indices(d = 1, j = 2, tree_info = info)
  expect_equal(cidx_l2_node2, c(8, 9, 10))
  # now for level=1, local j=3 => global idx=4 => children => global indices 11,12,13
  cidx_l2_node3 <- children_indices(d = 1, j = 3, tree_info = info)
  expect_equal(cidx_l2_node3, c(11, 12, 13))

  cidx_l2s <- lapply(cidx_l1, function(global_idx) {
    local_idx <- global_idx - 1
    children_indices(d = 1, j = local_idx, tree_info = info)
  })
  expect_equal(cidx_l2s[[1]], c(5, 6, 7))
  expect_equal(cidx_l2s[[2]], c(8, 9, 10))
  expect_equal(cidx_l2s[[3]], c(11, 12, 13))

  ## Now the next level in the tree. Should have 27 nodes in total. Three per each of the 9 nodes at level 1
  cidx_leaves <- lapply(unlist(cidx_l2s), function(global_idx) {
    local_idx <- global_idx - info$level_offsets[3]
    children_indices(d = 2, j = local_idx, tree_info = info)
  })
  expect_equal(length(unlist(cidx_leaves)), 27)
  expect_equal(unlist(cidx_leaves), seq(14, 40))

  ## No leaves should have children
  cidx_l3s <- sapply(seq(1, info$level_sizes[4]), function(local_indx) {
    tmp <- children_indices(d = 3, j = local_idx, tree_info = info)
    all.equal(tmp, integer(0))
  })
  expect_true(all(cidx_l3s))
})

test_that("The get_ancestors function works as it should", {
  info <- get_level_info(k = 2, l = 3)
  info_k3_l3 <- get_level_info(k = 3, l = 3)
  ## 1-based
  ##            1
  ##            2                      3       4
  ## 5         6             7     8 9 10     11 12 13
  ## 14 15 16  17 18 19 ...
  expect_equal(get_ancestors(node_id = 8, k = 3), c(3, 1))
  expect_equal(get_ancestors(node_id = 17, k = 3), c(6, 2, 1))
  expect_equal(get_ancestors(node_id = 14, k = 3), c(5, 2, 1))
  expect_equal(get_ancestors(node_id = 8, k = 2), c(4, 2, 1))
  expect_equal(get_ancestors(node_id = 15, k = 2), c(7, 4, 2, 1))
})

test_that("assign_alt() sets leaf alt with prob t and propagates up", {
  set.seed(123)
  info <- get_level_info(k = 2, l = 3)
  #  1
  # 2 3
  #  4    5    6      7
  # 8 9 10 11 12 13 14 15

  ## In a 3 level tree with k=2, we know that we have k^3=2^3=8 leaves (level 3 are levels, k=2 we ave k^l nodes at each level)
  ## And in our vector index format, the last 8 entries are the leavesa
  ## The general formaul is level_offset[l+1]+1 ... level_offsets[l+1]+k^d, where d is the depth
  global_leaf_indices <- seq.int(info$level_offsets[3 + 1] + 1, info$level_offsets[3 + 1] + 2^3)
  expect_equal(global_leaf_indices, seq.int(8, 15))

  # this is 4 out of 8
  alt_vec <- assign_alt(t = 0.5, tree_info = info)
  ## make sure that the correct number is assigned
  expect_number_assigned <- floor(.5 * length(global_leaf_indices))
  expect_equal(sum(alt_vec[global_leaf_indices]), expect_number_assigned)

  ## This rounds down (this is fine. We want to round down for the purposes of the simulation)
  alt_vec <- assign_alt(t = 0.3, tree_info = info)
  ## make sure that the correct number is assigned
  expect_number_assigned <- floor(.3 * length(global_leaf_indices))
  expect_equal(sum(alt_vec[global_leaf_indices]), expect_number_assigned)

  # Make sure that the vectors are the same length
  expect_equal(length(alt_vec), info$n_tot)

  # this is random, but we can check that if a leaf is alt => its ancestors are alt
  # we'll do multiple replicates in a loop to ensure no errors
  for (i in 1:10) {
    alt_vec <- assign_alt(t = 0.5, tree_info = info)
    # _check that each leaf's alt => parent's alt => root alt
    for (leaf_idx in global_leaf_indices) {
      if (alt_vec[leaf_idx]) { ## if a leaf is TRUE
        leaf_ancestors_idxs <- get_ancestors(node_id = leaf_idx, k = info$k)
        expect_true(all(alt_vec[leaf_ancestors_idxs]))
        # root is 1
        expect_true(alt_vec[1])
      }
    }
  }
  ##
})

test_that("draw_node_p_value() yields monotonic p-values in [parent_p,1]", {
  ## recall that when is_alt=TRUE we use the beta distribution but when is_alt=FALSE we use the uniform distribution
  parent_p <- 0.01
  for (i in 1:100) {
    p_null <- draw_node_p_value(is_alt = FALSE, parent_p, beta_params = c(0.1, 1))
    p_alt <- draw_node_p_value(is_alt = TRUE, parent_p, beta_params = c(0.1, 1))
    expect_true(p_null >= parent_p && p_null <= 1)
    expect_true(p_alt >= parent_p && p_alt <= 1)
  }
})

test_that("local_simes() example check", {
  x <- c(0.01, 0.05, 0.10)
  # sorted => 0.01,0.05,0.10; k=3
  # simes => min( (3/1)*0.01, (3/2)*0.05, (3/3)*0.10 ) = min(0.03, 0.075, 0.10 )=0.03
  psimes <- local_simes(x)
  expect_equal(psimes, 0.03, tolerance = 1e-8)
})

test_that("simulate_single_run() runs without error, returns T/F", {
  info <- get_level_info(k = 2, l = 3)
  alt_vec <- assign_alt(t = 0.3, tree_info = info)
  out <- simulate_single_run(
    tree_info = info,
    alpha = 0.05,
    alt = alt_vec,
    beta_params = c(0.1, 1)
  )
  expect_true(is.logical(out))
  expect_length(out, 1)
})

test_that("simulate_hier() returns plausible FWER", {
  skip_on_cran()
  set.seed(999)
  est <- simulate_hier_simes_local_modular(n_sim = 1000, k = 2, l = 2, t = 0.5, alpha = 0.05)
  # _just check that it's between 0 and 1
  expect_true(est >= 0 && est <= 1)
})
