library(testthat)

source("../../R/tree_test_sims.R") # Source the file with the code

test_that("get_level_info() returns correct sizes and offsets", {
  k <- 2
  l <- 3
  info <- get_level_info(k, l)
  expect_equal(info$level_sizes, c(1, 2, 4, 8)) # 2^(0:3)
  expect_equal(info$level_offsets, c(0, 1, 3, 7))
  expect_equal(info$n_tot, 15)
  ## We know that this is the formula for the total nodes in a k-ary tree
  total_nodes <- (k^(l+1) - 1) / (k - 1)
  expect_equal(info$n_tot,total_nodes)
})

test_that("children_indices() works for k=2, l=3", {
  # levels: d=0 => 1 node, d=1 => 2 nodes, d=2 => 4 leaves
  info <- get_level_info(k = 2, l = 3)
  # The root is level 0 => global idx=1 => local j=1 at d=0
  # children of root => level 1 => 2 nodes => global indices 2,3
  cidx <- children_indices(d = 0, j = 1, k = 2, level_offsets = info$level_offsets)
  expect_equal(cidx, c(2, 3))

  # now for level=1, local j=1 => global idx=2 => children => global indices 4,5
  cidx <- children_indices(d = 1, j = 1, k = 2, level_offsets = info$level_offsets)
  expect_equal(cidx, c(4, 5))
  cidx <- children_indices(d = 1, j = 2, k = 2, level_offsets = info$level_offsets)
  expect_equal(cidx, c(6, 7))
  # for level 2 (d=2) and the first node in that level (j=1), which is really global idx 4
  cidx <- children_indices(d = 2, j = 1, k = 2, level_offsets = info$level_offsets)
  expect_equal(cidx, c(8, 9))
  cidx <- children_indices(d = 2, j = 2, k = 2, level_offsets = info$level_offsets)
  expect_equal(cidx, c(10, 11))
  ## Look at the leaves --- should have no children
  cidx <- children_indices(d = 3, j = 1, k = 2, level_offsets = info$level_offsets)
  expect_equal(cidx, integer(0))
})

test_that("assign_alt() sets leaf alt with prob t and propagates up", {
  set.seed(123)
  info <- get_level_info(k = 2, l = 2)
  alt_vec <- assign_alt(
    k = 2, l = 2, t = 0.5,
    level_offsets = info$level_offsets,
    level_sizes = info$level_sizes
  )
  # _this is random, but we can check that if a leaf is alt => its ancestors are alt
  # _we'll do multiple replicates in a loop to ensure no errors
  for (i in 1:10) {
    alt_vec <- assign_alt(
      k = 2, l = 2, t = 0.5,
      level_offsets = info$level_offsets,
      level_sizes = info$level_sizes
    )
    # _check that each leaf's alt => parent's alt => root alt
    # level_2: global indices = 4,5,6,7
    for (leaf_idx in 4:7) {
      if (alt_vec[leaf_idx]) {
        # parent of leaf_idx is either 2 or 3
        parent_idx <- if (leaf_idx %in% c(4, 5)) 2 else 3
        expect_true(alt_vec[parent_idx])
        # root is 1
        expect_true(alt_vec[1])
      }
    }
  }
})

test_that("draw_node_p_value() yields monotonic p-values in [parent_p,1]", {
  parent_p <- 0.3
  for (i in 1:100) {
    p_null <- draw_node_p_value(FALSE, parent_p, c(0.1, 1))
    p_alt <- draw_node_p_value(TRUE, parent_p, c(0.1, 1))
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
  info <- get_level_info(k = 2, l = 2)
  alt_vec <- assign_alt(k = 2, l = 2, t = 0.5, info$level_offsets, info$level_sizes)
  out <- simulate_single_run(
    k = 2, l = 2, alpha = 0.05,
    alt = alt_vec,
    level_offsets = info$level_offsets,
    level_sizes = info$level_sizes,
    beta_params = c(0.1, 1)
  )
  expect_true(is.logical(out))
  expect_length(out, 1)
})

test_that("simulate_hier_simes_local_modular() returns plausible FWER", {
  skip_on_cran()
  set.seed(999)
  est <- simulate_hier_simes_local_modular(n_sim = 1000, k = 2, l = 2, t = 0.5, alpha = 0.05)
  # _just check that it's between 0 and 1
  expect_true(est >= 0 && est <= 1)
})
