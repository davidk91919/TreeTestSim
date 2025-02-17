context("Data.table-based Tree Simulation Functions")

devtools::load_all()

test_that("generate_tree_DT creates a complete tree with correct columns and levels", {
  max_level <- 2
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)

  # Expected number of nodes is sum(k^l) for l=0:max_level.
  expected_rows <- sum(k^(0:max_level))
  expect_equal(nrow(tree_dt), expected_rows)

  # The data.table should have these columns.
  expect_true(all(c("node", "level", "parent", "nonnull") %in% names(tree_dt)))

  # Check that the level assignments are correct:
  expect_equal(tree_dt[level == 0, .N], 1)
  expect_equal(tree_dt[level == 1, .N], k)
  expect_equal(tree_dt[level == 2, .N], k^2)

  # Check that the parent indices are as expected:
  # For node i > 1, parent = floor((i-2)/k) + 1.
  tree_dt[, computed_parent := ifelse(node == 1, NA_integer_, floor((node - 2) / k) + 1)]
  expect_equal(tree_dt$parent, tree_dt$computed_parent)
})

test_that("generate_tree_DT propagates nonnull flags correctly for extreme t", {
  max_level <- 3
  k <- 3

  # t = 0: all leaves false => all nodes should be false.
  dt0 <- generate_tree_DT(max_level, k, t = 0)
  expect_true(all(dt0[level == max_level, nonnull] == FALSE))
  expect_true(all(dt0$nonnull == FALSE))

  # t = 1: all leaves true => all nodes should be true.
  dt1 <- generate_tree_DT(max_level, k, t = 1)
  expect_true(all(dt1[level == max_level, nonnull] == TRUE))
  expect_true(all(dt1$nonnull == TRUE))

  # t=.5, about .5 of the leaves should be TRUE
  dt_half <- generate_tree_DT(max_level, k, t = .5)
  expect_equal(dt_half[level == max_level, sum(nonnull)], floor(.5 * nrow(dt_half[level == max_level])))
})

test_that("simulate_test_DT produces monotonic p-values", {
  # We simulate a small tree and verify that each child’s p-value is at least that of its parent.
  set.seed(1234)
  max_level <- 3
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)
  res <- simulate_test_DT(tree_dt, alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1)

  dt_sim <- res$treeDT

  # For each node with a parent, we require:
  #    child p_val >= parent's p_val.
  # Merge dt_sim with itself to obtain parent's p-values.
  parent_p_vals <- dt_sim[, .(node, parent_p = p_val)]
  children <- dt_sim[!is.na(parent)]
  children <- merge(children, parent_p_vals, by.x = "parent", by.y = "node", all.x = TRUE, sort = FALSE)

  # Test that for each child, p_val is at least as high as parent's p_val.
  expect_true(all(children$p_val >= children$parent_p, na.rm = TRUE))
})

test_that("simulate_test_DT gates branches when the Simes p-value exceeds alpha", {
  # Create a tree and force the children p-values to be high so that the Simes test fails.
  # One way to do this is to temporarily override runif() so that it returns 1.
  max_level <- 3
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)

  # Override runif locally.
  local_runif <- function(n, min, max) {
    rep(max, n)
  }
  old_runif <- stats::runif
  unlockBinding("runif", as.environment("package:stats"))
  assign("runif", local_runif, envir = as.environment("package:stats"))
  on.exit({
    assign("runif", old_runif, envir = as.environment("package:stats"))
    lockBinding("runif", as.environment("package:stats"))
  })

  ## recall that the mean of a beta distribution is a/(a+1)
  set.seed(12357)
  res <- simulate_test_DT(tree_dt, alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1)

  dt_sim <- res$treeDT

  # For any parent that is active, if the children branch was gated the
  # children p_val should be NA. We check that for each parent where Simes
  # test fails, none of its children have been assigned a p-value.

  gated_children <- dt_sim[!is.na(parent) & is.na(p_val)]
  expect_true(nrow(gated_children) > 0)

  ## check using the hommel adjustment too
  set.seed(12357)
  res <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, effN = 1000,
    N_total = 1000, beta_base = 0.1, local_adj_p_fn = local_hommel_all_ps
  )

  dt_sim <- res$treeDT

  # For any parent that is active, if the children branch was gated the
  # children p_val should be NA. We check that for each parent where the
  # hommel adjustment is > alpha, , none of its children have been assigned a
  # p-value.

  gated_children <- dt_sim[!is.na(parent) & is.na(p_val)]
  expect_true(nrow(gated_children) > 0)
})

test_that("simulating many p-values returns a value between 0 and 1", {
  set.seed(123456)

  ## We don't even need to adjust power as we "split" in order to control the
  ## FWER when t=0 or t=1 (which has no errors anyway and is only shown here
  ## for completeness)

  res1 <- simulate_many_runs_DT(
    n_sim = 1000, t = 0, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE
  )
  expect_true(is.numeric(res1))
  expect_true(res1 >= 0 && res1 <= 1)
  ## Should control res1 within simulation error
  expect_lt(res1, .05 + sqrt(.05 * (1 - .05) / 1000))

  res2 <- simulate_many_runs_DT(
    n_sim = 1000, t = 1, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE
  )
  ## No false positives possible here
  expect_lt(res2, .05 + sqrt(.05 * (1 - .05) / 1000))

  ## But we need to more closely simulate the power loss of splitting
  ## to control the FWER when we have a mix of nodes with effects

  res3 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = TRUE
  )
  res3
  expect_lt(res3, .05 + sqrt(.05 * (1 - .05) / 1000))

  ## Notice that this does not hold with adj_effN=FALSE. So, we need
  ## monotonicity and local gate but also need to reduce power when effects
  ## are mixed.
  res4 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE
  )
  res4
  expect_gt(res4, .05 + sqrt(.05 * (1 - .05) / 1000))

  ## Notice that the reduction in power is doing a lot of work here even
  ## without the local adjustment. For example, even with a large k and large l

  n_nodes <- sum(3^(0:10))

  ## the algorithm to reduce N is nonlinear -- starts by equal splitting by k
  ## but never goes to 0 (it goes no further than 1/100 of the original).

  res5 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 5, max_level = 5,
    alpha = 0.05, N_total = 10000000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_unadj_all_ps
  )
  res5
  expect_lt(res5, .05 + sqrt(.05 * (1 - .05) / 1000))

  ## Even local hommel doesn't control things, although it helps
  res6 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE,
    local_adj_p_fn = local_hommel_all_ps
  )
  res6
  expect_gt(res6, .05 + sqrt(.05 * (1 - .05) / 1000))

  res7 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_hommel_all_ps
  )
  res7
  expect_lt(res7, .05 + sqrt(.05 * (1 - .05) / 1000))
})

##
