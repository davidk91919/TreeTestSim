context("Tree Simulation Functions")

## If interactive need to load the library
devtools::load_all()

test_that("generate_tree returns a valid tree structure", {
  tree <- generate_tree(level = 0, max_level = 2, k = 2, t = 0.5)
  expect_true(is.list(tree))
  expect_true("level" %in% names(tree))
  expect_true("nonnull" %in% names(tree))
  expect_true("children" %in% names(tree))
  if (!is.null(tree$children)) {
    expect_equal(length(tree$children), 2)
  }
})

test_that("simulate_test returns a list with p_val and false_error", {
  tree <- generate_tree(level = 0, max_level = 1, k = 2, t = 0.5)
  res <- simulate_test(tree, parent_p = 0, alpha = 0.05, k = 2, effN = 1000, N_total = 1000, beta_base = 0.1)
  expect_true(is.list(res))
  expect_true("p_val" %in% names(res))
  expect_true("false_error" %in% names(res))
  expect_true(is.numeric(res$p_val))
  expect_true(is.logical(res$false_error))
})

test_that("simulate_FWER returns a value between 0 and 1", {
  tree <- generate_tree(level = 0, max_level = 3, k = 3, t = 0.5)
  fwer <- simulate_FWER(n_sim = 100, t = 0.5, k = 3, max_level = 3, alpha = 0.05, N_total = 1000, beta_base = 0.1)
  expect_true(is.numeric(fwer))
  expect_true(fwer >= 0 && fwer <= 1)
  ## We expect weak control of the FWER (i.e. control when all nodes have no effects, when the null of no effects is true for all nodes)
  fwer_null_true <- simulate_FWER(n_sim = 10000, t = 0, k = 3, max_level = 3, alpha = 0.05, N_total = 1000, beta_base = 0.1)
  expect_lt(fwer_null_true, .06)
  ## When all nodes have an effect, then there should be no errors
  fwer_null_false <- simulate_FWER(n_sim = 10000, t = 1, k = 3, max_level = 3, alpha = 0.05, N_total = 1000, beta_base = 0.1)
  expect_equal(fwer_null_false, 0)

  fwer_mixed <- simulate_FWER(n_sim = 10000, t = 0.5, k = 3, max_level = 3, alpha = 0.05, N_total = 1000, beta_base = 0.1)
  print(fwer_mixed)
  fwer_mixed <- simulate_FWER(n_sim = 10000, t = 0.5, k = 10, max_level = 2, alpha = 0.05, N_total = 1000, beta_base = 0.1)
  print(fwer_mixed)
  fwer_mixed <- simulate_FWER(n_sim = 10000, t = 0.5, k = 3, max_level = 5, alpha = 0.05, N_total = 1000, beta_base = 0.1)
  print(fwer_mixed)
})

# TODO: check power, check and stop if effN becomes too low (like if effN ==1 stop)
## I think that hommel()$adjusted or p.adjust(ps,method="hommel") will be more conservative than max(simes_p,thep).
