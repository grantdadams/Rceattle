# Integration test (fits a CEATTLE model): skipped on CRAN.
#
# sim_mod() returns the model's EXPECTED data (simulate = FALSE) or SAMPLED data
# (simulate = TRUE). The composition sampling branch (rmultinom / Dirichlet-
# multinomial draws, keyed by the string Comp_distribution) was previously
# untested -- it is where a comp/CAAL code-vs-string mismatch survived (PR 0.5).
#
# Scope: this covers the comp_data Multinomial path (the primary previously-
# untested branch). The CAAL sampling branch and the DirichletMultinomial branch
# are not yet exercised here -- a worthwhile follow-up.
testthat::skip_on_cran()

testthat::test_that("sim_mod returns expected comp proportions (simulate = FALSE)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
  fit <- suppressMessages(fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  testthat::expect_gt(nrow(fit$data_list$comp_data), 0)

  exp_dat <- Rceattle::sim_mod(fit, simulate = FALSE)
  comp_cols <- grep("^Comp_", names(exp_dat$comp_data))
  got <- unname(as.matrix(exp_dat$comp_data[, comp_cols]))
  # Expected data are exactly the model's predicted composition (comp_hat).
  testthat::expect_equal(got, unname(fit$quantities$comp_hat), tolerance = 1e-10)
})

testthat::test_that("sim_mod multinomial draws sum to sample size (simulate = TRUE)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
  fit <- suppressMessages(fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))

  set.seed(42)
  sim_dat <- Rceattle::sim_mod(fit, simulate = TRUE)
  comp_cols <- grep("^Comp_", names(sim_dat$comp_data))
  counts <- as.matrix(sim_dat$comp_data[, comp_cols])

  # rmultinom draws: non-negative integers summing to the row's Sample_size.
  testthat::expect_true(all(counts >= 0, na.rm = TRUE))
  testthat::expect_equal(counts, round(counts))
  testthat::expect_equal(unname(rowSums(counts, na.rm = TRUE)),
                         unname(sim_dat$comp_data$Sample_size),
                         tolerance = 1e-8)
  # A draw differs from the deterministic expectation (it actually sampled).
  exp_dat <- Rceattle::sim_mod(fit, simulate = FALSE)
  testthat::expect_false(isTRUE(all.equal(
    counts, unname(as.matrix(exp_dat$comp_data[, comp_cols])))))
})
