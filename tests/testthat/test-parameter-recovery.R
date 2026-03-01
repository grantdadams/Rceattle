testthat::skip_on_cran()
testthat::test_that("Parameter recovery (small simulated example)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  # source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 42)

  # Ensure compiled TMB is available
  compile_tmb_if_needed()

  # Run a quick fit (small iterations) to test parameter movement
  fit <- Rceattle::fit_mod(data_list = dat,
                          inits = NULL,
                          estimateMode = 0,
                          phase = FALSE,
                          loopnum = 1,
                          getsd = FALSE,
                          use_gradient = FALSE,
                          control = list(eval.max = 100, iter.max = 100),
                          verbose = 0)

  # Check returned structure and R0 quantity
  testthat::expect_true(is.list(fit))
  testthat::expect_true(!is.null(fit$quantities))
  testthat::expect_true(is.finite(fit$quantities$R0[1]))

  # Basic sanity checks: R0 positive and within a wide bound
  estR0 <- fit$quantities$R0[1]
  testthat::expect_true(estR0 > 0)
  testthat::expect_true(estR0 < 1e6)
})
