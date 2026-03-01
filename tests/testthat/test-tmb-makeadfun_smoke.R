testthat::skip_on_cran()
testthat::test_that("MakeADFun smoke test: obj$fn and obj$gr are finite", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  # source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)

  # Ensure compiled TMB is available
  compile_tmb_if_needed()

  # Run model in debug mode that only constructs the TMB object
  res <- Rceattle::fit_mod(data_list = dat,
                           inits = NULL,
                           estimateMode = 3,
                           phase = FALSE,
                           verbose = 0)

  testthat::expect_true(is.list(res))
  testthat::expect_true(!is.null(res$obj))

  # Evaluate objective and gradient
  fval <- res$obj$fn()
  gval <- res$obj$gr()

  testthat::expect_true(is.finite(fval))
  testthat::expect_true(all(is.finite(gval)))
})
