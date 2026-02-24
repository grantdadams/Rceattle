testthat::skip_on_cran()
testthat::test_that("Fit sanity regression: key quantities match baseline", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 1)

  # Ensure compiled TMB exists
  compile_tmb_if_needed()

  # Run a short/safe fit
  fit <- Rceattle::fit_mod(data_list = dat,
                          inits = NULL,
                          estimateMode = 0,
                          phase = FALSE,
                          loopnum = 1,
                          getsd = FALSE,
                          use_gradient = FALSE,
                          control = list(eval.max = 200, iter.max = 200),
                          verbose = 0)

  testthat::expect_s3_class(fit, "list")
  testthat::expect_true(!is.null(fit$quantities))

  # Collect a small set of key quantities to compare
  got <- list(
    R0 = as.numeric(fit$quantities$R0),
    ssb1 = as.numeric(fit$quantities$ssb[1, 1:min(ncol(fit$quantities$ssb), dat$endyr)]),
    biomass1 = as.numeric(fit$quantities$biomass[1, 1:min(ncol(fit$quantities$biomass), dat$endyr)])
  )

  baseline_path <- file.path("tests", "testthat", "fixtures", "fit_baseline.rds")

  # Regenerate baseline explicitly when requested by env var
  if (identical(tolower(Sys.getenv("REGEN_TEST_BASELINES")), "true")) {
    dir.create(dirname(baseline_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(got, baseline_path)
    testthat::skip("Baseline regenerated; rerun test without REGEN_TEST_BASELINES to validate.")
  }

  testthat::expect_true(file.exists(baseline_path), info = paste0("Baseline file missing: ", baseline_path, ". To create it run with REGEN_TEST_BASELINES=true."))
  expected <- readRDS(baseline_path)

  # Compare shapes
  testthat::expect_equal(length(got$R0), length(expected$R0))
  testthat::expect_equal(length(got$ssb1), length(expected$ssb1))
  testthat::expect_equal(length(got$biomass1), length(expected$biomass1))

  # Compare numeric values with tolerant thresholds
  testthat::expect_equal(got$R0, expected$R0, tolerance = 1e-6)
  testthat::expect_equal(got$ssb1, expected$ssb1, tolerance = 1e-6)
  testthat::expect_equal(got$biomass1, expected$biomass1, tolerance = 1e-6)
})
