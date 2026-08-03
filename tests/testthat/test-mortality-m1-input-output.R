# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test update M1 from data", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS") # ?BS2017SS for more information on the data
  data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower


  ################################################
  # Estimation
  ################################################
  # - Single species
  # M1_at_age at the input parameters equals the M1_base input, so build-only
  # (estimateMode = 3) is sufficient -- no optimization is exercised.
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Build only
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(verbose = 0))


  # - Multi species
  # -- Intialized from single-species parameter values
  ms_run <- fit_mod(data_list = BS2017MS,
                    inits = ss_run$estimated_params, # Initial parameters from single species
                    file = NULL, # Don't save
                    estimateMode = 3, # Build only
                    niter = 1, # 1 iterations around population and predation dynamics
                    random_rec = FALSE, # No random recruitment
                    M1Fun = build_M1(updateM1 = TRUE),
                    msmMode = 1, # MSVPA based
                    suitMode = 0, # empirical suitability
                    fit_control = fit_control(verbose = 0))

  testthat::expect_equal(as.numeric(ms_run$quantities$M1_at_age[1,1,1:12,1]), as.numeric(ms_run$data_list$M1_base[1, 3:(12+2)]))
})
