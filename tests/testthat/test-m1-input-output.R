test_that("M1 input and output the same", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower


  ################################################
  # Estimation
  ################################################
  # - Single species
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)


  # - Multi species
  # -- Intialized from single-species MLE parameter values
  ms_run <- fit_mod(data_list = BS2017MS,
                    inits = ss_run$estimated_params, # Initial parameters from single species ests
                    file = NULL, # Don't save
                    estimateMode = 0, # Estimate
                    niter = 3, # 3 iterations around population and predation dynamics
                    random_rec = FALSE, # No random recruitment
                    M1Fun = build_M1(updateM1 = TRUE),
                    msmMode = 1, # MSVPA based
                    suitMode = 0, # empirical suitability
                    verbose = 1)

  testthat::expect_equal(as.numeric(ms_run$quantities$M1[1,1,1:12]), as.numeric(ms_run$data_list$M1_base[1, 3:(12+2)]))
})
