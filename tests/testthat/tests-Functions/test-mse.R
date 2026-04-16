testthat::test_that("Test MSE - Tier 3 w no uncertainty", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip()

  library(Rceattle)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Data ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$projyr <- 2040
  BS2017SS$fleet_control$proj_F_prop <-rep(1,7)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Operating Model ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Single-species with fixed M
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # NPFMC Tier 3 ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  ss_run_Tier3 <- suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS,
                      inits = NULL, # Initial parameters = 0
                      estimateMode = 0, # Run projection only
                      HCR = build_hcr(HCR = 5, # Tier3 HCR
                                      Ftarget = 0.4, # F40%
                                      Flimit = 0.35, # F35%
                                      Plimit = 0.2, # No fishing when SB<SB20
                                      Alpha = 0.05),
                      msmMode = 0, # Single species mode
                      phase = TRUE,
                      verbose = 1)
  )


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # MSE ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  mse <- run_mse(om = ss_run, em = ss_run_Tier3, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = FALSE, sample_rec = FALSE)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Tests ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  nyrs <- length(1979:2040)
  testthat::expect_equal(24, length(mse$Sim_1$EM))
  testthat::expect_equal(rep(0.4, 3), as.numeric(mse$Sim_1$OM$quantities$ssb_depletion[,nyrs]), tolerance = 0.005)
})
