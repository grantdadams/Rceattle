test_that("M1 input and output the same", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$proj_F_prop <- rep(1, 7)

  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = NULL, # Initial parameters from ss_run
                                    estimateMode = 0, # Run projection only
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    FsprTarget = 0.4, # F40%
                                                    FsprLimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    msmMode = 0, # Single species mode
                                    phase = "default",
                                    verbose = 1)

  expect_equal((ss_run_Tier3$quantities$SBF/ss_run_Tier3$quantities$SB0)[,82], rep(0.4, 3), tolerance = 0.0001)
})
