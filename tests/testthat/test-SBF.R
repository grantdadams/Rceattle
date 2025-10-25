test_that("test SB40 under mean recruitment", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$proj_F_prop <- rep(1, 7)

  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = NULL, # Initial parameters from ss_run
                                    estimateMode = 0, # Run projection only
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    Ftarget = 0.4, # F40%
                                                    Flimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    msmMode = 0, # Single species mode
                                    phase = TRUE,
                                    verbose = 1)

  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$SBF/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001)
})


test_that("test SPR0", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify recruitment
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Set params
  inits <- build_params(GOA2018SS)
  inits$rec_pars[,1] <- R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 2,
                              verbose = 1)

  # Calculate SPR
  M <- 0.2
  wt <- 1
  sex_ratio <- 0.5
  pmature <- 1

  natage = rep(1, 10)
  for(age in 2:9){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[10] <-  natage[9] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  ssb_at_age <- natage * wt * pmature * sex_ratio
  sum(ssb_at_age)

  # Check SPRO
  testthat::expect_equal(as.numeric(ss_run$quantities$SPR0),  rep(sum(ssb_at_age), 3), tolerance = 0.0001)
})
