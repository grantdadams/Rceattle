testthat::test_that("Test SB40 under mean recruitment", {
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
  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$ssb/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001) # Default to mean R
})


testthat::test_that("Test SPR0 calculation", {
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
  # inits <- build_params(GOA2018SS)
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 2,
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  inits$rec_pars[,1] <- R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters at input
                              file = NULL, # Don't save
                              map = map,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 2,
                              verbose = 0)

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


testthat::test_that("Test hindcast the same across different HCRs/BRPs", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$proj_F_prop <- rep(1, 7)

  # -- NPFMC Tier 3
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                                    estimateMode = 0,
                                    msmMode = 0, # Single species mode
                                    phase = TRUE,
                                    verbose = 1)

  # -- Constant F as a percentage of SB0
  ss_run_fb0 <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 2, # Run projection only
                                  HCR = build_hcr(HCR = 3, # Constant F HCR
                                                  DynamicHCR = FALSE, # Use dynamic reference points
                                                  Ftarget = 0.4), # F that achieves 40% SB0
                                  msmMode = 0, # Single species mode
                                  verbose = 1)


  ss_run_dynamicfb0 <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 2, # Run projection only
                                         HCR = build_hcr(HCR = 3, # Constant F HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         Ftarget = 0.4), # F that achieves 40% SB0
                                         msmMode = 0, # Single species mode
                                         verbose = 1)


  # -- Constant Fspr
  ss_run_Fspr <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = ss_run$estimated_params, # Initial parameters from ss_run
                                   estimateMode = 2, # Run projection only
                                   HCR = build_hcr(HCR = 4, # Tier3 HCR
                                                   Ftarget = 0.4 # F40%
                                   ),
                                   msmMode = 0, # Single species mode
                                   verbose = 1)


  ss_run_dynamicFspr <- Rceattle::fit_mod(data_list = BS2017SS,
                                          inits = ss_run$estimated_params, # Initial parameters from ss_run
                                          estimateMode = 2, # Run projection only
                                          HCR = build_hcr(HCR = 4, # Tier3 HCR
                                                          DynamicHCR = TRUE, # Use dynamic reference points
                                                          Ftarget = 0.4 # F40%
                                          ),
                                          msmMode = 0, # Single species mode
                                          verbose = 1)


  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = ss_run$estimated_params, # Initial parameters from ss_run
                                    estimateMode = 2, # Run projection only
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    Ftarget = 0.4, # F40%
                                                    Flimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    msmMode = 0, # Single species mode
                                    verbose = 1)


  ss_run_dynamicTier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                           inits = ss_run$estimated_params, # Initial parameters from ss_run
                                           estimateMode = 2, # Run projection only
                                           HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                           DynamicHCR = TRUE, # Use dynamic reference points
                                                           Ftarget = 0.4, # F40%
                                                           Flimit = 0.35, # F35%
                                                           Plimit = 0.2, # No fishing when SB<SB20
                                                           Alpha = 0.05),
                                           msmMode = 0, # Single species mode
                                           verbose = 1)

  # -- PFMC Category 1
  ss_run_Cat1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = ss_run$estimated_params, # Initial parameters from ss_run
                                   estimateMode = 2, # Run projection only
                                   HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                   Flimit = 0.45, # F45%
                                                   Ptarget = 0.4, # Target is 40% B0
                                                   Plimit = 0.1, # No fishing when SB<SB10
                                                   Pstar = 0.45,
                                                   Sigma = 0.5),
                                   msmMode = 0, # Single species mode
                                   verbose = 1)

  ss_run_dynamicCat1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                          inits = ss_run$estimated_params, # Initial parameters from ss_run
                                          estimateMode = 2, # Run projection only
                                          HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                          DynamicHCR = TRUE, # Use dynamic reference points
                                                          Flimit = 0.45, # F45%
                                                          Ptarget = 0.4, # Target is 40% SB0
                                                          Plimit = 0.1, # No fishing when SB<SB10
                                                          Pstar = 0.45,
                                                          Sigma = 0.5),
                                          msmMode = 0, # Single species mode
                                          verbose = 1)

  # -- SESSF Tier 1
  ss_run_Tier1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = ss_run$estimated_params, # Initial parameters from ss_run
                                    estimateMode = 2, # Run projection only
                                    HCR = build_hcr(HCR = 7, # Tier 1 HCR
                                                    Ftarget = 0.48, # F40%
                                                    Flimit = 0.20, # F20%
                                                    Ptarget = 0.35, # Target is 35% SSB0
                                                    Plimit = 0.20, # No fishing when B<B20
                                    ),
                                    msmMode = 0, # Single species mode
                                    verbose = 1)


  ss_run_dynamicTier1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                           inits = ss_run$estimated_params, # Initial parameters from ss_run
                                           estimateMode = 2, # Run projection only
                                           HCR = build_hcr(HCR = 7, # Tier 1 HCR
                                                           DynamicHCR = TRUE,
                                                           Ftarget = 0.48, # F40%
                                                           Flimit = 0.20, # F20%
                                                           Ptarget = 0.35, # Target is 35% SSB0
                                                           Plimit = 0.20, # No fishing when B<B20
                                           ),
                                           msmMode = 0, # Single species mode
                                           verbose = 1)

  # -- Plot
  mod_list <- list(ss_run, ss_run_fb0, ss_run_Fspr, ss_run_Tier3, ss_run_Cat1, ss_run_Tier1, ss_run_dynamicfb0, ss_run_dynamicFspr, ss_run_dynamicTier3, ss_run_dynamicCat1, ss_run_dynamicTier1)
  nyrs <- length(BS2017SS$styr:BS2017SS$endyr)
  nyrs_proj <- length(BS2017SS$styr:BS2017SS$projyr)
  proj_yrs <- (nyrs+2):nyrs_proj

  # Test hindcast is the same
  for(i in 2:length(mod_list)){
    testthat::expect_equal(ss_run$quantities$ssb[,1:nyrs], mod_list[[i]]$quantities$ssb[,1:nyrs])
  }

  # Test forecast is different
  for(i in 2:length(mod_list)){
    testthat::expect_all_true(c(ss_run$quantities$ssb[,proj_yrs] != mod_list[[i]]$quantities$ssb[,proj_yrs]))
  }
})
