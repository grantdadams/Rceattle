testthat::test_that("Test SB0 under mean recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nprojyrs = 20, nages = nages, seed = 42)

  # Set params
  inits <- suppressMessages(build_params(dat))
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- 0
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1

  # Set logistic params
  inits$ln_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females
  inits$rec_dev[] <- 0.02 # Adding .1 to R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters from starting
                              file = NULL, # Don't save
                              recFun = build_srr(
                                proj_mean_rec = TRUE, # Project using mean rec over hindcast
                              ),
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  # Calculate SB0
  M <- dat$M1_base$Age1
  wt <- dat$weight$Age1
  sex_ratio <- dat$sex_ratio$Age1
  pmature <- dat$maturity$Age1

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0 + 0.02) * sex_ratio

  # Check SB0
  testthat::expect_equal(as.numeric(ss_run$quantities$SB0[1,nyrs + 20]), SB0, tolerance = 0.0001)
})

testthat::test_that("Test SB0 under R0", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nprojyrs = 20, nages = nages, seed = 42)

  # Set params
  inits <- suppressMessages(build_params(dat))
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- 0
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1

  # Set logistic params
  inits$ln_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females
  inits$rec_dev[] <- 0.02 # Adding .1 to R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters from starting
                              file = NULL, # Don't save
                              recFun = build_srr(
                                proj_mean_rec = FALSE, # Project using R0 or SRR over hindcast
                              ),
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  # Calculate SB0
  M <- dat$M1_base$Age1
  wt <- dat$weight$Age1
  sex_ratio <- dat$sex_ratio$Age1
  pmature <- dat$maturity$Age1

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0) * sex_ratio

  # Check SB0
  testthat::expect_equal(as.numeric(ss_run$quantities$SB0[1,nyrs + 20]), SB0, tolerance = 0.0001)
})


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

testthat::test_that("Test SB40 under R0", {
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
                                    recFun = build_srr(proj_mean_rec = FALSE),
                                    msmMode = 0, # Single species mode
                                    phase = TRUE,
                                    verbose = 1)

  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$SBF/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001)
  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$ssb/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001)
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
                                  estimateMode = 0, # Run projection only
                                  HCR = build_hcr(HCR = 3, # Constant F HCR
                                                  DynamicHCR = FALSE, # Use dynamic reference points
                                                  Ftarget = 0.4), # F that achieves 40% SB0
                                  msmMode = 0, # Single species mode
                                  verbose = 1)


  ss_run_dynamicfb0 <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 0, # Run projection only
                                         HCR = build_hcr(HCR = 3, # Constant F HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         Ftarget = 0.4), # F that achieves 40% SB0
                                         msmMode = 0, # Single species mode
                                         verbose = 1)


  # -- Constant Fspr
  ss_run_Fspr <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = ss_run$estimated_params, # Initial parameters from ss_run
                                   estimateMode = 0, # Run projection only
                                   HCR = build_hcr(HCR = 4, # Tier3 HCR
                                                   Ftarget = 0.4 # F40%
                                   ),
                                   msmMode = 0, # Single species mode
                                   verbose = 1)


  ss_run_dynamicFspr <- Rceattle::fit_mod(data_list = BS2017SS,
                                          inits = ss_run$estimated_params, # Initial parameters from ss_run
                                          estimateMode = 0, # Run projection only
                                          HCR = build_hcr(HCR = 4, # Tier3 HCR
                                                          DynamicHCR = TRUE, # Use dynamic reference points
                                                          Ftarget = 0.4 # F40%
                                          ),
                                          msmMode = 0, # Single species mode
                                          verbose = 1)


  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = ss_run$estimated_params, # Initial parameters from ss_run
                                    estimateMode = 0, # Run projection only
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    Ftarget = 0.4, # F40%
                                                    Flimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    msmMode = 0, # Single species mode
                                    verbose = 1)


  ss_run_dynamicTier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                           inits = ss_run$estimated_params, # Initial parameters from ss_run
                                           estimateMode = 0, # Run projection only
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
                                   estimateMode = 0, # Run projection only
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
                                          estimateMode = 0, # Run projection only
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
                                    estimateMode = 0, # Run projection only
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
                                           estimateMode = 0, # Run projection only
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


testthat::test_that("Test mean recruitment calculation", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    log_phi = log_phi # set to -Inf so no predation
  )

  # sum(sim$model_quantities$B_eaten_as_prey)

  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  inits <- suppressMessages( build_params(simData) )
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = inits, # Initial parameters from inits
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              recFun = build_srr(
                                proj_mean_rec = TRUE, # Project using mean rec over hindcast
                              ),
                              random_rec = FALSE, # No random recruitment
                              phase = FALSE,
                              msmMode = 0,
                              suitMode = 0,
                              niter = 5,
                              initMode = 2,
                              verbose = 0)

  # Recruitment
  nyrs <- length(simData$styr:simData$endyr)
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ss_run$quantities$R[,1:nyrs]))

  # Mean R
  testthat::expect_equal(as.numeric(ss_run$quantities$avg_R), as.numeric(rowMeans(ss_run$quantities$R[,1:nyrs])), tolerance = 1e-6)
 })
