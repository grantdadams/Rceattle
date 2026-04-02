testthat::skip_on_cran()
testthat::test_that("Fit sanity model: key quantities match baseline", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip()

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 1)

  # Ensure compiled TMB exists
  compile_tmb_if_needed()

  # Run a short/safe fit
  fit <- suppressMessages(
    Rceattle::fit_mod(data_list = dat,
                      inits = NULL,
                      estimateMode = 0,
                      phase = FALSE,
                      loopnum = 1,
                      getsd = FALSE,
                      use_gradient = FALSE,
                      control = list(eval.max = 200, iter.max = 200),
                      verbose = 0)
  )

  testthat::expect_s3_class(fit, "Rceattle")
  testthat::expect_true(!is.null(fit$quantities))

  # Collect a small set of key quantities to compare
  got <- list(
    R0 = as.numeric(fit$quantities$R0),
    ssb1 = as.numeric(fit$quantities$ssb[1,]),
    biomass1 = as.numeric(fit$quantities$biomass[1,])
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


testthat::test_that("Index, biomass, and catch = 0 match expected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  # Set params
  # inits <- suppressMessages(build_params(dat))
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- 0
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1


  # Specify logistic selectivity
  inf = 0; alpha = 0
  ages <- 1:5
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)


  # Set logistic params
  inits$ln_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              map = map,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  # Calculate SPR
  M <- dat$M1_base$Age1
  wt <- dat$weight$Age1
  sex_ratio <- dat$sex_ratio$Age1
  pmature <- dat$maturity$Age1

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0)

  # Check biomass
  testthat::expect_equal(as.numeric(ss_run$quantities$biomass[1,1:nyrs]), rep(SB0, nyrs), tolerance = 0.0001)

  # Check derived index
  testthat::expect_equal(ss_run$quantities$index_hat, 0.5 * SB0[ss_run$data_list$index_data$Species])

  # Check catch
  testthat::expect_equal(ss_run$quantities$catch_hat[1:nyrs], rep(0, nyrs))
})


testthat::test_that("Dynamics match CEATTLE single-species classic", {

  # Load old model
  load(system.file("extdata/CEATTLE_classic_ss.Rdata",package="Rceattle"))

  # Load data and set up inits
  data(BS2017SS) # ?BS2017SS for more information on the data

  BS2017SS$srr_prior_mean <- 9
  BS2017SS$initMode  <- 1
  BS2017SS$fleet_control$Sel_norm_bin1 <- -999 # Normalize by max
  # inits <- suppressMessages(build_params(BS2017SS))
  ss_run_old_params <- Rceattle::fit_mod(data_list = BS2017SS,
                                         estimateMode = 3,   # Dont' estimate
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0,        # Single species mode
                                         verbose = 0)
  inits <- ss_run_old_params$estimated_params
  map <- ss_run_old_params$map

  # - Update population dynamics from previous parameters
  inits$init_dev[,1:20] <- CEATTLE_classic_SS$estimated_params$init_dev
  inits$x_tj[1:39, 1:3] <- t(CEATTLE_classic_SS$estimated_params$rec_dev)
  inits$rec_pars[,1] <- CEATTLE_classic_SS$estimated_params$ln_mn_rec
  inits$ln_F[1:3, 1:39] <- CEATTLE_classic_SS$estimated_params$F_dev[,1:39] + CEATTLE_classic_SS$estimated_params$ln_mean_F
  inits$sel_coff[1:3,1,] <- CEATTLE_classic_SS$estimated_params$fsh_sel_coff
  BS2017SS$M1_base[,3:23] <- (BS2017SS$M1_base[,3:23] + 1e-4)


  ss_run_old_params <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = inits,      # Initial parameters
                                         estimateMode = 3,   # Dont' estimate
                                         map = map,
                                         M1Fun = build_M1(updateM1 = TRUE), # Update M1 from data
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0,        # Single species mode
                                         phase = TRUE,
                                         verbose = 0)

  # R
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$R[,1:39]), c(ss_run_old_params$quantities$R[,1:39]))

  # M
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$M[,,1:39]), c(ss_run_old_params$quantities$M_at_age[,1,,1:39]))


  # Fishery selectivity
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$fsh_sel), c(ss_run_old_params$quantities$sel_at_age[1:3,1,1:21,1]))

  # F
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$F[,,1:39]), c(ss_run_old_params$quantities$F_at_age[1:3,1,,1:39]))

  # Biomass
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$biomass[,1:39]), c(ss_run_old_params$quantities$biomass[,1:39]))

  # SSB
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$biomassSSB[,1:39]), c(ss_run_old_params$quantities$ssb[,1:39]))

  # N
  testthat::expect_equal(c(CEATTLE_classic_SS$quantities$NByage[,,1:39]), c(ss_run_old_params$quantities$N_at_age[,1,,1:39]))
}
)


testthat::test_that("Dynamics match multi-species CEATTLE classic", {

  # Load old model
  library(Rceattle)
  load(system.file("extdata/CEATTLE_classic_MS.Rdata",package="Rceattle"))

  # Load data and set up inits
  data("BS2017SS")
  data("BS2017MS")

  BS2017MS$srr_prior_mean <- 9
  BS2017MS$initMode  <- 2
  BS2017MS$fleet_control$Sel_norm_bin1 <- -999 # Normalize by max
  # BS2017MS$M1_base[,-c(1:2)] <- CEATTLE_classic_MS$quantities$M1 # + 0.0001 is in the old one

  # - Update diet data
  BS2017MS_new <- BS2017MS
  BS2017MS_new$diet_data <- as.data.frame(BS2017MS_new$diet_data)
  for(i in 1:nrow(BS2017MS_new$diet_data)){
    BS2017MS_new$diet_data$Stomach_proportion_by_weight[i] <- CEATTLE_classic_MS$data_list$UobsAge[
      BS2017MS_new$diet_data$Pred[i], BS2017MS_new$diet_data$Prey[i],
      BS2017MS_new$diet_data$Pred_age[i], BS2017MS_new$diet_data$Prey_age[i]
    ]
  }

  # - Update foraging data
  BS2017MS_new$ration_data <- as.data.frame(BS2017MS_new$ration_data)
  BS2017MS_new$ration_data[,4:ncol(BS2017MS_new$ration_data)] <- rbind(
    CEATTLE_classic_MS$data_list$Pyrs[-40,,1],
    CEATTLE_classic_MS$data_list$Pyrs[-40,,2],
    CEATTLE_classic_MS$data_list$Pyrs[-40,,3]
  )
  # inits <- suppressMessages(build_params(BS2017MS))

  # - Inits
  ms_run_old_params <- Rceattle::fit_mod(
    data_list = BS2017MS_new,
    M1Fun = build_M1(M1_model = 0,
                     updateM1 = TRUE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE),
    estimateMode = 3, # Fix at parameter values
    niter = CEATTLE_classic_MS$data_list$niter, # iterations around population and predation dynamics
    random_rec = FALSE, # No random recruitment
    msmMode = 1,  # MSVPA based
    suitMode = 0, # empirical suitability
    verbose = 0)
  inits <- ms_run_old_params$estimated_params
  map <- ms_run_old_params$map

  # - Update population dynamics from old parameters
  inits$init_dev[,1:20] <- CEATTLE_classic_MS$estimated_params$init_dev
  inits$x_tj[1:39, 1:3] <- t(CEATTLE_classic_MS$estimated_params$rec_dev)
  inits$rec_pars[,1] <- CEATTLE_classic_MS$estimated_params$ln_mn_rec
  inits$ln_F[1:3, 1:39] <- CEATTLE_classic_MS$estimated_params$F_dev[,1:39] + CEATTLE_classic_MS$estimated_params$ln_mean_F

  # -- Sel
  inits$sel_coff[1:3,1,] <- CEATTLE_classic_MS$estimated_params$fsh_sel_coff
  inits$sel_inf[1,4:6,1] <- CEATTLE_classic_MS$estimated_params$srv_sel_inf[1,]
  inits$ln_sel_slp[1,4:6,1] <- log(CEATTLE_classic_MS$estimated_params$srv_sel_slp[1,])
  inits$index_ln_q[4:6] <- CEATTLE_classic_MS$estimated_params$log_srv_q
  inits$index_ln_q[7] <- CEATTLE_classic_MS$estimated_params$log_eit_q # Need to scale by max sel because selectivity is rescaled to max = 1 in the CPP

  ms_run_old_params <- Rceattle::fit_mod(
    data_list = BS2017MS_new,
    inits = inits, # Initial parameters from old model
    map = map,
    M1Fun = build_M1(M1_model = 0,
                     updateM1 = TRUE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE),
    file = NULL, # Don't save
    estimateMode = 3, # Fix at parameter values
    niter = CEATTLE_classic_MS$data_list$niter, # iterations around population and predation dynamics
    random_rec = FALSE, # No random recruitment
    msmMode = 1,  # MSVPA based
    suitMode = 0, # empirical suitability
    verbose = 0)


  # R
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$R[,1:39]), c(ms_run_old_params$quantities$R[,1:39]))

  # M
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$M[,,1:39]), c(ms_run_old_params$quantities$M_at_age[,1,,1:39]))


  # Fishery selectivity
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$fsh_sel), c(ms_run_old_params$quantities$sel_at_age[1:3,1,1:21,1]))

  # F
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$F[,,1:39]), c(ms_run_old_params$quantities$F_at_age[1:3,1,,1:39]))

  # Biomass
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$biomass[,1:39]), c(ms_run_old_params$quantities$biomass[,1:39]))

  # SSB
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$biomassSSB[,1:39]), c(ms_run_old_params$quantities$ssb[,1:39]))

  # N
  testthat::expect_equal(c(CEATTLE_classic_MS$quantities$NByage[,,1:39]), c(ms_run_old_params$quantities$N_at_age[,1,,1:39]))
}
)
