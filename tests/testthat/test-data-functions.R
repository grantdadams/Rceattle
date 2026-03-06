
testthat::test_that("Test combine single-species data", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  # Make copy of data
  dat <- Rceattle::combine_data(dat, dat)

  # Set params
  inits <- suppressMessages(build_params(dat))
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd[] <- 0
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1


  # Set logistic params
  inits$ln_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  # Calculate SPR
  M <- dat$M1_base$Age1[1]
  wt <- dat$weight$Age1[1]
  sex_ratio <- dat$sex_ratio$Age1[1]
  pmature <- dat$maturity$Age1[1]

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0)

  # Check biomass
  testthat::expect_equal(as.numeric(ss_run$quantities$biomass[1,1:nyrs]), rep(SB0, nyrs), tolerance = 0.0001)
  testthat::expect_equal(as.numeric(ss_run$quantities$biomass[2,1:nyrs]), rep(SB0, nyrs), tolerance = 0.0001)
})



testthat::test_that("Test combine multi-species data", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )

  # Combine data
  simData <- Rceattle::combine_data(sim$data_list, sim$data_list)


  # * Fix parameters -----
  inits <- suppressMessages( build_params(simData) )
  inits$log_gam_a <- rep(log(gam_a), 2)
  inits$log_gam_b <- rep(log(gam_b), 2)
  inits$log_phi <- diag(2) %x% log_phi
  inits$log_phi[inits$log_phi == 0] <- -Inf
  inits$sel_inf[1,,1] <- rep(c(3,6,2.5,4), 2)
  inits$ln_sel_slp[1,,1] <- rep(log(c(2,2.5,2,2.5)), 2)
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$ln_F[6,] <- log(Fmort)
  inits$ln_F[8,] <- log(Fmort2)
  inits$rec_pars[,1] <- rep(log(c(1e2, 1e3)), 2)
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$rec_dev[1:2,1:30] <- sim$model_quantities$rec_devs
  inits$rec_dev[3:4,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[1:2,1:14] <- sim$model_quantities$init_devs
  inits$init_dev[3:4,1:14] <- sim$model_quantities$init_devs

  # Fit multi-species
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits,
                               estimateMode = 3, # Don't estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)


  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ms_run1$quantities$R[1:2,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ms_run1$quantities$R[3:4,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ms_run1$quantities$biomass[1:2,1:nyrs]), tolerance = 1e-6)
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ms_run1$quantities$biomass[3:4,1:nyrs]), tolerance = 1e-6)

  # Ration
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ms_run1$quantities$consumption_at_age[1:2,1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ms_run1$quantities$consumption_at_age[3:4,1,,1]))

  # Selectivity
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]), as.numeric(ms_run1$quantities$sel_at_age[c(1,3),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]), as.numeric(ms_run1$quantities$sel_at_age[c(5,7),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]), as.numeric(ms_run1$quantities$sel_at_age[c(2,4),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]), as.numeric(ms_run1$quantities$sel_at_age[c(6,8),,,1]))

  # F
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA), as.numeric(ms_run1$quantities$F_flt_age[c(2,4),1,,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA), as.numeric(ms_run1$quantities$F_flt_age[c(6,8),1,,1:nyrs]))

  # Q
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q), as.numeric(ms_run1$quantities$index_q[c(1,3),1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q), as.numeric(ms_run1$quantities$index_q[c(5,7),1]))

  # Suitability
  testthat::expect_equal(exp(ms_run1$estimated_params$log_gam_a), rep(gam_a, 2))
  testthat::expect_equal(exp(ms_run1$estimated_params$log_gam_b), rep(gam_b, 2))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run1$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run1$quantities$diet_hat[,2]))
})
