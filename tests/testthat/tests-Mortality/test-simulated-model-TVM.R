testthat::test_that("Test IID year time-varying M", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  M1rho = 0
  M1sd = 0.1
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Time-varying M
    rhoM = M1rho,
    sigmaM = M1sd,

    # Multispecies bits
    log_phi = log_phi # set to -Inf so no predation
  )

  # sum(sim$model_quantities$B_eaten_as_prey)

  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  inits <- suppressMessages(build_params(simData))
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$ln_M1[,1,] <- log(sim$model_quantities$M)
  for(sp in 1:2){
    for(age in 1:15){
      inits$ln_M1_dev[sp,1,age,1:nyrs] <- sim$model_quantities$M_vec[sp,]
    }
  }

  inits$M1_rho[] <- atanh(M1rho)
  inits$M1_dev_ln_sd[] <- log(M1sd)
  simData$M1_model = 1
  simData$M1_re = 2

  ss_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               M1Fun = build_M1(M1_model = 1, # Estimable M
                                                M1_re = 2), # IID year
                               phase = FALSE,
                               msmMode = 0,
                               initMode = 2,
                               verbose = 0)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ss_run1$quantities$R[,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ss_run1$quantities$biomass[,1:nyrs]), tolerance = 1e-6)

  # M
  testthat::expect_equal(as.numeric(sim$model_quantities$MAA), as.numeric(ss_run1$quantities$M_at_age[,1,,1:nyrs]), tolerance = 1e-6)

  # M map
  # - Devs
  testthat::expect_equal(matrix(ss_run1$map$mapList$ln_M1_dev[1,1,,], 15, 40), matrix(c(1:30, rep(NA, 10)), 15, 40, byrow = TRUE))
  testthat::expect_equal(matrix(ss_run1$map$mapList$ln_M1_dev[2,1,,], 15, 40), matrix(c(31:60, rep(NA, 10)), 15, 40, byrow = TRUE))

  # - SD and rho
  testthat::expect_equal(as.numeric(ss_run1$map$mapList$M1_rho[,1,]), as.numeric(rep(NA, 4)))
  testthat::expect_equal(matrix(ss_run1$map$mapList$M1_dev_ln_sd, 2, 1), matrix(1:2, 2, 1))

  # SD and rho
  ss_run1$quantities$rho_M_y
  ss_run1$quantities$Sigma_M

  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ss_run1$quantities$N_at_age[,1,,1:nyrs]))

  # Likelihood
  expected_nll <- calc_nll_ar1_1d(x = sim$model_quantities$M_vec[1,], sd = M1sd, rho = M1rho)
  obs_nll <- ss_run1$quantities$jnll_comp[16,1]

  testthat::expect_equal(expected_nll, obs_nll)

})



testthat::test_that("Test AR1 year time-varying M", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  M1rho = 0.9
  M1sd = 0.1
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Time-varying M
    rhoM = M1rho,
    sigmaM = M1sd,

    # Multispecies bits
    log_phi = log_phi # set to -Inf so no predation
  )

  # sum(sim$model_quantities$B_eaten_as_prey)

  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  inits <- suppressMessages(build_params(simData))
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$ln_M1[,1,] <- log(sim$model_quantities$M)
  for(sp in 1:2){
    for(age in 1:15){
      inits$ln_M1_dev[sp,1,age,1:nyrs] <- sim$model_quantities$M_vec[sp,]
    }
  }

  inits$M1_rho[] <- atanh(M1rho)
  inits$M1_dev_ln_sd[] <- log(M1sd)
  simData$M1_model = 1
  simData$M1_re = 5

  ss_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               M1Fun = build_M1(M1_model = 1, # Estimable M
                                                M1_re = 5), # AR1 year
                               phase = FALSE,
                               msmMode = 0,
                               initMode = 2,
                               verbose = 0)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ss_run1$quantities$R[,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ss_run1$quantities$biomass[,1:nyrs]), tolerance = 1e-6)

  # M
  testthat::expect_equal(as.numeric(sim$model_quantities$MAA), as.numeric(ss_run1$quantities$M_at_age[,1,,1:nyrs]), tolerance = 1e-6)

  # M map
  # - Devs
  testthat::expect_equal(matrix(ss_run1$map$mapList$ln_M1_dev[1,1,,], 15, 40), matrix(c(1:30, rep(NA, 10)), 15, 40, byrow = TRUE))
  testthat::expect_equal(matrix(ss_run1$map$mapList$ln_M1_dev[2,1,,], 15, 40), matrix(c(31:60, rep(NA, 10)), 15, 40, byrow = TRUE))

  # - SD and rho
  testthat::expect_equal(matrix(ss_run1$map$mapList$M1_rho, 2, 2), matrix(c(rep(NA, 2), 1:2), 2, 2))
  testthat::expect_equal(matrix(ss_run1$map$mapList$M1_dev_ln_sd, 2, 1), matrix(1:2, 2, 1))


  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ss_run1$quantities$N_at_age[,1,,1:nyrs]))

  # Likelihood
  expected_nll <- calc_nll_ar1_1d(x = sim$model_quantities$M_vec[1,], sd = M1sd, rho = M1rho)
  obs_nll <- ss_run1$quantities$jnll_comp[16,1]

  testthat::expect_equal(expected_nll, obs_nll)
})
