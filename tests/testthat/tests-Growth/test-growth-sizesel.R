testthat::test_that("Test internal VB growth and length-based logistic selectivity", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

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

    # Growth
    use_size_sel = TRUE,
    fish_CAAL_ISS = 1e6,
    srv_CAAL_ISS = 1e6,

    # Multispecies bits
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list

  # Fit multi-species
  # * Fix parameters -----
  mod0 <- suppressMessages( fit_mod(data_list = simData,
                                    inits = NULL,
                                    estimateMode = 3,
                                    growthFun = build_growth(fun = "vonBertalanffy"),
                                    random_rec = FALSE, msmMode = 0,
                                    fit_control = fit_control(phase = FALSE, verbose = 0)) )
  inits <- mod0$estimated_params
  inits$sel_inf[1,,1] <- c(20,35,15,30)
  inits$log_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$log_F[2,] <- log(Fmort)
  inits$log_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_log_q[] <- log(1)
  inits$x_tj[1:30, 1:simData$nspp] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$growth_log_sd[] <- log(3)
  inits$log_growth_pars[,1,] = matrix(log(c(0.3, 4.5, 90, 1.0,
                                           0.35, 4.5, 50, 1.0)), # K, L1, Linf, M
                                     nrow = nspp, ncol = 4, byrow = TRUE)

  # Fit Rceattle -------------------------------------------------------------
  ss_run_init <- Rceattle::fit_mod(data_list = simData,
                                   inits = inits, # Initial parameters from sim
                                   map = mod0$map,
                                   estimateMode = 3, # Don't estimate
                                   growthFun = build_growth(fun = "vonBertalanffy"), # Von bertalanffy growth
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(
                                     phase = FALSE,
                                     verbose = 1))


  # 1. Check size-selectivity ----
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[c(1,3),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[c(2,4),1,,1]))


  # 2. Check growth ----
  # - Biomass weight
  testthat::expect_equal(as.numeric(sim$model_quantities$length_at_age), as.numeric(ss_run_init$quantities$length_hat[c(1,3),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$growth_matrix), as.numeric(ss_run_init$quantities$growth_matrix[c(1,3),1,,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$WAA), as.numeric(ss_run_init$quantities$weight_hat[c(1,3),1,,1]))

  # - Spawn weight
  testthat::expect_equal(as.numeric(sim$model_quantities$length_at_age), as.numeric(ss_run_init$quantities$length_hat[c(2,4),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$growth_matrix), as.numeric(ss_run_init$quantities$growth_matrix[c(2,4),1,,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$WAA), as.numeric(ss_run_init$quantities$weight_hat[c(2,4),1,,1]))


  # 3. Check age-selectivity ----
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel), as.numeric(ss_run_init$quantities$sel_at_age[c(1,3),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel), as.numeric(ss_run_init$quantities$sel_at_age[c(2,4),1,,1]))


  # 4. Check biomass and N ----
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ss_run_init$quantities$biomass[,1:30]), tolerance = 1e-6)
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ss_run_init$quantities$N_at_age[,1,,1:30]))
})


testthat::test_that("Test Richard's growth and length-based logistic selectivity", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

  growth_params = matrix(c(0.3, 4.5, 90, 1.2,
                           0.35, 4.5, 50, 1.3), # K, L1, Linf, M
                         nrow = nspp, ncol = 4, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Growth
    growth_params = growth_params,
    growth_model = 2,
    use_size_sel = TRUE,
    fish_CAAL_ISS = 1e6,
    srv_CAAL_ISS = 1e6,

    # Multispecies bits
    log_phi = log_phi
  )

  # Set up Rceattle data
  simData <- sim$data_list

  # Fit multi-species
  # * Fix parameters -----
  mod0 <- suppressMessages( fit_mod(data_list = simData, inits = NULL, estimateMode = 3, growthFun = build_growth(fun = "Richards"), random_rec = FALSE, msmMode = 0, fit_control = fit_control(phase = FALSE, verbose = 0)) )
  inits <- mod0$estimated_params
  inits$sel_inf[1,,1] <- c(20,35,15,30)
  inits$log_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$log_F[2,] <- log(Fmort)
  inits$log_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_log_q[] <- log(1)
  inits$x_tj[1:30, 1:simData$nspp] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$growth_log_sd[] <- log(3)
  inits$log_growth_pars[,1,] = log(growth_params)

  # Fit Rceattle -------------------------------------------------------------
  ss_run_init <- Rceattle::fit_mod(data_list = simData,
                                   inits = inits, # Initial parameters from sim
                                   map = mod0$map,
                                   estimateMode = 3, # Don't estimate
                                   growthFun = build_growth(fun = "Richards"), # Richard's growth
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(
                                     phase = FALSE,
                                     verbose = 1))


  # 1. Check size-selectivity ----
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[c(1,3),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[c(2,4),1,,1]))


  # 2. Check growth ----
  # - Biomass weight
  testthat::expect_equal(as.numeric(sim$model_quantities$length_at_age), as.numeric(ss_run_init$quantities$length_hat[c(1,3),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$growth_matrix), as.numeric(ss_run_init$quantities$growth_matrix[c(1,3),1,,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$WAA), as.numeric(ss_run_init$quantities$weight_hat[c(1,3),1,,1]))

  # - Spawn weight
  testthat::expect_equal(as.numeric(sim$model_quantities$length_at_age), as.numeric(ss_run_init$quantities$length_hat[c(2,4),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$growth_matrix), as.numeric(ss_run_init$quantities$growth_matrix[c(2,4),1,,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$WAA), as.numeric(ss_run_init$quantities$weight_hat[c(2,4),1,,1]))


  # 3. Check age-selectivity ----
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel), as.numeric(ss_run_init$quantities$sel_at_age[c(1,3),1,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel), as.numeric(ss_run_init$quantities$sel_at_age[c(2,4),1,,1]))


  # 4. Check biomass and N ----
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ss_run_init$quantities$biomass[,1:30]), tolerance = 1e-6)
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ss_run_init$quantities$N_at_age[,1,,1:30]))
})
