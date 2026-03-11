

testthat::test_that("Simulated simple model with growth curve and size-based logistic selectivity", {
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
  inits <- suppressMessages( build_params(simData) )
  inits$sel_inf[1,,1] <- c(20,35,15,30)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$growth_ln_sd[] <- log(3)
  inits$ln_growth_pars[,1,] = matrix(c(0.3, 4.5, 90, 1.0,
                                       0.35, 4.5, 50, 1.0), # K, L1, Linf, M
                                     nrow = nspp, ncol = 4, byrow = TRUE)

  # Fit Rceattle -------------------------------------------------------------
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Estimate
                              growthFun = build_growth(growth_model = 0),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = FALSE,
                              verbose = 1)
  ss_run$quantities$jnll_comp

  # Set up inits
  inits <- ss_run$estimated_params
  map = ss_run$map
  map$mapList$ln_growth_pars[1,1,] <- c(1, NA, 2, NA)
  map$mapFactor$ln_growth_pars <- factor(map$mapList$ln_growth_pars)
  inits$ln_growth_pars[1,1,] <- log(c(0.3, 4.5, 90, 1.0))
  inits$growth_ln_sd[1,1,] <- log(3)
  ss_run_init1 <- Rceattle::fit_mod(data_list = simData,
                                    inits = inits, # Initial parameters = 0
                                    map = map,
                                    file = NULL, # Don't save
                                    estimateMode = 0, # Estimate
                                    growthFun = build_growth(growth_model = 1),
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 0, # Single species mode
                                    phase = TRUE,
                                    verbose = 1)


  inits$ln_F[2,] <- log(c(seq(0.02, 0.3, length.out = 10), seq(0.3, 0.05, length.out = 10)))
  inits$rec_pars[1,1] <- log(50)
  inits$rec_dev[1,1:nyrs] <- sim$rec_devs
  inits$init_dev[1,1:(nages-2)] <- sim$init_devs
  inits$ln_sel_slp[1,,1] <- log(c(1, 0.5))
  inits$sel_inf[1,,1] <- c(30, 60)
  ss_run_init2 <- Rceattle::fit_mod(data_list = simData,
                                    inits = inits, # Initial parameters = 0
                                    map = ss_run$map,
                                    file = NULL, # Don't save
                                    estimateMode = 3, # Estimate
                                    growthFun = build_growth(growth_model = 1),
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 0, # Single species mode
                                    phase = TRUE,
                                    verbose = 1)

  ss_run_init1$quantities$jnll_comp
  ss_run_init2$quantities$jnll_comp
  sum(ss_run_init1$quantities$jnll_comp)
  sum(ss_run_init2$quantities$jnll_comp)


  plot(x = sim$SSB, y = ss_run_init1$quantities$ssb[1,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom, y = ss_run_init$quantities$biomass[1,1:nyrs], ylab = "Rceattle biomass", xlab = "True biomass"); abline(1,1)
  plot(x = sim$NAA[1:nyrs,1], y = ss_run_init$quantities$R[1,1:nyrs]); abline(1,1)

  ss_run_init$quantities$weight_hat[4,1,,1]
  sim$WAA

  # 1. Check size-selectivity ----
  fish_size_sel = 1 / (1 + exp(-0.5 * (lengths - 60)))
  srv_size_sel = 1 / (1 + exp(-1 * (lengths - 30)))

  testthat::expect_equal(as.numeric(srv_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[1,1,,1]))
  testthat::expect_equal(as.numeric(fish_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[2,1,,1]))


  # 2. Check age-selectivity ----
  testthat::expect_equal(as.numeric(sim$srv_sel), as.numeric(ss_run_init$quantities$sel_at_age[1,1,,1]))
  testthat::expect_equal(as.numeric(sim$fish_sel), as.numeric(ss_run_init$quantities$sel_at_age[2,1,,1]))

})
