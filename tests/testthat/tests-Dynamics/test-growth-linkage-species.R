testthat::test_that("species-specific priors fire only on the matching rows", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(11)
  nyrs <- 30
  nspp <- 2
  fmort <- matrix(seq(0.05, 0.25, length.out = nyrs * nspp),
                  nspp, nyrs, byrow = TRUE)
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = fmort,
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  yrs      <- sim_data$styr:sim_data$endyr
  sim_data$env_data <- data.frame(
    Year = yrs,
    temp = seq(-1, 1, length.out = length(yrs))
  )

  # Different prior strength per species: tight on species 1, weak on
  # species 2.
  sd_sp1 <- 0.2
  sd_sp2 <- 1.0
  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species,
        priors  = list(temp = list(`1` = normal(0, sd_sp1),
                                   `2` = normal(0, sd_sp2)))
      )
    )
  )

  base_fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  tbl <- base_fit$data_list$linkage_table
  temp_col <- match("temp", colnames(base_fit$data_list$linkage_X))
  is_temp  <- tbl$X_col == temp_col
  testthat::expect_equal(
    tbl$prior_p2[is_temp & tbl$species == 1L], sd_sp1
  )
  testthat::expect_equal(
    tbl$prior_p2[is_temp & tbl$species == 2L], sd_sp2
  )
  # Intercept rows still have no prior.
  testthat::expect_true(
    all(tbl$prior_family[!is_temp] == "none")
  )

  # Inject a beta of 0.4 on every temp row and verify the per-species
  # prior NLL contribution matches -dnorm(0.4, 0, sd_sp_x, log = TRUE).
  inits <- base_fit$estimated_params
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  beta_temp <- 0.4
  inits$beta_linkage[is_temp] <- beta_temp

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    inits       = inits,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  expected <- -stats::dnorm(beta_temp, 0, sd_sp1, log = TRUE) +
    -stats::dnorm(beta_temp, 0, sd_sp2, log = TRUE)
  testthat::expect_equal(
    sum(fit$quantities$jnll_comp["Linkage-table priors", ]),
    expected,
    tolerance = 1e-8
  )
})


testthat::test_that("per-species formulas fit through fit_mod end-to-end", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(13)
  nyrs <- 30
  nspp <- 2
  fmort <- matrix(seq(0.05, 0.25, length.out = nyrs * nspp),
                  nspp, nyrs, byrow = TRUE)
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = fmort,
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  yrs <- sim_data$styr:sim_data$endyr
  yr_idx <- seq_along(yrs)
  sim_data$env_data <- data.frame(
    Year = yrs,
    temp = seq(-1, 1, length.out = length(yrs)),
    PDO  = stats::rnorm(length(yrs))
  )

  # Species 1: linear in temp only. Species 2: linear in temp + PDO.
  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = list(
        Rceattle::linkage_spec(formula = ~ temp,
                               by      = ~ species,
                               species = 1L),
        Rceattle::linkage_spec(formula = ~ temp + PDO,
                               by      = ~ species,
                               species = 2L)
      )
    )
  )

  base_fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  tbl <- base_fit$data_list$linkage_table
  pdo_col <- match("PDO", colnames(base_fit$data_list$linkage_X))

  # Species 1 has 2 rows (Intercept, temp); never references PDO.
  sp1 <- tbl[tbl$species == 1L, ]
  testthat::expect_equal(nrow(sp1), 2L)
  testthat::expect_false(any(sp1$X_col == pdo_col))

  # Species 2 has 3 rows (Intercept, temp, PDO).
  sp2 <- tbl[tbl$species == 2L, ]
  testthat::expect_equal(nrow(sp2), 3L)
  testthat::expect_true(any(sp2$X_col == pdo_col))

  # Inject coefficients: temp = 0.3 for both species; PDO = -0.2 for sp 2.
  inits <- base_fit$estimated_params
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  temp_col_global <- match("temp", colnames(base_fit$data_list$linkage_X))
  is_temp <- tbl$X_col == temp_col_global
  is_pdo  <- tbl$X_col == pdo_col
  inits$beta_linkage[is_temp] <- 0.3
  inits$beta_linkage[is_pdo]  <- -0.2

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    inits       = inits,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # Species 1: log_K offset = 0.3 * temp_yr (no PDO contribution).
  # Species 2: log_K offset = 0.3 * temp_yr - 0.2 * PDO_yr.
  k_idx <- 1L
  sex   <- 1L
  off <- fit$quantities$growth_linkage_offset
  obs_sp1 <- as.numeric(off[1, sex, yr_idx, k_idx])
  obs_sp2 <- as.numeric(off[2, sex, yr_idx, k_idx])
  testthat::expect_equal(obs_sp1, 0.3 * sim_data$env_data$temp,
                         tolerance = 1e-10)
  testthat::expect_equal(
    obs_sp2,
    0.3 * sim_data$env_data$temp - 0.2 * sim_data$env_data$PDO,
    tolerance = 1e-10
  )
})


testthat::test_that("per-species formulas + species-specific priors compose", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(17)
  nyrs <- 25
  nspp <- 2
  fmort <- matrix(seq(0.05, 0.25, length.out = nyrs * nspp),
                  nspp, nyrs, byrow = TRUE)

  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = fmort,
    use_size_sel = TRUE,
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )

  sim_data <- sim$data_list
  yrs <- sim_data$styr:sim_data$projyr
  sim_data$env_data <- data.frame(
    Year = yrs,
    temp = seq(0, 1, length.out = length(yrs)),
    PDO  = stats::rnorm(length(yrs))
  )

  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = list(
        Rceattle::linkage_spec(
          formula = ~ temp,
          by      = ~ species,
          species = 1L,
          init  = list(
            "(Intercept)" = log(0.3),
            temp = 1),
          priors  = list(temp = normal(0, 0.3))
        ),
        Rceattle::linkage_spec(
          formula = ~ temp + PDO,
          by      = ~ species,
          species = 2L,
          init  = list(
            "(Intercept)" = log(0.3),
            temp = 2,
            PDO = 3),
          priors  = list(
            temp = normal(0, 0.7),
            PDO  = normal(0, 1.0)
          )
        )
      )
    )
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  tbl <- fit$data_list$linkage_table
  temp_col <- match("temp", colnames(fit$data_list$linkage_X))
  pdo_col  <- match("PDO",  colnames(fit$data_list$linkage_X))

  sp1_temp <- tbl[tbl$species == 1L & tbl$X_col == temp_col, ]
  sp2_temp <- tbl[tbl$species == 2L & tbl$X_col == temp_col, ]
  sp2_pdo  <- tbl[tbl$species == 2L & tbl$X_col == pdo_col, ]

  testthat::expect_equal(sp1_temp$prior_family, "normal")
  testthat::expect_equal(sp1_temp$prior_p2, 0.3)
  testthat::expect_equal(sp2_temp$prior_family, "normal")
  testthat::expect_equal(sp2_temp$prior_p2, 0.7)
  testthat::expect_equal(sp2_pdo$prior_family, "normal")
  testthat::expect_equal(sp2_pdo$prior_p2, 1.0)


  # Model
  # - Params
  testthat::expect_equal(fit$estimated_params$beta_linkage, c(log(0.3),1,log(0.3),2,3))

  # - Prior
  priornll <- dnorm(1, 0, 0.3, log = TRUE) +  dnorm(2, 0, 0.7, log = TRUE) + dnorm(3, 0, 1, log = TRUE)
  testthat::expect_equal(sum(fit$quantities$jnll_comp[20,]), -priornll)
})



testthat::test_that("Test internal K-linked growth", {
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
  yrs <- simData$styr:simData$projyr
  simData$env_data <- data.frame(
    Year = yrs,
    temp = seq(0, 1, length.out = length(yrs)),
    PDO  = stats::rnorm(length(yrs))
  )

  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = list(
        Rceattle::linkage_spec(
          formula = ~ temp,
          by      = ~ species,
          species = 1L,
          init  = list(
            "(Intercept)" = log(0.3),
            temp = 1),
          priors  = list(temp = normal(0, 0.3))
        ),
        Rceattle::linkage_spec(
          formula = ~ temp + PDO,
          by      = ~ species,
          species = 2L,
          init  = list(
            "(Intercept)" = log(0.35),
            temp = 2,
            PDO = 3),
          priors  = list(
            temp = normal(0, 0.7),
            PDO  = normal(0, 1.0)
          )
        )
      )
    )
  )

  # Fit multi-species
  # * Fix parameters -----
  mod0 <- suppressMessages( fit_mod(data_list = simData,
                                    inits = NULL,
                                    estimateMode = 3,
                                    growthFun = growth_spec,
                                    random_rec = FALSE,
                                    msmMode = 0,
                                    fit_control = fit_control(phase = FALSE, verbose = 0)) )

  # - Params
  testthat::expect_equal(as.numeric(mod0$estimated_params$log_growth_pars[,1,1]), c(0, 0))
  testthat::expect_equal(mod0$estimated_params$beta_linkage, c(log(0.3),1,log(0.35),2,3))

  # - Map
  testthat::expect_all_true(is.na(c(mod0$map$mapList$log_growth_pars[,1,c(1,4)]))) # base k and m are off
  testthat::expect_all_true(!is.na(c(mod0$map$mapList$log_growth_pars[,1,-c(1,4)]))) # l1 and linf are on

  # - Set other inits
  inits <- mod0$estimated_params
  inits$sel_inf[1,,1] <- c(20,35,15,30)
  inits$log_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$log_F[2,] <- log(Fmort)
  inits$log_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_log_q[] <- log(1)
  inits$R_log_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$growth_log_sd[] <- log(3)
  inits$log_growth_pars[,1,] = matrix(log(c(1, 4.5, 90, 1.0,
                                           1, 4.5, 50, 1.0)), # K set to 0 b/c linkage, L1, Linf, M
                                     nrow = nspp, ncol = 4, byrow = TRUE)

  # Fit Rceattle -------------------------------------------------------------
  ss_run_init <- Rceattle::fit_mod(data_list = simData,
                                   inits = inits, # Initial parameters from sim
                                   map = mod0$map,
                                   estimateMode = 3, # Don't estimate
                                   growthFun = growth_spec, # Von bertalanffy growth
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(
                                     phase = FALSE,
                                     verbose = 1))


  #  Check growth parameters ----
  # - L1
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,2]), rep(4.5, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,2]))))

  # - Linf
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[2,,,3]), rep(50, length(as.numeric(ss_run_init$quantities$growth_parameters[2,,,3]))))
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[1,,,3]), rep(90, length(as.numeric(ss_run_init$quantities$growth_parameters[1,,,3]))))

  # - m
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,4]), rep(1, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,4]))))

  # - K (has linkage)
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_linkage_offset[1,1,,1]), log(0.3) + (simData$env_data$temp)) # Species 1
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[1,1,,1]), (0.3) * exp(simData$env_data$temp))

  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_linkage_offset[2,1,,1]), log(0.35) + (2 * simData$env_data$temp + 3 * simData$env_data$PDO)) # Species2
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[2,1,,1]), (0.35) * exp(2 * simData$env_data$temp + 3 * simData$env_data$PDO))


  # - Prior
  priornll <- dnorm(1, 0, 0.3, log = TRUE) +  dnorm(2, 0, 0.7, log = TRUE) + dnorm(3, 0, 1, log = TRUE)
  testthat::expect_equal(sum(ss_run_init$quantities$jnll_comp[20,]), -priornll)
})


testthat::test_that("Test internal L1-linked growth", {
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
  yrs <- simData$styr:simData$projyr
  simData$env_data <- data.frame(
    Year = yrs,
    temp = seq(0, 1, length.out = length(yrs)),
    PDO  = stats::rnorm(length(yrs))
  )

  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_L1 = list(
        Rceattle::linkage_spec(
          formula = ~ temp,
          by      = ~ species,
          species = 1L,
          init  = list(
            "(Intercept)" = log(4.5),
            temp = 1),
          priors  = list(temp = normal(0, 0.3))
        ),
        Rceattle::linkage_spec(
          formula = ~ temp + PDO,
          by      = ~ species,
          species = 2L,
          init  = list(
            "(Intercept)" = log(5),
            temp = 2,
            PDO = 3),
          priors  = list(
            temp = normal(0, 0.7),
            PDO  = normal(0, 1.0)
          )
        )
      )
    )
  )

  # Fit multi-species
  # * Fix parameters -----
  mod0 <- suppressMessages( fit_mod(data_list = simData,
                                    inits = NULL,
                                    estimateMode = 3,
                                    growthFun = growth_spec,
                                    random_rec = FALSE,
                                    msmMode = 0,
                                    fit_control = fit_control(phase = FALSE, verbose = 0)) )

  # - Params
  testthat::expect_equal(as.numeric(mod0$estimated_params$log_growth_pars[,1,2]), c(0, 0))
  testthat::expect_equal(mod0$estimated_params$beta_linkage, c(log(4.5),1,log(5),2,3))

  # - Map
  testthat::expect_all_true(is.na(c(mod0$map$mapList$log_growth_pars[,1,c(2,4)]))) # base k and m are off
  testthat::expect_all_true(!is.na(c(mod0$map$mapList$log_growth_pars[,1,-c(2,4)]))) # l1 and linf are on

  # - Set other inits
  inits <- mod0$estimated_params
  inits$sel_inf[1,,1] <- c(20,35,15,30)
  inits$log_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$log_F[2,] <- log(Fmort)
  inits$log_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_log_q[] <- log(1)
  inits$R_log_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$growth_log_sd[] <- log(3)
  inits$log_growth_pars[,1,] = matrix(log(c(0.3, 1, 90, 1.0,
                                           0.3, 1, 50, 1.0)), # K set to 0 b/c linkage, L1, Linf, M
                                     nrow = nspp, ncol = 4, byrow = TRUE)

  # Fit Rceattle -------------------------------------------------------------
  ss_run_init <- Rceattle::fit_mod(data_list = simData,
                                   inits = inits, # Initial parameters from sim
                                   map = mod0$map,
                                   estimateMode = 3, # Don't estimate
                                   growthFun = growth_spec, # Von bertalanffy growth
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(
                                     phase = FALSE,
                                     verbose = 1))



  # Check growth parameters ----
  # - L1
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_linkage_offset[1,1,,2]), log(4.5) + (simData$env_data$temp)) # Species1
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[1,1,,2]), (4.5) * exp(simData$env_data$temp))

  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_linkage_offset[2,1,,2]), log(5) + (2 * simData$env_data$temp + 3 * simData$env_data$PDO)) # Species2
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[2,1,,2]), (5) * exp(2 * simData$env_data$temp + 3 * simData$env_data$PDO))

  # - Linf
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[2,,,3]), rep(50, length(as.numeric(ss_run_init$quantities$growth_parameters[2,,,3]))))
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[1,,,3]), rep(90, length(as.numeric(ss_run_init$quantities$growth_parameters[1,,,3]))))

  # - m
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,4]), rep(1, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,4]))))

  # - K
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,1]), rep(0.3, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,1]))))


  # - Prior
  priornll <- dnorm(1, 0, 0.3, log = TRUE) +  dnorm(2, 0, 0.7, log = TRUE) + dnorm(3, 0, 1, log = TRUE)
  testthat::expect_equal(sum(ss_run_init$quantities$jnll_comp[20,]), -priornll)
})


testthat::test_that("Test internal Linf-linked growth", {
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
  yrs <- simData$styr:simData$projyr
  simData$env_data <- data.frame(
    Year = yrs,
    temp = seq(0, 1, length.out = length(yrs)),
    PDO  = stats::rnorm(length(yrs))
  )

  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_Linf = list(
        Rceattle::linkage_spec(
          formula = ~ temp,
          by      = ~ species,
          species = 1L,
          init  = list(
            "(Intercept)" = log(90),
            temp = 1),
          priors  = list(temp = normal(0, 0.3))
        ),
        Rceattle::linkage_spec(
          formula = ~ temp + PDO,
          by      = ~ species,
          species = 2L,
          init  = list(
            "(Intercept)" = log(50),
            temp = 2,
            PDO = 3),
          priors  = list(
            temp = normal(0, 0.7),
            PDO  = normal(0, 1.0)
          )
        )
      )
    )
  )

  # Fit multi-species
  # * Fix parameters -----
  mod0 <- suppressMessages( fit_mod(data_list = simData,
                                    inits = NULL,
                                    estimateMode = 3,
                                    growthFun = growth_spec,
                                    random_rec = FALSE,
                                    msmMode = 0,
                                    fit_control = fit_control(phase = FALSE, verbose = 0)) )

  # - Params
  testthat::expect_equal(as.numeric(mod0$estimated_params$log_growth_pars[,1,3]), c(0, 0))
  testthat::expect_equal(mod0$estimated_params$beta_linkage, c(log(90),1,log(50),2,3))

  # - Map
  testthat::expect_all_true(is.na(c(mod0$map$mapList$log_growth_pars[,1,c(3,4)]))) # base k and m are off
  testthat::expect_all_true(!is.na(c(mod0$map$mapList$log_growth_pars[,1,-c(3,4)]))) # l1 and linf are on

  # - Set other inits
  inits <- mod0$estimated_params
  inits$sel_inf[1,,1] <- c(20,35,15,30)
  inits$log_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$log_F[2,] <- log(Fmort)
  inits$log_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_log_q[] <- log(1)
  inits$R_log_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs
  inits$growth_log_sd[] <- log(3)
  inits$log_growth_pars[,1,] = matrix(log(c(0.3, 4.5, 1, 1.0,
                                           0.3, 4.5, 1, 1.0)), # K set to 0 b/c linkage, L1, Linf, M
                                     nrow = nspp, ncol = 4, byrow = TRUE)

  # Fit Rceattle -------------------------------------------------------------
  ss_run_init <- Rceattle::fit_mod(data_list = simData,
                                   inits = inits, # Initial parameters from sim
                                   map = mod0$map,
                                   estimateMode = 3, # Don't estimate
                                   growthFun = growth_spec, # Von bertalanffy growth
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(
                                     phase = FALSE,
                                     verbose = 1))

  # Check growth parameters ----
  # - L1
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,2]), rep(4.5, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,2]))))

  # - Linf
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_linkage_offset[1,1,,3]), log(90) + (simData$env_data$temp)) # Species1
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[1,1,,3]), (90) * exp(simData$env_data$temp))

  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_linkage_offset[2,1,,3]), log(50) + (2 * simData$env_data$temp + 3 * simData$env_data$PDO)) # Species2
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[2,1,,3]), (50) * exp(2 * simData$env_data$temp + 3 * simData$env_data$PDO))

  # - m
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,4]), rep(1, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,4]))))

  # - K
  testthat::expect_equal(as.numeric(ss_run_init$quantities$growth_parameters[,,,1]), rep(0.3, length(as.numeric(ss_run_init$quantities$growth_parameters[,,,1]))))

  # - Prior
  priornll <- dnorm(1, 0, 0.3, log = TRUE) +  dnorm(2, 0, 0.7, log = TRUE) + dnorm(3, 0, 1, log = TRUE)
  testthat::expect_equal(sum(ss_run_init$quantities$jnll_comp[20,]), -priornll)
})
