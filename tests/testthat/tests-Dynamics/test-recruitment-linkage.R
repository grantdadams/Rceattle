testthat::test_that("recruitment_linkage_offset propagates into log_R0", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Single-species fixture; mean-recruitment srr_fun = 0 lets us
  # verify R = R0 * exp(rec_dev + log_R0_offset) cleanly.
  set.seed(7)
  nyrs <- 20
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 7)
  yrs    <- dat$styr:dat$endyr
  yr_idx <- seq_along(yrs)
  temp_v <- seq(-1, 1, length.out = length(yrs))
  dat$env_data <- data.frame(Year = yrs, temp = temp_v)

  # Linkage on log_R0 with srr_fun = 0 (mean only).
  rec_spec <- Rceattle::build_srr(
    srr_fun  = 0,
    linkages = list(
      log_R0 = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )

  base_run <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  testthat::expect_s3_class(base_run, "Rceattle")

  tbl <- base_run$data_list$linkage_table
  testthat::expect_setequal(tbl$process, "recruitment")
  testthat::expect_setequal(tbl$param,   "log_R0")
  # 2 design cols (Intercept + temp) x 1 species = 2 rows.
  testthat::expect_equal(nrow(tbl), 2L)

  # With beta_linkage = 0, the offset must be exactly zero.
  off0 <- base_run$quantities$recruitment_linkage_offset
  testthat::expect_true(all(off0 == 0))

  # Inject a known beta on temp; rerun in estimateMode = 3.
  inits <- base_run$estimated_params
  temp_col <- match("temp", colnames(base_run$data_list$linkage_X))
  is_temp  <- tbl$X_col == temp_col
  beta_temp <- 0.3
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  inits$beta_linkage[is_temp] <- beta_temp
  # Set rec_devs to zero to isolate the linkage effect.
  inits$rec_dev[] <- 0

  pert <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = inits,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # The offset tensor for log_R0 must equal beta * temp[yr] for the
  # only species; entries for log_alpha and log_beta remain zero.
  off1 <- pert$quantities$recruitment_linkage_offset
  obs_R0 <- as.numeric(off1[1, 1, yr_idx])  # 1 = log_R0 (cpp index 0 -> R 1)
  testthat::expect_equal(obs_R0, beta_temp * temp_v, tolerance = 1e-10)
  testthat::expect_true(all(off1[, 2, ] == 0))   # log_alpha
  testthat::expect_true(all(off1[, 3, ] == 0))   # log_beta

  # And recruitment R(sp, yr) = R0 * exp(beta * temp[yr]) for the
  # hindcast years (with rec_dev set to 0).
  base_R <- as.numeric(base_run$quantities$R[1, yr_idx])
  pert_R <- as.numeric(pert$quantities$R[1, yr_idx])
  # base_run also had rec_dev != 0 (estimated zeros). Re-fit base
  # with rec_dev = 0 so the comparison is clean.
  inits_base <- base_run$estimated_params
  inits_base$rec_dev[] <- 0
  base_dev0 <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = inits_base,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  base_R <- as.numeric(base_dev0$quantities$R[1, yr_idx])

  ratio <- pert_R / base_R
  testthat::expect_equal(ratio, exp(beta_temp * temp_v), tolerance = 1e-8)
})

testthat::test_that("linked log_R0 intercept is reflected in R_init at year 0", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 20, nages = 5, seed = 7)

  rec_spec <- Rceattle::build_srr(
    srr_fun  = 0,
    linkages = list(
      log_R0 = Rceattle::linkage_spec(formula = ~ 1, by = ~ species)
    )
  )

  base_run <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  testthat::expect_s3_class(base_run, "Rceattle")

  base_inits <- base_run$estimated_params
  base_inits$rec_dev[] <- 0
  base_inits$beta_linkage[] <- 0
  base_dev0 <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = base_inits,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  pert_inits <- base_run$estimated_params
  pert_inits$rec_dev[] <- 0
  intercept_row <- which(base_run$data_list$linkage_table$design_col == "(Intercept)")
  pert_inits$beta_linkage[intercept_row] <- 0.75
  pert <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = pert_inits,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  base_R_init <- as.numeric(base_dev0$quantities$R_init[1])
  pert_R_init <- as.numeric(pert$quantities$R_init[1])
  testthat::expect_equal(pert_R_init / base_R_init, exp(0.75), tolerance = 1e-8)
})

testthat::test_that("linked log_R0 offsets propagate into recruitment for later years", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(7)
  nyrs <- 20
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 7)
  yrs    <- dat$styr:dat$endyr
  yr_idx <- seq_along(yrs)
  temp_v <- seq(-1, 1, length.out = length(yrs))
  dat$env_data <- data.frame(Year = yrs, temp = temp_v)

  rec_spec <- Rceattle::build_srr(
    srr_fun  = 0,
    linkages = list(
      log_R0 = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )

  base_run <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  testthat::expect_s3_class(base_run, "Rceattle")

  inits <- base_run$estimated_params
  temp_col <- match("temp", colnames(base_run$data_list$linkage_X))
  is_temp  <- base_run$data_list$linkage_table$X_col == temp_col
  beta_temp <- 0.3
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  inits$beta_linkage[is_temp] <- beta_temp
  inits$rec_dev[] <- 0

  pert <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = inits,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  inits_base <- base_run$estimated_params
  inits_base$rec_dev[] <- 0
  base_dev0 <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = inits_base,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  year_to_check <- 10L
  base_R <- as.numeric(base_dev0$quantities$R[1, year_to_check])
  pert_R <- as.numeric(pert$quantities$R[1, year_to_check])
  testthat::expect_equal(pert_R / base_R, exp(beta_temp * temp_v[year_to_check]), tolerance = 1e-8)
})

testthat::test_that("recruitment linkage on log_alpha works for Beverton-Holt", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(13)
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 13)
  yrs    <- dat$styr:dat$endyr
  yr_idx <- seq_along(yrs)
  temp_v <- seq(-1, 1, length.out = length(yrs))
  dat$env_data <- data.frame(Year = yrs, temp = temp_v)

  # BH SRR with linkage on log_alpha.
  rec_spec <- Rceattle::build_srr(
    srr_fun  = 2,
    linkages = list(
      log_alpha = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )

  base_run <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  testthat::expect_s3_class(base_run, "Rceattle")

  tbl <- base_run$data_list$linkage_table
  testthat::expect_setequal(tbl$param, "log_alpha")

  # Perturb the temp coefficient and verify the offset tensor.
  inits <- base_run$estimated_params
  temp_col <- match("temp", colnames(base_run$data_list$linkage_X))
  is_temp  <- tbl$X_col == temp_col
  beta_temp <- 0.25
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  inits$beta_linkage[is_temp] <- beta_temp

  pert <- suppressMessages(Rceattle::fit_mod(
    data_list   = dat,
    inits       = inits,
    recFun      = rec_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # Offset on log_alpha matches beta * temp; log_R0 / log_beta zero.
  off <- pert$quantities$recruitment_linkage_offset
  obs_alpha <- as.numeric(off[1, 2, yr_idx])  # log_alpha is cpp index 1
  testthat::expect_equal(obs_alpha, beta_temp * temp_v, tolerance = 1e-10)
  testthat::expect_true(all(off[, 1, ] == 0))   # log_R0
  testthat::expect_true(all(off[, 3, ] == 0))   # log_beta
})


testthat::test_that("growth + M + recruitment linkages compose in one fit", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(59)
  nyrs <- 20
  nspp <- 2
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = matrix(seq(0.05, 0.25, length.out = nyrs * nspp),
                     nspp, nyrs, byrow = TRUE),
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  sim_data$env_data <- data.frame(
    Year = sim_data$styr:sim_data$endyr,
    temp = stats::rnorm(nyrs)
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    growthFun   = Rceattle::build_growth(
      fun      = "vonBertalanffy",
      linkages = list(
        log_K = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    ),
    M1Fun       = Rceattle::build_M1(
      M1_model = 1,
      linkages = list(
        log_M1 = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    ),
    recFun      = Rceattle::build_srr(
      srr_fun  = 0,
      linkages = list(
        log_R0 = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    ),
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  tbl <- fit$data_list$linkage_table
  testthat::expect_setequal(tbl$process, c("growth", "M", "recruitment"))
  # 4 growth (2 cols x 2 sp) + 4 M (2 x 2) + 4 recruitment (2 x 2) = 12
  testthat::expect_equal(nrow(tbl), 12L)

  # All three offset tensors are present and zero by default.
  testthat::expect_true(all(fit$quantities$growth_linkage_offset      == 0))
  testthat::expect_true(all(fit$quantities$M_linkage_offset           == 0))
  testthat::expect_true(all(fit$quantities$recruitment_linkage_offset == 0))
})



testthat::test_that("Test species-specific recruitment linkeages with mean rec (hindcast and projection)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  projyrs <- GOA2018SS$styr:GOA2018SS$projyr
  nyrs <- length(yrs)
  nyrsproj <- length(projyrs)
  R0 = 10:12

  # Env data
  GOA2018SS$env_data <- data.frame(Year = projyrs,
                                   EnvIndex = seq(0,1, length.out = nyrsproj),
                                   EnvIndex2 = seq(0,1, length.out = nyrsproj),
                                   EnvIndex3 = seq(0,1, length.out = nyrsproj))

  # Set params
  rec_spec <- build_srr(
    srr_fun = "mean",
    proj_mean_rec = FALSE,
    srr_est_mode = 1,
    linkages = list(
      log_R0 = Rceattle::linkage_spec(
        formula = ~ EnvIndex + EnvIndex2 + EnvIndex3,
        init = list("(Intercept)" = 10,
                    "EnvIndex" =  1,
                    "EnvIndex3" = 1),
        by = ~ species,
        priors = list("EnvIndex2" = normal(2, 0.5)),
        species = 1
      ),
      log_R0 = Rceattle::linkage_spec(
        formula = ~ EnvIndex,
        init = list("(Intercept)" = 11,
                    "EnvIndex" = 1),
        by = ~ species,
        species = 2
      ),
      log_R0 = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = 12),
        priors = list("(Intercept)" = normal(12, 0.5)),
        by = ~ species,
        species = 3
      )
    )
  )

  ss_run <- fit_mod(data_list = GOA2018SS,
                    estimateMode = 3,
                    recFun = rec_spec,
                    initMode = "Equilibrium",
                    fit_control = fit_control(verbose = 0))

  # Check coefs
  testthat::expect_equal(ss_run$estimated_params$beta_linkage, c(10, 1, 0, 1, 11, 1, 12))

  # Check ssb
  # - Species 1 (multiple)
  testthat::expect_equal(as.numeric(ss_run$quantities$R[1,]),
                         as.numeric(exp(R0[1] + as.matrix(GOA2018SS$env_data[,-1]) %*% c(1,0,1))),
                         tolerance = 0.0001)

  # - Species 2 (single)
  testthat::expect_equal(as.numeric(ss_run$quantities$R[2,]),
                         as.numeric(exp(R0[2] + as.matrix(GOA2018SS$env_data[,2]))),
                         tolerance = 0.0001)

  # - Species 3 (none)
  testthat::expect_equal(as.numeric(ss_run$quantities$R[3,]),
                         as.numeric(rep(exp(R0[3]), nyrsproj)),
                         tolerance = 0.0001)

  testthat::expect_equal(sum(ss_run$quantities$jnll_comp[20,]),
                         -dnorm(12,12,0.5, log = TRUE) - dnorm(0, 2, 0.5, log = TRUE),
                         tolerance = 0.0001)


})



testthat::test_that("Test environmental linkeage with mean rec", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- stats::rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs))

  # Set params
  rec_spec <- build_srr(
    srr_fun = "mean",
    proj_mean_rec = FALSE,
    srr_est_mode = 1,
    linkages = list(
      log_R0 = Rceattle::linkage_spec(
        formula = ~ 0 + EnvIndex,
        by = ~ species
      )
    )
  )

  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      recFun = rec_spec,
                      fit_control = fit_control(
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$rec_pars[,1] <- R0
  tbl <- ss_run$data_list$linkage_table
  env_col <- match("EnvIndex", colnames(ss_run$data_list$linkage_X))
  is_env <- tbl$param == "log_R0" & tbl$X_col == env_col
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  inits$beta_linkage[is_env] <- 1:3
  inits$log_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              recFun = rec_spec,
                              initMode = "NonEquilibrium",
                              fit_control = fit_control(
                                verbose = 0))


  # Check ssb
  for(sp in 1:3){
    testthat::expect_equal(as.numeric(ss_run$quantities$R[sp,1:nyrs]),  exp(R0[sp] + GOA2018SS$env_data$EnvIndex * sp), tolerance = 0.0001)
  }

})


testthat::test_that("Test multiple recruitment linkeages with mean rec", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- stats::rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs), EnvIndex2 = seq(0,1, length.out = nyrs), EnvIndex3 = seq(0,1, length.out = nyrs))

  # Set params
  rec_spec <- build_srr(
    srr_fun = "mean",
    proj_mean_rec = FALSE,
    srr_est_mode = 1,
    linkages = list(
      log_R0 = Rceattle::linkage_spec(
        formula = ~ 0 + EnvIndex + EnvIndex2 + EnvIndex3,
        by = ~ species
      )
    )
  )
  mod0 <- suppressMessages(
    fit_mod(data_list = GOA2018SS,
            inits = NULL,
            estimateMode = 3,
            random_rec = FALSE,
            msmMode = 0,
            recFun = rec_spec,
            initMode = 1,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$rec_pars[,1] <- R0
  tbl <- mod0$data_list$linkage_table
  env_cols <- c("EnvIndex", "EnvIndex2", "EnvIndex3")
  is_log_r0 <- tbl$param == "log_R0"
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  for (sp in 1:3) {
    for (j in seq_along(env_cols)) {
      idx <- which(is_log_r0 & tbl$species == sp & tbl$design_col == env_cols[j])
      inits$beta_linkage[idx] <- (j - 1) * 3 + sp
    }
  }
  inits$log_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              recFun = rec_spec,
                              initMode = 1,
                              fit_control = fit_control(
                                verbose = 0))


  expected_beta <- matrix(c(1,4,7,
                            2,5,8,
                            3,6,9), nrow = 3, byrow = TRUE)

  # Check ssb
  for(sp in 1:3){
    testthat::expect_equal(as.numeric(ss_run$quantities$R[sp,1:nyrs]),
                           as.numeric(exp(R0[sp] + as.matrix(GOA2018SS$env_data[,-1]) %*% expected_beta[sp,])),
                           tolerance = 0.0001)
  }
})



testthat::test_that("Beverton alpha-linked recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  projyrs <- GOA2018SS$styr:GOA2018SS$projyr
  nyrs <- length(yrs)
  nyrsproj <- length(projyrs)
  R0 = 10:12

  # Env data
  GOA2018SS$env_data <- data.frame(Year = projyrs,
                                   EnvIndex = seq(0,1, length.out = nyrsproj),
                                   EnvIndex2 = seq(0,1, length.out = nyrsproj),
                                   EnvIndex3 = seq(0,1, length.out = nyrsproj))

  # Set params
  alpha = 0.4
  beta = 1e-6

  rec_spec <- build_srr(
    srr_fun = "BevertonHolt",
    proj_mean_rec = FALSE,
    srr_est_mode = 1,
    linkages = list(
      log_beta = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = log(beta)),
        by = ~ species
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ EnvIndex + EnvIndex2 + EnvIndex3,
        init = list("(Intercept)" = log(0.3),
                    "EnvIndex" =  1,
                    "EnvIndex3" = 1),
        by = ~ species,
        priors = list("EnvIndex2" = normal(2, 0.5)),
        species = 1
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ EnvIndex,
        init = list("(Intercept)" = log(0.35),
                    "EnvIndex" = 2),
        by = ~ species,
        species = 2
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = log(0.4)),
        priors = list("(Intercept)" = normal(12, 0.5)),
        by = ~ species,
        species = 3
      )
    )
  )

  ss_run <- fit_mod(data_list = GOA2018SS,
                    estimateMode = 3,
                    recFun = rec_spec,
                    initMode = "Equilibrium",
                    fit_control = fit_control(verbose = 0))

  # Check linkage coefs
  testthat::expect_equal(ss_run$estimated_params$beta_linkage, c(log(beta), log(beta), log(beta), log(0.3), 1, 0, 1, log(0.35), 2, log(0.4)))

  # - Base params
  testthat::expect_equal(as.numeric(ss_run$estimated_params$rec_pars[,2:3]), rep(0, 6))

  # - Map
  testthat::expect_all_true(is.na(c(ss_run$map$mapList$rec_pars))) # all base pars are off
  testthat::expect_all_true(!is.na(c(ss_run$map$mapList$beta_linkage))) # linkages are on


  # Check offsets
  sp1_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[1,2,])
  sp2_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[2,2,])
  sp3_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[3,2,])

  testthat::expect_equal(sp1_offset, log(0.3) + (GOA2018SS$env_data$EnvIndex) + GOA2018SS$env_data$EnvIndex3) # Species1
  testthat::expect_equal(sp2_offset, log(0.35) + (2 * GOA2018SS$env_data$EnvIndex)) # Species2
  testthat::expect_all_true(sp3_offset == log(0.4)) # Species3
})



testthat::test_that("Ricker alpha-linked recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  projyrs <- GOA2018SS$styr:GOA2018SS$projyr
  nyrs <- length(yrs)
  nyrsproj <- length(projyrs)
  R0 = 10:12

  # Env data
  GOA2018SS$env_data <- data.frame(Year = projyrs,
                                   EnvIndex = seq(0,1, length.out = nyrsproj),
                                   EnvIndex2 = seq(0,1, length.out = nyrsproj),
                                   EnvIndex3 = seq(0,1, length.out = nyrsproj))

  # Set params
  alpha = 0.4
  beta = 1e-6

  rec_spec <- build_srr(
    srr_fun = "Ricker",
    proj_mean_rec = FALSE,
    srr_est_mode = 1,
    linkages = list(
      log_beta = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = log(beta)),
        by = ~ species
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ EnvIndex + EnvIndex2 + EnvIndex3,
        init = list("(Intercept)" = log(0.3),
                    "EnvIndex" =  1,
                    "EnvIndex3" = 1),
        by = ~ species,
        priors = list("EnvIndex2" = normal(2, 0.5)),
        species = 1
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ EnvIndex,
        init = list("(Intercept)" = log(0.35),
                    "EnvIndex" = 2),
        by = ~ species,
        species = 2
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = log(0.4)),
        priors = list("(Intercept)" = normal(12, 0.5)),
        by = ~ species,
        species = 3
      )
    )
  )

  ss_run <- fit_mod(data_list = GOA2018SS,
                    estimateMode = 3,
                    recFun = rec_spec,
                    initMode = "Equilibrium",
                    fit_control = fit_control(verbose = 0))

  # Check linkage coefs
  testthat::expect_equal(ss_run$estimated_params$beta_linkage, c(log(beta), log(beta), log(beta), log(0.3), 1, 0, 1, log(0.35), 2, log(0.4)))

  # - Base params
  testthat::expect_equal(as.numeric(ss_run$estimated_params$rec_pars[,2:3]), rep(0, 6))

  # - Map
  testthat::expect_all_true(is.na(c(ss_run$map$mapList$rec_pars))) # all base pars on off
  testthat::expect_all_true(!is.na(c(ss_run$map$mapList$beta_linkage))) # linkages are on


  # Check offsets
  sp1_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[1,2,])
  sp2_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[2,2,])
  sp3_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[3,2,])

  testthat::expect_equal(sp1_offset, log(0.3) + (GOA2018SS$env_data$EnvIndex) + GOA2018SS$env_data$EnvIndex3) # Species1
  testthat::expect_equal(sp2_offset, log(0.35) + (2 * GOA2018SS$env_data$EnvIndex)) # Species2
  testthat::expect_all_true(sp3_offset == log(0.4)) # Species3
})



testthat::test_that("Beverton alpha-linked recruitment (penalty approach)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  projyrs <- GOA2018SS$styr:GOA2018SS$projyr
  nyrs <- length(yrs)
  nyrsproj <- length(projyrs)
  R0 = 10:12

  # Env data
  GOA2018SS$env_data <- data.frame(Year = projyrs,
                                   EnvIndex = seq(0,1, length.out = nyrsproj),
                                   EnvIndex2 = seq(0,1, length.out = nyrsproj),
                                   EnvIndex3 = seq(0,1, length.out = nyrsproj))

  # Set params
  alpha = 0.4
  beta = 1e-6

  rec_spec <- build_srr(
    srr_fun = "mean",
    srr_pred_fun = "BevertonHolt",
    proj_mean_rec = FALSE,
    srr_est_mode = 1,
    linkages = list(
      log_beta = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = log(beta)),
        by = ~ species
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ EnvIndex + EnvIndex2 + EnvIndex3,
        init = list("(Intercept)" = log(0.3),
                    "EnvIndex" =  1,
                    "EnvIndex3" = 1),
        by = ~ species,
        priors = list("EnvIndex2" = normal(2, 0.5)),
        species = 1
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ EnvIndex,
        init = list("(Intercept)" = log(0.35),
                    "EnvIndex" = 2),
        by = ~ species,
        species = 2
      ),
      log_alpha = Rceattle::linkage_spec(
        formula = ~ 1,
        init = list("(Intercept)" = log(0.4)),
        priors = list("(Intercept)" = normal(12, 0.5)),
        by = ~ species,
        species = 3
      )
    )
  )

  ss_run <- fit_mod(data_list = GOA2018SS,
                    estimateMode = 3,
                    recFun = rec_spec,
                    initMode = "Equilibrium",
                    fit_control = fit_control(verbose = 0))

  # Check linkage coefs
  testthat::expect_equal(ss_run$estimated_params$beta_linkage, c(log(beta), log(beta), log(beta), log(0.3), 1, 0, 1, log(0.35), 2, log(0.4)))

  # - Base params
  testthat::expect_equal(as.numeric(ss_run$estimated_params$rec_pars[,2:3]), rep(0, 6))

  # - Map
  testthat::expect_all_true(is.na(c(ss_run$map$mapList$rec_pars[,2:3]))) # base alpha and beta are off
  testthat::expect_all_true(!is.na(c(ss_run$map$mapList$rec_pars[,1]))) # base R0 is on
  testthat::expect_all_true(!is.na(c(ss_run$map$mapList$beta_linkage))) # linkages are on


  # Check offsets
  sp1_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[1,2,])
  sp2_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[2,2,])
  sp3_offset <-as.numeric(ss_run$quantities$recruitment_linkage_offset[3,2,])

  testthat::expect_equal(sp1_offset, log(0.3) + (GOA2018SS$env_data$EnvIndex) + GOA2018SS$env_data$EnvIndex3) # Species1
  testthat::expect_equal(sp2_offset, log(0.35) + (2 * GOA2018SS$env_data$EnvIndex)) # Species2
  testthat::expect_all_true(sp3_offset == log(0.4)) # Species3

  # Check there is a likelihood penalty
  testthat::expect_true(sum(ss_run$quantities$jnll_comp[12,]) != 0)
})


