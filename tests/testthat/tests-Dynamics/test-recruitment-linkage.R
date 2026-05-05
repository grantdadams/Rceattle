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

  # With ln_beta_linkage = 0, the offset must be exactly zero.
  off0 <- base_run$quantities$recruitment_linkage_offset
  testthat::expect_true(all(off0 == 0))

  # Inject a known beta on temp; rerun in estimateMode = 3.
  inits <- base_run$estimated_params
  temp_col <- match("temp", colnames(base_run$data_list$linkage_X))
  is_temp  <- tbl$X_col == temp_col
  beta_temp <- 0.3
  inits$ln_beta_linkage <- as.numeric(inits$ln_beta_linkage)
  inits$ln_beta_linkage[is_temp] <- beta_temp
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
  inits$ln_beta_linkage <- as.numeric(inits$ln_beta_linkage)
  inits$ln_beta_linkage[is_temp] <- beta_temp

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
