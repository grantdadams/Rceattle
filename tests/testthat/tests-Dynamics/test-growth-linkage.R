testthat::test_that("linkage offset propagates to growth_parameters[K]", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Synthetic 2-species, 30-yr von Bertalanffy fixture (single sex).
  set.seed(123)
  nyrs <- 30
  nspp <- 2
  fmort1 <- c(seq(0.02, 0.3, length.out = nyrs / 2),
              seq(0.3,  0.05, length.out = nyrs / 2))
  fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = matrix(c(fmort1, fmort2), nspp, nyrs, byrow = TRUE),
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list

  # Linear temperature index, mean zero so the intercept stays at 0.
  yrs       <- sim_data$styr:sim_data$endyr
  yr_idx    <- seq_along(yrs)
  temp_vec  <- seq(-1, 1, length.out = length(yrs))
  sim_data$env_data <- data.frame(Year = yrs, temp = temp_vec)

  # Growth spec: vB with a temperature linkage on log_K, by species.
  growth_spec <- Rceattle::build_growth(
    fun      = "vonBertalanffy",
    linkages = list(
      log_K = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )

  # Baseline run with default beta = 0.
  base_run <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  testthat::expect_s3_class(base_run, "Rceattle")

  # The linkage table should have:
  #   2 design cols (Intercept + temp) x 2 species x 1 param = 4 rows.
  tbl <- base_run$data_list$linkage_table
  testthat::expect_equal(nrow(tbl), 4L)
  testthat::expect_setequal(tbl$param,   "log_K")
  testthat::expect_setequal(tbl$species, 1:2)

  # With ln_beta_linkage = 0 the offset tensor must be exactly zero.
  off0 <- base_run$quantities$growth_linkage_offset
  testthat::expect_true(all(off0 == 0))

  # Now perturb only the `temp` coefficient (X_col == 2) and rerun.
  inits <- base_run$estimated_params
  temp_col <- match("temp", colnames(base_run$data_list$linkage_X))
  is_temp  <- tbl$X_col == temp_col
  beta_temp <- 0.5
  inits$ln_beta_linkage <- as.numeric(inits$ln_beta_linkage)
  inits$ln_beta_linkage[is_temp] <- beta_temp

  perturbed <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    inits       = inits,
    growthFun   = growth_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # 1. The offset tensor for K must equal beta_temp * temp[yr] for each
  #    species; entries for L1 / Linf / m must remain zero.
  #    cpp dim order on the offset is [sp, sex, par, yr].
  off1 <- perturbed$quantities$growth_linkage_offset
  k_idx <- 1L
  for (sp in seq_len(nspp)) {
    sex <- 1L
    obs <- as.numeric(off1[sp, sex, k_idx, yr_idx])
    testthat::expect_equal(obs, beta_temp * temp_vec, tolerance = 1e-10)
  }
  for (par_idx in c(2L, 3L, 4L)) {
    slice <- off1[, , par_idx, , drop = FALSE]
    testthat::expect_true(all(slice == 0))
  }

  # 2. growth_parameters[K] must equal K_base * exp(beta_temp * temp)
  #    over the env-data years. cpp dim order is [sp, sex, yr, par].
  gp_base <- base_run$quantities$growth_parameters
  gp_pert <- perturbed$quantities$growth_parameters
  for (sp in seq_len(nspp)) {
    sex <- 1L
    base_k <- as.numeric(gp_base[sp, sex, yr_idx, k_idx])
    pert_k <- as.numeric(gp_pert[sp, sex, yr_idx, k_idx])
    testthat::expect_equal(pert_k / base_k, exp(beta_temp * temp_vec),
                           tolerance = 1e-8)
  }

  # 3. L1, Linf, m must be unchanged from the baseline.
  for (par_idx in c(2L, 3L, 4L)) {
    base_par <- gp_base[, , , par_idx, drop = FALSE]
    pert_par <- gp_pert[, , , par_idx, drop = FALSE]
    testthat::expect_equal(as.numeric(pert_par), as.numeric(base_par),
                           tolerance = 1e-12)
  }
})
