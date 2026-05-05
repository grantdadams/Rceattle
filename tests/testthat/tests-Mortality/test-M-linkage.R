testthat::test_that("M_linkage_offset propagates into M1_at_age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(123)
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
  yr_idx   <- seq_along(yrs)
  temp_vec <- seq(-1, 1, length.out = length(yrs))
  sim_data$env_data <- data.frame(Year = yrs, temp = temp_vec)

  # M1 spec: estimate sex/age-invariant ln_M1, plus a linkage on log_M
  # with `~ temp` by species. Use the legacy structural switch
  # (M1_model = 1) so ln_M1 has degrees of freedom; the linkage adds
  # an environmental offset on top of that.
  m1_spec <- Rceattle::build_M1(
    M1_model = 1,
    linkages = list(
      log_M = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )

  # Baseline run with default ln_beta_linkage = 0.
  base_run <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    M1Fun       = m1_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  testthat::expect_s3_class(base_run, "Rceattle")

  tbl <- base_run$data_list$linkage_table
  testthat::expect_setequal(tbl$process, "M")
  testthat::expect_setequal(tbl$param, "log_M")
  # 2 design cols (Intercept + temp) x 2 species = 4 rows; age_bin
  # stays NA so the offset broadcasts across ages.
  testthat::expect_equal(nrow(tbl), 4L)
  testthat::expect_true(all(is.na(tbl$age_bin)))

  # Step-3B verification: with no betas estimated the offset is zero.
  off0 <- base_run$quantities$M_linkage_offset
  testthat::expect_true(all(off0 == 0))

  # Now perturb only the `temp` coefficient and rerun.
  inits <- base_run$estimated_params
  temp_col <- match("temp", colnames(base_run$data_list$linkage_X))
  is_temp  <- tbl$X_col == temp_col
  beta_temp <- 0.3
  inits$ln_beta_linkage <- as.numeric(inits$ln_beta_linkage)
  inits$ln_beta_linkage[is_temp] <- beta_temp

  pert <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    inits       = inits,
    M1Fun       = m1_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # 1. The offset tensor must equal beta * temp[yr] across every
  #    (sex, age) for each species, since the row's age_bin is NA.
  off1 <- pert$quantities$M_linkage_offset
  for (sp in seq_len(nspp)) {
    for (sex in seq_len(sim_data$nsex[sp])) {
      for (age in seq_len(sim_data$nages[sp])) {
        obs <- as.numeric(off1[sp, sex, age, yr_idx])
        testthat::expect_equal(obs, beta_temp * temp_vec,
                               tolerance = 1e-10)
      }
    }
  }

  # 2. M1_at_age perturbed = M1_at_age baseline * exp(beta * temp)
  m_base <- base_run$quantities$M1_at_age
  m_pert <- pert$quantities$M1_at_age
  for (sp in seq_len(nspp)) {
    for (sex in seq_len(sim_data$nsex[sp])) {
      for (age in seq_len(sim_data$nages[sp])) {
        ratio <- as.numeric(m_pert[sp, sex, age, yr_idx]) /
          as.numeric(m_base[sp, sex, age, yr_idx])
        testthat::expect_equal(ratio, exp(beta_temp * temp_vec),
                               tolerance = 1e-8)
      }
    }
  }
})


testthat::test_that("M linkage with species-keyed priors fires per species", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(31)
  nyrs <- 25
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

  # Per-species priors on the temp coefficient. The sim helper
  # produces single-sex data, so this exercises species-only
  # keying; the nested species->sex form is unit-tested in
  # tests/testthat/tests-Linkage/test-priors.R.
  m1_spec <- Rceattle::build_M1(
    M1_model = 1,
    linkages = list(
      log_M = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species,
        priors  = list(temp = list(`1` = normal(0, 0.1),
                                   `2` = normal(0, 0.5)))
      )
    )
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    M1Fun       = m1_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  tbl <- fit$data_list$linkage_table
  temp_col <- match("temp", colnames(fit$data_list$linkage_X))
  temp_rows <- tbl[tbl$X_col == temp_col, ]

  testthat::expect_equal(temp_rows$prior_p2[temp_rows$species == 1L], 0.1)
  testthat::expect_equal(temp_rows$prior_p2[temp_rows$species == 2L], 0.5)
})


testthat::test_that("growth and M linkages compose in the same fit", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(53)
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

  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )
  m1_spec <- Rceattle::build_M1(
    M1_model = 1,
    linkages = list(
      log_M = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list   = sim_data,
    growthFun   = growth_spec,
    M1Fun       = m1_spec,
    estimateMode = 3,
    msmMode     = 0,
    random_rec  = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  tbl <- fit$data_list$linkage_table
  testthat::expect_setequal(tbl$process, c("growth", "M"))
  # 4 growth rows (2 cols x 2 sp) + 4 M rows (2 cols x 2 sp) = 8.
  testthat::expect_equal(nrow(tbl), 8L)
  testthat::expect_equal(sum(tbl$process == "growth"), 4L)
  testthat::expect_equal(sum(tbl$process == "M"),      4L)

  # Both per-process REPORT() tensors are present and zero by default.
  testthat::expect_true(all(fit$quantities$growth_linkage_offset == 0))
  testthat::expect_true(all(fit$quantities$M_linkage_offset      == 0))
})
