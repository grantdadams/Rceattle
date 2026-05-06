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
  obs_sp1 <- as.numeric(off[1, sex, k_idx, yr_idx])
  obs_sp2 <- as.numeric(off[2, sex, k_idx, yr_idx])
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
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  yrs <- sim_data$styr:sim_data$endyr
  sim_data$env_data <- data.frame(
    Year = yrs,
    temp = seq(-1, 1, length.out = length(yrs)),
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
          priors  = list(temp = normal(0, 0.3))
        ),
        Rceattle::linkage_spec(
          formula = ~ temp + PDO,
          by      = ~ species,
          species = 2L,
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
})
