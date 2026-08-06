# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("build_bounds honors lower/upper from linkage_table", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(42)
  nyrs <- 20
  nspp <- 2
  fmort1 <- seq(0.05, 0.25, length.out = nyrs)
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = matrix(rep(fmort1, nspp), nspp, nyrs, byrow = TRUE),
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  yrs <- sim_data$styr:sim_data$endyr
  sim_data$env_data <- data.frame(
    Year = yrs,
    temp = seq(-1, 1, length.out = length(yrs))
  )

  spec <- Rceattle::linkage_spec(
    formula = ~ temp,
    by      = ~ species,
    bounds  = list(temp = c(-1.5, 1.5))
  )
  growth_spec <- Rceattle::build_growth(
    fun      = "vonBertalanffy",
    linkages = list(K = spec)
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
  is_temp <- tbl$X_col == temp_col

  testthat::expect_equal(unique(fit$bounds$lower$beta_linkage[is_temp]),
                         -1.5)
  testthat::expect_equal(unique(fit$bounds$upper$beta_linkage[is_temp]),
                         1.5)
  # Intercept rows have no user-supplied bound so they stay at the
  # finite values pulled from the linkage_table's own defaults.
  intercept_lower <- fit$bounds$lower$beta_linkage[!is_temp]
  testthat::expect_true(all(is.finite(intercept_lower) |
                            is.infinite(intercept_lower)))
})


testthat::test_that("linkage_spec validates init and bounds as named lists", {
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, init = 10),
    "`init` must be a named list keyed by design-matrix column name"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, bounds = c(-1, 1)),
    "`bounds` must be a named list keyed by design-matrix column name"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, init = list(temp = c(1, 2))),
    "init\\$temp must be a numeric scalar"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, bounds = list(temp = 1)),
    "bounds\\$temp must be a numeric vector of length 2"
  )
})


testthat::test_that("normal prior on beta_linkage adds expected NLL", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(7)
  nyrs <- 20
  nspp <- 2
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = matrix(seq(0.05, 0.25, length.out = nyrs * nspp),
                     nspp, nyrs, byrow = TRUE),
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  yrs <- sim_data$styr:sim_data$endyr
  sim_data$env_data <- data.frame(
    Year = yrs,
    temp = seq(-1, 1, length.out = length(yrs))
  )

  # Two builds: one with no prior, one with a normal(0, 0.5) prior on temp.
  growth_no_prior <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      K = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )
  growth_with_prior <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      K = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species,
        priors  = list(temp = normal(0, 0.5))
      )
    )
  )

  # Fit both with the same nonzero beta_linkage so the prior fires.
  base <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, growthFun = growth_no_prior,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  inits <- base$estimated_params
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  tbl <- base$data_list$linkage_table
  temp_col <- match("temp", colnames(base$data_list$linkage_X))
  is_temp <- tbl$X_col == temp_col
  beta_temp <- 0.4
  inits$beta_linkage[is_temp] <- beta_temp

  no_prior <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, growthFun = growth_no_prior,
    inits = inits,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))
  with_prior <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, growthFun = growth_with_prior,
    inits = inits,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # Slot 20 (index 20 in R, 19 in cpp) holds linkage prior contributions.
  jnll_no_prior   <- no_prior$quantities$jnll_comp
  jnll_with_prior <- with_prior$quantities$jnll_comp
  testthat::expect_equal(sum(jnll_no_prior[20, ]), 0)
  testthat::expect_gt(sum(jnll_with_prior[20, ]), 0)

  # Total NLL must match: no_prior_total + sum(prior_contrib) ==
  # with_prior_total (the only other change is the prior block).
  delta_nll <- sum(jnll_with_prior) - sum(jnll_no_prior)
  testthat::expect_equal(delta_nll, sum(jnll_with_prior[20, ]),
                         tolerance = 1e-6)

  # Analytical check: each species' temp-coefficient row contributes
  # -dnorm(0.4, 0, 0.5, log = TRUE) to slot 19 (cpp slot is 0-based).
  expected_per_row <- -stats::dnorm(beta_temp, 0, 0.5, log = TRUE)
  expected_total   <- sum(is_temp) * expected_per_row
  testthat::expect_equal(sum(jnll_with_prior[20, ]), expected_total,
                         tolerance = 1e-8)
})


testthat::test_that("sd_L1 intercept init lands on growth_log_sd, not log_growth_pars", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  set.seed(11)
  nyrs <- 15
  nspp <- 2
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = matrix(seq(0.05, 0.25, length.out = nyrs * nspp),
                     nspp, nyrs, byrow = TRUE),
    log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)
  )
  sim_data <- sim$data_list
  yrs <- sim_data$styr:sim_data$endyr
  # Need a `temp` column even though the SD spec uses ~ 1 -- the test
  # data uses Ceq = 4 (temperature-dependent consumption) which reads
  # env_data$temp downstream of the linkage system.
  sim_data$env_data <- data.frame(Year = yrs, temp = 0)

  # init is supplied on the natural scale of the SD (4.4.0 contract);
  # build_params() applies log() before writing into growth_log_sd.
  init_sd <- 5
  growth_spec <- Rceattle::build_growth(
    fun      = "vonBertalanffy",
    linkages = list(
      sd_L1 = Rceattle::linkage_spec(
        formula = ~ 1,
        init    = list(`(Intercept)` = init_sd)
      )
    )
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list    = sim_data,
    growthFun    = growth_spec,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # Init must have landed in growth_log_sd[, , 1] (SD at L1) on the log
  # scale, NOT in log_growth_pars (which would corrupt K).
  testthat::expect_equal(
    as.numeric(fit$estimated_params$growth_log_sd[, , 1]),
    rep(log(init_sd), nspp * dim(fit$estimated_params$growth_log_sd)[2])
  )
  # K (first slot of log_growth_pars) keeps its build_params default
  # (log(0.3)), unaffected by the SD init.
  testthat::expect_equal(
    as.numeric(fit$estimated_params$log_growth_pars[, , 1]),
    rep(log(0.3), nspp * dim(fit$estimated_params$log_growth_pars)[2])
  )

  # And the intercept row in the linkage table is still mapped out at 0
  # (the level is carried by the base parameter, not by beta_linkage).
  tbl <- fit$data_list$linkage_table
  sd_rows <- tbl$param == "sd_L1"
  testthat::expect_true(all(tbl$design_col[sd_rows] == "(Intercept)"))
  testthat::expect_true(all(is.na(fit$map$beta_linkage[sd_rows])))
})


testthat::test_that("a linkage coefficient's bounds are applied to that coefficient", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Regression test for a silent misalignment: fit_mod assembles lower/upper by
  # walking build_params()'s parameter order, but TMB orders obj$par by the
  # sequence the PARAMETER_* macros appear in the template -- and build_params()
  # lists the linkage coefficients AFTER log_F while ceattle.cpp declares them
  # BEFORE it. The two vectors are the same LENGTH either way, so nothing
  # downstream can notice; the constraints simply land on the wrong parameters,
  # putting the user's coefficient bounds on log_F and log_F's [-1000, 10] rail
  # on the coefficient. Any model carrying a linkage was affected.
  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$endyr
  d$env_data <- data.frame(Year = yrs, temp = as.numeric(scale(seq_along(yrs))))

  cap <- 0.02
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    qFun = Rceattle::build_catchability(linkages = list(
      q = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 7L,
                                 bounds = list(temp = c(-cap, cap))))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  beta <- as.numeric(fit$estimated_params$beta_linkage)
  beta <- beta[fit$data_list$linkage_table$design_col == "temp"]
  testthat::expect_length(beta, 1L)
  testthat::expect_lte(beta, cap + 1e-8)
  testthat::expect_gte(beta, -cap - 1e-8)
})
