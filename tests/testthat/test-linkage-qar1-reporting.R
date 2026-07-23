# =============================================================================
# Formula-linkage effect sizes are reported. The estimated linkage coefficients
# (beta_linkage), the random-effect deviations (beta_linkage_re), and the
# Rogers-2024 QAR1 effect size (beta_linkage_obs) are REPORT'd (so they appear in
# fit$quantities) and, for the effect sizes, ADREPORT'd (so they carry a standard
# error in the sdreport). Previously the effect size -- the whole point of the
# QAR1 model -- was buried in the raw parameter vector with no readable exposure.
# =============================================================================

testthat::test_that("a QAR1 q-linkage reports its effect size and deviations", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  d <- Rceattle::BS2017SS
  nyr <- d$endyr - d$styr + 1
  d$env_data <- data.frame(Year = d$styr:d$endyr,
                           qcov = as.numeric(scale(seq_len(nyr))))
  qflt <- which(d$fleet_control$Catchability %in% c(1L, 2L))[1]

  cf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = qflt,
                               observe = "qcov", obs_sd = 0.1)))
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, qFun = cf, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  q <- fit$quantities
  testthat::expect_true(all(c("beta_linkage", "beta_linkage_re", "beta_linkage_obs")
                            %in% names(q)))
  testthat::expect_equal(length(q$beta_linkage_obs), 1L)   # one effect size (one observed group)
  testthat::expect_equal(length(q$beta_linkage_re), nyr)   # one AR1 deviate per hindcast year
  testthat::expect_true(is.finite(q$beta_linkage_obs))
})

testthat::test_that("QAR1 observation SD is held fixed at the input obs_sd", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  d <- Rceattle::BS2017SS
  nyr <- d$endyr - d$styr + 1
  set.seed(42)
  d$env_data <- data.frame(Year = d$styr:d$endyr, qcov = stats::rnorm(nyr))
  qflt <- which(d$fleet_control$Catchability %in% c(1L, 2L))[1]
  cf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = qflt,
                               observe = "qcov", obs_sd = 0.5)))
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, qFun = cf, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  # One per observed group, held FIXED at the input obs_sd (mapped out, so its
  # value stays put and it is not among the free parameters).
  testthat::expect_equal(length(fit$estimated_params$log_obs_sd_linkage), 1L)
  testthat::expect_equal(exp(fit$estimated_params$log_obs_sd_linkage), 0.5,
                         tolerance = 1e-8)
  testthat::expect_false("log_obs_sd_linkage" %in% names(fit$obj$par))
})

testthat::test_that(".extend_env_data prepends to styr with NA, leaves complete/malformed data alone", {
  testthat::skip_if_not_installed("Rceattle")
  # Late-starting env_data: prepend styr..first-1 with NA.
  out <- Rceattle:::.extend_env_data(data.frame(Year = 1985:1990, x = 1:6), 1979L)
  testthat::expect_equal(out$Year, 1979:1990)
  testthat::expect_true(all(is.na(out$x[out$Year < 1985])))
  testthat::expect_equal(out$x[out$Year >= 1985], as.numeric(1:6))
  # Already contiguous from styr: unchanged.
  ok <- data.frame(Year = 1979:1985, x = 1:7)
  testthat::expect_equal(Rceattle:::.extend_env_data(ok, 1979L), ok)
  # Starts BEFORE styr / unsorted: left for .check_env_data_years to reject.
  before <- data.frame(Year = 1975:1985)
  testthat::expect_equal(Rceattle:::.extend_env_data(before, 1979L), before)
})

testthat::test_that("QAR1 observation skips years absent from env_data (mask)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  # covariate present only for the second half of the hindcast
  yrs <- d$styr:d$endyr
  keep <- yrs[yrs >= stats::median(yrs)]
  d$env_data <- data.frame(Year = keep, qcov = as.numeric(scale(seq_along(keep))))
  d$env_data <- Rceattle:::.extend_env_data(d$env_data, d$styr)
  qflt <- which(d$fleet_control$Catchability %in% c(1L, 2L))[1]
  cf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = qflt,
                               observe = "qcov", obs_sd = 0.5)))
  pool <- Rceattle:::pool_linkages(
    spec_groups = list(q = cf$linkages), env_data = d$env_data,
    strata = list(fleet = seq_len(nrow(d$fleet_control)), species = 1L,
                  sex = 1L, age_bin = seq_len(d$nages[1])))
  enc <- Rceattle:::encode_linkage_for_tmb(pool$table, pool$X)
  # one latent per hindcast year; observations only on the present years.
  testthat::expect_equal(length(enc$linkage_re_obs_mask), length(yrs))
  testthat::expect_equal(sum(enc$linkage_re_obs_mask), length(keep))
})

testthat::test_that("a fixed-effect covariate with missing years errors (not silent NaN)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  # a fixed-effect covariate present only late; auto-extend NA-fills it, which
  # would NaN the design matrix -> must error, not fit.
  keep <- (d$styr:d$endyr); keep <- keep[keep >= stats::median(keep)]
  d$env_data <- data.frame(Year = keep, temp = as.numeric(scale(seq_along(keep))))
  qflt <- which(d$fleet_control$Catchability %in% c(1L, 2L))[1]
  cf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = qflt)))
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(d, qFun = cf, estimateMode = 3, msmMode = 0,
      random_rec = FALSE, fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "missing \\(NA\\) values")
})

testthat::test_that("obs_sd_est = TRUE makes the QAR1 observation SD an estimated parameter", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  d <- Rceattle::BS2017SS
  nyr <- d$endyr - d$styr + 1
  set.seed(7)
  d$env_data <- data.frame(Year = d$styr:d$endyr, qcov = stats::rnorm(nyr))
  qflt <- which(d$fleet_control$Catchability %in% c(1L, 2L))[1]
  cf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = qflt,
                               observe = "qcov", obs_sd = 0.5, obs_sd_est = TRUE)))
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, qFun = cf, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  # now a FREE parameter (contrast with the fixed default above)
  testthat::expect_true("log_obs_sd_linkage" %in% names(fit$obj$par))
})

testthat::test_that("obs_sd_est is rejected without observe", {
  testthat::skip_if_not_installed("Rceattle")
  testthat::expect_error(
    Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, obs_sd_est = TRUE),
    "only used with")
})

testthat::test_that("a model with no linkages reports empty linkage effect-size vectors", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  testthat::expect_equal(length(fit$quantities$beta_linkage_obs), 0L)
})
