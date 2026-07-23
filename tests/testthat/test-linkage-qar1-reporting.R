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

testthat::test_that("QAR1 observation SD is an estimated parameter", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  d <- Rceattle::BS2017SS
  nyr <- d$endyr - d$styr + 1
  # A noisy covariate (not a smooth series the AR1 can perfectly track) so the
  # observation SD is informative/identified. (On a smooth series obs_sd is only
  # weakly identified and can collapse toward 0 -- the documented caveat.)
  set.seed(42)
  d$env_data <- data.frame(Year = d$styr:d$endyr, qcov = stats::rnorm(nyr))
  qflt <- which(d$fleet_control$Catchability %in% c(1L, 2L))[1]
  cf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = qflt,
                               observe = "qcov", obs_sd = 0.5)))
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, qFun = cf, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  # It is a free parameter (one per observed group) and estimates to a finite SD.
  testthat::expect_equal(length(fit$estimated_params$log_obs_sd_linkage), 1L)
  osd <- exp(fit$estimated_params$log_obs_sd_linkage)
  testthat::expect_true(is.finite(osd) && osd > 0)
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
