# =============================================================================
# Back-compat: the fleet_control columns `Time_varying_q_sd_prior` /
# `Time_varying_sel_sd_prior` were renamed to `Time_varying_q_sd` /
# `Time_varying_sel_sd` (they are the input VALUE of the time-varying deviate
# SD, not a prior on it). switch_check() accepts the deprecated spellings and
# upgrades them in place, so existing data lists / .rda / scripts keep fitting
# identically. These tests pin that alias -- both the rename mechanics and
# bit-identical fit equivalence -- so the released-package contract holds.
# =============================================================================

# Rename the modern SD columns back to their deprecated spellings.
.deprecate_sd_names <- function(data_list) {
  fc <- data_list$fleet_control
  names(fc)[names(fc) == "Time_varying_sel_sd"] <- "Time_varying_sel_sd_prior"
  names(fc)[names(fc) == "Time_varying_q_sd"]   <- "Time_varying_q_sd_prior"
  data_list$fleet_control <- fc
  data_list
}

testthat::test_that("switch_check() renames Q_prior/Index_sd_prior/Catch_sd_prior", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  fc <- d$fleet_control
  names(fc)[names(fc) == "Catchability_init"] <- "Q_prior"
  names(fc)[names(fc) == "Index_sd"] <- "Index_sd_prior"
  names(fc)[names(fc) == "Catch_sd"] <- "Catch_sd_prior"
  d$fleet_control <- fc

  out <- suppressMessages(suppressWarnings(Rceattle::switch_check(d)))
  fc_new <- out$fleet_control
  # `Q_prior` now upgrades to the renamed canonical `Catchability_init`.
  testthat::expect_true(all(c("Catchability_init", "Index_sd", "Catch_sd") %in% names(fc_new)))
  testthat::expect_false(any(c("Q_prior", "Index_sd_prior", "Catch_sd_prior") %in% names(fc_new)))
  testthat::expect_equal(fc_new$Catchability_init, fc$Q_prior)
  testthat::expect_equal(fc_new$Index_sd, fc$Index_sd_prior)
  testthat::expect_equal(fc_new$Catch_sd, fc$Catch_sd_prior)
})

testthat::test_that("switch_check() renames deprecated Time_varying_*_sd_prior columns", {
  testthat::skip_if_not_installed("Rceattle")
  dat_old <- .deprecate_sd_names(make_test_data())
  fc_old  <- dat_old$fleet_control

  testthat::expect_message(
    out <- suppressWarnings(Rceattle::switch_check(dat_old)),
    "Time_varying_q_sd_prior' is deprecated")

  fc_new <- out$fleet_control
  # New names present, old names gone.
  testthat::expect_true(all(c("Time_varying_q_sd", "Time_varying_sel_sd") %in% names(fc_new)))
  testthat::expect_false(any(c("Time_varying_q_sd_prior", "Time_varying_sel_sd_prior") %in% names(fc_new)))
  # Values carried over unchanged.
  testthat::expect_equal(fc_new$Time_varying_q_sd,   fc_old$Time_varying_q_sd_prior)
  testthat::expect_equal(fc_new$Time_varying_sel_sd, fc_old$Time_varying_sel_sd_prior)
})

testthat::test_that("a deprecated-SD-column data list fits identically to the renamed one", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  ctl     <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  dat_new <- make_test_data(seed = 7)
  dat_old <- .deprecate_sd_names(dat_new)

  # estimateMode = 3 returns the real objective without optimizing (fast, and
  # deterministic for a bit-identity check).
  fit_new <- suppressMessages(Rceattle::fit_mod(
    data_list = dat_new, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = ctl))
  fit_old <- suppressMessages(Rceattle::fit_mod(
    data_list = dat_old, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = ctl))

  obj_new <- sum(fit_new$quantities$jnll_comp)
  obj_old <- sum(fit_old$quantities$jnll_comp)
  # Finite (the old failure mode was log(numeric(0)) -> NaN) and bit-identical.
  testthat::expect_true(is.finite(obj_new))
  testthat::expect_equal(obj_old, obj_new, tolerance = 0)
})

# HCR is a fit_mod() argument, not a data field, so a data_list read from a
# workbook carries none. NULL %in% ... is logical(0), which propagated into the
# HCR checks and left them holding logical(0)/NA -- validating a workbook before
# fitting then died with "missing value where TRUE/FALSE needed" or "argument is
# of length zero" instead of reporting anything.
testthat::test_that("data_check tolerates a data_list with no HCR", {
  testthat::skip_if_not_installed("TMB")

  dat <- make_test_data(nyrs = 10, nages = 5, seed = 123)
  dat$HCR <- NULL
  testthat::expect_null(dat$HCR)
  testthat::expect_no_error(suppressWarnings(Rceattle:::data_check(dat)))

  # ... and the check it guards still fires when an HCR IS supplied.
  dat$HCR <- 5
  dat$fleet_control$Proj_F_proportion <- 0
  err <- tryCatch({ suppressWarnings(Rceattle:::data_check(dat)); NULL },
                  error = function(e) conditionMessage(e))
  testthat::expect_true(!is.null(err))
  testthat::expect_match(err, "Proj_F_proportion")
})
