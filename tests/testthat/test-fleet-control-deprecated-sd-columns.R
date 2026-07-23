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
