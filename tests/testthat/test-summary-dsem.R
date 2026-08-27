# Constructive tests for summary.Rceattle()'s DSEM branch.
#
# These build a synthetic Rceattle object carrying a dsem$sem_full shaped exactly
# as dsem::dsem() returns it (a data.frame with columns path, lag, name, start,
# parameter, first, second, direction) and drive summary() directly. Working from
# a fixture rather than a real fit keeps the tests fast and lets them cover the
# map and sdreport combinations a single fit cannot reach at once. The contract
# the downstream analysis scripts depend on is summary(fit)$coefficients and
# $recruitment_sd$R_sd.

# One estimated path, one variance, one path FIXED at its start value.
fake_sem_full <- function() {
  data.frame(
    path      = c("BT -> recdevs1", "recdevs1 <-> recdevs1", "X -> recdevs1"),
    lag       = c(1, 0, 0),
    name      = c("BT_to_R", "sigmaR1", "X_fixed"),
    start     = c(0, 1, 0.25),
    parameter = c(1L, 2L, 0L),          # 0 => fixed at `start`
    first     = c("BT", "recdevs1", "X"),
    second    = c("recdevs1", "recdevs1", "recdevs1"),
    direction = c(1, 2, 1),
    stringsAsFactors = FALSE
  )
}

fake_fit <- function(beta_z = c(0.5, 0.8), sem_full = fake_sem_full(),
                     bmap = factor(c(1, 2)), sdrep = NULL, nspp = 1L) {
  structure(
    list(
      data_list = list(nspp = nspp, spnames = paste0("spp", seq_len(nspp)),
                       random_rec = TRUE),
      estimated_params = list(dsem_beta_z = beta_z),
      quantities = list(R_sd = rep(0.8, nspp)),
      map   = list(mapList = list(dsem_beta_z = bmap)),
      sdrep = sdrep,
      dsem  = list(sem_full = sem_full,
                   tmb_inputs = list(data = list(rec_sd_idx = rep(2L, nspp))))
    ),
    class = "Rceattle"
  )
}

testthat::test_that("summary() returns the column contract the analysis scripts read", {
  out <- suppressMessages(summary(fake_fit()))

  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("coefficients", "recruitment_sd") %in% names(out)))

  # The exact columns Results/Initial_DSEM.csv carries, in order.
  testthat::expect_equal(
    colnames(out$coefficients),
    c("path", "lag", "name", "start", "parameter", "first", "second",
      "direction", "Estimate", "Std_Error", "z_value", "p_value")
  )
  testthat::expect_true(is.character(out$coefficients$direction) ||
                          is.numeric(out$coefficients$direction))
  testthat::expect_true(is.numeric(out$coefficients$Estimate))

  testthat::expect_true(all(c("Species", "R_sd", "Estimated", "Std_Error") %in%
                              names(out$recruitment_sd)))
  testthat::expect_equal(nrow(out$recruitment_sd), 1L)
})

testthat::test_that("estimated paths take beta_z; only parameter == 0 falls back to start", {
  out <- suppressMessages(summary(fake_fit(beta_z = c(0.5, 0.8),
                                           bmap = factor(c(1, 2, 3)))))
  co <- out$coefficients

  testthat::expect_equal(co$Estimate[co$name == "BT_to_R"], 0.5)
  testthat::expect_equal(co$Estimate[co$name == "sigmaR1"], 0.8)
  # parameter == 0 -> the start value, not an estimate
  testthat::expect_equal(co$Estimate[co$name == "X_fixed"], 0.25)
})

testthat::test_that("a NaN estimate is reported as NaN, not silently rewritten to its start value", {
  # Regression: keying the fixed-path fallback on is.na(Estimate) rewrote a
  # diverged path (NaN MLE) into its start value, reporting a failed estimate as
  # a fixed input.
  out <- suppressMessages(summary(fake_fit(beta_z = c(NaN, 0.8),
                                           bmap = factor(c(1, 2, 3)))))
  co <- out$coefficients
  testthat::expect_true(is.nan(co$Estimate[co$name == "BT_to_R"]))
  testthat::expect_false(isTRUE(co$Estimate[co$name == "BT_to_R"] == 0))
})

testthat::test_that("a missing or short beta_z errors instead of returning start values", {
  # Regression: a NULL beta_z produced a full table of start values with no
  # error -- a plausible-looking sheet of zeros written straight to Results/.
  testthat::expect_error(
    suppressMessages(summary(fake_fit(beta_z = NULL))),
    "no estimated_params\\$dsem_beta_z"
  )
  testthat::expect_error(
    suppressMessages(summary(fake_fit(beta_z = 0.5))),   # SEM references index 2
    "references parameter index"
  )
})

testthat::test_that("Std_Error/z_value/p_value are always present, NA without an sdreport", {
  # Regression: making them conditional changed the column count between fits,
  # breaking the rbind() of several models' tables in the analysis scripts.
  out <- suppressMessages(summary(fake_fit(bmap = factor(c(1, 2, 3)))))
  testthat::expect_true(all(c("Std_Error", "z_value", "p_value") %in%
                              colnames(out$coefficients)))
  testthat::expect_true(all(is.na(out$coefficients$Std_Error)))

  a <- suppressMessages(summary(fake_fit(bmap = factor(c(1, 2, 3)))))$coefficients
  b <- suppressMessages(summary(fake_fit(beta_z = c(0.1, 0.2),
                                         bmap = factor(c(1, 2, 3)))))$coefficients
  testthat::expect_equal(ncol(a), ncol(b))
  testthat::expect_silent(rbind(a, b))
})

testthat::test_that("recruitment SD reports the right species and masks fixed SDs", {
  f <- fake_fit(beta_z = c(0.5, 0.8), bmap = factor(c(1, NA, 3)), nspp = 1L)
  rs <- suppressMessages(summary(f))$recruitment_sd
  # rec_sd_idx is 1-based; entry 2 is mapped off => fixed, not estimated.
  testthat::expect_false(rs$Estimated[1])
  testthat::expect_true(is.na(rs$Std_Error[1]))
})

testthat::test_that("a fit with no DSEM carries no DSEM tables", {
  plain <- structure(list(data_list = list(nspp = 1L, spnames = "spp1")),
                     class = "Rceattle")
  out <- suppressMessages(summary(plain))
  testthat::expect_s3_class(out, "summary.Rceattle")
  testthat::expect_identical(out$spec, plain)
  testthat::expect_null(out$dsem_coefficients)
  testthat::expect_null(out$recruitment_sd)
  # On a non-DSEM fit $coefficients is the fixed effects, and this fixture has none.
  testthat::expect_null(out$coefficients)
})

# $coefficients means the SEM path table on a DSEM fit and the fixed effects on
# every other fit. That is deliberate -- the DSEM analysis scripts read the path
# table from that name -- and it is the one place in summary() where a slot's
# meaning depends on the model, so pin both halves.
testthat::test_that("$coefficients is the path table on a DSEM fit, fixed effects elsewhere", {
  out <- suppressMessages(summary(fake_fit()))
  testthat::expect_identical(out$coefficients, out$dsem_coefficients)
  testthat::expect_true("path" %in% colnames(out$coefficients))
  # The fixed effects are still reachable, under a name that means one thing.
  testthat::expect_true("fixed_coefficients" %in% names(out))
})

testthat::test_that("summary() rejects a non-Rceattle object", {
  testthat::expect_error(summary.Rceattle(structure(list(), class = "nope")),
                         "not an Rceattle model")
})

# The map filter is bridged by `parameter`, not by row position. Sem rows
# routinely outnumber beta_z entries -- the fixture has three paths and two
# parameters, because one path is fixed at its start value -- so a filter keyed
# on length(bmap) == nrow(coefs) never fires and warns on every ordinary fit.

testthat::test_that("an ordinary fit reports every path and does not warn", {
  # 3 sem rows, 2 beta_z entries, both mapped on: nothing is dropped.
  testthat::expect_no_warning(out <- suppressMessages(summary(fake_fit())))
  testthat::expect_equal(nrow(out$coefficients), 3L)
  testthat::expect_true("X_fixed" %in% out$coefficients$name)
})

testthat::test_that("a mapped-off path is dropped, and a fixed path is kept", {
  # beta_z entry 1 (BT_to_R) mapped off; entry 2 (sigmaR1) still estimated.
  # X_fixed owns no beta_z entry (parameter == 0) so it cannot be mapped off.
  f <- fake_fit(bmap = factor(c(NA, 2)))
  testthat::expect_no_warning(out <- suppressMessages(summary(f)))
  testthat::expect_setequal(out$coefficients$name, c("sigmaR1", "X_fixed"))
  testthat::expect_false("BT_to_R" %in% out$coefficients$name)
})

# Everything above drives summary() from a hand-built fixture, which cannot catch
# the fixture naming a slot the fit does not actually have -- the accessor then
# only ever validates against itself. summary() read `estimated_params$beta_z`
# for exactly that reason while a fit stores `dsem_beta_z` (every DSEM parameter
# is prefixed), so it errored on every real DSEM fit and every fixture test
# passed. This runs a real one.

testthat::test_that("summary() reads the slot names a real DSEM fit carries", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # The contract summary() depends on, asserted against the fit rather than a
  # fixture: name these wrong and summary() cannot report anything.
  testthat::expect_false(is.null(fit$estimated_params$dsem_beta_z))
  testthat::expect_false(is.null(fit$map$mapList$dsem_beta_z))
  testthat::expect_false(is.null(fit$dsem$sem_full))

  out <- suppressMessages(summary(fit))
  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("coefficients", "recruitment_sd") %in% names(out)))
  testthat::expect_true(is.numeric(out$recruitment_sd$R_sd))
})

testthat::test_that("a map too short for the SEM warns and filters nothing", {
  # The genuine inconsistency: the SEM references parameter 2 but the map has
  # only one entry. Report everything rather than drop rows on a bad index.
  f <- fake_fit(bmap = factor(1))
  testthat::expect_warning(out <- suppressMessages(summary(f)),
                           "references parameter index")
  testthat::expect_equal(nrow(out$coefficients), 3L)
})
