# Constructive tests for summary.Rceattle()'s DSEM branch.
#
# DSEM is not yet wired into fit_mod() on this branch, so the branch is
# unreachable from a real fit. These build a synthetic Rceattle object carrying a
# dsem$sem_full shaped exactly as dsem::dsem() returns it (a data.frame with
# columns path, lag, name, start, parameter, first, second, direction) and drive
# summary() directly. That keeps the contract the downstream analysis scripts
# depend on -- summary(fit)$coefficients and $recruitment_sd$R_sd -- under test
# before Tier 2 can regress it.

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
      estimated_params = list(beta_z = beta_z),
      quantities = list(R_sd = rep(0.8, nspp)),
      map   = list(mapList = list(beta_z = bmap)),
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
    "no estimated_params\\$beta_z"
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

testthat::test_that("mapped-off paths are dropped, and a length mismatch warns instead of recycling", {
  # bmap length == nrow(sem_full): the filter applies.
  out <- suppressMessages(summary(fake_fit(bmap = factor(c(1, NA, 3)))))
  testthat::expect_equal(nrow(out$coefficients), 2L)
  testthat::expect_false("sigmaR1" %in% out$coefficients$name)

  # Length mismatch: report unfiltered with a warning rather than let the
  # shorter logical recycle and silently select the wrong rows.
  testthat::expect_warning(
    out2 <- summary(fake_fit(bmap = factor(c(1, NA)))),
    "differ"
  )
  testthat::expect_equal(nrow(out2$coefficients), 3L)
})

testthat::test_that("recruitment SD reports the right species and masks fixed SDs", {
  f <- fake_fit(beta_z = c(0.5, 0.8), bmap = factor(c(1, NA, 3)), nspp = 1L)
  rs <- suppressMessages(summary(f))$recruitment_sd
  # rec_sd_idx is 1-based; entry 2 is mapped off => fixed, not estimated.
  testthat::expect_false(rs$Estimated[1])
  testthat::expect_true(is.na(rs$Std_Error[1]))
})

testthat::test_that("a fit with no DSEM is untouched by the restored method", {
  plain <- structure(list(data_list = list(nspp = 1L, spnames = "spp1")),
                     class = "Rceattle")
  out <- suppressMessages(summary(plain))
  testthat::expect_identical(out, plain)
  testthat::expect_null(out$coefficients)
})

testthat::test_that("summary() rejects a non-Rceattle object", {
  testthat::expect_error(summary.Rceattle(structure(list(), class = "nope")),
                         "not an Rceattle model")
})
