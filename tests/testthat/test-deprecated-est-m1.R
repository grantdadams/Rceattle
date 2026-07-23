# =============================================================================
# `est_M1` was renamed to `M1_model`. A data list carrying only the old name
# must still be honoured: switch_check() aliases it (with a deprecation message)
# BEFORE the M1_model default fills in 0 (fixed M). Without the alias the M1
# estimation setting is silently dropped -- M1 is fixed and the fit is wrong.
# =============================================================================

testthat::test_that("est_M1 is aliased to M1_model (not silently dropped to fixed M)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$M1_model <- NULL
  d$est_M1   <- 1L   # estimate sex- and age-invariant M1

  testthat::expect_message(
    out <- suppressWarnings(Rceattle::switch_check(d)),
    "'est_M1' is deprecated")
  testthat::expect_equal(out$M1_model, 1L)          # carried over, NOT defaulted to 0
})

testthat::test_that("M1_model takes precedence when both names are present", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$M1_model <- 2L
  d$est_M1   <- 1L
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$M1_model, 2L)
})

testthat::test_that("an all-NA est_M1 is treated as unset (falls through to the default)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$M1_model <- NULL
  d$est_M1   <- rep(NA_integer_, d$nspp)   # the GOA2018SS-style placeholder
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$M1_model, rep(0L, d$nspp))   # default, NOT NA
})

testthat::test_that("est_M1 is folded into M1_model in the fit path (warns, not silent)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  # est_M1 differs from the default build_M1() (M1_model = 0): fit_mod aliases it
  # early, so the existing "M1_model in data is different" warning fires instead
  # of the setting being silently dropped.
  d <- make_test_data()
  d$M1_model <- NULL
  d$est_M1   <- 1L
  testthat::expect_warning(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))),
    "M1_model in data is different")
})

testthat::test_that("a data list with only M1_model is unaffected (no est_M1 message)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$M1_model <- 0L
  d$est_M1   <- NULL
  msgs <- character(0)
  out <- withCallingHandlers(
    suppressWarnings(Rceattle::switch_check(d)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  testthat::expect_false(any(grepl("est_M1", msgs)))
  testthat::expect_equal(out$M1_model, 0L)
})
