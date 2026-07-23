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
