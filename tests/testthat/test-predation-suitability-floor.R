# Under empirical suitability a prey age whose average numbers approach zero
# while its diet proportion stays positive sends stom_div_bio to infinity, so
# that one cell takes the whole normalized share. "Other food" is
# 1 - sum(prey suitability), so it collapses, available food collapses with it,
# and M2 = ration * N / avail_food diverges -- finite, no NaN, and before this
# the only warning raised was an unrelated one about missing diet rows.
#
# Measured, and attributed to the model each number came from:
#   two-species fixture, window narrowed : other_food 1e+05, avail_food 3.6e-06,
#                                          suit_other 3.6e-11, max M2 1.3e+11
#   BS2017MS, its own default cold start : suit_other 4.5e-10, max M2 2.0e+04
#   BS2017MS fitted                      : suit_other 2.9e-03, max M2 1.785
#   goa_ms  fitted                       : suit_other 5.3e-03, max M2 0.449
#
# predation.hpp holds other food at a fixed FRACTION (1e-6) of the prey term, so
# suit_other >= 1e-6/(1 + 1e-6). A ratio, so it carries no units. It binds at
# both multispecies reference models' COLD STARTS and at neither's optimum,
# which is why golden -- evaluated at the stored MLE -- is bit-identical and
# could never have caught this.

testthat::skip_on_cran()

valid_scan <- function(fit) {
  dl <- fit$data_list
  q <- fit$quantities
  nyh <- dl$endyr - dl$styr + 1L
  out <- list(suit_other = numeric(0), m2 = numeric(0), avail = numeric(0))
  for (sp in seq_len(dl$nspp)) {
    for (sx in seq_len(dl$nsex[sp])) {
      for (ag in seq_len(dl$nages[sp])) {
        out$suit_other <- c(out$suit_other,
                            q$suit_other[sp, sx, ag, seq_len(nyh)])
        out$m2 <- c(out$m2, q$M2_at_age[sp, sx, ag, seq_len(nyh)])
        out$avail <- c(out$avail, q$avail_food[sp, sx, ag, seq_len(nyh)])
      }
    }
  }
  out
}

fit_ms <- function(...) {
  data("BS2017MS", package = "Rceattle", envir = environment())
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017MS, inits = NULL, random_rec = FALSE, msmMode = 1,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0), ...)))
}


testthat::test_that("the floor bounds a collapsing other-food share", {
  testthat::skip_if_not_installed("TMB")
  data("BS2017MS", package = "Rceattle", envir = environment())
  # Narrowing the suitability window to the last five years is a cheap way into
  # the pathology: the window no longer matches the years the diet describes.
  ey <- BS2017MS$endyr
  v <- valid_scan(fit_ms(estimateMode = 3, suit_styr = ey - 4, suit_endyr = ey))

  # The floor, not merely "small": 1e-6/(1 + 1e-6).
  testthat::expect_gte(min(v$suit_other), 1e-6 / (1 + 1e-6) * (1 - 1e-9))
  testthat::expect_lt(max(v$m2), 1e4)
  testthat::expect_true(all(is.finite(v$m2)))
})


testthat::test_that("a fitted model does not reach the floor", {
  # The floor must be inert on anything anyone fits, or it is changing advice.
  testthat::skip_if_not_installed("TMB")
  v <- valid_scan(fit_ms(estimateMode = 1))
  testthat::expect_gt(min(v$suit_other), 1e-3)
  testthat::expect_lt(max(v$m2), 5)
  # Pinned: a change that moved these would be moving a real multispecies fit.
  testthat::expect_equal(min(v$suit_other), 0.002901, tolerance = 1e-3)
  testthat::expect_equal(max(v$m2), 1.785, tolerance = 1e-3)
})


testthat::test_that("hitting the floor warns, and a healthy fit does not", {
  # A bounded number is still wrong, so the message is the deliverable.
  testthat::skip_if_not_installed("TMB")
  data("BS2017MS", package = "Rceattle", envir = environment())
  ey <- BS2017MS$endyr

  testthat::expect_warning(
    suppressMessages(Rceattle::fit_mod(
      data_list = BS2017MS, inits = NULL, random_rec = FALSE, msmMode = 1,
      estimateMode = 3, suit_styr = ey - 4, suit_endyr = ey,
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "other food' suitability floor bound")

  # Negative control: the same model over its own window must be silent, or the
  # warning says nothing.
  testthat::expect_no_warning(
    suppressMessages(Rceattle::fit_mod(
      data_list = BS2017MS, inits = NULL, random_rec = FALSE, msmMode = 1,
      estimateMode = 1,
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))))
})


testthat::test_that("the diagnostic reads valid index ranges only", {
  # BS2017MS has nages 12, 12, 21, so 18 of its 63 age slots are padding at 0
  # and a whole-array scan would read every one of them as floored.
  testthat::skip_if_not_installed("TMB")
  fit <- fit_ms(estimateMode = 1)
  dl <- fit$data_list
  testthat::expect_true(any(dl$nages < max(dl$nages)))
  testthat::expect_equal(min(fit$quantities$suit_other), 0)
  testthat::expect_gt(min(valid_scan(fit)$suit_other), 1e-3)
  # So the helper must stay silent here despite the zeros in the padding.
  testthat::expect_null(
    Rceattle:::.warn_suitability_floor(dl, fit$quantities))
})


testthat::test_that("a single-species fit is silent", {
  # suit_other is setZero() and written only under predation, so on msmMode = 0
  # every cell is 0. Without a guard that reads as floored: measured before the
  # guard, a 6-year single-species fit warned about 24 cells and BS2017SS about
  # 1755, on models with no predation at all.
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 6, nages = 4, seed = 1)
  testthat::expect_no_warning(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))))
  testthat::expect_null(Rceattle:::.warn_suitability_floor(
    list(msmMode = 0L, nspp = 1L, nsex = 1L, nages = 2L, styr = 1L, endyr = 2L),
    list(suit_other = array(0, dim = c(1, 1, 2, 2)))))
})


testthat::test_that("the warning counts every floored cell, not just one", {
  # Measured before the fix: five floored cells in one species plus one lower
  # cell in another were reported as "1 cell(s)" and attributed to the species
  # holding the lower one, because the accumulator only ran in the else branch.
  dl <- list(msmMode = 1L, nspp = 2L, nsex = c(1L, 1L), nages = c(3L, 3L),
             styr = 1980L, endyr = 1989L, spnames = c("Pollock", "Cod"))
  so <- array(0.5, dim = c(2, 1, 3, 10))
  so[1, 1, 1, 1:5] <- 1e-7
  so[2, 1, 2, 3] <- 1e-9
  h <- suppressWarnings(
    Rceattle:::.warn_suitability_floor(dl, list(suit_other = so)))
  testthat::expect_equal(nrow(h), 6L)
  testthat::expect_warning(
    Rceattle:::.warn_suitability_floor(dl, list(suit_other = so)),
    "bound 6 predator-age-year")
  # The lowest cell names its own species, not the one with the most cells.
  testthat::expect_warning(
    Rceattle:::.warn_suitability_floor(dl, list(suit_other = so)), "Cod age 2")
})
