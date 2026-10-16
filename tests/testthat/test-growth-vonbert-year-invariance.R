# Regression baseline for minage > 0 (the historical default).
#
# Replacing `age_L1` with `age_L1_safe` in last_linear is algebraically a
# no-op when minage > 0, but a future refactor could break it. This test
# asserts year-invariance of length_hat across all ages under constant
# growth params -- catches any regression in either the cohort recursion
# or the year-0 boundary branch.

testthat::skip_on_cran()
testthat::test_that("vonBertalanffy + minage = 1: length-at-age is year-invariant for every age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  dat <- make_test_data(nyrs = 6, nages = nages, minage = 1, seed = 42,
                        growth = "vonBertalanffy")
  testthat::expect_equal(dat$minage, 1)
  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  nyrs <- dat$endyr - dat$styr + 1L
  # length_hat: [wt_slot, sex, age, yr]. Pop slot = 1 for sp = 1.
  lhat <- res$quantities$length_hat[1, 1, , seq_len(nyrs), drop = FALSE]
  testthat::expect_true(all(is.finite(lhat)))

  # For every age, the across-year range must collapse under constant params.
  for (a in seq_len(nages)) {
    rng <- diff(range(lhat[1, 1, a, ]))
    testthat::expect_lt(rng, 1e-8,
                        label = paste0("year-range at age index ", a))
  }

  # Age 0 (index 1) sits at current_age = 1 = age_L1, so the linear ramp
  # endpoint equals l1. Pull l1 from the second age's closed-form: at
  # age = 1 (current_age = 2), L = linf + (l1 - linf) * exp(-K).
  # We only assert that the youngest age is strictly below the next age
  # (monotone growth) -- a coarse but param-agnostic sanity check.
  testthat::expect_lt(lhat[1, 1, 1, 1], lhat[1, 1, 2, 1])
})


# A NON-INTEGER anchor above max(minage, 1) is what the block above cannot see,
# because every default anchor is an integer or 0.5. The cohort recursion reads
# the slot one year younger, so the first slot above the anchor has to be read
# off the anchor instead -- its predecessor is still on the linear ramp. Fixing
# that moves no anchor anyone uses (5.57.3; `inst/dev/TRAPS.md` has the measured
# error and the reroute rule).
#
# `make_test_data` gives Lmin == l1, so its ramp is FLAT and these blocks see
# only the "advanced a whole year instead of a fraction" half of the defect, not
# the "recursed from a mid-ramp length" half. Real stocks have l1 > Lmin.

.anchor_fit <- function(fun, minage, age_L1, nages = 6, nyrs = 6) {
  dat <- make_test_data(nyrs = nyrs, nages = nages, minage = minage, seed = 42,
                        growth = fun)
  compile_tmb_if_needed()
  res <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    growthFun = Rceattle::build_growth(fun = fun, growth_age_L1 = age_L1),
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))))
  testthat::expect_equal(unname(res$obj$env$data$growth_age_L1[1]), age_L1)
  res
}

# The closed form off the anchor, both growth forms:
# L = linf + (l1 - linf) * exp(-K(age - a_L1)), and its Richards power mean.
.anchor_closed_form <- function(g, fun, age, age_L1) {
  kappa <- g[1, 1, 1, 1]
  l1 <- g[1, 1, 1, 2]
  linf <- g[1, 1, 1, 3]
  d <- exp(-kappa * (age - age_L1))
  if (fun == "vonBertalanffy") return(linf + (l1 - linf) * d)
  m <- g[1, 1, 1, 4]
  (linf^m + (l1^m - linf^m) * d)^(1 / m)
}

testthat::test_that("a non-integer growth anchor stays year-invariant", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # minage 0 is in the grid because it is the SS3 converter's own default, and
  # it moves by 3.9 cm at this anchor.
  for (fun in c("vonBertalanffy", "Richards")) {
    for (minage in c(0L, 1L)) {
      # Anchor 2.5, so the cohort at Jan-1 age 3 is the first above it and the
      # one whose predecessor is still on the ramp. It sits in 0-based slot
      # `3 - minage`, hence array index `4 - minage` (rule 11).
      res <- .anchor_fit(fun, minage = minage, age_L1 = 2.5)
      nyrs <- res$data_list$endyr - res$data_list$styr + 1L
      lhat <- res$quantities$length_hat[1, 1, , seq_len(nyrs)]
      lab <- paste(fun, "minage", minage)
      testthat::expect_true(all(is.finite(lhat)), label = paste(lab, "finite"))
      for (a in seq_len(nrow(lhat))) {
        testthat::expect_lt(diff(range(lhat[a, ])), 1e-8,
                            label = paste(lab, "year-range at age index", a))
      }
      # And the value, so a uniformly-wrong trajectory cannot pass.
      testthat::expect_equal(
        lhat[4L - minage, nyrs],
        .anchor_closed_form(res$quantities$growth_parameters, fun, 3, 2.5),
        tolerance = 1e-8, label = paste(lab, "slot above the anchor"))
    }
  }
})


# Two structural guards on the gate's second clause and its floored anchor.
# NEITHER reproduces a defect -- both configurations are bit-identical on the
# 5.57.3 parent, where `max(minage, 1)` happened to name the right slot. They
# exist because the two halves of the new test are otherwise unpinned: drop
# `|| age == 0` and the first block still passes while slot 0 recurses from
# index -1, and swap `age_L1_safe` for `age_L1` and nothing fails at all.
testthat::test_that("slot 0 stays off the recursion, anchor floored", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # (a) minage 2 with the SS3-style 0.5 anchor: slot 0's Jan-1 age is already
  # past the anchor, and it has no one-year-younger predecessor to recurse from.
  # Without the `age == 0` clause this reads index -1, which is an out-of-bounds
  # load this build (safebounds = FALSE) does not report.
  for (fun in c("vonBertalanffy", "Richards")) {
    res <- .anchor_fit(fun, minage = 2, age_L1 = 0.5)
    nyrs <- res$data_list$endyr - res$data_list$styr + 1L
    lhat <- res$quantities$length_hat[1, 1, , seq_len(nyrs)]
    testthat::expect_true(all(is.finite(lhat)), label = paste(fun, "finite"))
    testthat::expect_lt(diff(range(lhat[1, ])), 1e-8)
    testthat::expect_equal(
      lhat[1, nyrs],
      .anchor_closed_form(res$quantities$growth_parameters, fun, 2, 0.5),
      tolerance = 1e-8, label = paste(fun, "slot 0 above the anchor"))
  }

  # (b) An anchor of 0 is the one configuration that tells `age_L1_safe` apart
  # from `age_L1` in the gate, and it is SS3's own default (`AFIX = 0.`). Slot 1
  # must be read off the anchor: with the unfloored test it recurses from slot 0
  # instead, 6.23 cm low at realistic parameters.
  # Slot 0 itself is NOT asserted -- at this anchor it returns the lowest
  # population bin edge where SS3 returns `l1`, which is an open defect in
  # the developer backlog, not something to pin here.
  res <- .anchor_fit("vonBertalanffy", minage = 0, age_L1 = 0)
  nyrs <- res$data_list$endyr - res$data_list$styr + 1L
  lhat <- res$quantities$length_hat[1, 1, , seq_len(nyrs)]
  testthat::expect_true(all(is.finite(lhat)))
  testthat::expect_equal(
    lhat[2, nyrs],
    .anchor_closed_form(res$quantities$growth_parameters,
                        "vonBertalanffy", 1, 0),
    tolerance = 1e-8)
})
