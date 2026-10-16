# The within-year growth routine chose its branch on `age + 1.0` -- the 0-based
# slot index plus one -- where the quantity it meant is the cohort's JAN-1 age,
# `age + minage(sp)`. The two coincide only at `minage == 1`. Fixed, and the
# plus-group pin moved above the anchor test, which the fix would otherwise
# preempt. `inst/dev/TRAPS.md` carries the measured numbers.
#
# Two things this fixture has to get right, both of which a first attempt did
# not:
#
#   * `fracyr = spawn_month / 12`, and at `fracyr == 0` the first branch
#     (`current_age <= age_L1`) takes everything. So `spawn_month > 0`.
#   * the anchor must sit ABOVE `minage`. At the default `max(0.5, minage)` the
#     diverging slot is the one whose Jan-1 age IS the anchor, where its Jan-1
#     length is exactly `l1` and the two arms give the same number -- the
#     branches differ and the values do not. `growth_age_L1 = 3.5` puts
#     the divergence on slot 2, whose Jan-1 length is on the curve -- true only
#     since 5.57.3; before it that Jan-1 value was itself a mid-ramp recursion.
#     These assertions survived either way because each is taken relative to
#     the Jan-1 value it reads back, or sits on the mixed branch.

testthat::skip_on_cran()

# length_hat[wt_slot, sex, age, yr]; slot 1 Jan-1, slot 2 the mid-year
# (spawning) evaluation, as test-growth-plus-group-within-year.R uses it.
.m2_fit <- function(minage, spawn_month = 6, age_L1 = NA, nages = 5,
                    fun = "vonBertalanffy") {
  dat <- make_test_data(nyrs = 6, nages = nages, minage = minage, seed = 42,
                        growth = fun)
  dat$spawn_month <- spawn_month
  compile_tmb_if_needed()
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    growthFun = Rceattle::build_growth(fun = fun, growth_age_L1 = age_L1),
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))))
}

# The advance arm, rearranged: L + (L - linf)(e - 1) == linf + (L - linf)e.
.m2_advance <- function(jan1, linf, kappa, fracyr) {
  linf + (jan1 - linf) * exp(-kappa * fracyr)
}


testthat::test_that("a cohort above the anchor is advanced, not re-read from it", {
  testthat::skip_if_not_installed("TMB")
  # minage 2, anchor 3.5, fracyr 0.5. Slot 2 holds the age-4 cohort: Jan-1 age
  # 4, mid-year 4.5, both above the anchor, so it must be advanced from its
  # Jan-1 length. The old test `2 + 1 = 3 < 3.5` sent it to the anchor arm
  # instead, discarding the cohort history the Jan-1 value carries.
  fit <- .m2_fit(minage = 2, age_L1 = 3.5)
  q <- fit$quantities
  testthat::expect_equal(unname(fit$obj$env$data$growth_age_L1[1]), 3.5)

  kappa <- q$growth_parameters[1, 1, 1, 1]
  linf  <- q$growth_parameters[1, 1, 1, 3]
  jan1  <- q$length_hat[1, 1, 3, 2]   # slot 2, year 2 (so a recursion exists)
  mid   <- q$length_hat[2, 1, 3, 2]
  testthat::expect_true(all(is.finite(c(kappa, linf, jan1))) && jan1 > 0)
  testthat::expect_equal(mid, .m2_advance(jan1, linf, kappa, 0.5),
                         tolerance = 1e-8)
})


testthat::test_that("minage = 1 is bit-identical, so nothing shipping moves", {
  testthat::skip_if_not_installed("TMB")
  # `age + 1.0` and `age + minage(sp)` are the same expression at minage = 1,
  # so the fix is a no-op on the configuration every bundled dataset uses. To
  # mean anything this has to REACH the anchor arm, which needs
  # fracyr > age_L1 - floor(age_L1): anchor 3.5 with fracyr 0.75, where slot 2
  # (Jan-1 age 3, mid-year 3.75) crosses it. At fracyr 0.5 every slot below the
  # anchor is still on the ramp and the arm never fires.
  fit <- .m2_fit(minage = 1, spawn_month = 9, age_L1 = 3.5)
  q <- fit$quantities
  kappa <- q$growth_parameters[1, 1, 1, 1]
  l1    <- q$growth_parameters[1, 1, 1, 2]
  linf  <- q$growth_parameters[1, 1, 1, 3]
  anchor <- linf + (l1 - linf) * exp(-kappa * ((2 + 1 + 0.75) - 3.5))
  testthat::expect_equal(q$length_hat[2, 1, 3, 2], anchor, tolerance = 1e-8)
  # And the slot above it is advanced from its Jan-1 length, not re-read.
  testthat::expect_equal(q$length_hat[2, 1, 4, 2],
                         .m2_advance(q$length_hat[1, 1, 4, 2], linf, kappa, 0.75),
                         tolerance = 1e-8)
})


testthat::test_that("the plus group stays pinned when the anchor is above it", {
  testthat::skip_if_not_installed("TMB")
  # The regression the reorder prevents. With the whole age range under the
  # anchor, the anchor test also matches the plus group; placed first it would
  # discard the Jan-1 mean length estimate_growth() wrote for that slot.
  # minage 0, nages 4, anchor 3.5: the plus slot's Jan-1 age is 3 < 3.5.
  fit <- .m2_fit(minage = 0, spawn_month = 9, age_L1 = 3.5, nages = 4)
  q <- fit$quantities
  n <- dim(q$length_hat)[3]
  testthat::expect_equal(q$length_hat[2, 1, n, 2], q$length_hat[1, 1, n, 2],
                         tolerance = 0)
})


testthat::test_that("minage = 0 past mid-year gets the anchor arm it is owed", {
  testthat::skip_if_not_installed("TMB")
  # The configuration a routine SS3 update could walk into: the converter's
  # default is minage = 0, the anchor then defaults to 0.5, and a survey month
  # after June puts fracyr above it. Slot 0's Jan-1 age is 0, below the anchor,
  # and its mid-year age 0.75 is above -- so it crossed during the year and is
  # read off the curve from the anchor. `age + 1.0 < 0.5` was false, so the old
  # code advanced it from a ramp length instead, and the length fell as the
  # season crossed July 1.
  fit <- .m2_fit(minage = 0, spawn_month = 9, nages = 5)   # anchor -> 0.5
  q <- fit$quantities
  testthat::expect_equal(unname(fit$obj$env$data$growth_age_L1[1]), 0.5)
  kappa <- q$growth_parameters[1, 1, 1, 1]
  l1    <- q$growth_parameters[1, 1, 1, 2]
  linf  <- q$growth_parameters[1, 1, 1, 3]
  # The anchor arm: linf + (l1 - linf) * exp(-kappa * (current_age - age_L1)).
  expect <- linf + (l1 - linf) * exp(-kappa * ((0 + 0 + 0.75) - 0.5))
  testthat::expect_equal(q$length_hat[2, 1, 1, 2], expect, tolerance = 1e-8)
})
