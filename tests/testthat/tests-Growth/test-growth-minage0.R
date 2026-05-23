# Regression test for a TMB::MakeADFun segfault that occurred when
# combining `minage = 0` with a parametric growth function (vonBertalanffy
# or Richards). The cause was a negative-index access in growth.hpp:
#
#   length_hat(wtind, sex, age - 1, yr - 1)   // age = 0  =>  -1  =>  segfault
#
# The bug fired only at minage = 0 because the boundary branch
# `if (current_age == age_L1_ceil)` evaluated to `1 == 0` and was never
# taken, causing the unguarded recursive branch to run at age = 0.
#
# The fix sets `age_L1_ceil = max(minage, 1)` so the youngest age (with
# current_age = 1.0) always hits the boundary branch regardless of minage.
#
# This test must construct + evaluate the TMB AD function. If the bug
# returns the test will crash R (and CI) rather than fail gracefully, which
# is the expected signal -- a segfault is a regression we cannot ignore.

testthat::skip_on_cran()
testthat::test_that("minage = 0 + vonBertalanffy growth builds without segfault", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 6, nages = 5, minage = 0, seed = 42,
                        growth = "vonBertalanffy")
  testthat::expect_equal(dat$minage, 0)

  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  testthat::expect_true(is.list(res))
  testthat::expect_true(!is.null(res$obj))

  # Evaluating the AD function exercises the templated growth code; a
  # negative-index access would crash here.
  fval <- res$obj$fn()
  gval <- res$obj$gr()
  testthat::expect_true(is.finite(fval))
  testthat::expect_true(all(is.finite(gval)))

  # Slot 1 (= age 0 at minage = 0) must have a finite, non-zero weight.
  wt_slot1_yr1 <- res$quantities$weight_hat[1, 1, 1, 1]
  testthat::expect_true(is.finite(wt_slot1_yr1))
  testthat::expect_true(wt_slot1_yr1 > 0)

  # --- Numerical invariant for the minage = 0 anchor bug ---
  # With constant growth parameters across years, length at the end of age 0
  # (stored as length_hat[slot, sex, age = 1, yr]) must be the same in every
  # year. The pre-fix bug used `last_linear = Lmin_sp + b_len * age_L1`, which
  # collapses to Lmin_sp when minage = 0 (age_L1 = 0); this made year 1 use
  # `l1` as the cohort anchor (via the yr == 0 closed form) and years > 1 use
  # `Lmin_sp` (via the boundary branch). The fix uses `age_L1_safe`, restoring
  # `last_linear = l1` for both cases. The check below would fail under the
  # old code whenever l1 != Lmin_sp.
  lhat_age0 <- res$quantities$length_hat[1, 1, 1, seq_len(dat$endyr - dat$styr + 1L)]
  testthat::expect_true(all(is.finite(lhat_age0)))
  testthat::expect_lt(diff(range(lhat_age0)), 1e-8)
})

testthat::test_that("minage = 0 + Richards growth builds without segfault", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 6, nages = 5, minage = 0, seed = 42,
                        growth = "Richards")
  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = "Richards"),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  testthat::expect_true(is.list(res))
  fval <- res$obj$fn()
  gval <- res$obj$gr()
  testthat::expect_true(is.finite(fval))
  testthat::expect_true(all(is.finite(gval)))
})
