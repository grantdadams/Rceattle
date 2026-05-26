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
