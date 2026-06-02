# Richards counterpart to test-growth-minage0.R's VBGF invariant.
#
# Same `last_linear = Lmin_sp + b_len * age_L1` -> `* age_L1_safe` fix applies
# to the Richards branch in estimate_growth(). Under constant growth params
# across years, length at the end of age 0 (length_hat[slot, sex, 1, yr])
# must be identical for every hindcast year. The pre-fix code at minage = 0
# would have year 1 use `l1` and years > 1 use `Lmin_sp` as the cohort
# anchor (in L^m space), so this fails whenever l1 != Lmin_sp.

testthat::skip_on_cran()
testthat::test_that("Richards + minage = 0: length at age 0 is year-invariant", {
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

  nyrs <- dat$endyr - dat$styr + 1L
  lhat_age0 <- res$quantities$length_hat[1, 1, 1, seq_len(nyrs)]
  testthat::expect_true(all(is.finite(lhat_age0)))
  testthat::expect_lt(diff(range(lhat_age0)), 1e-8)
})
