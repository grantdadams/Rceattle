# Sanity check on the size-transition matrix:
# P(length | age) must integrate to 1 across length bins for every age.
#
# Catches future edits to the pnorm bin construction in section 4 of
# estimate_growth() (and estimate_growth_within_yr()), where the first
# bin is a length minus-group, interior bins are pnorm differences, and
# the last bin is a length plus-group. Drift in any of those branches
# would leave rows summing to != 1.

testthat::skip_on_cran()
testthat::test_that("growth_matrix rows sum to 1 for every age and year (vonBertalanffy)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  dat <- make_test_data(nyrs = 4, nages = nages, minage = 1, seed = 42,
                        growth = "vonBertalanffy")
  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  # growth_matrix: [wt_slot, sex, age, length_bin, yr]
  gm <- res$quantities$growth_matrix
  nyrs <- dat$endyr - dat$styr + 1L
  nlen <- dat$nlengths

  for (yr in seq_len(nyrs)) {
    for (a in seq_len(nages)) {
      row <- gm[1, 1, a, seq_len(nlen), yr]
      testthat::expect_true(all(is.finite(row)),
                            label = paste0("finite row at age=", a, " yr=", yr))
      testthat::expect_lt(abs(sum(row) - 1), 1e-10,
                          label = paste0("row sum at age=", a, " yr=", yr))
    }
  }
})

testthat::test_that("growth_matrix rows sum to 1 with minage = 0 (vonBertalanffy)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  dat <- make_test_data(nyrs = 4, nages = nages, minage = 0, seed = 42,
                        growth = "vonBertalanffy")
  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  gm <- res$quantities$growth_matrix
  nlen <- dat$nlengths
  for (a in seq_len(nages)) {
    row <- gm[1, 1, a, seq_len(nlen), 1]
    testthat::expect_lt(abs(sum(row) - 1), 1e-10,
                        label = paste0("minage=0 row sum at age=", a))
  }
})
