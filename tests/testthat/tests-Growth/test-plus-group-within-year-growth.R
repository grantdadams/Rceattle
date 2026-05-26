# Verify the plus group ADVANCES through within-year growth (SS3
# convention) rather than being pinned at its Jan-1 value (the old
# WHAM-style behavior the previous code implemented).
#
# In growth.hpp's estimate_growth_within_yr() the plus-group branch was:
#   else if (age + 1 == nages(sp)) length_hat = length_hat(id_pop, ...)
# which held the plus group at its Jan-1 length. That branch is gone;
# the plus group now falls through to the standard VBGF advance:
#   L_within = linf + (L_jan1 - linf) * exp(-K * fracyr)
#
# Setting spawn_month > 0 makes the SSB slot evaluate at fracyr > 0, so
# its plus-group length must exceed the pop (Jan-1) plus-group length.

testthat::skip_on_cran()
testthat::test_that("plus group advances under within-year growth (SSB > pop at fracyr > 0)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  dat <- make_test_data(nyrs = 4, nages = nages, minage = 1, seed = 42,
                        growth = "vonBertalanffy")
  dat$spawn_month <- 6  # mid-year SSB -> fracyr = 0.5
  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  # length_hat: [wt_slot, sex, age, yr]. Pop slot = 1, SSB slot = 2.
  pop_plus_jan1 <- res$quantities$length_hat[1, 1, nages, 1]
  ssb_plus_mid  <- res$quantities$length_hat[2, 1, nages, 1]

  testthat::expect_true(is.finite(pop_plus_jan1))
  testthat::expect_true(is.finite(ssb_plus_mid))

  # The plus-group mean length at start of year sits below linf (because
  # of the static exp(-0.2*a) correction in estimate_growth's section 2).
  # Advancing by fracyr > 0 of VBGF must move it strictly toward linf.
  testthat::expect_gt(ssb_plus_mid, pop_plus_jan1)

  # Sanity check non-plus ages also advance, to confirm fracyr is active.
  pop_a2_jan1 <- res$quantities$length_hat[1, 1, 2, 1]
  ssb_a2_mid  <- res$quantities$length_hat[2, 1, 2, 1]
  testthat::expect_gt(ssb_a2_mid, pop_a2_jan1)
})
