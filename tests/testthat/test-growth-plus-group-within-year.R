# Verify the plus group is PINNED at its Jan-1 length under within-year growth
# (WHAM-style), consistently for von Bertalanffy and Richards. Interior ages
# advance via within-year growth; only the plus group is held at its Jan-1
# (id_pop) value.
#
# In growth.hpp's estimate_growth_within_yr() the plus-group branch is:
#   else if (age + 1.0 == nages(sp)) length_hat = length_hat(id_pop, ...)
# for both growth families. The recruitment-weighted-mean correction for fish
# promoted into the plus group is reapplied at every year boundary by
# estimate_growth(), so id_pop already carries the corrected Jan-1 length.
#
# Setting spawn_month > 0 makes the SSB slot evaluate at fracyr > 0: the
# plus-group length must EQUAL the pop (Jan-1) plus-group length (pinned), while
# interior ages must exceed their Jan-1 length (advanced).

testthat::skip_on_cran()

# Shared check for one growth family.
.expect_plus_group_pinned <- function(growth_fun) {
  nages <- 5
  dat <- make_test_data(nyrs = 4, nages = nages, minage = 1, seed = 42,
                        growth = growth_fun)
  dat$spawn_month <- 6  # mid-year SSB -> fracyr = 0.5
  compile_tmb_if_needed()

  res <- Rceattle::fit_mod(
    data_list    = dat,
    inits        = NULL,
    estimateMode = 3,
    growthFun    = Rceattle::build_growth(fun = growth_fun),
    fit_control  = Rceattle::fit_control(phase = FALSE, verbose = 0)
  )

  # length_hat: [wt_slot, sex, age, yr]. Pop slot = 1, SSB slot = 2.
  pop_plus_jan1 <- res$quantities$length_hat[1, 1, nages, 1]
  ssb_plus_mid  <- res$quantities$length_hat[2, 1, nages, 1]
  testthat::expect_true(is.finite(pop_plus_jan1))
  testthat::expect_true(is.finite(ssb_plus_mid))

  # Plus group pinned at Jan-1: the mid-year SSB-slot plus-group length equals
  # the Jan-1 pop-slot plus-group length.
  testthat::expect_equal(ssb_plus_mid, pop_plus_jan1)

  # Interior (non-plus) ages still advance within-year, confirming fracyr is active.
  pop_a2_jan1 <- res$quantities$length_hat[1, 1, 2, 1]
  ssb_a2_mid  <- res$quantities$length_hat[2, 1, 2, 1]
  testthat::expect_gt(ssb_a2_mid, pop_a2_jan1)
}

testthat::test_that("plus group is pinned under within-year growth (von Bertalanffy)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  .expect_plus_group_pinned("vonBertalanffy")
})

testthat::test_that("plus group is pinned under within-year growth (Richards)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  .expect_plus_group_pinned("Richards")
})
