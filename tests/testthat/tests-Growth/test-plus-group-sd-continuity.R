# Verify the plus group's SD-at-length comes from the same length-based
# linear interpolation as every other age above age_L1, with no
# discontinuity at the boundary.
#
# The pre-fix branch in growth.hpp used `else if (age == nages - 1)
# length_sd = exp(growth_log_sd(sp, sex, 1))`, which pinned the plus
# group at the upper SD anchor regardless of where its mean length sits.
# After the fix the plus group is computed via:
#   sd = sd0 + (sd1 - sd0) / (linf - l1) * (length_hat - l1)
# so its SD should land "near" the SD at the previous age, scaled by
# the small change in length.
#
# We can't read length_sd directly (it's a local in estimate_growth),
# so we recover it from each row of growth_matrix as the second-moment
# spread of P(length | age). Bin midpoints are unknown in physical
# units, so we use the bin index as a unitless proxy -- the comparison
# only relies on the RATIO of two implied SDs, which is dimensionless.

testthat::skip_on_cran()
testthat::test_that("plus-group SD is continuous with the prior age (no exp(sd1) jump)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 6  # need enough ages that the plus group is well above age_L1
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

  gm <- res$quantities$growth_matrix
  nlen <- dat$nlengths
  bin_idx <- seq_len(nlen)  # unitless proxy for bin midpoints

  implied_sd <- function(a, yr = 1) {
    p <- gm[1, 1, a, bin_idx, yr]
    mu <- sum(p * bin_idx)
    sqrt(sum(p * (bin_idx - mu)^2))
  }

  sd_prev <- implied_sd(nages - 1L)
  sd_plus <- implied_sd(nages)

  testthat::expect_true(is.finite(sd_prev))
  testthat::expect_true(is.finite(sd_plus))
  testthat::expect_gt(sd_prev, 0)
  testthat::expect_gt(sd_plus, 0)

  # Under linear length-based interpolation with constant growth params,
  # the SD difference between consecutive ages near the plus group should
  # be small compared to the SD itself. A factor-of-2 jump would indicate
  # the old `length_sd = exp(sd1)` plus-group override is back.
  ratio <- sd_plus / sd_prev
  testthat::expect_gt(ratio, 0.5)
  testthat::expect_lt(ratio, 2.0)
})
