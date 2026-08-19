# The natural-scale index likelihood (Index_distribution = "Normal") is a normal
# LEFT-TRUNCATED AT ZERO, not a plain normal.
#
# An index cannot be negative and data_check() will not accept one, so the
# support is (0, inf) and the density has to be renormalized over it:
#   log f(x) = log phi(x; mu, sd) - log(1 - Phi(-mu/sd)),  1 - Phi(-mu/sd) = Phi(mu/sd)
# Without that term the density does not integrate to one over the values the
# data can take, and the simulator could not draw from the same distribution the
# likelihood scores.
#
# No golden model uses this family -- all three are lognormal -- so /golden-check
# is silent here and this is the regression net.

.tn_fixture <- function(sd, nyrs = 20) {
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 123)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- "Normal"
  dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
  dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
}

testthat::test_that("the Normal index likelihood carries the truncation constant", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # sd large relative to the index, so the constant is far from zero and the
  # check has power: at mu/sd ~ 1 it is worth ~0.17 nll per observation.
  fit <- .tn_fixture(sd = 80)
  q   <- fit$quantities
  idx <- fit$data_list$index_data
  srv <- idx$Fleet_name == "Survey" & idx$Year > 0 &
    idx$Year <= fit$data_list$endyr & idx$Observation > 0

  hat <- as.numeric(q$index_hat)[srv]
  obs <- idx$Observation[srv]
  sd_ <- as.numeric(exp(q$index_ln_sd %||% log(80)))[1]
  if (!is.finite(sd_)) sd_ <- 80

  plain <- -sum(stats::dnorm(obs, hat, sd_, log = TRUE))
  trunc <- plain + sum(log(stats::pnorm(hat / sd_)))
  got   <- sum(q$jnll_comp[1, ])          # JNLL_INDEX

  # The reported index nll must match the TRUNCATED form, not the plain one.
  testthat::expect_equal(got, trunc, tolerance = 1e-6)
  # ...and the two must actually differ, or the test proves nothing.
  testthat::expect_gt(abs(trunc - plain), 1)
})

testthat::test_that("a truncated Normal index draw is always positive", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Inverse-CDF sampling on (0, inf): positivity holds by construction, with no
  # rejection loop and no retry budget that could run out. sd is deliberately
  # larger than the index, where an untruncated normal would go negative often.
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    local({
      dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
      srv <- dat$fleet_control$Fleet_name == "Survey"
      dat$fleet_control$Index_distribution[srv] <- "Normal"
      dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
      dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- 150
      dat
    }),
    file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  set.seed(4)
  for (k in 1:15) {
    sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))
    srv <- sim$index_data$Fleet_name == "Survey"
    testthat::expect_true(all(sim$index_data$Observation[srv] > 0))
  }
})
