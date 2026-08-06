# =============================================================================
# The simulator draws from the mean the estimator's likelihood assumes.
#
# The index/catch lognormal likelihood fits to
#   mean = log(hat) - bias_adjust_obs * sigma^2 / 2
# (ceattle.cpp, JNLL_INDEX / JNLL_CATCH). sim_mod() used to hardcode the
# -sigma^2/2 offset, which only agreed with the estimator because .refit_like()
# happened to force bias adjustment back on in every refit. Once a refit keeps
# the source model's setting, a `bias_adjust_obs = FALSE` model would simulate
# data offset by exp(-sigma^2/2) and then fit a likelihood expecting no offset:
# a systematic bias in scale and catchability that no number of simulations
# averages away, silently corrupting self_test() and run_mse(simulate_data).
# =============================================================================

.sma_fit <- function(ba) {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0, bias_adjust_obs = ba))))
}

testthat::test_that("sim_mod() applies the model's own observation bias adjustment", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  for (ba in list(TRUE, FALSE)) {
    fit <- .sma_fit(ba)
    stored <- as.numeric(fit$data_list$bias_adjust_obs)[1]
    testthat::expect_equal(stored, as.numeric(isTRUE(ba)))

    set.seed(42)
    sim <- suppressWarnings(suppressMessages(
      Rceattle::sim_mod(fit, simulate = TRUE)))

    # Reproduce the draw the cpp's likelihood implies for this model.
    ih <- fit$quantities$index_hat
    sd <- fit$quantities$log_index_sd
    set.seed(42)
    want <- exp(stats::rnorm(length(ih), mean = log(ih) - stored * sd^2 / 2, sd = sd))

    testthat::expect_equal(as.numeric(sim$index_data$Observation), as.numeric(want),
                           info = paste("bias_adjust_obs =", ba))
  }
})

testthat::test_that("turning bias adjustment off actually moves the simulated data", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Fixture guard: if the two draws coincided the test above would pass
  # vacuously. The offset is exp(-sigma^2/2), so it must be a real shift.
  on  <- .sma_fit(TRUE)
  off <- .sma_fit(FALSE)
  set.seed(1); a <- suppressWarnings(suppressMessages(Rceattle::sim_mod(on,  TRUE)))
  set.seed(1); b <- suppressWarnings(suppressMessages(Rceattle::sim_mod(off, TRUE)))

  testthat::expect_false(isTRUE(all.equal(as.numeric(a$index_data$Observation),
                                          as.numeric(b$index_data$Observation))))
  # And in the expected direction: no bias adjustment means a larger mean.
  testthat::expect_gt(mean(b$index_data$Observation), mean(a$index_data$Observation))
})
