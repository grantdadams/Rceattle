# =============================================================================
# A refit keeps the source model's likelihood definition.
#
# fit_control() carries two kinds of setting: optimizer knobs, which a
# diagnostic is free to choose for itself (a retrospective phases, an MSE does
# not), and settings that define the OBJECTIVE. Bias adjustment is the second
# kind -- a model fitted with it off and refitted with it on is a different
# likelihood, so Mohn's rho, a jitter spread or an MSE trajectory would be
# comparing two models rather than one.
#
# It cannot be left to fit_control()'s defaults. fit_mod() resets the data_list
# copies to 1 and then re-applies fit_control's values, whose defaults are TRUE,
# so a control rebuilt from scratch silently switches bias adjustment back on
# (worth ~880 jnll units on BS2017SS). .refit_like() recovers the resolved
# values from the data_list instead, which covers every caller --
# retrospective(), jitter(), self_test(), profile(), run_mse(), remove_F(),
# sample_rec(), reweight_comps() -- without each having to pass them.
# =============================================================================

.rlc_fit <- function(obs, proc) {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3,
    fit_control = Rceattle::fit_control(
      phase = FALSE, getsd = FALSE, verbose = 0,
      bias_adjust_obs = obs, bias_adjust_proc = proc))))
}

.rlc_refit <- function(fit) {
  suppressWarnings(suppressMessages(Rceattle:::.refit_like(
    data_list = fit$data_list, inits = fit$estimated_params,
    estimateMode = 3)))
}

testthat::test_that("a refit preserves the source model's bias adjustment", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  for (case in list(c(FALSE, FALSE), c(TRUE, TRUE),
                    c(TRUE, FALSE), c(FALSE, TRUE))) {
    fit   <- .rlc_fit(case[1], case[2])
    refit <- .rlc_refit(fit)
    label <- paste0("obs=", case[1], " proc=", case[2])

    testthat::expect_equal(refit$data_list$bias_adjust_obs,
                           fit$data_list$bias_adjust_obs, info = label)
    testthat::expect_equal(refit$data_list$bias_adjust_proc,
                           fit$data_list$bias_adjust_proc, info = label)
    # Mode 3 returns the real objective, so this is a genuine numeric check.
    testthat::expect_equal(refit$obj$fn(refit$obj$par),
                           fit$obj$fn(fit$obj$par),
                           tolerance = 1e-10, info = label)
  }
})

testthat::test_that("bias adjustment actually moves the objective", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Fixture guard: if these ever coincide the test above passes vacuously.
  off <- .rlc_fit(FALSE, FALSE)
  on  <- .rlc_fit(TRUE, TRUE)
  testthat::expect_false(isTRUE(all.equal(off$obj$fn(off$obj$par),
                                          on$obj$fn(on$obj$par))))
})

testthat::test_that("a data_list that never went through a fit gets the default", {
  # .refit_like() reads the multipliers off the data_list; a data_list with no
  # resolved values must fall back to fit_control()'s own default rather than
  # NA-poisoning the control.
  testthat::expect_true(Rceattle:::.refit_bias_adjust(NULL))
  testthat::expect_true(Rceattle:::.refit_bias_adjust(numeric(0)))
  testthat::expect_true(Rceattle:::.refit_bias_adjust(NA))
})

testthat::test_that("a fractional bias adjustment is carried, not rounded", {
  # These are DATA_SCALARs the cpp uses as a plain multiplier on the correction
  # (bias_adjust_obs * sigma^2 / 2), so a partial ramp is a real setting.
  # Coercing to a logical would quantize 0.5 up to 1 and reintroduce a smaller
  # version of the bug this recovery exists to fix -- so the value comes back
  # numeric and unchanged.
  testthat::expect_identical(Rceattle:::.refit_bias_adjust(1), 1)
  testthat::expect_identical(Rceattle:::.refit_bias_adjust(0), 0)
  testthat::expect_identical(Rceattle:::.refit_bias_adjust(0.5), 0.5)
  testthat::expect_identical(Rceattle:::.refit_bias_adjust(c(0.25, 9)), 0.25)
})

testthat::test_that("a fractional bias adjustment survives a refit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0, bias_adjust_proc = 0.5))))
  testthat::expect_equal(as.numeric(fit$data_list$bias_adjust_proc), 0.5)

  refit <- .rlc_refit(fit)
  testthat::expect_equal(as.numeric(refit$data_list$bias_adjust_proc), 0.5)
  testthat::expect_equal(refit$obj$fn(refit$obj$par), fit$obj$fn(fit$obj$par),
                         tolerance = 1e-10)
})

testthat::test_that("comp_offset also survives a refit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # comp_offset is the other data-side likelihood constant on fit_control. It
  # already survived (its default is NULL, so fit_mod() keeps the data_list
  # value) -- pinned so it stays that way.
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0, comp_offset = 1e-3))))
  testthat::expect_equal(fit$data_list$comp_offset, 1e-3)

  refit <- .rlc_refit(fit)
  testthat::expect_equal(refit$data_list$comp_offset, 1e-3)
  testthat::expect_equal(refit$obj$fn(refit$obj$par), fit$obj$fn(fit$obj$par),
                         tolerance = 1e-10)
})
