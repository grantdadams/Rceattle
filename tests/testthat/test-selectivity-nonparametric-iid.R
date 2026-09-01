# Time-varying non-parametric (Ianelli, type 2) selectivity.
#
# Provenance. `Time_varying_sel = "IID"` on a non-parametric fleet was refused by
# build_map() and data_check() even though the bundled `Atka2022` -- and the ADMB
# bridge it comes from (`Rceattle-models/BSAI atka mackerel/2022 atka
# bridging.R`) -- is configured exactly that way, so that dataset could not be
# fitted as shipped. 5.25.0 scores the mode instead of refusing it:
# `dnorm(sel_coff_dev, 0, sel_dev_sd)` on each estimated coefficient.
#
# Why the two modes are not interchangeable, and why only one may be integrated:
# "RandomWalk" scores the year-to-year change in the REALIZED log-selectivity,
# which is renormalized to mean 1 within each year. The level of a year's
# coefficients does not enter that density, so `sel_dev_sd` is not identified
# from it and Laplace-integrating the deviates drives it to zero -- measured at
# 2.7e-8 on Atka2022, which is a time-invariant selectivity reported as a
# time-varying one. "IID" scores the deviations themselves and does identify the
# sd. But the AMAK shape penalty beside it is one-sided, so the Laplace objective
# is only piecewise smooth: fit_mod() refuses `random_sel = TRUE` for the walk
# because the density is improper, and for IID while `Sel_curve_pen1` is non-zero
# because the optimizer stops at a kink rather than an optimum.
testthat::skip_on_cran()

testthat::test_that("NonParametric accepts IID and the bundled Atka2022 fits", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  # As shipped: fishery is Selectivity 2 with Time_varying_sel 1 (IID), sd 0.35.
  testthat::expect_equal(as.numeric(Atka2022$fleet_control$Selectivity[2]), 2)
  testthat::expect_equal(as.numeric(Atka2022$fleet_control$Time_varying_sel[2]), 1)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = Atka2022, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  testthat::expect_true(is.finite(fit$opt$objective))

  # The deviates are estimated and the curve actually moves between years.
  sel <- fit$quantities$sel_at_age[2, 1, , ]
  testthat::expect_gt(max(apply(sel, 1, function(x) diff(range(x)))), 0)
})


testthat::test_that("random_sel = TRUE is refused for both non-parametric modes", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  d <- Atka2022
  d$fleet_control$Time_varying_sel <- as.character(d$fleet_control$Time_varying_sel)

  # The walk cannot be integrated: its density is blind to the level of each
  # year's coefficients, so the sd collapses. Refused by name, not by a bare TMB
  # "NA/NaN gradient evaluation" that mentions neither selectivity nor random_sel.
  d$fleet_control$Time_varying_sel[2] <- "RandomWalk"
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
      random_sel = TRUE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "random_sel"
  )

  # ... but it still fits as the penalized effects the AMAK formulation intends.
  testthat::expect_true(is.finite(suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, inits = NULL, msmMode = 0,
                      estimateMode = "Hindcast",
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))
  ))$opt$objective))

  # IID is refused too while the one-sided shape penalty is on, for a different
  # reason: that penalty is not twice differentiable, so the Laplace objective is
  # only piecewise smooth. Atka2022 stops at a kink with a maximum gradient of
  # 6.8 and an sd 27% from the value it reaches with the penalty off.
  d$fleet_control$Time_varying_sel[2] <- "IID"
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
      random_sel = TRUE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "Sel_curve_pen1"
  )

  # With that penalty off the deviates integrate cleanly, which is the control
  # showing the refusal is about the penalty and not about the IID density.
  d$fleet_control$Sel_curve_pen1[2] <- 0
  re <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
    random_sel = TRUE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  testthat::expect_gt(exp(re$estimated_params$sel_dev_log_sd[2]), 1e-3)
  testthat::expect_lt(
    max(abs(re$obj$gr(re$obj$env$last.par.best[re$obj$env$lfixed()]))), 1e-2)
})


testthat::test_that("NonParametricPM still refuses IID, naming the alternative", {
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  d <- Atka2022
  d$fleet_control$Selectivity <- as.character(d$fleet_control$Selectivity)
  d$fleet_control$Selectivity[2] <- "NonParametricPM"
  # Its sel_coff_dev ARE walk increments (selectivity.hpp case 9 builds each
  # year from the previous one), so an independent-deviate reading of them would
  # describe a different curve than the model draws.
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "NonParametricPM"
  )
})
