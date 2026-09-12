# Dynamic B0 is the fitted history with no fishing: the stock-recruit curve at
# the dynamic SSB times each year's realized deviation from the curve. Under the
# Ianelli form the hindcast is on mean recruitment until srr_mse_switchyr, so
# rec_dev there is measured from R0, not from the curve; the dynamic runs take
# log R - log(curve at the hindcast SSB) instead. With fishing near zero the
# dynamic SSB equals the hindcast's, so dynamic B0 equals hindcast biomass in
# every year, while numbers-at-age stay above the hindcast's 0.001 floor (the
# dynamic runs have none). Taking rec_dev directly (up to 5.33.0) put the curve over a
# mean-based deviation, which on the Pacific hake operating model ran dynamic
# recruitment about 1.5 times too high in every year before the switch.

.ianelli_no_f_pars <- function(m, seed = 1) {
  set.seed(seed)
  pars <- m$obj$env$last.par
  pars[names(pars) == "log_F"] <- -30
  rid <- which(names(pars) == "rec_dev")
  testthat::expect_gt(length(rid), 0)
  pars[rid] <- stats::rnorm(length(rid), 0, 0.3)
  pars
}

test_that("under the Ianelli form dynamic B0 follows the realized deviation from the curve", {
  testthat::skip_on_cran()
  d  <- make_test_data()
  sw <- d$endyr - 5
  for (curve in c("BevertonHolt", "Ricker")) {
    m <- suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = NULL, estimateMode = 3, random_rec = FALSE,
      recFun = build_srr(srr_fun = "mean", srr_pred_fun = curve,
                         srr_mse_switchyr = sw),
      fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
    pars <- .ianelli_no_f_pars(m)
    r    <- m$obj$report(pars)
    nh   <- length(d$styr:d$endyr)
    expect_lt(max(abs(r$DynamicB0[1, 1:nh] / r$biomass[1, 1:nh] - 1)), 1e-8, label = curve)

    # Not vacuous: before the switch the curve is not R0, so rec_dev is not the
    # deviation from the curve there.
    pre <- which((d$styr:d$endyr) <= sw)[-1]
    expect_gt(max(abs(r$R_hat[1, pre] / r$R0[1] - 1)), 1e-3, label = curve)
  }
})

test_that("the same holds under predation, with the fitted suitability", {
  testthat::skip_on_cran()
  set.seed(1)
  d  <- make_msm_test_data()$data_list
  sw <- d$endyr - 5
  m  <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    recFun = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                       srr_est_mode = "Estimated", srr_mse_switchyr = sw),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  pars <- .ianelli_no_f_pars(m)
  r    <- m$obj$report(pars)
  nh   <- length(d$styr:d$endyr)

  # The hindcast floors numbers-at-age at 0.001 and the dynamic runs do not, so
  # compare only species-years whose every age has stayed well above it. At these
  # start values the prey reaches the floor in year 2, and the predator only after
  # the switch year, so the predator's pre-switch years -- where the change acts --
  # are all compared.
  pre   <- which((d$styr:d$endyr) <= sw)[-1]
  min_n <- apply(r$N_at_age[, , , 1:nh, drop = FALSE], c(1, 4), min)
  clear <- t(apply(min_n > 1, 1, cumprod)) == 1
  expect_true(all(clear[1, pre]))
  rel <- r$DynamicB0[, 1:nh] / r$biomass[, 1:nh] - 1
  expect_lt(max(abs(rel[clear])), 1e-8)

  expect_gt(max(abs(r$R_hat[1, pre] / exp(m$obj$env$parList(pars)$rec_pars[1, 1]) - 1)), 1e-3)
})
