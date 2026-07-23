# =============================================================================
# Non-parametric selectivity penalties can be specified as standard deviations
# instead of the cryptic penalty WEIGHTS. Every such penalty is a Gaussian SSQ
# (weight * x^2 = x^2 / (2*sd^2)), so switch_check() converts the SD columns
# (Sel_shape_sd [+ Sel_shape_dir], Sel_curvature_sd, Sel_devmag_sd) into
# Sel_curve_pen1/2/3 via weight = 1/(2*sd^2). Legacy Sel_curve_pen columns are
# left untouched (so existing models are bit-identical); a fleet supplying the SD
# columns fits equivalently to the same model expressed as weights.
# =============================================================================

testthat::test_that("switch_check converts penalty SDs to Sel_curve_pen weights", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  testthat::skip_if(length(np) == 0)
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_   # drop legacy so conversion fires
  d$fleet_control$Sel_curve_pen2[np] <- NA_real_
  d$fleet_control$Sel_shape_sd      <- NA_real_; d$fleet_control$Sel_shape_sd[np]      <- 1 / sqrt(2 * 20)
  d$fleet_control$Sel_shape_dir     <- NA;       d$fleet_control$Sel_shape_dir[np]     <- "Decreasing"
  d$fleet_control$Sel_curvature_sd  <- NA_real_; d$fleet_control$Sel_curvature_sd[np]  <- 1 / sqrt(2 * 12.5)

  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$fleet_control$Sel_curve_pen1[np], rep(20, length(np)),  tolerance = 1e-8)
  testthat::expect_equal(out$fleet_control$Sel_curve_pen2[np], rep(12.5, length(np)), tolerance = 1e-8)

  # "Increasing" flips the shape-penalty sign; legacy weight is never overwritten.
  d2 <- d; d2$fleet_control$Sel_shape_dir[np] <- "Increasing"
  out2 <- suppressMessages(Rceattle::switch_check(d2))
  testthat::expect_equal(out2$fleet_control$Sel_curve_pen1[np], rep(-20, length(np)), tolerance = 1e-8)
})

testthat::test_that("penalty SD columns are rejected on non-non-parametric forms", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  # fleet 4 is logistic (Selectivity == 1): a penalty SD there is meaningless and
  # (for AR1 forms) would corrupt the logit-rho reuse of Sel_curve_pen.
  logi <- which(d$fleet_control$Selectivity == 1)
  testthat::skip_if(length(logi) == 0)
  d$fleet_control$Sel_shape_sd <- NA_real_; d$fleet_control$Sel_shape_sd[logi[1]] <- 0.2
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "NonParametric")
})

testthat::test_that("a non-positive penalty SD is rejected", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  testthat::skip_if(length(np) == 0)
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_
  d$fleet_control$Sel_shape_sd <- NA_real_; d$fleet_control$Sel_shape_sd[np[1]] <- 0
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "positive standard deviation")
})

testthat::test_that("penalty SD columns fit equivalently to the legacy weights", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  ctl <- Rceattle::fit_control(phase = TRUE, verbose = 0)

  base <- suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 0, random_rec = FALSE,
    msmMode = 0, fit_control = ctl))

  d <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  d$fleet_control$Sel_shape_sd     <- NA_real_; d$fleet_control$Sel_shape_sd[np]     <- 1 / sqrt(2 * d$fleet_control$Sel_curve_pen1[np])
  d$fleet_control$Sel_curvature_sd <- NA_real_; d$fleet_control$Sel_curvature_sd[np] <- 1 / sqrt(2 * d$fleet_control$Sel_curve_pen2[np])
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_
  d$fleet_control$Sel_curve_pen2[np] <- NA_real_
  sdfit <- suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 0, random_rec = FALSE, msmMode = 0, fit_control = ctl))

  # Equivalent up to the floating-point of the sd<->weight round-trip amplified
  # by the optimisation (~1e-7 on an objective ~1e4, i.e. relative ~1e-11).
  testthat::expect_equal(sdfit$opt$objective, base$opt$objective, tolerance = 1e-5)
})
