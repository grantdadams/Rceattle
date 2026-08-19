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

testthat::test_that("LogisticPM accepts Sel_shape_sd / Sel_devmag_sd (two-sided, positive)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  f <- which(d$fleet_control$Selectivity == 1)[1]   # repurpose a logistic fleet
  testthat::skip_if(length(f) == 0)
  d$fleet_control$Selectivity[f] <- 11              # LogisticPM
  d$fleet_control$Sel_shape_sd  <- NA_real_; d$fleet_control$Sel_shape_sd[f]  <- 1 / sqrt(2 * 20)
  d$fleet_control$Sel_devmag_sd <- NA_real_; d$fleet_control$Sel_devmag_sd[f] <- 1 / sqrt(2 * 8)
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$fleet_control$Sel_curve_pen1[f], 20, tolerance = 1e-8)  # pen1 (shape)
  testthat::expect_equal(out$fleet_control$Sel_curve_pen3[f], 8,  tolerance = 1e-8)  # pen3 (dev-mag)
})

testthat::test_that("Sel_curvature_sd is rejected on LogisticPM (pen2 unused there)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  f <- which(d$fleet_control$Selectivity == 1)[1]
  testthat::skip_if(length(f) == 0)
  d$fleet_control$Selectivity[f] <- 11
  d$fleet_control$Sel_curvature_sd <- NA_real_; d$fleet_control$Sel_curvature_sd[f] <- 0.2
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)), "NonParametric")
})

testthat::test_that("Sel_shape_dir = 'Increasing' is rejected on LogisticPM (two-sided)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  f <- which(d$fleet_control$Selectivity == 1)[1]
  testthat::skip_if(length(f) == 0)
  d$fleet_control$Selectivity[f] <- 11
  d$fleet_control$Sel_shape_sd  <- NA_real_; d$fleet_control$Sel_shape_sd[f]  <- 0.2
  d$fleet_control$Sel_shape_dir <- NA;       d$fleet_control$Sel_shape_dir[f] <- "Increasing"
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)), "two-sided")
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


testthat::test_that("mode 5 does not feed sel_dev_log_sd from unestimated deviates", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Time_varying_sel = "RandomWalkAscending" varies the ascending limb only, so
  # build_map() estimates no descending deviate. Those deviates sit at 0, and the
  # descending random-walk penalty on them used to be accumulated anyway. With
  # random_sel = FALSE that is a constant (sel_dev_log_sd is mapped out); with
  # random_sel = TRUE the SD IS estimated, so the spurious term is
  # 2 * nyrs * nsex * log(sigma) and it biases the SD downward. Pin the gradient
  # so that term cannot come back unnoticed for either setting.
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  flt <- which(d$fleet_control$Time_varying_sel == 5)
  testthat::skip_if(length(flt) == 0)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, random_sel = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # The SD is estimable under random_sel = TRUE -- otherwise this test is vacuous.
  testthat::expect_true("sel_dev_log_sd" %in% names(fit$obj$par))

  # Its gradient must come only from deviates the model actually estimates. The
  # descending deviates are mapped out at 0, so a descending penalty would add a
  # term with no data behind it; the ascending pair is what legitimately informs
  # the SD.
  g <- fit$obj$gr(fit$obj$par)
  gsd <- g[names(fit$obj$par) == "sel_dev_log_sd"]
  testthat::expect_true(all(is.finite(gsd)))
  testthat::expect_equal(length(gsd), sum(names(fit$obj$par) == "sel_dev_log_sd"))
})
