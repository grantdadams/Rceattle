# Lognormal priors and the Ianelli stock-recruit penalty follow bias_adjust_proc:
# centred at -sd^2/2 when TRUE, so the prior value (or R_hat) is a mean, and
# uncentred when FALSE, so it is a median. Each check writes the density out by
# hand against a built (not optimized) model and reads the reported jnll row.

.ba_build <- function(bias, data = Rceattle::BS2017SS, ...) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = data, estimateMode = 3, msmMode = 0, random_rec = FALSE, ...,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0, getsd = FALSE,
                                        bias_adjust_proc = bias))))
}
.ba_row <- function(fit, row) fit$quantities$jnll_comp[row, ]

.ricker_prior <- function() {   # deprecated form, still checked against the linkage form
  suppressWarnings(Rceattle::build_srr(srr_fun = "Ricker", srr_est_mode = "LognormalPrior",
                                       srr_prior = 5, srr_prior_sd = 0.5))
}

# At the starting value a centre shifted the wrong way gives the same density, so
# a sign error would pass; each check evaluates 0.3 away from the start.
.ba_report_at <- function(fit, par_name, shift = 0.3) {
  p <- fit$obj$par
  i <- names(p) == par_name
  p[i] <- p[i] + shift
  list(rep = fit$obj$report(p), pars = fit$obj$env$parList(p))
}

testthat::test_that("the Ricker alpha prior is mean-centred under bias_adjust_proc", {
  testthat::skip_on_cran()
  for (bias in c(TRUE, FALSE)) {
    fit <- .ba_build(bias, recFun = .ricker_prior())
    at  <- .ba_report_at(fit, "rec_pars")
    la  <- at$pars$rec_pars[, 2]
    testthat::expect_false(isTRUE(all.equal(unname(la), rep(log(5), length(la)))))
    testthat::expect_equal(
      unname(at$rep$jnll_comp[9, seq_along(la)]),   # row 9: stock-recruit prior
      unname(-dnorm(la, log(5) - bias * 0.5^2 / 2, 0.5, log = TRUE)),
      tolerance = 1e-10, info = as.character(bias))
  }
})

testthat::test_that("the linkage lognormal prior on alpha gives the same objective as LognormalPrior", {
  testthat::skip_on_cran()
  for (bias in c(TRUE, FALSE)) {
    a <- .ba_build(bias, recFun = .ricker_prior())
    b <- .ba_build(bias, recFun = Rceattle::build_srr(
      srr_fun = "Ricker",
      linkages = list(alpha = Rceattle::linkage_spec(
        ~ 1, init = list(`(Intercept)` = 5),
        priors = list(`(Intercept)` = lognormal(log(5), 0.5))))))
    testthat::expect_equal(sum(.ba_row(b, "Linkage-table priors")),
                           sum(.ba_row(a, "Stock-recruit prior")),
                           tolerance = 1e-10, info = as.character(bias))
    testthat::expect_equal(b$obj$fn(), a$obj$fn(), tolerance = 1e-8, info = as.character(bias))
  }
})

testthat::test_that("the M prior is mean-centred under bias_adjust_proc, the median without", {
  testthat::skip_on_cran()
  for (bias in c(TRUE, FALSE)) {
    fit <- .ba_build(bias, M1Fun = Rceattle::build_M1(
      M1_model = 1, M1_use_prior = TRUE, M_prior = 0.2, M_prior_sd = 0.1))
    # M1_model 1 estimates one M1 for all ages, stored in log_M1[sp, 1, 1] of the
    # parameter list TMB uses (the build's starting values differ by age).
    lm  <- fit$obj$env$parList()$log_M1[, 1, 1]
    testthat::expect_equal(
      unname(.ba_row(fit, "M prior")[seq_along(lm)]),
      unname(-dnorm(lm, log(0.2) - bias * 0.1^2 / 2, 0.1, log = TRUE)),
      tolerance = 1e-10, info = as.character(bias))
  }
})

testthat::test_that("the catchability prior is mean-centred under bias_adjust_proc", {
  testthat::skip_on_cran()
  d   <- Rceattle::BS2017SS
  srv <- which(d$fleet_control$Fleet_type == 2)[1]
  d$fleet_control$Catchability[srv]          <- "Estimated-with-prior"
  d$fleet_control$Catchability_init[srv]     <- 1
  d$fleet_control$Catchability_prior_sd[srv] <- 0.2
  flt <- d$fleet_control$Fleet_code[srv]
  for (bias in c(TRUE, FALSE)) {
    fit <- .ba_build(bias, data = d)
    at  <- .ba_report_at(fit, "index_log_q")
    lq  <- unname(at$pars$index_log_q[flt])
    testthat::expect_false(isTRUE(all.equal(lq, 0)))
    testthat::expect_equal(
      unname(at$rep$jnll_comp[7, flt]),   # row 7: catchability prior
      -dnorm(lq, log(1) - bias * 0.2^2 / 2, 0.2, log = TRUE),
      tolerance = 1e-10, info = as.character(bias))
  }
})

testthat::test_that(".srr_hat_cols() matches the template's penalty years", {
  dl <- list(styr = 1980, srr_hat_styr = 1981, srr_hat_endyr = 2000)
  testthat::expect_identical(Rceattle:::.srr_hat_cols(dl, 30), 2:21)
  # styr is never in the penalty, and the window ends with the hindcast or peel.
  dl$srr_hat_styr <- 1980
  testthat::expect_identical(Rceattle:::.srr_hat_cols(dl, 10), 2:10)
  # No penalty years in range: empty, so each caller decides what to do.
  dl$srr_hat_styr <- 2005
  testthat::expect_identical(Rceattle:::.srr_hat_cols(dl, 30), integer(0))
  dl$srr_hat_styr <- NA
  testthat::expect_error(Rceattle:::.srr_hat_cols(dl, 30), "must be set")
})

testthat::test_that("a lognormal prior on a random-effect SD is mean-centred under bias_adjust_proc", {
  testthat::skip_on_cran()
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L,
                               priors = list(sigma = lognormal(log(0.05), 0.3)))))
  for (bias in c(TRUE, FALSE)) {
    fit <- .ba_build(bias, qFun = qfun)
    ls  <- fit$estimated_params$log_sigma_linkage
    testthat::expect_equal(
      sum(.ba_row(fit, "Linkage-table priors")),
      sum(-dnorm(ls, log(0.05) - bias * 0.3^2 / 2, 0.3, log = TRUE)),
      tolerance = 1e-10, info = as.character(bias))
  }
})

testthat::test_that("the Ianelli penalty centres R_hat as the mean under bias_adjust_proc", {
  testthat::skip_on_cran()
  for (bias in c(TRUE, FALSE)) {
    fit <- .ba_build(bias, recFun = Rceattle::build_srr(
      srr_fun = "mean", srr_pred_fun = "BevertonHolt"))
    q   <- fit$quantities
    dl  <- fit$data_list
    yrs <- (dl$srr_hat_styr:dl$srr_hat_endyr) - dl$styr + 1
    sd  <- exp(fit$estimated_params$R_log_sd)
    expected <- vapply(seq_len(dl$nspp), function(sp) -sum(dnorm(
      log(q$R[sp, yrs]), log(q$R_hat[sp, yrs]) - bias * sd[sp]^2 / 2, sd[sp],
      log = TRUE)), numeric(1))
    testthat::expect_equal(
      unname(.ba_row(fit, "Stock-recruit penalty")[seq_len(dl$nspp)]),
      expected, tolerance = 1e-8, info = as.character(bias))
  }
})

testthat::test_that("sample_rec() projects the Ianelli form at the mean ratio over the penalty years", {
  testthat::skip_on_cran()
  fit <- .ba_build(TRUE, recFun = Rceattle::build_srr(
    srr_fun = "mean", srr_pred_fun = "BevertonHolt"))
  q   <- fit$quantities
  dl  <- fit$data_list
  nh  <- dl$endyr - dl$styr + 1
  yrs <- (dl$srr_hat_styr:dl$srr_hat_endyr) - dl$styr + 1
  # The first year, never in the penalty, is outside the average.
  testthat::expect_false(1 %in% yrs)
  s <- sample_rec(fit, sample_rec = FALSE, update_model = FALSE)
  testthat::expect_equal(unname(s$estimated_params$rec_dev[, nh + 1]),
                         unname(log(rowMeans((q$R / q$R_hat)[, yrs]))),
                         tolerance = 1e-12)
})
