# Constructive tests for the stock-recruit starting values and for which
# parameters are declared random. These paths are NOT covered by the four golden
# models (all of which use mean recruitment with random_rec = FALSE), so they are
# the regression net for:
#   * srr_prior being an alpha for Ricker but a STEEPNESS for Beverton-Holt
#   * the srr_alpha_init / srr_beta_init starting values from build_srr()
#   * init_dev being declared random only where the cpp gives it a density
# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true).
testthat::skip_on_cran()

# Minimal data_list carrying the stock-recruit switches build_params() reads.
srr_dat <- function(dat, srr_fun = 2, srr_est_mode = 1, srr_prior = 4,
                    srr_prior_sd = 1,
                    srr_alpha_init = NULL, srr_beta_init = NULL) {
  spec <- Rceattle::build_srr(srr_fun = srr_fun, srr_est_mode = srr_est_mode,
                              srr_prior = srr_prior,
                              srr_prior_sd = srr_prior_sd,
                              srr_alpha_init = srr_alpha_init,
                              srr_beta_init = srr_beta_init)
  dat$srr_fun        <- spec$srr_fun
  dat$srr_pred_fun   <- spec$srr_pred_fun
  dat$srr_est_mode   <- spec$srr_est_mode
  dat$srr_prior      <- spec$srr_prior
  dat$srr_alpha_init <- spec$srr_alpha_init
  dat$srr_beta_init  <- spec$srr_beta_init
  dat
}

testthat::test_that("srr_prior seeds alpha only where it is an alpha", {
  testthat::skip_if_not_installed("TMB")
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  # Beverton-Holt with a steepness prior: srr_prior is a steepness in (0, 1),
  # NOT an alpha. Seeding alpha with log(steepness) starts the optimiser at a
  # near-zero recruits-per-spawner and yields an NA/NaN gradient.
  for (mode in c(2, 3)) {
    p <- Rceattle::build_params(srr_dat(dat, srr_fun = 2, srr_est_mode = mode,
                                        srr_prior = 0.777,
                                        srr_prior_sd = 0.2))
    testthat::expect_false(isTRUE(all.equal(unname(p$rec_pars[1, 2]), log(0.777))))
  }

  # Ricker: srr_prior IS a prior on alpha, so it remains the alpha seed.
  for (mode in c(0, 2)) {
    p <- Rceattle::build_params(srr_dat(dat, srr_fun = 4, srr_est_mode = mode,
                                        srr_prior = 4))
    testthat::expect_equal(unname(p$rec_pars[1, 2]), log(4))
  }

  # Beverton-Holt with no prior applied (srr_est_mode 1 = "Estimated"):
  # srr_prior is inert as a prior and is retained as the alpha seed, which is
  # the long-standing behaviour existing models depend on.
  p <- Rceattle::build_params(srr_dat(dat, srr_fun = 2, srr_est_mode = 1,
                                      srr_prior = 4))
  testthat::expect_equal(unname(p$rec_pars[1, 2]), log(4))
})

testthat::test_that("build_srr() starting values reach rec_pars", {
  testthat::skip_if_not_installed("TMB")
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  alpha <- 1170.2
  beta <- 1.0725e-3
  p <- Rceattle::build_params(srr_dat(dat, srr_fun = 2, srr_alpha_init = alpha,
                                      srr_beta_init = beta))
  testthat::expect_equal(unname(p$rec_pars[1, 2]), log(alpha))
  testthat::expect_equal(unname(p$rec_pars[1, 3]), log(beta))

  # Explicit starting values win over srr_prior.
  p2 <- Rceattle::build_params(srr_dat(dat, srr_fun = 4, srr_est_mode = 2,
                                       srr_prior = 4, srr_alpha_init = alpha))
  testthat::expect_equal(unname(p2$rec_pars[1, 2]), log(alpha))

  # Absent, the defaults are untouched.
  p3 <- Rceattle::build_params(srr_dat(dat, srr_fun = 2, srr_est_mode = 3,
                                       srr_prior = 0.777, srr_prior_sd = 0.2))
  testthat::expect_equal(unname(p3$rec_pars[1, 3]), log(3))
})

testthat::test_that("srr_est_mode = 0 fixes alpha to the prior mean for every curve", {
  testthat::skip_if_not_installed("TMB")
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
  alpha <- 1170.2

  # The gate used to be Ricker-only, so a Beverton-Holt fit with supplied `inits`
  # had alpha mapped out but never set to the prior mean -- it stayed pinned at
  # build_params()' placeholder. Supplying `inits` is the path that regressed.
  for (fun in c(2, 4)) {
    inits <- suppressMessages(Rceattle::build_params(dat))
    mod <- Rceattle::fit_mod(
      data_list = dat, inits = inits, msmMode = 0,
      recFun = Rceattle::build_srr(srr_fun = fun, srr_est_mode = 0,
                                   srr_prior = alpha),
      estimateMode = 3,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))
    # alpha is fixed ...
    testthat::expect_true(is.na(mod$map$mapList$rec_pars[1, 2]),
                          info = paste0("srr_fun = ", fun))
    # ... and fixed AT the prior mean, not at the placeholder. Read the value TMB
    # actually holds: `initial_params` is snapshotted before this override runs.
    testthat::expect_equal(
      unname(exp(mod$obj$env$parList()$rec_pars[1, 2])), alpha,
      info = paste0("srr_fun = ", fun))
  }

  # Under the Ianelli configuration the curve is estimated as a recruitment
  # penalty rather than driving the hindcast, but alpha is a free parameter of it
  # just the same, so srr_est_mode = 0 must fix it there too. R0 stays estimated:
  # recruitment is R0 * exp(rec_dev) in that configuration.
  mod <- Rceattle::fit_mod(
    data_list = dat, msmMode = 0,
    recFun = Rceattle::build_srr(srr_fun = 0, srr_pred_fun = 2,
                                 srr_est_mode = 0, srr_prior = alpha),
    estimateMode = 3,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))
  map_rec <- mod$map$mapList$rec_pars[1, ]
  testthat::expect_true(is.na(map_rec[2]))    # alpha fixed
  testthat::expect_false(is.na(map_rec[1]))   # R0 still estimated
  testthat::expect_false(is.na(map_rec[3]))   # beta still estimated
  testthat::expect_equal(
    unname(exp(mod$obj$env$parList()$rec_pars[1, 2])), alpha)
})

testthat::test_that("build_srr() rejects an invalid steepness prior", {
  # For Beverton-Holt, srr_est_mode 2/3 put the prior on steepness, so srr_prior
  # must be in (0, 1).
  testthat::expect_error(
    Rceattle::build_srr(srr_fun = 2, srr_est_mode = 3, srr_prior = 4),
    "must be in \\(0, 1\\)")
  testthat::expect_error(
    Rceattle::build_srr(srr_fun = 2, srr_est_mode = 2, srr_prior = 1.5),
    "must be in \\(0, 1\\)")

  # The beta prior converts (mean, sd) to shapes by moments, which are positive
  # only when sd^2 < mu(1 - mu). Outside that the prior fails SILENTLY, because
  # TMB's lgamma-based dbeta returns a finite value for negative shapes. The
  # package default srr_prior_sd = 1 is never valid here.
  testthat::expect_error(
    Rceattle::build_srr(srr_fun = 2, srr_est_mode = 3, srr_prior = 0.75),
    "srr_prior_sd")
  testthat::expect_silent(
    Rceattle::build_srr(srr_fun = 2, srr_est_mode = 3, srr_prior = 0.75,
                        srr_prior_sd = 0.2))

  # Ricker takes an alpha-valued prior, so none of this applies.
  testthat::expect_silent(
    Rceattle::build_srr(srr_fun = 4, srr_est_mode = 2, srr_prior = 4))
  # Neither does it for a curve that applies no prior at all.
  testthat::expect_silent(
    Rceattle::build_srr(srr_fun = 2, srr_est_mode = 1, srr_prior = 4))
})

testthat::test_that("init_dev is random only where the cpp gives it a density", {
  testthat::skip_if_not_installed("TMB")
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  # ceattle applies dnorm(init_dev, -sigma^2/2, R_sd) only when initMode > 1.
  # Declaring init_dev random anywhere else makes the Laplace approximation
  # integrate over an improper (flat) prior instead of estimating it as a fixed
  # effect.
  expected <- c("0" = FALSE, "1" = FALSE, "2" = TRUE, "3" = TRUE, "4" = TRUE)

  for (mode in names(expected)) {
    mod <- Rceattle::fit_mod(data_list = dat, initMode = as.integer(mode),
                             random_rec = TRUE, estimateMode = 3,
                             fit_control = Rceattle::fit_control(getsd = FALSE,
                                                                 verbose = 0))
    is_random <- "init_dev" %in% unique(names(mod$obj$env$par)[mod$obj$env$random])
    testthat::expect_identical(
      is_random, unname(expected[[mode]]),
      info = paste0("initMode = ", mode)
    )

    # And the declaration must agree with whether a density was actually
    # accumulated, which is the invariant that matters.
    jnll <- mod$quantities$jnll_comp
    density <- sum(jnll[grep("Initial abundance", rownames(jnll)), ])
    testthat::expect_identical(is_random, density != 0,
                               info = paste0("initMode = ", mode))
  }

  # rec_dev stays random throughout -- this is about init_dev only.
  mod <- Rceattle::fit_mod(data_list = dat, initMode = 0, random_rec = TRUE,
                           estimateMode = 3,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  testthat::expect_true(
    "rec_dev" %in% unique(names(mod$obj$env$par)[mod$obj$env$random]))
})
