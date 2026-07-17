# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
#
# Covers the covariance (MVN) survey-index likelihood (fleet_control
# Index_loglike == "MVN") and its companion arithmetic-mean analytical q
# (Catchability == "AnalyticalArith"), which together reproduce the AMAK/ebswp
# DoCovBTS covariance survey likelihood 0.5 * r' Sigma^-1 r.
testthat::skip_on_cran()

testthat::test_that("MVN survey likelihood equals 0.5 * r' Sigma^-1 r on the model state", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  nyrs  <- 8
  dat   <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  # A non-diagonal, positive-definite covariance for the Survey's nyrs biomass
  # observations (equicorrelated, rho = 0.3), so the off-diagonal Sigma^-1 terms
  # are actually exercised.
  sds   <- rep(20, nyrs)
  Rho   <- matrix(0.3, nyrs, nyrs); diag(Rho) <- 1
  Sigma <- diag(sds) %*% Rho %*% diag(sds)

  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_loglike[srv] <- "MVN"
  dat$fleet_control$Catchability[srv]  <- "AnalyticalArith"
  dat$index_cov <- list(Survey = Sigma)

  fit <- Rceattle::fit_mod(data_list = dat, estimateMode = 3,
                           fit_control = fit_control(phase = FALSE, verbose = 0))

  rep  <- fit$obj$report()
  td   <- fit$obj$env$data
  ictl <- td$index_ctl
  iobs <- td$index_obs
  ihat <- rep$index_hat

  # Fitted Survey (Fleet_code 1) rows, in index_obs order (== Sigma row order)
  sel <- which(ictl[, 1] == 1 & ictl[, 3] > 0 & ictl[, 3] <= td$endyr & iobs[, 1] > 0)
  testthat::expect_equal(length(sel), nyrs)

  # (a) AnalyticalArith q: total predicted index == total observed index
  testthat::expect_equal(sum(ihat[sel]), sum(iobs[sel, 1]), tolerance = 1e-6)

  # (b) the TMB covariance likelihood equals the plain-R quadratic form on the
  #     same predicted state (this is the core correctness guarantee)
  r    <- iobs[sel, 1] - ihat[sel]
  ll_R <- as.numeric(0.5 * t(r) %*% solve(Sigma) %*% r)
  testthat::expect_equal(rep$jnll_comp[1, 1], ll_R, tolerance = 1e-8)
})

testthat::test_that("MVNORM reports the full density = MVN (bare) + normalizing constant", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  nyrs  <- 8
  dat   <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)
  sds   <- rep(20, nyrs)
  Rho   <- matrix(0.3, nyrs, nyrs); diag(Rho) <- 1
  Sigma <- diag(sds) %*% Rho %*% diag(sds)
  const <- 0.5 * (as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) + nyrs * log(2 * pi))

  bts_jnll <- function(mode) {
    d <- dat
    srv <- d$fleet_control$Fleet_name == "Survey"
    d$fleet_control$Index_loglike[srv] <- mode
    d$fleet_control$Catchability[srv]  <- "AnalyticalArith"
    d$index_cov <- list(Survey = Sigma)
    fit <- Rceattle::fit_mod(data_list = d, estimateMode = 3,
                             fit_control = fit_control(phase = FALSE, verbose = 0))
    fit$obj$report()$jnll_comp[1, 1]
  }
  # Same fixed state, so the two families differ only by the fixed constant
  testthat::expect_equal(bts_jnll("MVNORM") - bts_jnll("MVN"), const, tolerance = 1e-6)
})

testthat::test_that("MVN and lognormal index likelihoods differ, and lognormal is the default", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  nyrs  <- 8
  base  <- make_test_data(nyrs = nyrs, nages = nages, seed = 7)

  # Default (no Index_loglike override) is Lognormal
  testthat::expect_true(all(base$fleet_control$Index_loglike == "Lognormal"))

  ss_ln <- Rceattle::fit_mod(data_list = base, estimateMode = 3,
                             fit_control = fit_control(phase = FALSE, verbose = 0))

  mvn <- base
  sds   <- rep(15, nyrs)
  Sigma <- diag(sds^2)                      # diagonal Sigma (independent, but MVN form)
  srv <- mvn$fleet_control$Fleet_name == "Survey"
  mvn$fleet_control$Index_loglike[srv] <- "MVN"
  mvn$index_cov <- list(Survey = Sigma)
  ss_mvn <- Rceattle::fit_mod(data_list = mvn, estimateMode = 3,
                              fit_control = fit_control(phase = FALSE, verbose = 0))

  # Different likelihood family -> different survey-biomass jnll component
  testthat::expect_false(isTRUE(all.equal(ss_ln$obj$report()$jnll_comp[1, 1],
                                          ss_mvn$obj$report()$jnll_comp[1, 1])))
})

testthat::test_that("data_check rejects invalid MVN covariance input", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nyrs <- 8
  mk_mvn <- function() {
    d <- make_test_data(nyrs = nyrs, nages = 5, seed = 1)
    d$fleet_control$Index_loglike[d$fleet_control$Fleet_name == "Survey"] <- "MVN"
    d
  }
  fc <- fit_control(phase = FALSE, verbose = 0)

  # No covariance matrix supplied for an MVN fleet
  d1 <- mk_mvn()
  testthat::expect_error(
    Rceattle::fit_mod(data_list = d1, estimateMode = 3, fit_control = fc),
    "no covariance matrix"
  )

  # Wrong-dimension covariance matrix
  d2 <- mk_mvn(); d2$index_cov <- list(Survey = diag(3))
  testthat::expect_error(
    Rceattle::fit_mod(data_list = d2, estimateMode = 3, fit_control = fc),
    "fitted survey observations"
  )

  # Non-symmetric covariance matrix
  d3 <- mk_mvn()
  S  <- diag(nyrs); S[1, 2] <- 5
  d3$index_cov <- list(Survey = S)
  testthat::expect_error(
    Rceattle::fit_mod(data_list = d3, estimateMode = 3, fit_control = fc),
    "not symmetric"
  )
})
