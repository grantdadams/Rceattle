# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
#
# Covers the covariance (MVN) survey-index likelihood (fleet_control
# Index_distribution == "MVN") and its companion arithmetic-mean analytical q
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
  dat$fleet_control$Index_distribution[srv] <- "MVN"
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
    d$fleet_control$Index_distribution[srv] <- mode
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

  # Default (no Index_distribution override) is Lognormal
  testthat::expect_true(all(base$fleet_control$Index_distribution == "Lognormal"))

  ss_ln <- Rceattle::fit_mod(data_list = base, estimateMode = 3,
                             fit_control = fit_control(phase = FALSE, verbose = 0))

  mvn <- base
  sds   <- rep(15, nyrs)
  Sigma <- diag(sds^2)                      # diagonal Sigma (independent, but MVN form)
  srv <- mvn$fleet_control$Fleet_name == "Survey"
  mvn$fleet_control$Index_distribution[srv] <- "MVN"
  mvn$index_cov <- list(Survey = Sigma)
  ss_mvn <- Rceattle::fit_mod(data_list = mvn, estimateMode = 3,
                              fit_control = fit_control(phase = FALSE, verbose = 0))

  # Different likelihood family -> different survey-biomass jnll component
  testthat::expect_false(isTRUE(all.equal(ss_ln$obj$report()$jnll_comp[1, 1],
                                          ss_mvn$obj$report()$jnll_comp[1, 1])))
})

testthat::test_that("index_cov is re-aligned when the fitted year range changes (retro / subset / MSE growth)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages <- 5
  nyrs  <- 8
  dat   <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)
  sds   <- rep(20, nyrs)
  Rho   <- matrix(0.3, nyrs, nyrs); diag(Rho) <- 1
  Sigma <- diag(sds) %*% Rho %*% diag(sds)
  srv   <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- "MVN"
  dat$fleet_control$Catchability[srv]  <- "AnalyticalArith"
  dat$index_cov <- list(Survey = Sigma)

  # A fresh fit tags the stored Sigma with its fitted years, so it is
  # self-describing for every subsequent re-fit (retrospective / MSE / jitter).
  fit <- Rceattle::fit_mod(data_list = dat, estimateMode = 3,
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  testthat::expect_equal(rownames(fit$data_list$index_cov$Survey), as.character(1:nyrs))

  # (a) Shrink: a retrospective peel drops the fitted set to 1:5; the Sigma is
  #     subset to the retained years with the covariance block preserved. A real
  #     peel also clamps suit_endyr to the new endyr (see retrospective()), so
  #     mirror that here rather than leaving suit_endyr past the fitted range.
  dl2 <- fit$data_list; dl2$endyr <- 5
  dl2$suit_endyr <- pmin(dl2$suit_endyr, dl2$endyr)
  fit2 <- Rceattle::fit_mod(data_list = dl2, estimateMode = 3,
                            fit_control = fit_control(phase = FALSE, verbose = 0))
  S2 <- fit2$data_list$index_cov$Survey
  testthat::expect_equal(dim(S2), c(5L, 5L))
  testthat::expect_equal(unname(S2), unname(Sigma[1:5, 1:5]))

  # ...and the covariance likelihood equals 0.5 r' Sigma_sub^-1 r on the subset.
  rep  <- fit2$obj$report(); td <- fit2$obj$env$data
  sel  <- which(td$index_ctl[, 1] == 1 & td$index_ctl[, 3] > 0 &
                  td$index_ctl[, 3] <= td$endyr & td$index_obs[, 1] > 0)
  testthat::expect_equal(length(sel), 5L)
  r    <- td$index_obs[sel, 1] - rep$index_hat[sel]
  ll_R <- as.numeric(0.5 * t(r) %*% solve(Sigma[1:5, 1:5]) %*% r)
  testthat::expect_equal(rep$jnll_comp[1, 1], ll_R, tolerance = 1e-8)

  # (b) Grow: an MSE assessment step appends a future survey observation. The new
  #     year is added as an independent diagonal block with variance
  #     (Observation * Log_sd)^2, retained years keep their covariance.
  ic <- Rceattle:::.align_index_cov  # exercise the aligner directly for the grow path
  dl3 <- fit$data_list
  new_obs <- data.frame(Fleet_code = 1L, Year = nyrs + 1L, Observation = 123,
                        Log_sd = 0.25)
  # match existing index_data columns
  base_row <- dl3$index_data[dl3$index_data$Fleet_code == 1L, ][1, ]
  new_row  <- base_row
  new_row$Year <- nyrs + 1L; new_row$Observation <- 123; new_row$Log_sd <- 0.25
  dl3$index_data <- rbind(dl3$index_data, new_row)
  dl3$endyr <- nyrs + 1L
  dl3 <- ic(dl3)
  S3  <- dl3$index_cov$Survey
  testthat::expect_equal(dim(S3), c(nyrs + 1L, nyrs + 1L))
  testthat::expect_equal(unname(S3[1:nyrs, 1:nyrs]), unname(Sigma))          # old block intact
  testthat::expect_equal(S3[nyrs + 1L, nyrs + 1L], (123 * 0.25)^2)           # new indep. variance
  testthat::expect_equal(S3[nyrs + 1L, 1L], 0)                               # no spurious cross term
})

# The two natural-scale families, fitted on the same data, differ by exactly the
# truncation constant. "Normal" is the plain density -- it is what the ADMB
# bridges compare against term for term, so it must stay untruncated --
# and "TruncatedNormal" renormalizes it over (0, Inf).
.natural_scale_fixture <- function(family, sd = 20) {
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 42)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- family
  dat$fleet_control$Catchability[srv]  <- "AnalyticalArith"
  dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd   # ABSOLUTE sd
  Rceattle::fit_mod(data_list = dat, estimateMode = 3,
                    fit_control = fit_control(phase = FALSE, verbose = 0))
}

testthat::test_that("Normal index likelihood is the natural-scale -dnorm with an absolute sd", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  fit <- .natural_scale_fixture("Normal")
  rep <- fit$obj$report(); td <- fit$obj$env$data
  sel <- which(td$index_ctl[, 1] == 1 & td$index_ctl[, 3] > 0 &
                 td$index_ctl[, 3] <= td$endyr & td$index_obs[, 1] > 0)
  # 0.5*(obs - q*pred)^2 / sd^2 + normalizing constant = -sum dnorm(obs, pred, sd).
  # No truncation term: that is "TruncatedNormal", below.
  ll_R <- -sum(dnorm(td$index_obs[sel, 1], rep$index_hat[sel], 20, log = TRUE))
  testthat::expect_equal(rep$jnll_comp[1, 1], ll_R, tolerance = 1e-8)

  # Regression: a non-MVN family (0/3/4) carries a 1x1 dummy Sigma and must NOT
  # enter the MVN block (Index_distribution >= 1 there previously segfaulted for
  # "Normal" by applying MVNORM(1x1) to a length-nyrs residual).
  testthat::expect_false(is.null(fit$obj))
})

testthat::test_that("TruncatedNormal adds the log Phi(mu/sd) renormalization and nothing else", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # sd large relative to the index, so the constant is far from zero and the
  # comparison has power. mu > 0 always, so log Phi(mu/sd) lies in [log(0.5), 0]
  # and the term is bounded however extreme the fixture.
  fit_t <- .natural_scale_fixture("TruncatedNormal", sd = 80)
  fit_n <- .natural_scale_fixture("Normal",          sd = 80)

  rep_t <- fit_t$obj$report(); td <- fit_t$obj$env$data
  rep_n <- fit_n$obj$report()
  sel <- which(td$index_ctl[, 1] == 1 & td$index_ctl[, 3] > 0 &
                 td$index_ctl[, 3] <= td$endyr & td$index_obs[, 1] > 0)

  ll_R <- -sum(dnorm(td$index_obs[sel, 1], rep_t$index_hat[sel], 80, log = TRUE)) +
    sum(log(pnorm(rep_t$index_hat[sel] / 80)))
  testthat::expect_equal(rep_t$jnll_comp[1, 1], ll_R, tolerance = 1e-8)

  # The two families sit on the same predicted index -- the truncation constant
  # is the ONLY difference between them...
  testthat::expect_equal(rep_t$index_hat, rep_n$index_hat, tolerance = 1e-10)
  testthat::expect_equal(rep_t$jnll_comp[1, 1] - rep_n$jnll_comp[1, 1],
                         sum(log(pnorm(rep_t$index_hat[sel] / 80))), tolerance = 1e-8)
  # ...and it is large enough here that the test would have caught a sign error.
  testthat::expect_gt(abs(rep_t$jnll_comp[1, 1] - rep_n$jnll_comp[1, 1]), 1)
})

testthat::test_that("data_check rejects invalid MVN covariance input", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nyrs <- 8
  mk_mvn <- function() {
    d <- make_test_data(nyrs = nyrs, nages = 5, seed = 1)
    d$fleet_control$Index_distribution[d$fleet_control$Fleet_name == "Survey"] <- "MVN"
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
