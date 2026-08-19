# Integration test (fits a CEATTLE model): skipped on CRAN.
#
# sim_mod() returns the model's EXPECTED data (simulate = FALSE) or SAMPLED data
# (simulate = TRUE). The composition sampling branch (rmultinom / Dirichlet-
# multinomial draws, keyed by the string Comp_distribution) was previously
# untested -- it is where a comp/CAAL code-vs-string mismatch survived (PR 0.5).
#
# Scope: this covers the comp_data Multinomial path (the primary previously-
# untested branch). The CAAL sampling branch and the DirichletMultinomial branch
# are not yet exercised here -- a worthwhile follow-up.
testthat::skip_on_cran()

testthat::test_that("sim_mod returns expected comp proportions (simulate = FALSE)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
  fit <- suppressMessages(fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  testthat::expect_gt(nrow(fit$data_list$comp_data), 0)

  exp_dat <- Rceattle::sim_mod(fit, simulate = FALSE)
  comp_cols <- grep("^Comp_", names(exp_dat$comp_data))
  got <- unname(as.matrix(exp_dat$comp_data[, comp_cols]))
  # Expected data are exactly the model's predicted composition (comp_hat).
  testthat::expect_equal(got, unname(fit$quantities$comp_hat), tolerance = 1e-10)
})

testthat::test_that("sim_mod multinomial draws sum to sample size (simulate = TRUE)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
  fit <- suppressMessages(fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))

  set.seed(42)
  sim_dat <- Rceattle::sim_mod(fit, simulate = TRUE)
  comp_cols <- grep("^Comp_", names(sim_dat$comp_data))
  counts <- as.matrix(sim_dat$comp_data[, comp_cols])

  # rmultinom draws: non-negative integers summing to the row's Sample_size.
  testthat::expect_true(all(counts >= 0, na.rm = TRUE))
  testthat::expect_equal(counts, round(counts))
  testthat::expect_equal(unname(rowSums(counts, na.rm = TRUE)),
                         unname(sim_dat$comp_data$Sample_size),
                         tolerance = 1e-8)
  # A draw differs from the deterministic expectation (it actually sampled).
  exp_dat <- Rceattle::sim_mod(fit, simulate = FALSE)
  testthat::expect_false(isTRUE(all.equal(
    counts, unname(as.matrix(exp_dat$comp_data[, comp_cols])))))
})


# --- Survey index: one draw per Index_distribution -------------------------
# sim_mod() drew every fleet as an independent lognormal regardless of
# Index_distribution, so an MVN/MVNORM fleet was simulated on the wrong scale,
# with the wrong spread, and with none of the correlation its Sigma carries.
# That does not error: self_test() runs and reports recovery against a
# data-generating process the likelihood never assumed.
#
# Build a survey whose predicted index is flat, so the realised spread and
# correlation of the draws can be read off directly and compared with the
# fleet's own likelihood.
.sim_index_fixture <- function(dist, nyrs = 20, sd = 20, rho = 0.6, cv = 0.2) {
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 123)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- dist
  dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
  if (dist %in% c("MVN", "MVNORM")) {
    Rho <- matrix(rho, nyrs, nyrs)
    diag(Rho) <- 1
    dat$index_cov <- list(Survey = diag(rep(sd, nyrs)) %*% Rho %*%
                            diag(rep(sd, nyrs)))
  }
  if (dist == "Normal") {
    dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd
  }
  if (dist == "Lognormal") {
    dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- cv
  }
  suppressMessages(suppressWarnings(fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
}

# Realised sd and mean pairwise correlation of the survey draws. Correlation is
# measured ACROSS replicates: within a single replicate the sample mean absorbs
# the shared factor of a compound-symmetry Sigma, returning about -1/(n-1) no
# matter what rho is.
.sim_index_moments <- function(fit, nrep = 400, log_scale = FALSE) {
  srv <- fit$data_list$index_data$Fleet_name == "Survey"
  hat <- fit$quantities$index_hat[srv]
  reps <- replicate(nrep, {
    suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))$
      index_data$Observation[srv]
  })
  dev <- if (log_scale) log(reps) - log(hat) else reps - hat
  C <- stats::cor(t(dev))
  list(sd = mean(apply(dev, 1, stats::sd)), cor = mean(C[upper.tri(C)]))
}

testthat::test_that("sim_mod draws the survey index under each Index_distribution", {
  testthat::skip_if_not_installed("TMB")
  set.seed(42)

  # Lognormal: log scale, independent.
  m <- .sim_index_moments(.sim_index_fixture("Lognormal"), log_scale = TRUE)
  testthat::expect_equal(m$sd, 0.2, tolerance = 0.05)
  testthat::expect_equal(m$cor, 0, tolerance = 0.05)

  # Normal: NATURAL scale with an absolute sd, independent -- not lognormal.
  m <- .sim_index_moments(.sim_index_fixture("Normal"))
  testthat::expect_equal(m$sd, 20, tolerance = 0.1)
  testthat::expect_equal(m$cor, 0, tolerance = 0.05)

  # MVN / MVNORM: natural scale AND correlated. The two differ only by a
  # constant in the likelihood, so they must simulate identically.
  for (dist in c("MVN", "MVNORM")) {
    m <- .sim_index_moments(.sim_index_fixture(dist))
    testthat::expect_equal(m$sd, 20, tolerance = 0.1, info = dist)
    testthat::expect_equal(m$cor, 0.6, tolerance = 0.15, info = dist)
  }
})

testthat::test_that("sim_mod keeps natural-scale index draws positive", {
  testthat::skip_if_not_installed("TMB")

  # A natural-scale draw can go negative, which no index can be, and
  # data_check() rejects Observation <= 0 -- so the data set would not refit at
  # all and self_test() would count the run as not converged, reading as a
  # convergence problem rather than a simulation one. Non-positive draws are
  # redrawn instead (normal truncated at zero).
  for (dist in c("Normal", "TruncatedNormal", "MVN")) {
    fit <- .sim_index_fixture(dist, sd = 60)      # index_hat is 100
    srv <- fit$data_list$index_data$Fleet_name == "Survey"
    set.seed(1)
    for (k in 1:20) {
      sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))
      testthat::expect_true(all(sim$index_data$Observation[srv] > 0), info = dist)
    }
  }
})

testthat::test_that("sim_mod warns when truncation is doing the work, on both branches", {
  testthat::skip_if_not_installed("TMB")

  # Positive draws are not enough on their own: a row that keeps being redrawn
  # follows a normal truncated at zero, not the normal the likelihood assumes,
  # so a self_test() built on it tests a different data-generating process.
  # The rate is per ROW and read off the worst one -- a fleet mean would hide a
  # single marginal row, which is how truncation actually presents.
  # Both untruncated natural-scale families. "Normal" is fitted as a plain normal
  # (it reproduces the AMAK avo_like / cpue_like term for term, which is what the
  # ADMB bridges compare against), and MVN cannot be truncated at all -- a
  # multivariate normal truncated to the positive orthant has no closed-form
  # sampler. Either way the draw is redrawn until positive while the likelihood
  # scores the untruncated normal, and that gap is what the warning is about.
  for (dist in c("Normal", "MVN")) {
    fit <- .sim_index_fixture(dist, sd = 115)    # index_hat is 100
    set.seed(1)
    testthat::expect_warning(Rceattle::sim_mod(fit, simulate = TRUE),
                             "truncated at zero", info = dist)
  }

  # The warning on a univariate fleet names the family that removes the gap.
  fit <- .sim_index_fixture("Normal", sd = 115)
  set.seed(1)
  w <- testthat::capture_warnings(Rceattle::sim_mod(fit, simulate = TRUE))
  testthat::expect_true(any(grepl("TruncatedNormal", w)))

  # TruncatedNormal has no gap: it is FITTED as a normal left-truncated at zero
  # and drawn from that same distribution by inverse CDF, so however hard the
  # truncation bites the draw still follows its own likelihood. Positive by
  # construction, and nothing to warn about.
  fit <- .sim_index_fixture("TruncatedNormal", sd = 115)
  set.seed(1)
  sim <- testthat::expect_no_warning(Rceattle::sim_mod(fit, simulate = TRUE))
  srv <- sim$index_data$Fleet_name == "Survey"
  testthat::expect_true(all(sim$index_data$Observation[srv] > 0))

  # A fleet with a small sd must stay silent on every branch.
  for (dist in c("Normal", "TruncatedNormal", "MVN")) {
    fit <- .sim_index_fixture(dist, sd = 5)
    set.seed(1)
    testthat::expect_no_warning(Rceattle::sim_mod(fit, simulate = TRUE))
  }
})

testthat::test_that("sim_mod keeps a wide correlated fleet positive", {
  testthat::skip_if_not_installed("TMB")

  # An MVN vector is rejected if ANY row is non-positive, so the joint rejection
  # probability climbs with the number of rows. A flat retry budget failed here
  # while every row on its own was fine; the budget scales with the row count.
  fit <- .sim_index_fixture("MVN", nyrs = 42, sd = 70, rho = 0)  # index_hat 100
  srv <- fit$data_list$index_data$Fleet_name == "Survey"
  set.seed(1)
  for (k in 1:10) {
    sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))
    testthat::expect_true(all(sim$index_data$Observation[srv] > 0))
  }
})

testthat::test_that("an MVN fleet simulates from the covariance it was fitted with", {
  testthat::skip_if_not_installed("TMB")

  # The correlated draw is taken by the model, from the Sigma the model holds, so
  # editing index_cov after the fit does not change it -- simulating reflects the
  # model as fitted, not the data_list as since edited. (A model fitted WITHOUT a
  # covariance for an MVN fleet is rejected by data_check() before it can get
  # here, which is where that error belongs.)
  fit <- .sim_index_fixture("MVN")
  set.seed(5); want <- Rceattle::sim_mod(fit, simulate = TRUE)$index_data$Observation

  edited <- fit
  edited$data_list$index_cov <- list()   # Sigma lost (e.g. a lossy round-trip)
  set.seed(5); got <- Rceattle::sim_mod(edited, simulate = TRUE)$index_data$Observation
  testthat::expect_identical(got, want)
})


testthat::test_that("sim_mod handles a model with no diet data", {
  testthat::skip_if_not_installed("TMB")

  # Stomach contents are drawn by the TMB model now, and only for a predator
  # whose suitability is estimated. A single-species fixture with no diet table
  # must simulate without complaint. What diet does when it IS fitted -- and when
  # it is present but not fitted -- needs a multispecies model, and lives in
  # test-functions-sim-mod-diet.R.
  fit <- .sim_index_fixture("Lognormal")
  testthat::expect_equal(nrow(fit$data_list$diet_data), 0)
  testthat::expect_no_warning(Rceattle::sim_mod(fit, simulate = TRUE))
  testthat::expect_no_error(Rceattle::sim_mod(fit, simulate = FALSE))
})


testthat::test_that("a covariance survey fleet warns about rows it cannot simulate", {
  testthat::skip_if_not_installed("TMB")

  # The MVN/MVNORM covariance covers the fitted years only, so rows outside that
  # window keep their original values. run_mse() reveals exactly those rows to
  # the estimation model as each assessment's new survey data, so silence here
  # would mean an MSE evaluating against a survey that never varied.
  fit <- .sim_index_fixture("MVNORM")

  # As fitted, every survey row is inside the window: nothing to warn about.
  testthat::expect_no_warning(Rceattle::sim_mod(fit, simulate = TRUE))

  # Hide one year from the estimation model, as run_mse() does.
  hidden <- fit
  srv <- which(hidden$data_list$index_data$Fleet_name == "Survey")
  hidden$data_list$index_data$Year[utils::tail(srv, 1)] <-
    -hidden$data_list$index_data$Year[utils::tail(srv, 1)]
  testthat::expect_warning(Rceattle::sim_mod(hidden, simulate = TRUE),
                           "were NOT simulated")

  # An independent family has no such limit -- it redraws every row.
  ln <- .sim_index_fixture("Lognormal")
  srv2 <- which(ln$data_list$index_data$Fleet_name == "Survey")
  ln$data_list$index_data$Year[utils::tail(srv2, 1)] <-
    -ln$data_list$index_data$Year[utils::tail(srv2, 1)]
  testthat::expect_no_warning(Rceattle::sim_mod(ln, simulate = TRUE))
})
