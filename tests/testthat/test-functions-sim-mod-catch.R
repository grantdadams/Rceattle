# Total catch is drawn by the TMB model's SIMULATE block, beside the catch
# density (ceattle.cpp, slot 1). These tests pin what would otherwise drift
# silently: the draw's location and spread against what the likelihood assumes,
# the standard deviation being indexed per row, and simulating leaving nothing
# behind on the fitted object.

.sim_catch_fixture <- function(log_sd = 0.2, per_row = FALSE, seed = 123) {
  d <- make_test_data(nyrs = 20, nages = 5, seed = seed)
  if (per_row) {
    # A single shared Log_sd cannot tell a per-row sd from a per-fleet one, or
    # from one off by a row -- every candidate gives the same answer.
    set.seed(seed)
    d$catch_data$Log_sd <- stats::runif(nrow(d$catch_data), 0.05, 0.9)
  } else {
    d$catch_data$Log_sd <- log_sd
  }
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
}

# Rows the model actually predicts. clean_data() appends projection rows whose
# catch_hat is 0; their draw is exactly 0 by construction and carries no
# information about location or spread.
.sim_catch_fitted <- function(fit) fit$quantities$catch_hat > 0

.sim_catch_reps <- function(fit, nrep, seed = 11) {
  fitted <- .sim_catch_fitted(fit)
  set.seed(seed)
  replicate(nrep, suppressWarnings(
    Rceattle::sim_mod(fit, simulate = TRUE)$catch_data$Catch[fitted]))
}


testthat::test_that("the catch draw is centred and spread as the likelihood assumes", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .sim_catch_fixture(log_sd = 0.3)
  fitted <- .sim_catch_fitted(fit)
  hat <- as.numeric(fit$quantities$catch_hat[fitted])
  sd_row <- as.numeric(fit$quantities$log_catch_sd[fitted])

  ba <- fit$data_list$bias_adjust_obs
  if (is.null(ba)) ba <- fit$obj$env$data$bias_adjust_obs
  if (is.null(ba)) ba <- 1
  ba <- as.numeric(ba)[1]

  reps <- .sim_catch_reps(fit, nrep = 1500)

  # Location, checked separately from spread: an sd check passes straight over a
  # biased mean. ceattle.cpp fits log(obs) ~ N(log(hat) - ba*sd^2/2, sd).
  # Tolerance is set from the Monte Carlo error of the offset (~0.002 here).
  offset <- mean(rowMeans(log(reps)) - log(hat))
  testthat::expect_equal(offset, -ba * sd_row[1]^2 / 2, tolerance = 0.15)

  # Spread on the log scale.
  testthat::expect_equal(mean(apply(log(reps), 1, stats::sd)),
                         sd_row[1], tolerance = 0.03)

  # And the draw is not simply the expectation.
  testthat::expect_gt(mean(apply(log(reps), 1, stats::sd)), 0)
})


testthat::test_that("the catch sd is applied per row, not per fleet", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .sim_catch_fixture(per_row = TRUE)
  fitted <- .sim_catch_fitted(fit)
  nominal <- as.numeric(fit$quantities$log_catch_sd[fitted])
  testthat::expect_gt(stats::sd(nominal), 0.1)   # the fixture really does vary

  reps <- .sim_catch_reps(fit, nrep = 1500)
  realised <- apply(log(reps), 1, stats::sd)

  # Each row's realised spread must track its OWN nominal sd. A per-fleet or
  # off-by-one sd would break this while leaving the pooled mean sd intact.
  testthat::expect_gt(stats::cor(realised, nominal), 0.99)
  testthat::expect_equal(realised, nominal, tolerance = 0.15)
})


testthat::test_that("projection rows are drawn, not left NA", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .sim_catch_fixture()
  # clean_data() appends projection rows with Catch = NA. The SIMULATE block
  # sits outside the hindcast gate so they get filled; left NA, data_check()
  # rejects the refit and self_test() reports a convergence failure instead.
  set.seed(4)
  sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))

  testthat::expect_equal(nrow(sim$catch_data), nrow(fit$data_list$catch_data))
  testthat::expect_false(anyNA(sim$catch_data$Catch))
  # catch_hat is 0 there, so log(0) = -inf, rnorm returns without drawing, and
  # exp(-inf) = 0.
  proj <- !.sim_catch_fitted(fit)
  testthat::expect_true(all(sim$catch_data$Catch[proj] == 0))
})


testthat::test_that("simulate = FALSE returns the expected values and draws nothing", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .sim_catch_fixture()

  set.seed(3); before <- .Random.seed
  sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = FALSE))

  # run_mse(simulate_data = FALSE) and tools/verify/verify-mse-om-horizon.R both
  # rest on this path consuming no random numbers at all.
  testthat::expect_identical(.Random.seed, before)
  testthat::expect_equal(as.numeric(sim$catch_data$Catch),
                         as.numeric(fit$quantities$catch_hat))
})


testthat::test_that("simulating leaves the fitted object usable", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .sim_catch_fixture()
  before_obj <- fit$obj$fn(fit$obj$env$last.par.best)
  before_dat <- fit$obj$env$data
  before_rep <- names(fit$obj$report(fit$obj$env$last.par.best))

  set.seed(5); invisible(suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE)))

  testthat::expect_identical(fit$obj$fn(fit$obj$env$last.par.best), before_obj)
  testthat::expect_identical(fit$obj$env$data, before_dat)

  # TMB never clears the report environment, so anything a SIMULATE block
  # reports stays visible on the object. Hence the `_sim` names: a leftover entry
  # must never be readable as the observed data it replaced.
  after_rep <- names(fit$obj$report(fit$obj$env$last.par.best))
  leaked <- setdiff(after_rep, before_rep)
  testthat::expect_true(all(grepl("_sim$", leaked)))
  testthat::expect_false(any(c("catch_obs", "obsvec") %in% leaked))
})


testthat::test_that("sim_mod() refuses to guess when it cannot simulate", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .sim_catch_fixture()

  # mse_summary() and model_average() drop $obj; simulating needs it.
  no_obj <- fit; no_obj$obj <- NULL
  testthat::expect_error(suppressWarnings(Rceattle::sim_mod(no_obj, TRUE)),
                         "needs the model's TMB object")
  # ... but the expected-value path reads only $quantities.
  testthat::expect_no_error(suppressWarnings(Rceattle::sim_mod(no_obj, FALSE)))

  # A data_list edited after the fit no longer lines up row for row. The old R
  # draw recycled the shorter vector silently and returned a wrong answer.
  short <- fit
  short$data_list$catch_data <- short$data_list$catch_data[1:5, ]
  testthat::expect_error(suppressWarnings(Rceattle::sim_mod(short, TRUE)),
                         "line up row for row")
})
