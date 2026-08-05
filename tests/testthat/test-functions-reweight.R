# =============================================================================
# Iterative McAllister-Ianelli composition reweighting.
#
# The composition weight is a PARAMETER, so it warm-starts from `inits` like
# any other. That is what makes an iterative tuning loop possible at all: each
# pass hands the previous fit's implied weight to the next through `inits`,
# not through the fleet_control column (which is only the starting value used
# when a model is built from scratch).
# =============================================================================

.rw_fit <- function() {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE,
                                        verbose = 0))))
}


testthat::test_that("a composition weight warm-starts from inits", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")

  # The premise the loop rests on: a weight handed in through `inits` survives.
  # It used to be overwritten from fleet_control on every fit, which made an
  # iterative loop impossible -- each pass silently reset to the column.
  d <- Rceattle::BS2017SS
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                        verbose = 0))))
  testthat::expect_equal(unname(fit$estimated_params$comp_weights),
                         unname(d$fleet_control$Comp_weights))

  p <- fit$estimated_params
  p$comp_weights[1] <- 0.37
  warm <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = p, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                        verbose = 0))))
  testthat::expect_equal(unname(warm$estimated_params$comp_weights[1]), 0.37)
})


testthat::test_that("reweight_comps() converges and returns the tuned fit", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")

  tuned <- suppressWarnings(
    Rceattle::reweight_comps(.rw_fit(), n_iter = 8, tol = 0.02, verbose = FALSE))

  testthat::expect_true(tuned$reweight$converged)
  testthat::expect_gt(tuned$reweight$iterations, 1)

  # The model returned must be the one fitted WITH the final weights, not the
  # one that implied them -- otherwise the reported history and the model
  # disagree by one iteration.
  h <- tuned$reweight$history
  last <- h[h$iteration == max(h$iteration), ]
  testthat::expect_equal(
    unname(tuned$estimated_params$comp_weights[last$Fleet_code]),
    last$weight)

  # Tuning has to actually move the weights off their starting value of 1.
  testthat::expect_false(isTRUE(all.equal(last$weight,
                                          rep(1, nrow(last)))))
  # And the relative change must shrink -- a loop that is not converging is a
  # different failure from one that has not finished.
  worst <- tapply(h$rel_change, h$iteration, max)
  testthat::expect_lt(worst[[length(worst)]], worst[[1]])
})


testthat::test_that("reweight_comps() reports rather than silently stopping short", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")

  # One iteration cannot converge from a standing start, so it must warn and
  # still hand back the partially tuned fit with its history.
  testthat::expect_warning(
    tuned <- Rceattle::reweight_comps(.rw_fit(), n_iter = 1, tol = 1e-8,
                                      verbose = FALSE),
    "not converged")
  testthat::expect_false(tuned$reweight$converged)
  testthat::expect_equal(tuned$reweight$iterations, 1)
  testthat::expect_s3_class(tuned$reweight$history, "data.frame")
})


testthat::test_that("Dirichlet-multinomial fleets are left alone", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")

  # A DM fleet estimates its own weight inside the likelihood, so tuning it
  # externally would fight the estimate.
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) < 2)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts[1]] <- "DirichletMultinomial"

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1, msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE,
                                        verbose = 0))))
  tuned <- suppressWarnings(
    Rceattle::reweight_comps(fit, n_iter = 2, tol = 0.05, verbose = FALSE))

  testthat::expect_false(comp_flts[1] %in% tuned$reweight$fleets)
  testthat::expect_true(all(comp_flts[-1] %in% tuned$reweight$fleets))
})


testthat::test_that("a requested fleet that cannot be tuned is named", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")

  # Asking for a fleet and quietly getting a subset would report a weighting
  # the caller never asked for -- the failure mode this loop exists to avoid.
  # A Dirichlet-multinomial fleet is ineligible by construction, so requesting
  # one alongside a tunable fleet is the case to catch.
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) < 2)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts[1]] <- "DirichletMultinomial"

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1, msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE,
                                        verbose = 0))))

  testthat::expect_warning(
    tuned <- Rceattle::reweight_comps(fit, fleets = comp_flts[1:2],
                                      n_iter = 1, tol = 1, verbose = FALSE),
    "Not tuning requested fleet")
  testthat::expect_equal(tuned$reweight$fleets, comp_flts[2])
})
