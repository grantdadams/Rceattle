# =============================================================================
# Editing a composition-weight column and re-fitting from an existing fit is a
# no-op, and says so.
#
# The weights warm-start from `inits` like every other parameter, so
# fleet_control$Comp_weights (and CAAL_weights / Diet_comp_weights) is read only
# on a build from scratch. Tuning them by editing the column and re-fitting was
# the long-standing workflow, so the case where the edit cannot take effect
# warns rather than passing silently.
#
# A Dirichlet-multinomial weight is deliberately exempt: the likelihood
# estimates it, so a refit's `inits` differing from the column is the normal
# state and warning on it would fire on every diagnostic refit.
# =============================================================================

.cwarn_fit <- function() {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
}

.cwarn_refit <- function(fit, edit_fleet, value, distribution = NULL) {
  d <- fit$data_list
  d$fleet_control$Comp_weights[edit_fleet] <- value
  if (!is.null(distribution)) {
    d$fleet_control$Comp_distribution[edit_fleet] <- distribution
  }
  Rceattle::fit_mod(
    data_list = d, inits = fit$estimated_params, estimateMode = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))
}

testthat::test_that("a Comp_weights edit that cannot take effect warns", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .cwarn_fit()
  flt <- fit$data_list$fleet_control$Fleet_code[
    which(!is.na(fit$estimated_params$comp_weights))[1]]

  testthat::expect_warning(
    suppressMessages(.cwarn_refit(fit, flt, 4.2)),
    "has no effect")

  # The message has to name the column and the fleet, or it cannot be acted on.
  w <- tryCatch(suppressMessages(.cwarn_refit(fit, flt, 4.2)),
                warning = conditionMessage)
  testthat::expect_match(w, "Comp_weights")
  testthat::expect_match(w, paste0("fleet ", flt))

  # And the fit really does keep the parameter, not the column.
  refit <- suppressWarnings(suppressMessages(.cwarn_refit(fit, flt, 4.2)))
  testthat::expect_equal(unname(refit$estimated_params$comp_weights[flt]),
                         unname(fit$estimated_params$comp_weights[flt]))
})

testthat::test_that("an unedited column does not warn", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .cwarn_fit()
  flt <- fit$data_list$fleet_control$Fleet_code[
    which(!is.na(fit$estimated_params$comp_weights))[1]]

  # Same value as the fit already carries: nothing is being discarded.
  w <- tryCatch({
    suppressMessages(.cwarn_refit(
      fit, flt, unname(fit$estimated_params$comp_weights[flt])))
    NA_character_
  }, warning = conditionMessage)
  testthat::expect_false(isTRUE(grepl("has no effect", w)))
})

testthat::test_that("a Dirichlet-multinomial weight is exempt under a fixed map", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # This has to be exercised through a map that FIXES the weight. Under the
  # ordinary build_map() path a DM weight is left free, so `is.na(.fix)` already
  # excludes it and the DM clause is never consulted -- a test built that way
  # passes whether or not the clause exists.
  #
  # The case that matters is the one run_mse() / retrospective() create: a
  # debug map fixes every weight, and for a DM model the fitted parameter has
  # legitimately moved away from its column. Warning there would fire once per
  # fleet per assessment year of every simulation.
  d <- Rceattle::BS2017SS
  d$fleet_control$Comp_distribution[1] <- "DirichletMultinomial"
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  testthat::expect_equal(
    as.character(fit$data_list$fleet_control$Comp_distribution[1]),
    "DirichletMultinomial")

  map <- suppressWarnings(Rceattle::build_map(
    fit$data_list, fit$estimated_params, debug = TRUE, random_rec = FALSE))
  map$mapFactor$dummy <- as.factor(NA); map$mapList$dummy <- NA
  map$mapList$log_F[] <- seq_along(map$mapList$log_F)
  map$mapFactor$log_F <- factor(map$mapList$log_F)
  # Fixture guard: both weights must be fixed, or neither case is under test.
  testthat::expect_true(all(is.na(map$mapList$comp_weights[1:2])))

  # Fleet 1 (DM) and fleet 2 (multinomial) both diverge from their columns.
  pars <- fit$estimated_params
  pars$comp_weights[1] <- 2.5
  pars$comp_weights[2] <- 4.2

  w <- character(0)
  invisible(withCallingHandlers(
    suppressMessages(Rceattle::fit_mod(
      data_list = fit$data_list, inits = pars, map = map, estimateMode = 3,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0))),
    warning = function(x) { w <<- c(w, conditionMessage(x))
                            invokeRestart("muffleWarning") }))

  hits <- grep("has no effect", w, value = TRUE)
  # Exactly one warning, naming the multinomial fleet and not the DM one.
  testthat::expect_length(hits, 1)
  testthat::expect_match(hits, "fleet 2")
  testthat::expect_false(grepl("fleet 1", hits))
})
