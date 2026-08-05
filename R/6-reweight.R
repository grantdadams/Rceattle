#' Iteratively reweight composition data (McAllister-Ianelli)
#'
#' @description
#' Refits a model, each time setting each fleet's composition weight to the
#' McAllister & Ianelli (1997) weight implied by the previous fit, until the
#' weights stop moving. This is the standard iterative tuning loop for
#' multinomial composition likelihoods: the weight multiplies the input sample
#' size, so fitting once and reading `Comp_weights_mcallister` gives a better
#' weight, which changes the fit, which changes the weight.
#'
#' Every fit already reports the implied weights in
#' `fit$data_list$fleet_control$Comp_weights_mcallister`. This function closes
#' the loop around them.
#'
#' The weight is a **parameter** (`comp_weights`), not a data column, so each
#' iteration carries it forward through `inits`. The `Comp_weights` column is
#' the starting value used to build a model from scratch; editing it and
#' re-fitting from an existing fit would have no effect.
#'
#' @details
#' Only fleets whose `Comp_distribution` is a multinomial family are tuned.
#' A Dirichlet-multinomial fleet estimates its own weight inside the
#' likelihood -- that is the point of the DM -- so an external tuning loop
#' would fight the estimate; those fleets are reported and left alone.
#'
#' Convergence is on the largest relative change in any tuned weight between
#' iterations. Reaching `n_iter` without meeting `tol` is a warning, not an
#' error: the partially tuned fit is still returned, and the weight history
#' shows whether it was converging.
#'
#' @param fit A fitted `Rceattle` model.
#' @param n_iter Maximum number of reweighting iterations (default 10).
#' @param tol Relative change in the weights below which to stop (default 0.01).
#' @param fleets Fleet codes to tune. Default `NULL` tunes every eligible fleet.
#' @param verbose Print the weights each iteration (default `TRUE`).
#'
#' @return The refitted `Rceattle` model, fitted with the final weights, and
#'   carrying a `reweight` element holding the weight history (one row per
#'   iteration), the number of iterations run, and whether it converged. The
#'   last row of the history is the weight the returned model was fitted with.
#'
#' @references
#' McAllister, M.K. and Ianelli, J.N. (1997) Bayesian stock assessment using
#' catch-age data and the sampling-importance resampling algorithm.
#' \emph{Canadian Journal of Fisheries and Aquatic Sciences} 54:284-300.
#'
#' @examples
#' \dontrun{
#' fit  <- fit_mod(BS2017SS, estimateMode = "Hindcast")
#' tuned <- reweight_comps(fit)
#' tuned$reweight$history      # weight per fleet, per iteration
#' }
#' @export
reweight_comps <- function(fit, n_iter = 10, tol = 0.01, fleets = NULL,
                           verbose = TRUE) {

  if (!inherits(fit, "Rceattle")) {
    stop("`fit` must be a fitted Rceattle model.", call. = FALSE)
  }
  fc <- fit$data_list$fleet_control

  # Fleets the loop can tune: they have composition data, and their likelihood
  # takes the weight as a multiplier rather than estimating it.
  # Negative years mark rows carried but not fitted, and calc_mcall_ianelli()
  # averages over fitted rows only. A fleet with no fitted rows yields a NaN
  # weight, so it must not be eligible in the first place.
  cd <- fit$data_list$comp_data
  has_comp <- unique(cd$Fleet_code[cd$Fleet_code > 0 & cd$Year > 0])
  dist <- as.character(fc$Comp_distribution)
  multinom <- fc$Fleet_code[dist %in% c("Multinomial", "MultinomialAFSC")]
  eligible <- intersect(has_comp, multinom)

  skipped <- setdiff(has_comp, eligible)
  if (length(skipped) && verbose) {
    message("Not tuning fleet(s) ", paste(skipped, collapse = ", "),
            ": their composition likelihood estimates its own weight.")
  }
  if (!is.null(fleets)) eligible <- intersect(eligible, fleets)
  if (!length(eligible)) {
    stop("No fleets to reweight: none have composition data fitted with a ",
         "multinomial likelihood.", call. = FALSE)
  }

  estimate_mode <- fit$data_list$estimateMode
  if (!is.null(estimate_mode) && estimate_mode > 1) {
    stop("`fit` was run at estimateMode ", estimate_mode, ", which does not ",
         "optimise the hindcast, so the weights cannot be tuned: every pass ",
         "would return the same value and report convergence. Refit at ",
         "\"Estimate\" or \"Hindcast\" first.", call. = FALSE)
  }
  # Refit the way the input fit was produced. .refit_like() otherwise defaults
  # to map = NULL, phase = FALSE, getsd = FALSE, which would rebuild a
  # user-supplied map, drop phasing the model may need to reach its optimum,
  # and return the tuned fit with no standard errors.
  # `phase` may be a named list of per-parameter phases, so it is carried
  # through verbatim rather than coerced to a flag.
  fc_cfg     <- fit$run_config$fit_control
  fit_phase  <- if (is.null(fc_cfg$phase)) FALSE else fc_cfg$phase
  fit_getsd  <- if (is.null(fc_cfg$getsd)) FALSE else fc_cfg$getsd
  fit_map    <- fit$map

  history <- NULL
  converged <- FALSE
  iter_run <- 0L

  for (iter in seq_len(n_iter)) {
    iter_run <- iter

    # Weight implied by the current fit, per fleet.
    implied <- fit$data_list$fleet_control$Comp_weights_mcallister[
      match(eligible, fit$data_list$fleet_control$Fleet_code)]
    if (all(is.na(implied))) {
      stop("The fit carries no McAllister-Ianelli weights; it must be fitted ",
           "with composition data before it can be reweighted.", call. = FALSE)
    }
    # A non-finite or non-positive weight cannot be applied: zero would delete
    # the fleet's composition data from the likelihood, and NaN would poison
    # every later warm start while making convergence unreachable.
    bad <- !is.finite(implied) | implied <= 0
    if (any(bad)) {
      stop("McAllister-Ianelli weight is ",
           paste(unique(implied[bad]), collapse = ", "), " for fleet(s) ",
           paste(eligible[bad], collapse = ", "),
           "; it cannot be applied. Check that those fleets have fitted ",
           "composition data with non-zero sample sizes.", call. = FALSE)
    }

    current <- fit$estimated_params$comp_weights[eligible]
    # Relative move, guarding a zero denominator.
    rel <- abs(implied - current) / pmax(abs(current), .Machine$double.eps)
    history <- rbind(history, data.frame(
      iteration = iter, Fleet_code = eligible, weight = implied,
      rel_change = rel, row.names = NULL))

    if (verbose) {
      message("Reweight ", iter, ": ",
              paste0("fleet ", eligible, " = ", signif(implied, 4),
                     collapse = "; "))
    }

    converged <- all(is.finite(rel)) && max(rel) < tol

    # Apply and refit even on the final pass, so the model returned is the one
    # fitted with the converged weights rather than the one that implied them.
    # The weight is a parameter, so it travels in `inits`, not the data column.
    inits <- fit$estimated_params
    inits$comp_weights[eligible] <- implied

    fit <- .refit_like(data_list = fit$data_list, inits = inits,
                       map = fit_map, estimateMode = estimate_mode,
                       phase = fit_phase, getsd = fit_getsd)

    if (converged) break
  }

  if (!converged) {
    warning("Composition weights had not converged to within tol = ", tol,
            " after ", iter_run, " iteration(s); returning the last fit. ",
            "Inspect `$reweight$history`.", call. = FALSE)
  }

  # The tuned weights are parameters, but the column is what a model rebuilt
  # from this data_list starts from -- write it back so saving the data and
  # re-fitting reproduces the tuned model rather than the untuned one.
  fit$data_list$fleet_control$Comp_weights[eligible] <-
    fit$estimated_params$comp_weights[eligible]

  fit$reweight <- list(history = history, iterations = iter_run,
                       converged = converged, tol = tol, fleets = eligible)
  fit
}
