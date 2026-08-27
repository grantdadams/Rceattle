#' Specify the stock-recruit relationship (SRR) for Rceattle
#'
#' @param srr_fun Stock-recruit function used in the hindcast estimation (see the list below). Default = 0
#' @param srr_pred_fun Stock-recruit function used for projection, reference points, and penalties (see below). When \code{srr_fun == 0}, the stock-recruit curve is added as a penalty on the annually estimated hindcast recruitment (following AMAK and Jim Ianelli's pollock model). If \code{srr_fun > 0}, then \code{srr_pred_fun = srr_fun} and no extra penalty is added.
#' @param proj_mean_rec Recruitment used in the projection: `TRUE`/1 (default) = mean recruitment, the average R over the hindcast; `FALSE`/0 = the stock-recruit relationship given by `srr_pred_fun`. Equilibrium and dynamic reference points follow the curve whenever `srr_pred_fun` is a stock-recruit form, regardless of this switch. A DSEM built with `estimate_projection = TRUE` projects recruitment through the SEM and overrides this; `fit_mod()` says so.
#' @param srr_hat_styr Integer. First year used to estimate the recruitment-penalty function (the AMAK/Ianelli penalty, active when \code{srr_pred_fun > 0} and \code{srr_fun = 0}), starting at \code{styr + 1}. Defaults to \code{styr + 1} in \code{data_list}. Useful when the environmental data conditioning the stock-recruit relationship is not available until the terminal year but projections are still wanted.
#' @param srr_hat_endyr Integer. Last year used to estimate the recruitment-penalty function (the AMAK/Ianelli penalty, active when \code{srr_pred_fun > 0} and \code{srr_fun = 0}). Defaults to \code{endyr} in \code{data_list}. Useful when the environmental data conditioning the stock-recruit relationship does not span the full time series but projections are still wanted.
#' @param srr_est_mode Switch to determine estimation mode. Accepts integer codes or the equivalent readable strings: 0 / "Fixed" = fix alpha to prior mean; 1 / "Estimated" = freely estimate R0, alpha, and/or beta (default); 2 / "LognormalPrior" = lognormally distributed prior for alpha (Ricker) or steepness (Beverton); 3 / "BetaPrior" = beta distributed prior for steepness (Beverton) given mean and sd.
#' @param srr_prior mean for normally distributed prior for stock-recruit parameter.
#'   Note this is not the same quantity for every curve: for Ricker
#'   (\code{srr_pred_fun} 4 / 5) it is a prior on \eqn{\alpha} itself, while for
#'   Beverton-Holt (2 / 3) it is a prior on \strong{steepness} and so must lie in
#'   (0, 1). See \code{srr_est_mode}.
#' @param srr_prior_sd Prior standard deviation for stock-recruit parameter
#' @param srr_alpha_init,srr_beta_init Optional starting values for the
#'   stock-recruit \eqn{\alpha} and \eqn{\beta} parameters, on the natural
#'   (not log) scale, one value per species. Only used when the curve actually
#'   estimates them (\code{srr_fun} or \code{srr_pred_fun} above 1); ignored for
#'   mean recruitment, where both are mapped out.
#'
#'   The package defaults (\eqn{\alpha = e^3}, \eqn{\beta = 3}) are placeholders
#'   with no knowledge of the stock's scale. \eqn{\beta} in particular sets the
#'   density dependence in \eqn{R = \alpha S / (1 + \beta S)}, so it must be on
#'   the order of \eqn{(\alpha - 1/\phi_0) / R_0} -- typically \eqn{10^{-3}} or
#'   smaller for a stock measured in tonnes. Starting three orders of magnitude
#'   away drives predicted recruitment to near zero and the optimizer returns
#'   \code{NA/NaN gradient evaluation}. For a Beverton-Holt seeded from a
#'   steepness \eqn{h} and unfished spawning biomass per recruit \eqn{\phi_0}:
#'   \deqn{\alpha = \frac{4h}{\phi_0 (1 - h)}, \qquad
#'         \beta  = \frac{\alpha - 1/\phi_0}{R_0}.}
#' @param srr_indices Soft-deprecated. Use the `linkages` argument instead. See `vignette("environmental-linkages-and-priors")`.
#' @param Bmsy_lim Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a likelihood penalty if beta is estimated above this limit. Default `NA` is not used.
#' @param srr_mse_switchyr Year at which an MSE switches from the annual recruitment-penalty estimate to the stock-recruit function (the \code{srr_fun = 0}, \code{srr_pred_fun > 0} case).
#' @param linkages Optional named list of [linkage_spec()] objects keyed by recruitment parameter name (must be one of `"R0"`, `"alpha"`, `"beta"`). Each spec describes how that parameter depends on environmental covariates and on stratifying factors (species, sex). The offset enters additively (on the log scale) inside the recruitment compute. See `vignette("environmental-linkages-and-priors")` for details.
#'
#' @description
#'
#' **Stock recruitment relationships currently implemented in Rceattle:**
#'
#' - \code{srr_fun = 0} or \code{"mean"}: No stock recruit relationship. Recruitment is a function of \eqn{R0} (on the log scale) and annual deviates (i.e. steepness = 0.99).
#'  \deqn{R_y = exp(R0 + R_{dev,y})}
#'
#' - \code{srr_fun = 2} or \code{"BevertonHolt"}: Beverton-holt stock-recruitment relationship
#'   \deqn{R_y = \frac{\alpha_{srr} * SB_{y-minage}}{1+\beta_{srr} * SB_{y-minage}}}
#'
#' - \code{srr_fun = 4} or \code{"Ricker"}: Ricker stock-recruitment relationship
#'   \deqn{R_y = \alpha_{srr} * SB_{y-minage} * exp(-\beta_{srr} * SB_{y-minage})}
#'
#' The Beverton-Holt and Ricker curves above are the deterministic mean; realized
#' recruitment applies the annual log deviation, \eqn{R_y \cdot exp(R_{dev,y})}, as in the
#' mean form. For numerical stability the Ricker \eqn{\beta_{srr}} is estimated on a scale
#' divided by 1,000,000, so the fitted \code{beta} is 1e6 times the density-dependence
#' coefficient in the equation above; \code{Bmsy_lim} (\eqn{\approx 1/\beta_{srr}}) carries
#' the same scaling.
#'
#' When \code{srr_pred_fun > 0} and \code{srr_fun = 0} recruitment in the hindcast is estimated as in \code{srr_fun = 0} \deqn{R_y = exp(R0 + R_{dev,y})}, but an additional stock recruitment relationship defined by \code{srr_pred_fun} is estimated between \code{srr_hat_styr} and \code{srr_hat_endyr} and treated as an additional penalty. The stock recruitment relationship defined by \code{srr_pred_fun} is then used in the projection.
#'
#'
#' @return A \code{list} containing the stock recruitment relationship settings
#' @examples
#' # Beverton-Holt fitted to the hindcast, with a prior on steepness.
#' build_srr(srr_fun = "BevertonHolt", srr_pred_fun = "BevertonHolt",
#'           srr_est_mode = "Estimated",
#'           srr_prior = 0.8, srr_prior_sd = 0.15)
#'
#' # Mean recruitment: no stock-recruit relationship fitted.
#' build_srr(srr_fun = "mean")
#' @export
#'
build_srr <- function(srr_fun = 0,  #srr_model
                      srr_pred_fun = srr_fun, #srr_forecast_model
                      proj_mean_rec = TRUE,
                      srr_mse_switchyr = NULL,
                      srr_hat_styr = NULL,
                      srr_hat_endyr = NULL,
                      srr_est_mode = 1,
                      srr_prior = 4,
                      srr_prior_sd = 1,
                      srr_alpha_init = NULL,
                      srr_beta_init = NULL,
                      srr_indices = NA,
                      Bmsy_lim = NA,
                      linkages = NULL){

  srr_fun      <- .coerce_srr_fun(srr_fun,      "srr_fun")
  srr_pred_fun <- .coerce_srr_fun(srr_pred_fun, "srr_pred_fun")
  srr_est_mode <- .map_switch(srr_est_mode, srr_est_mode_map, "srr_est_mode")

  # Set pred/RP/penalty to same as SR curve if SR fun > 0
  if(srr_fun > 0){
    srr_pred_fun = srr_fun
  }

  if(!srr_pred_fun %in% c(4,5, "Ricker")){
    Bmsy_lim = -999
  }

  # For a Beverton-Holt curve, srr_est_mode 2 and 3 put the prior on steepness,
  # which must lie in (0, 1). The beta prior converts (mean, sd) to shape
  # parameters by moments,
  #   a = ((1 - mu)/sd^2 - 1/mu) * mu^2,   b = a * (1/mu - 1),
  # which are positive only when sd^2 < mu * (1 - mu). Outside that range the
  # prior is not a density; TMB's dbeta returns a finite value for negative
  # shapes rather than NaN, so it would otherwise pass unnoticed.
  if (srr_est_mode %in% c(2, 3) && srr_pred_fun %in% c(2, 3)) {
    bad_h <- !is.na(srr_prior) & (srr_prior <= 0 | srr_prior >= 1)
    if (any(bad_h)) {
      stop("For a Beverton-Holt curve, `srr_est_mode = ", srr_est_mode,
           "` puts the prior on steepness, so `srr_prior` must be in (0, 1). ",
           "Got: ", paste(srr_prior[bad_h], collapse = ", "),
           ".\n  (For Ricker, `srr_prior` is a prior on alpha instead.)",
           call. = FALSE)
    }
    if (srr_est_mode == 3) {
      max_sd <- sqrt(srr_prior * (1 - srr_prior))
      bad_sd <- !is.na(srr_prior) & !is.na(srr_prior_sd) & (srr_prior_sd >= max_sd)
      if (any(bad_sd)) {
        stop("`srr_est_mode = 3` (beta prior on steepness) needs ",
             "`srr_prior_sd` < sqrt(srr_prior * (1 - srr_prior)) = ",
             paste(signif(max_sd[bad_sd], 4), collapse = ", "),
             ", otherwise the beta shape parameters are negative and the prior ",
             "is silently meaningless. Got srr_prior_sd = ",
             paste(srr_prior_sd[bad_sd], collapse = ", "),
             ".\n  Note the default srr_prior_sd = 1 is never valid here; ",
             "supply a smaller value.", call. = FALSE)
      }
    }
  }

  # With srr_est_mode = 0 the alpha parameter is fixed AT srr_prior, so here
  # srr_prior is an alpha and not a steepness -- the opposite of modes 2 and 3
  # above. A Beverton-Holt alpha below 1/SPR0 puts the curve under the
  # replacement line, and a value in (0, 1) is nearly always a steepness passed
  # to the wrong mode. Warn rather than stop: 1/SPR0 is not known until the
  # model is built, and a genuinely small alpha is legitimate for a stock with
  # a large SPR0.
  if (isTRUE(srr_est_mode == 0) && srr_pred_fun %in% c(2, 3)) {
    looks_like_h <- !is.na(srr_prior) & srr_prior > 0 & srr_prior < 1
    if (any(looks_like_h)) {
      warning("`srr_est_mode = 0` fixes alpha at `srr_prior`, so `srr_prior` ",
              "is an alpha here, not a steepness. Got ",
              paste(srr_prior[looks_like_h], collapse = ", "),
              ", which is below the replacement line 1/SPR0 for most stocks ",
              "and would give a steepness under 0.2.\n  If you meant a ",
              "steepness h, either use `srr_est_mode = 2` / `3` (which do put ",
              "the prior on steepness) or convert it: ",
              "alpha = 4h / (SPR0 * (1 - h)).", call. = FALSE)
    }
  }

  linkages <- .validate_recruitment_linkages(linkages, srr_pred_fun)

  # `srr_indices` is soft-deprecated in favour of `linkages = list(R0
  # = ..., alpha = ..., beta = ...)`. NA is "not supplied".
  if (!(length(srr_indices) == 1L && is.na(srr_indices))) {
    .warn_srr_indices_deprecation()
  }

  list(srr_fun = srr_fun,
       srr_pred_fun = srr_pred_fun,
       proj_mean_rec = proj_mean_rec,
       srr_mse_switchyr = srr_mse_switchyr,
       srr_hat_styr = srr_hat_styr,
       srr_hat_endyr = srr_hat_endyr,
       srr_est_mode = srr_est_mode,
       srr_prior = srr_prior,
       srr_prior_sd = srr_prior_sd,
       srr_alpha_init = srr_alpha_init,
       srr_beta_init = srr_beta_init,
       srr_indices = srr_indices,
       Bmsy_lim = Bmsy_lim,
       linkages = linkages
  )
}


#' Is `srr_prior` an alpha, or a steepness?
#'
#' `srr_prior` is a prior on **steepness** where the model consumes it as one:
#' the lognormal (`srr_est_mode` 2) and beta (`srr_est_mode` 3) priors on a
#' Beverton-Holt curve (`srr_pred_fun` 2 or 3). Everywhere else -- Ricker at any
#' `srr_est_mode`, and `srr_est_mode` 0 ("fix alpha to prior mean") or 1
#' ("estimate") for any curve -- it is an alpha, and so is a valid starting
#' value for `rec_pars[, "Alpha"]`.
#'
#' `build_params()` and `fit_mod()` both seed alpha and share this rule.
#'
#' @param data_list A `data_list` carrying `srr_est_mode` / `srr_pred_fun`.
#' @return `TRUE` when `srr_prior` may be used as an alpha starting value.
#' @keywords internal
#' @noRd
.srr_prior_is_alpha <- function(data_list) {
  as_int <- function(x) {
    if (is.null(x)) return(NA_integer_)
    suppressWarnings(as.integer(x)[1])
  }
  steepness_case <-
    isTRUE(as_int(data_list$srr_est_mode) %in% c(2L, 3L)) &&
    isTRUE(as_int(data_list$srr_pred_fun) %in% c(2L, 3L))
  !steepness_case
}


#' String<->integer mapping for `srr_fun` / `srr_pred_fun` in
#' [build_srr()]
#'
#' Either form is accepted; the canonical integer code is what the
#' TMB template ultimately consumes. Only the structural codes (0,
#' 2, 4) get string aliases. The env-driven codes (1, 3,
#' 5) still work with a soft-deprecation warning -- their structural
#' part is identical to 0 / 2 / 4 respectively, and the env effect
#' is expressed via the `linkages` argument to [build_srr()].
#'
#' @keywords internal
.SRR_FUNS <- c(
  mean         = 0L,
  BevertonHolt = 2L,
  Ricker       = 4L
)


#' Deprecated env-driven `srr_fun` / `srr_pred_fun` integer codes.
#' @keywords internal
#' @noRd
.SRR_DEPRECATED_FUNS <- c(1L, 3L, 5L)


#' Coerce an `srr_fun` / `srr_pred_fun` value to canonical integer.
#'
#' Accepts either a string from [.SRR_FUNS] (length-1) or a length-1
#' integer in 0..5. Integer codes 1, 3, 5 emit a soft-deprecation
#' warning pointing users at the linkage table.
#'
#' @keywords internal
#' @noRd
.coerce_srr_fun <- function(x, what) {
  .coerce_switch_arg(
    x, map = .SRR_FUNS, what = what,
    deprecated = .SRR_DEPRECATED_FUNS,
    warn_fn = function(int) .warn_srr_fun_deprecation(int, what),
    length_exact_one = TRUE,
    legacy_note = ", plus 1/3/5 for legacy env modes")
}


#' @keywords internal
#' @noRd
.warn_srr_fun_deprecation <- function(int, what) {
  warning(
    sprintf("%s = %d is soft-deprecated: the structural part of ", what, int),
    "mode 1 / 3 / 5 is identical to mode 0 / 2 / 4 respectively, ",
    "and the environmental effect formerly driven by ",
    "`srr_indices` is now better expressed through the linkages ",
    "argument to build_srr():\n\n",
    "  build_srr(srr_fun = ", switch(as.character(int),
                                     "1" = 0, "3" = 2, "5" = 4), ",\n",
    "            linkages = list(",
    switch(as.character(int), "1" = "R0", "alpha"),
    " = linkage_spec(formula = ~ <env_col>)))\n\n",
    "See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.warn_srr_indices_deprecation <- function() {
  warning(
    "`srr_indices` is deprecated. Environmental ",
    "effects are now expressed through the linkages argument ",
    "to build_srr():\n\n",
    "  build_srr(srr_fun = ...,\n",
    "            linkages = list(R0 = linkage_spec(\n",
    "              formula = ~ <env_col>)))\n\n",
    " See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' Allowed recruitment-parameter names for `linkages` in [build_srr()]
#'
#' Natural-scale names of the underlying recruitment parameters
#' that the linkage system can address. Linkages on `R0` are
#' meaningful for any `srr_fun` (the offset is added to the log of
#' equilibrium / mean recruitment when the default log link is used);
#' linkages on `alpha` and `beta` only do work when the chosen
#' `srr_fun` actually uses alpha / beta (Beverton-Holt, Ricker).
#'
#' @keywords internal
RECRUITMENT_LINKAGE_PARAMS <- c("R0", "alpha", "beta")


#' Map recruitment linkage param names to columns of `rec_pars`.
#' @keywords internal
#' @noRd
.REC_PARAM_TO_INDEX <- c(R0 = 1L, alpha = 2L, beta = 3L)


#' Validate and canonicalize the `linkages` argument of [build_srr()]
#'
#' Returns either `NULL` (no linkages) or a named list of
#' `Rceattle_linkage_spec` objects (or lists thereof) with `param`
#' filled in from the list keys. Errors loudly on invalid param
#' names so the user catches typos at build time. Warns when a
#' linkage references a parameter that the chosen `srr_fun` does
#' not consume (e.g. `alpha` with the mean-only `srr_fun = 0`).
#'
#' @keywords internal
#' @noRd
.validate_recruitment_linkages <- function(linkages, srr_fun) {
  linkages <- .validate_process_linkages(
    linkages, RECRUITMENT_LINKAGE_PARAMS, "recruitment"
  )
  if (is.null(linkages)) return(NULL)
  # Soft consistency check: a BH/Ricker SRR uses alpha and
  # beta; the mean-only srr_fun (0) uses neither.
  uses_alpha_beta <- srr_fun %in% c(2L, 3L, 4L, 5L)
  if (!uses_alpha_beta) {
    flagged <- intersect(names(linkages), c("alpha", "beta"))
    if (length(flagged) > 0) {
      warning("linkages$", paste(flagged, collapse = " / "),
              " is supplied but srr_pred_fun = ", srr_fun, " does not ",
              "use alpha / beta; the offset will be retained on the ",
              "object but will not affect recruitment. Use `srr_pred_fun` ",
              "in 'Ricker' or 'BevertonHolt' for an SRR that consumes these ",
              "parameters. Note that `srr_pred_fun = srr_fun`, if not supplied", call. = FALSE)
    }
  }
  linkages
}
