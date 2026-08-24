#' Specify the stock-recruit relationship (SRR) for Rceattle
#'
#' @param srr_fun Stock-recruit function used in the hindcast estimation (see the list below). Default = 0
#' @param srr_pred_fun Stock-recruit function used for projection, reference points, and penalties (see below). When \code{srr_fun == 0}, the stock-recruit curve is added as a penalty on the annually estimated hindcast recruitment (following AMAK and Jim Ianelli's pollock model). If \code{srr_fun > 0}, then \code{srr_pred_fun = srr_fun} and no extra penalty is added.
#' @param proj_mean_rec Project the model at mean recruitment (the average R of the hindcast) rather than from the stock-recruit relationship and recruitment deviations. \code{TRUE} / 1 = mean recruitment (default); \code{FALSE} / 0 = SRR(omega, srr_devs). A DSEM built with \code{estimate_projection = TRUE} projects through the SEM and overrides this.
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


#' Coerce a string/integer switch argument to its canonical integer code.
#'
#' Shared validation for the integer-returning `build_*()` switch arguments
#' (`srr_fun`, `M1_model`). Accepts a canonical string key from `map` or a
#' legacy integer code (a value of `map`, or a `deprecated` code); errors on
#' anything else and calls `warn_fn(int)` when a deprecated code is supplied.
#'
#' @param x user input (length-1+ character or numeric).
#' @param map named integer vector mapping canonical string -> integer code.
#' @param what argument name, used in messages.
#' @param deprecated integer codes that are accepted but soft-deprecated.
#' @param warn_fn `function(int)` emitting the deprecation warning, or `NULL`.
#' @param length_exact_one require length 1 (`srr_fun`) vs length >= 1 (`M1_model`).
#' @param legacy_note text appended inside the integer-range error message.
#' @keywords internal
#' @noRd
.coerce_switch_arg <- function(x, map, what,
                               deprecated = integer(0),
                               warn_fn = NULL,
                               length_exact_one = FALSE,
                               legacy_note = "") {
  if (length_exact_one) {
    if (length(x) != 1L) {
      stop(sprintf("`%s` must be length 1", what), call. = FALSE)
    }
  } else if (length(x) == 0L) {
    stop(sprintf("`%s` must have length >= 1", what), call. = FALSE)
  }
  if (is.numeric(x)) {
    int <- as.integer(x)
    allowed <- c(unname(map), deprecated)
    if (anyNA(int) || any(!int %in% allowed)) {
      stop(sprintf(
        "integer `%s` must be one of: %s (= %s%s)",
        what,
        paste(map, collapse = ", "),
        paste(names(map), collapse = "/"),
        legacy_note), call. = FALSE)
    }
    if (length(deprecated) > 0L && any(int %in% deprecated) && !is.null(warn_fn)) {
      warn_fn(int)
    }
    return(int)
  }
  if (is.character(x)) {
    bad <- setdiff(x, names(map))
    if (length(bad) > 0L) {
      stop(sprintf(
        "unknown `%s` value(s): %s; allowed: %s",
        what,
        paste(unique(bad), collapse = ", "),
        paste(names(map), collapse = ", ")), call. = FALSE)
    }
    return(unname(map[x]))
  }
  stop(sprintf("`%s` must be a string or integer; got %s",
               what, class(x)[1]), call. = FALSE)
}


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



#' Allowed M-parameter names for `linkages` in [build_M1()]
#'
#' Natural-scale names of the underlying natural-mortality
#' parameters that the linkage system can address. Currently just
#' `M1` -- with the default log link the offset is added to M1 on the
#' log scale (applied across all ages unless the linkage row
#' pins a specific `age_bin`).
#'
#' Note: predation mortality `M2` is a derived quantity in CEATTLE
#' (a function of predator abundance, suitability, and ration), not
#' a parameter. There is no `M2` linkage target; environmental
#' effects on predation are mediated upstream via recruitment,
#' growth, suitability, or ration inputs.
#'
#' @keywords internal
M_LINKAGE_PARAMS <- c("M1")


#' String<->integer mapping for `M1_model` in [build_M1()]
#'
#' Either form is accepted by [build_M1()]; the canonical integer
#' code is what the TMB template ultimately consumes. The
#' env-driven integer codes 4 and 5 still work with a deprecation
#' warning -- their structural part is identical to 1 and 2
#' respectively, and the env effect is expressed via the
#' \code{linkages} argument to [build_M1()] (see
#' `vignette("environmental-linkages-and-priors")`). No string alias is offered
#' for 4 or 5 to discourage their use in new code.
#'
#' @keywords internal
.M1_MODELS <- c(
  fixed             = 0L,
  sex_age_invariant = 1L,
  sex_specific      = 2L,
  sex_age_specific  = 3L
)


#' Deprecated env-driven `M1_model` integer codes.
#' @keywords internal
#' @noRd
.M1_DEPRECATED_MODELS <- c(4L, 5L)


#' String<->integer mapping for `M1_re` in [build_M1()]
#'
#' Either form is accepted by [build_M1()]; the canonical integer
#' code is what the TMB template ultimately consumes. The
#' env-driven options (4, 5) on M1_model can also be expressed
#' through the `linkages` argument to [build_M1()]; the integer
#' codes continue to work for backwards compatibility.
#'
#' @keywords internal
.M1_RES <- c(
  none         = 0L,
  iid_age      = 1L,
  iid_year     = 2L,
  iid_age_year = 3L,
  ar1_age      = 4L,
  ar1_year     = 5L,
  ar1_age_year = 6L
)


#' Coerce an `M1_model` or `M1_re` value to canonical integer code(s).
#'
#' Accepts either a string from the corresponding map (length-1 or
#' length-N for per-species values) or a length-1/length-N integer
#' vector that is already in the allowed set. Errors loudly on
#' anything else.
#'
#' For `M1_model` only, the env-driven integer codes 4
#' and 5 (controlled by `M1_indices`) are accepted for backwards
#' compatibility but emit a soft-deprecation warning pointing users
#' at the linkage table -- see [build_M1()].
#'
#' @keywords internal
#' @noRd
.coerce_M1_arg <- function(x, map, what) {
  deprecated <- if (identical(what, "M1_model")) .M1_DEPRECATED_MODELS else integer(0)
  .coerce_switch_arg(
    x, map = map, what = what,
    deprecated = deprecated,
    warn_fn = function(int) .warn_M1_model_deprecation(int),
    length_exact_one = FALSE)
}


#' @keywords internal
#' @noRd
.warn_M1_model_deprecation <- function(int) {
  hits <- intersect(int, .M1_DEPRECATED_MODELS)
  warning(
    "M1_model %in% c(", paste(hits, collapse = ", "),
    ") is soft-deprecated: the structural part of mode 4 is ",
    "identical to mode 1, and mode 5 to mode 2 (or 1 in single-",
    "sex models). The environmental effect formerly driven by ",
    "`M1_indices` is now better expressed through the linkages ",
    "argument to build_M1():\n\n",
    "  build_M1(M1_model = c(1, 2, 1),\n",
    "           linkages = list(M1 = linkage_spec(\n",
    "             formula = ~ <env_col>, by = ~species,\n",
    "             species = <which species had the env effect>)))\n\n",
    "See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.warn_M1_indices_deprecation <- function() {
  warning(
    "`M1_indices` is soft-deprecated. The same environmental ",
    "effect is now better expressed through the linkages argument ",
    "to build_M1():\n\n",
    "  build_M1(M1_model = ...,\n",
    "           linkages = list(M1 = linkage_spec(\n",
    "             formula = ~ <env_col>, by = ~species)))\n\n",
    "Both paths add additively to log_M1 on the log scale, so do ",
    "NOT supply both for the same coefficient or you will ",
    "double-count. See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' Specify the residual natural mortality (M1) model for Rceattle
#'
#' @param M1_model Vector or scalar specifying the M1 structural fixed-
#'   effects model. Either an integer code or the equivalent string
#'   alias (both forms are accepted; the integer code is canonical):
#'   * `0` / `"fixed"` -- use the input `M1_base` (no estimation).
#'   * `1` / `"sex_age_invariant"` -- estimate one `M1_{spp}`.
#'   * `2` / `"sex_specific"` -- estimate `M1_{spp, sex}`.
#'   * `3` / `"sex_age_specific"` -- estimate `M1_{spp, sex, age}`.
#'   * `4`, `5` -- soft-deprecated env-driven codes; use the
#'     `linkages` argument instead. See `vignette("environmental-linkages-and-priors")`.
#' @param M1_re Vector or scalar specifying the M1 random-effects
#'   model. Either an integer code or the equivalent string alias:
#'   `0` / `"none"`, `1` / `"iid_age"`, `2` / `"iid_year"`,
#'   `3` / `"iid_age_year"`, `4` / `"ar1_age"`, `5` / `"ar1_year"`,
#'   `6` / `"ar1_age_year"`.
#' @param updateM1 If using initial parameters, use M1 fixed effects
#'   from data (`M1_base`) instead. Default `FALSE`.
#' @param M1_use_prior Vector or scalar; if `TRUE`, apply the
#'   lognormal `M_prior` / `M_prior_sd` to `M1` directly.
#' @param M2_use_prior Vector or scalar; if `TRUE`, apply the
#'   lognormal prior to `M1 + M2` in multi-species models.
#' @param M_prior Mean (natural-scale) of the lognormal prior on M.
#' @param M_prior_sd SD (log-scale) of the lognormal prior on M.
#' @param M1_indices Soft-deprecated. Vector of column indices into
#'   `env_data` (excluding `Year`) for environmentally linked M1 when
#'   `M1_model %in% c(4, 5)`. Use the `linkages` argument instead;
#'   see \code{vignette("environmental-linkages-and-priors")}.
#' @param linkages Optional named list of [linkage_spec()] objects
#'   keyed by M parameter name (currently the only valid key is
#'   `"M1"`). Each spec describes how `M1` depends on
#'   environmental covariates and on stratifying factors (species,
#'   sex, age). The offset enters additively (on the log scale)
#'   inside the `M1_at_age` compute. A row's `age_bin == NA`
#'   broadcasts the offset across ages; specific values pin it to
#'   that age slice.
#'
#' @return A list of switches for defining the M1 model.
#' @export
#'
#' @examples
#' \donttest{
#' # Sex/age-invariant M with a temperature linkage on M1
#' build_M1(
#'   M1_model = "sex_age_invariant",
#'   linkages = list(
#'     M1 = linkage_spec(
#'       formula = ~ temp,
#'       by      = ~ species,
#'       priors  = list(temp = normal(0, 0.5))
#'     )
#'   )
#' )
#' }
build_M1 <- function(M1_model = 0,
                     M1_re = 0,
                     updateM1 = FALSE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE,
                     M_prior = 0.40,
                     M_prior_sd = 0.35,
                     M1_indices = NA,
                     linkages = NULL){
  M1_model <- .coerce_M1_arg(M1_model, .M1_MODELS, "M1_model")
  M1_re    <- .coerce_M1_arg(M1_re,    .M1_RES,    "M1_re")
  linkages <- .validate_M_linkages(linkages)
  # `M1_indices` is soft-deprecated in favour of `linkages = list(log_M
  # = ...)`. NA / NA_real_ / NA_integer_ are all "not supplied". Any
  # other value triggers the deprecation note.
  if (!(length(M1_indices) == 1L && is.na(M1_indices))) {
    .warn_M1_indices_deprecation()
  }
  list(
    M1_model     = M1_model,
    M1_re        = M1_re,
    updateM1     = updateM1,
    M1_use_prior = M1_use_prior,
    M2_use_prior = M2_use_prior,
    M_prior      = M_prior,
    M_prior_sd   = M_prior_sd,
    M1_indices   = M1_indices,
    linkages     = linkages
  )
}


#' Validate and canonicalize the `linkages` argument of [build_M1()]
#'
#' Returns either `NULL` (no linkages) or a named list of
#' `Rceattle_linkage_spec` objects (or lists thereof) with `param`
#' filled in from the list keys. Errors loudly on invalid param
#' names so the user catches typos at build time.
#'
#' @keywords internal
#' @noRd
.validate_M_linkages <- function(linkages) {
  .validate_process_linkages(linkages, M_LINKAGE_PARAMS, "M")
}



#' Allowed growth functions for `fun` in [build_growth()]
#'
#' Currently:
#'   * `"empirical"`: use the empirical weight-at-age input directly,
#'     no parameters estimated.
#'   * `"vonBertalanffy"`: sex-specific von Bertalanffy (parameters
#'     `K`, `L1`, `Linf`).
#'   * `"Richards"`: sex-specific Richards (adds `m`).
#'
#' @keywords internal
GROWTH_FUNS <- c("empirical", "vonBertalanffy", "Richards")


#' Internal: integer codes consumed by the TMB template
#' @keywords internal
#' @noRd
.GROWTH_FUN_TO_INT <- c(empirical = 0L, vonBertalanffy = 1L, Richards = 2L)


#' Internal: plus-group SD-at-age treatment codes consumed by the TMB template
#'
#' `"WHAM"` (1) pins the oldest age class's SD-at-age to the upper anchor
#' `exp(sd_Linf)` (the WHAM SDAA convention). `"SS3"` (2) instead interpolates
#' it by length like any interior age. Only affects estimated growth
#' (`growth_model > 0`); the default is WHAM.
#' @keywords internal
#' @noRd
.GROWTH_SD_STYLE <- c(WHAM = 1L, SS3 = 2L)


#' Allowed growth-parameter names for `linkages` in [build_growth()]
#'
#' Natural-scale names of the underlying growth-function parameters.
#' Von Bertalanffy uses `K`, `L1`, `Linf`; Richards adds `m`.
#' `sd_L1` / `sd_Linf` are the standard deviations of length-at-age
#' anchored at `L1` and `Linf` (the SD-at-age interpolation endpoints).
#' Only intercept-only specs (`~ 1`) are honored
#' on the SD endpoints -- they thread through `init` / `bounds` /
#' `priors` onto the growth SD-at-age but do not vary by year.
#' The empirical weight-at-age model admits no linkages.
#'
#' @keywords internal
GROWTH_LINKAGE_PARAMS <- c("K", "L1", "Linf", "m", "sd_L1", "sd_Linf")


#' Mean-growth subset of [GROWTH_LINKAGE_PARAMS]
#'
#' Names that index into `log_growth_pars` / `growth_parameters` /
#' `growth_linkage_offset` along their last dim. The SD endpoints live
#' on a separate parameter (`growth_log_sd`) and are intentionally
#' excluded here.
#'
#' @keywords internal
.GROWTH_MEAN_PARAMS <- c("K", "L1", "Linf", "m")


#' Map mean-growth linkage param names to slices along the third dim of
#' `log_growth_pars` (`[nspp, nsex, n_growth_pars]`).
#' @keywords internal
#' @noRd
.GROWTH_PARAM_TO_INDEX <- c(K = 1L, L1 = 2L, Linf = 3L, m = 4L)


#' Map growth SD linkage param names to slices along the third dim of
#' `growth_log_sd` (`[nspp, nsex, 2]`).
#' @keywords internal
#' @noRd
.GROWTH_SD_PARAM_TO_INDEX <- c(sd_L1 = 1L, sd_Linf = 2L)


#' Specify the growth model for Rceattle
#'
#' @param fun Growth function. Either a string ([GROWTH_FUNS]:
#'   `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
#'   equivalent integer code (`0`, `1`, `2`). The canonical string form
#'   is stored on the returned object.
#' @param growth_age_L1 Von Bertalanffy / Richards anchor age (the age at
#'   which mean length equals `L1`). Matches SS3's `Growth_Age_for_L1`
#'   control input. Scalar (recycled across species) or a length-`nspp`
#'   vector for per-species values. Default `NA` inherits
#'   `data_list$growth_age_L1` if supplied (e.g. from the SS3 converter),
#'   otherwise falls back to `max(0.5, minage[sp])` so `minage >= 1`
#'   models stay backwards-compatible and `minage = 0` models pick up an
#'   SS3-consistent half-year anchor.
#' @param sd_plus_group How the oldest age class's SD-at-age is treated (only
#'   affects estimated growth, `fun != "empirical"`). `"WHAM"` pins the
#'   plus-group SD to the upper anchor `exp(sd_Linf)` (the WHAM SDAA
#'   convention); `"SS3"` instead interpolates it by length like any interior
#'   age. Accepts a string or the integer code (`1`/`2`), scalar or a
#'   length-`nspp` vector. Default `NA` inherits `data_list$growth_sd_style`
#'   if present (so a refit keeps the original choice), otherwise `"WHAM"`.
#' @param linkages Optional named list of [linkage_spec()] objects
#'   keyed by parameter name (must be one of [GROWTH_LINKAGE_PARAMS]).
#'   The mean-growth keys (`K`, `L1`, `Linf`, `m`)
#'   accept arbitrary one-sided formulas and make that growth parameter
#'   year-varying (a per-year offset around its mean). The SD-endpoint keys
#'   (`sd_L1`, `sd_Linf`) only honor intercept-bearing
#'   formulas (typically `~ 1`) -- they thread `init`, `bounds`, and
#'   `priors` onto the growth SD-at-age, giving
#'   the SDs the same prior/fix/initial-value contract as the mean
#'   parameters. Slope rows on SD specs raise a warning and have no
#'   effect; slope-only formulas (`~ 0 + temp`) error.
#'
#' @return A list of switches defining the growth model.
#' @export
#'
#' @examples
#' \donttest{
#' # Sex-specific von Bertalanffy with temperature on K, by species + sex
#' build_growth(
#'   fun = "vonBertalanffy",   # or fun = 1
#'   linkages = list(
#'     K = linkage_spec(
#'       formula = ~ temp,
#'       by      = ~ species + sex,
#'       priors  = list(temp = normal(0, 1))
#'     )
#'   )
#' )
#' }
build_growth <- function(fun = "empirical",
                         growth_age_L1 = NA,
                         sd_plus_group = NA,
                         linkages = NULL) {
  fun <- .coerce_growth_fun(fun)
  # NA (the default) means "inherit": fit_mod() resolves it from
  # data_list$growth_sd_style if present, else the WHAM fallback -- exactly like
  # growth_age_L1. This keeps a refit that rebuilds growth via build_growth(fun=)
  # from silently resetting an SS3 model to WHAM.
  sd_style_int <- if (all(is.na(sd_plus_group))) {
    NA_integer_
  } else {
    .coerce_switch_arg(sd_plus_group, .GROWTH_SD_STYLE, "sd_plus_group")
  }
  linkages <- .validate_growth_linkages(linkages, fun)
  list(
    fun            = fun,
    linkages       = linkages,
    # Internal integer code consumed by the TMB template until the
    # linkage-driven path replaces it. Vectorized so per-species
    # growth functions (e.g. fun = c("vonBertalanffy", "Richards")) each
    # propagate their own code downstream.
    growth_model   = unname(.GROWTH_FUN_TO_INT[fun]),
    # Plus-group SD-at-age treatment. Like `fun`/`growth_model`, the canonical
    # string is kept under the argument name (for save_config round-trip) and
    # the int code (1 = WHAM, 2 = SS3) feeds the TMB template. Per-species.
    # NA_integer_ = inherit (resolved in fit_mod); the string field is NA too and
    # is dropped from a saved config, so unspecified models re-derive on load.
    sd_plus_group   = names(.GROWTH_SD_STYLE)[match(sd_style_int, .GROWTH_SD_STYLE)],
    growth_sd_style = sd_style_int,
    # VB anchor age (= age at which `l1` is the length). Matches SS3's
    # `Growth_Age_for_L1` ctl input. Scalar or length-nspp; pass NA
    # (default) to inherit max(0.5, minage[sp]) downstream so old
    # configurations stay unchanged at minage >= 1 and minage = 0
    # models pick up an SS3-consistent half-year anchor.
    growth_age_L1  = growth_age_L1
  )
}


#' Coerce a `fun` argument (string or integer) to canonical strings
#'
#' Accepts either a [GROWTH_FUNS] string or its integer code from
#' `.GROWTH_FUN_TO_INT`. Length-N vectors are allowed for per-species
#' specifications. Returns a character vector of canonical names.
#'
#' @keywords internal
#' @noRd
.coerce_growth_fun <- function(fun) {
  if (length(fun) == 0L) {
    stop("`fun` must have length >= 1", call. = FALSE)
  }
  if (is.numeric(fun)) {
    int <- as.integer(fun)
    hit <- match(int, .GROWTH_FUN_TO_INT)
    if (anyNA(hit)) {
      stop("integer `fun` must be one of: ",
           paste(.GROWTH_FUN_TO_INT, collapse = ", "),
           " (= ",
           paste(names(.GROWTH_FUN_TO_INT), collapse = "/"),
           ")", call. = FALSE)
    }
    return(names(.GROWTH_FUN_TO_INT)[hit])
  }
  if (is.character(fun)) {
    bad <- setdiff(fun, GROWTH_FUNS)
    if (length(bad) > 0L) {
      stop("unknown growth fun(s): ",
           paste(unique(bad), collapse = ", "),
           "; allowed: ", paste(GROWTH_FUNS, collapse = ", "),
           call. = FALSE)
    }
    return(fun)
  }
  stop("`fun` must be a string or integer; got ", class(fun)[1],
       call. = FALSE)
}


#' Validate and canonicalize the `linkages` argument of [build_growth()]
#'
#' Returns either `NULL` (no linkages) or a named list of fully-typed
#' `Rceattle_linkage_spec` objects with `param` filled in from the list
#' keys. Errors loudly on invalid param names so the user catches typos
#' at build time rather than at materialization time.
#'
#' @keywords internal
#' @noRd
.validate_growth_linkages <- function(linkages, fun) {
  linkages <- .validate_process_linkages(
    linkages, GROWTH_LINKAGE_PARAMS, "growth"
  )
  if (is.null(linkages)) return(NULL)
  if (any(fun == "empirical")) {
    warning("at least one species uses fun = 'empirical', which does ",
            "not consume linkages; the supplied specs will be ",
            "retained on the object but ignored at fit time for those ",
            "species.", call. = FALSE)
  }
  if (!all(fun == "Richards") && "m" %in% names(linkages)) {
    stop("linkages$m is only valid when every species uses ",
         "fun = 'Richards'; von Bertalanffy has no `m` parameter.",
         call. = FALSE)
  }
  # SD endpoints (log_sd_L1, log_sd_Linf) plug into `growth_log_sd`,
  # which has no year dim. Only intercept-bearing formulas (~ 1 plus
  # any covariates that the user wants to anchor on the intercept) are
  # honored: the intercept init/prior/bounds threads onto
  # `growth_log_sd`, and slope rows are silently dropped at the TMB
  # encoder (no offset path exists for SD). Block slope-only formulas
  # (`~ 0 + temp`) outright -- they'd leave the SD pinned at its
  # `build_params()` default with no escape.
  sd_keys <- intersect(names(linkages), names(.GROWTH_SD_PARAM_TO_INDEX))
  for (nm in sd_keys) {
    specs <- linkages[[nm]]
    if (inherits(specs, "Rceattle_linkage_spec")) specs <- list(specs)
    for (sp_i in specs) {
      f <- sp_i$formula
      rhs <- attr(stats::terms(f), "intercept")
      has_slope <- length(attr(stats::terms(f), "term.labels")) > 0L
      if (!isTRUE(rhs == 1L)) {
        stop(sprintf(
          "linkages$%s must use an intercept-bearing formula (e.g. ~ 1); ",
          nm),
          "year-varying offsets on growth SD are not wired through ",
          "growth.hpp. Drop the `0 +` from the formula, or move the ",
          "env effect onto K / L1 / Linf / m.",
          call. = FALSE)
      }
      if (has_slope) {
        warning(sprintf(
          "linkages$%s carries slope terms but growth SD has no year ",
          "dim in growth.hpp; the intercept init/prior/bounds will be ",
          "honored on `growth_log_sd`, but slope coefficients have no ",
          "effect on the SD.", nm), call. = FALSE)
      }
    }
  }
  linkages
}


#' @keywords internal
#' @noRd
# Base stratum a process uses when the user omits `by` in linkage_spec(): the
# fleet-indexed processes key by fleet, everything else by species. Composition
# is split -- the fleet DM weights (theta_comp / theta_caal) key by fleet, the
# per-predator diet weight (theta_diet) by species.
.default_stratum <- function(process_label, param) {
  if (process_label %in% c("q", "sel")) return(~ fleet)
  if (process_label == "comp") {
    return(if (identical(as.character(param), "theta_diet")) ~ species else ~ fleet)
  }
  ~ species                                    # recruitment, M, growth
}

# Fill an omitted `by` with the process's base stratum; an explicitly-set `by`
# (a formula, or NULL for a single shared coefficient) is left as the user gave it.
.resolve_auto_by <- function(spec, process_label, param) {
  if (isTRUE(spec$by_auto)) spec$by <- .default_stratum(process_label, param)
  spec
}

.stamp_param <- function(val, param, process_label = NULL) {
  stamp1 <- function(s) {
    s <- .set_linkage_param(s, param)
    if (!is.null(process_label)) s <- .resolve_auto_by(s, process_label, param)
    s
  }
  if (inherits(val, "Rceattle_linkage_spec")) return(stamp1(val))
  lapply(val, stamp1)
}


#' Catchability parameters that accept a linkage
#' @keywords internal
#' @noRd
Q_LINKAGE_PARAMS <- c("q")


#' @keywords internal
#' @noRd
.validate_q_linkages <- function(linkages) {
  .validate_process_linkages(linkages, Q_LINKAGE_PARAMS, "q")
}


#' Catchability specification
#'
#' @description
#' Carry environmental linkages on survey/index catchability `q`. The effect of
#' an `env_data` covariate is written as a formula and can carry priors, bounds,
#' and an estimation phase like any other linkage.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   catchability parameter. The only parameter is `q`. Coefficients are per
#'   fleet by default (`by = ~ fleet`); use the `fleet` argument of
#'   [linkage_spec()] to restrict a spec to particular fleets.
#'
#' @return A list of catchability settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # One temperature effect on q per fleet. `by = ~ fleet` is the default, and a
#' # fleet without an estimated survey catchability is not linkable, so at fit time
#' # restrict the linkage to the estimated-q fleets (here fleet 7 of BS2017SS).
#' build_catchability(linkages = list(
#'   q = linkage_spec(~ temp, fleet = 7)))
#'
#' # Restrict it to fleets 1 and 3, with a prior on the slope
#' build_catchability(linkages = list(
#'   q = linkage_spec(~ temp, fleet = c(1, 3),
#'                    priors = list(temp = prior_normal(0, 1)))))
#' }
#'
#' @export
build_catchability <- function(linkages = NULL) {
  list(linkages = .validate_q_linkages(linkages))
}


#' Selectivity parameters that accept a linkage
#'
#' Slot names shared across the parametric forms, plus DoubleNormal aliases.
#' @keywords internal
#' @noRd
SEL_LINKAGE_PARAMS <- c("slp_asc", "slp_desc", "inf_asc", "inf_desc", "coff",
                        "sigma_asc", "sigma_desc", "peak", "right_floor")


#' @keywords internal
#' @noRd
.validate_sel_linkages <- function(linkages) {
  .validate_process_linkages(linkages, SEL_LINKAGE_PARAMS, "sel")
}


#' Selectivity specification
#'
#' @description
#' Carries environmental linkages on selectivity parameters. The effect on a
#' parameter is written as a formula and composes additively with any
#' `Time_varying_sel` process error on the same fleet (the two are separate
#' mechanisms: a covariate effect versus a deviation).
#'
#' The parameter names are the shape parameters of the parametric selectivity
#' forms:
#' \describe{
#'   \item{`slp_asc`, `slp_desc`}{ascending / descending logistic slope (log
#'     scale); for a double-normal the ascending / descending width, aliased
#'     `sigma_asc` / `sigma_desc`.}
#'   \item{`inf_asc`, `inf_desc`}{ascending / descending inflection age/length
#'     (natural scale); for a double-normal the peak and the logit right-floor,
#'     aliased `peak` / `right_floor`.}
#'   \item{`coff`}{non-parametric selectivity-at-bin coefficients.}
#' }
#'
#' Every parameter accepts `link = "log"` (multiplicative on the natural
#' parameter) or `link = "identity"` (additive), like the other processes.
#'
#' **Priors on a selectivity parameter.** An intercept-only formula (`~ 1`) with
#' a `priors` entry places a prior on the selectivity parameter itself (no
#' year-to-year offset is added). Read the prior on the parameter's own scale:
#' the slopes (`slp_asc` / `slp_desc`) are on the log scale (use `lognormal()`),
#' the inflections (`inf_asc` / `inf_desc`) on the natural scale (use `normal()`).
#' See Examples for a normal prior on the ascending inflection. This mirrors the
#' prior-only [build_composition()] path.
#'
#' A selectivity prior targets one parameter, so in a two-sex model an
#' unstratified `~ 1` prior constrains sex 1 only -- use `by = ~ sex` for a
#' per-sex prior. An `init` on a selectivity intercept has no effect (the
#' starting value comes from the data), and a prior on the double-normal
#' `right_floor` is not supported.
#' For a fleet that mirrors another fleet's selectivity (shared
#' `Selectivity_index`), place the prior on the lead fleet so the shared
#' parameter block is not penalized more than once.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   selectivity parameter. Coefficients are per fleet by default
#'   (`by = ~ fleet`); use the `fleet` argument of [linkage_spec()] to restrict
#'   a spec to particular fleets.
#'
#' @return A list of selectivity settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # A cold-pool effect on the ascending inflection of a logistic fleet
#' build_selectivity(linkages = list(
#'   inf_asc = linkage_spec(~ cold_pool, by = ~ fleet)))
#'
#' # A normal prior on the ascending inflection (intercept-only formula)
#' build_selectivity(linkages = list(
#'   inf_asc = linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 3)))))
#' }
#'
#' @export
build_selectivity <- function(linkages = NULL) {
  list(linkages = .validate_sel_linkages(linkages))
}


#' @keywords internal
#' @noRd
COMP_LINKAGE_PARAMS <- c("theta_comp", "theta_caal", "theta_diet")


#' @keywords internal
#' @noRd
.validate_comp_linkages <- function(linkages) {
  .validate_process_linkages(linkages, COMP_LINKAGE_PARAMS, "comp")
}


#' Composition-weighting specification
#'
#' @description
#' Carries **priors** on the Dirichlet-multinomial composition-weighting
#' overdispersion. The DM weight (the "theta" that scales the effective sample
#' size) is otherwise an unpenalized free parameter; a linkage lets you put a
#' prior on it through the same grammar as every other parameter. The three
#' parameters target the three DM likelihoods:
#' \describe{
#'   \item{`theta_comp`}{age / length composition weight (`comp_weights`), per
#'     fleet.}
#'   \item{`theta_caal`}{conditional-age-at-length weight (`caal_weights`), per
#'     fleet.}
#'   \item{`theta_diet`}{diet composition weight (`diet_comp_weights`), per
#'     predator.}
#' }
#'
#' This is a **prior-only** process: the DM weight is a scalar, not a
#' year-varying quantity, so only intercept formulas (`~ 1`) with a `priors`
#' entry are meaningful. The prior is placed on the natural-scale DM weight
#' `theta = exp(weight)`, so a `gamma()` prior reads naturally. A linkage on a
#' fleet whose `Comp_distribution` (or `CAAL_distribution`) is not
#' `"DirichletMultinomial"` errors, since the weight has no effect there.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   `theta_comp` / `theta_caal` (per fleet by default, `by = ~ fleet`) or
#'   `theta_diet` (per predator by default, `by = ~ species`).
#'
#' @return A list of composition-weighting settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # A weak gamma prior on the DM overdispersion of every fleet's age comps
#' build_composition(linkages = list(
#'   theta_comp = linkage_spec(~ 1, by = ~ fleet,
#'                             priors = list(`(Intercept)` = gamma(2, 0.5)))))
#' }
#'
#' @export
build_composition <- function(linkages = NULL) {
  list(linkages = .validate_comp_linkages(linkages))
}


# Map a selectivity linkage param name to (array, 1-based slot). log_sel_slp
# and sel_inf are [2, fleet, sex]; sel_coff is [fleet, sex, bin].
.SEL_PARAM_TO_SLOT <- list(
  slp_asc     = list(arr = "log_sel_slp", slot = 1L),
  slp_desc    = list(arr = "log_sel_slp", slot = 2L),
  sigma_asc   = list(arr = "log_sel_slp", slot = 1L),
  sigma_desc  = list(arr = "log_sel_slp", slot = 2L),
  inf_asc     = list(arr = "sel_inf",     slot = 1L),
  inf_desc    = list(arr = "sel_inf",     slot = 2L),
  peak        = list(arr = "sel_inf",     slot = 1L),
  right_floor = list(arr = "sel_inf",     slot = 2L),
  coff        = list(arr = "sel_coff",    slot = NA_integer_)
)


# Selectivity forms whose consume site reads the linkage offset: every
# PARAMETRIC form. DoubleNormal reuses the slp/inf slots (peak/sigma/floor
# aliases) and LogisticPM's multiplicative deviates carry the offset inside
# their exp. The non-parametric forms are excluded on purpose -- see the
# `coff` note in .check_sel_linkage_support().
.SEL_LINKAGE_WIRED_FORMS <- c("Logistic", "DoubleLogistic", "DescendingLogistic",
                              "DoubleNormal", "LogisticPM")
.SEL_LINKAGE_WIRED_PARAMS <- c("slp_asc", "slp_desc", "inf_asc", "inf_desc",
                               "sigma_asc", "sigma_desc", "peak", "right_floor")


#' Reject selectivity linkages the template does not yet consume
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param fleet_control the fleet control table.
#' @return invisibly NULL; errors on an unsupported sel linkage.
#' @keywords internal
#' @noRd
.check_sel_linkage_support <- function(linkage_table, fleet_control) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  sel <- linkage_table[linkage_table$process == "sel", , drop = FALSE]
  if (nrow(sel) == 0L) return(invisible())

  bad_param <- setdiff(unique(sel$param), .SEL_LINKAGE_WIRED_PARAMS)
  if (length(bad_param) > 0) {
    extra <- if ("coff" %in% bad_param) paste0(
      "\n  `coff` (non-parametric selectivity) cannot carry a linkage: those ",
      "forms mean-centre their coefficients each year, so a per-year offset ",
      "applied across all bins cancels exactly. A meaningful effect would need ",
      "a per-bin covariate, which the formula grammar does not express.") else ""
    stop(sprintf(
      paste0("selectivity linkage parameter(s) not supported: %s.\n",
             "  Supported: %s.%s"),
      paste(bad_param, collapse = ", "),
      paste(.SEL_LINKAGE_WIRED_PARAMS, collapse = ", "), extra), call. = FALSE)
  }

  # Every fleet a sel row targets must use a wired selectivity form.
  flts <- unique(sel$fleet)
  flts <- flts[!is.na(flts)]
  if (length(flts) == 0L) flts <- seq_len(nrow(fleet_control))  # NA = all fleets
  forms <- as.character(fleet_control$Selectivity[flts])
  bad_flt <- flts[!forms %in% .SEL_LINKAGE_WIRED_FORMS]
  if (length(bad_flt) > 0) {
    stop(sprintf(
      paste0("selectivity linkage on fleet(s) %s whose form (%s) is not yet ",
             "wired for linkages.\n  Wired forms: %s."),
      paste(fleet_control$Fleet_name[bad_flt], collapse = ", "),
      paste(unique(as.character(fleet_control$Selectivity[bad_flt])),
            collapse = ", "),
      paste(.SEL_LINKAGE_WIRED_FORMS, collapse = ", ")), call. = FALSE)
  }

  # A PRIOR on a selectivity intercept re-targets the base parameter, whose scale
  # depends on the fleet's form and whose ownership depends on Selectivity_index.
  # Two cases the re-target cannot yet express correctly are rejected up front
  # (a covariate linkage on the same slot is fine -- only priors are affected).
  prior_rows <- sel[!is.na(sel$prior_family) & sel$prior_family != "none", ,
                    drop = FALSE]
  if (nrow(prior_rows) > 0L) {
    row_flt <- function(f) if (is.na(f)) 1L else as.integer(f)   # NA fleet = cell 1 (cpp default)

    # (a) DoubleNormal stores sel_inf(1) as logit(right_floor), so a natural-scale
    # prior on `inf_desc` / `right_floor` would be evaluated on the logit scale.
    # Reject until the logit transform is wired -- the ascending peak (`inf_asc`)
    # and the sigmas/slopes (log scale) are unaffected.
    dn <- prior_rows[prior_rows$param %in% c("inf_desc", "right_floor"), , drop = FALSE]
    dn_flt <- unique(vapply(dn$fleet, row_flt, integer(1)))
    dn_flt <- dn_flt[as.character(fleet_control$Selectivity[dn_flt]) == "DoubleNormal"]
    if (length(dn_flt) > 0L) {
      stop(sprintf(paste0(
        "prior on `inf_desc` / `right_floor` for DoubleNormal fleet(s) %s is not ",
        "supported: that slot holds logit(right_floor), so a natural-scale prior ",
        "would be applied on the logit scale. Prior the ascending peak / sigmas ",
        "instead."),
        paste(fleet_control$Fleet_name[dn_flt], collapse = ", ")), call. = FALSE)
    }

    # (b) Fleets that mirror another fleet's selectivity (Selectivity_index != own
    # Fleet_code) share one parameter block; a prior on the mirror double-counts
    # the block (cf. the shared-block penalty trap). Require the prior on the lead
    # fleet (Selectivity_index == Fleet_code).
    sidx <- fleet_control$Selectivity_index
    mir_flt <- unique(vapply(prior_rows$fleet, row_flt, integer(1)))
    mir_flt <- mir_flt[!is.na(sidx[mir_flt]) & sidx[mir_flt] != mir_flt]
    if (length(mir_flt) > 0L) {
      stop(sprintf(paste0(
        "selectivity prior on fleet(s) %s that mirror another fleet's ",
        "selectivity (Selectivity_index != Fleet_code): the shared block would be ",
        "penalized once per sharing fleet. Place the prior on the lead fleet ",
        "(the one whose Selectivity_index equals its Fleet_code)."),
        paste(fleet_control$Fleet_name[mir_flt], collapse = ", ")), call. = FALSE)
    }

    # (c) A prior on a limb the fleet's own curve never uses. Logistic reads only
    # the ascending slots, DescendingLogistic only the descending ones; the other
    # pair stays at its build default and never enters selectivity-at-age. The
    # prior would still be added to the objective -- a constant that shifts the
    # reported likelihood and moves with an unrelated default, while doing
    # nothing to the fit. Silently accepting it is how a reconciliation against
    # another model picks up an unexplained offset.
    used <- list(Logistic           = c("slp_asc", "inf_asc"),
                 DescendingLogistic = c("slp_desc", "inf_desc"))
    for (form in names(used)) {
      f_rows <- prior_rows[vapply(prior_rows$fleet, function(f)
        as.character(fleet_control$Selectivity[row_flt(f)]) == form,
        logical(1)), , drop = FALSE]
      unused <- f_rows[!f_rows$param %in% used[[form]], , drop = FALSE]
      if (nrow(unused) > 0L) {
        stop(sprintf(paste0(
          "selectivity prior on `%s` for %s fleet(s) %s: that %s curve does not ",
          "use those parameters, so the prior would add a constant to the ",
          "objective without affecting the fit. Prior %s instead, or drop the ",
          "fleet from this prior's fleet list."),
          paste(unique(unused$param), collapse = "`, `"), form,
          paste(unique(fleet_control$Fleet_name[
            vapply(unused$fleet, row_flt, integer(1))]), collapse = ", "),
          form, paste(used[[form]], collapse = " / ")), call. = FALSE)
      }
    }
  }
  invisible()
}


# Catchability forms that do NOT estimate q, so a q linkage must not attach:
# "Fixed" holds index_log_q at its input (a linkage would silently turn a fixed
# q time-varying, contrary to the assessor's Fixed setting), and
# "Analytical"/"AnalyticalArith" solve q from the data (index_log_q is mapped
# out, so a linkage targets a non-free parameter and does nothing). Both are
# rejected up front rather than quietly changing q. See build_map_catchability.
.Q_LINKAGE_UNESTIMATED_FORMS <- c("Fixed", "Analytical", "AnalyticalArith")

# A second, different reason to refuse a q linkage. "Environmental" and "AR1" DO
# estimate q, but they rebuild index_q from their own formula in the cpp,
# assigning over the value that carries q_linkage_offset rather than adding to
# it. The linkage is then accepted, REPORTed in q_linkage_offset as a live
# covariate effect, and never enters the likelihood -- beta_linkage is left free
# with an identically-zero gradient. Refuse it, and point at the linkage
# equivalent: these two forms exist to be replaced by exactly that grammar, so a
# user reaching for both wants the linkage on its own.
.Q_LINKAGE_SELFBUILT_FORMS <- c("Environmental", "AR1")


#' Reject q linkages on fleets whose catchability is not estimated
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param fleet_control the fleet control table.
#' @return invisibly NULL; errors on a q linkage targeting a fleet whose q is
#'   fixed or analytically solved.
#' @keywords internal
#' @noRd
# Informational one-time message: a catchability/selectivity linkage whose `by` was
# auto-filled (the user omitted it) and that names no fleets applies to EVERY fleet
# of that process -- one coefficient each. Tell the user so the per-fleet expansion
# (and any resulting eligibility error) is not a surprise; naming fleets via
# linkage_spec(fleet = ) restricts it.
.message_auto_fleet_linkages <- function(spec_groups) {
  labs <- c(q = "catchability", sel = "selectivity")
  for (proc in names(labs)) {
    specs <- spec_groups[[proc]]
    if (is.null(specs)) next
    for (nm in names(specs)) {
      val <- specs[[nm]]
      lst <- if (inherits(val, "Rceattle_linkage_spec")) list(val) else val
      for (s in lst) {
        if (isTRUE(s$by_auto) && is.null(s$fleet) && "fleet" %in% all.vars(s$by)) {
          message(sprintf(
            paste0("A %s linkage (`%s`) did not name fleets, so it applies to every ",
                   "eligible fleet -- one coefficient each. Pass ",
                   "linkage_spec(fleet = ) to restrict it to specific fleets."),
            labs[[proc]], nm))
        }
      }
    }
  }
  invisible()
}

.check_q_linkage_support <- function(linkage_table, fleet_control) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  q <- linkage_table[linkage_table$process == "q", , drop = FALSE]
  if (nrow(q) == 0L) return(invisible())

  flts <- unique(q$fleet)
  flts <- flts[!is.na(flts)]
  if (length(flts) == 0L) flts <- seq_len(nrow(fleet_control))  # NA = all fleets
  forms <- as.character(fleet_control$Catchability[flts])
  # A fleet does not estimate q if its Catchability holds q fixed / solves it from
  # the data (Fixed / Analytical), or is absent (NA) -- a fleet with no survey index
  # has no catchability to link. A linkage on any of these cannot be estimated.
  self <- flts[!is.na(forms) & forms %in% .Q_LINKAGE_SELFBUILT_FORMS]
  if (length(self) > 0) {
    stop(sprintf(
      paste0("catchability linkage on fleet(s) %s whose Catchability (%s) builds ",
             "q from its own formula: the cpp assigns index_q from that formula, ",
             "overwriting the linkage offset instead of adding to it, so the ",
             "linkage would be reported in q_linkage_offset but never enter the ",
             "likelihood.\n",
             "  Express the whole relationship as the linkage and set ",
             "Catchability to \"Estimated\": a covariate effect is ",
             "build_catchability(linkages = list(q = linkage_spec(~ covariate, ",
             "by = ~ fleet))), an AR1 deviation is linkage_spec(ar1(1 | Year), ",
             "by = ~ fleet); see vignette('environmental-linkages-and-priors')."),
      paste(fleet_control$Fleet_name[self], collapse = ", "),
      paste(unique(as.character(fleet_control$Catchability[self])), collapse = ", ")),
      call. = FALSE)
  }

  bad <- flts[is.na(forms) | forms %in% .Q_LINKAGE_UNESTIMATED_FORMS]
  if (length(bad) > 0) {
    stop(sprintf(
      paste0("catchability linkage on fleet(s) %s whose Catchability (%s) does ",
             "not estimate q: index_log_q is held fixed (Fixed), solved from the ",
             "data (Analytical), or absent (NA -- the fleet has no survey index), ",
             "so a linkage would turn a fixed q time-varying or have no effect.\n",
             "  Give the fleet an estimated survey catchability (\"Estimated\" / ",
             "\"Estimated-with-prior\") with index data, or restrict the linkage to ",
             "estimated-q fleets via linkage_spec(fleet = ...)."),
      paste(fleet_control$Fleet_name[bad], collapse = ", "),
      paste(unique(as.character(fleet_control$Catchability[bad])),
            collapse = ", ")), call. = FALSE)
  }
  invisible()
}


#' Reject composition-weighting (comp) linkages that cannot take effect
#'
#' @description
#' `comp` linkages are prior-only (the Dirichlet-multinomial weight is a scalar,
#' not a year-varying quantity), so only intercept formulas are allowed; a
#' covariate slope would be estimated to no effect. And the weight is only a
#' free parameter under a `"DirichletMultinomial"` likelihood, so a prior on a
#' non-DM fleet/predator would target a fixed parameter. Both error up front.
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param data_list the data list (for `fleet_control` and `Diet_distribution`).
#' @return invisibly NULL; errors on an unsupported comp linkage.
#' @keywords internal
#' @noRd
.check_comp_linkage_support <- function(linkage_table, data_list) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  cmp <- linkage_table[linkage_table$process == "comp", , drop = FALSE]
  if (nrow(cmp) == 0L) return(invisible())
  fc <- data_list$fleet_control

  # (a) prior-only: only intercept rows (theta is a scalar; no accumulator).
  if (any(cmp$design_col != "(Intercept)")) {
    stop("composition-weighting (comp) linkages are prior-only: use an ",
         "intercept formula `~ 1` with a `priors` entry (the DM weight is a ",
         "scalar, not year-varying), not a covariate slope.", call. = FALSE)
  }

  # (a2) `est_phase` still does not apply: the intercept coefficient it would
  # control is mapped out at 0 for comp, so phasing it does nothing. `init` and
  # `bounds` DO apply -- they re-target the DM weight itself, the same contract
  # every other process gives the intercept (see build_params() "Push
  # (Intercept) inits to the base parameter").
  if (any(cmp$est_phase != 1L)) {
    stop("composition-weighting (comp) linkages are prior-only: `est_phase` on ",
         "the spec does not fix or phase the DM weight (the intercept ",
         "coefficient it controls is mapped out at 0); fix the weight via ",
         "`map` instead.", call. = FALSE)
  }

  # (b) each row must name the ONE weight it targets. The cpp prior loop runs
  # once per linkage row and addresses a single weight -- comp/caal weights are
  # fleet-indexed, diet weights predator(species)-indexed -- so a row must carry
  # a concrete fleet (theta_comp / theta_caal) or species (theta_diet). There is
  # NO "all fleets/species" broadcast: an unstratified sentinel (NA) collapses
  # to the first index in the cpp, silently attaching the prior to fleet 1 /
  # predator 1 (and dropping the rest). `theta_comp` / `theta_caal` therefore
  # need `by = ~ fleet` (comp weights are not species-indexed); `theta_diet`
  # needs `by = ~ species`. Reject the wrong / missing stratum up front.
  dm_str <- function(v) as.character(v) == "DirichletMultinomial"
  for (r in seq_len(nrow(cmp))) {
    prm <- cmp$param[r]
    fl  <- cmp$fleet[r]
    sp  <- cmp$species[r]
    if (prm %in% c("theta_comp", "theta_caal")) {
      if (is.na(fl)) {
        stop(sprintf(paste0(
          "%s is a fleet-indexed DM weight: stratify the linkage by fleet ",
          "(`by = ~ fleet`, naming fleets via `fleet = `). An unstratified ",
          "spec would attach the prior to a single fleet only."), prm),
          call. = FALSE)
      }
      col <- if (prm == "theta_comp") "Comp_distribution" else "CAAL_distribution"
      if (!dm_str(fc[[col]][fl])) {
        stop(sprintf(paste0(
          "%s linkage on fleet %s whose %s is not 'DirichletMultinomial'; ",
          "the DM weight is not estimated there, so the prior has no effect."),
          prm, fc$Fleet_name[fl], col), call. = FALSE)
      }
    } else if (prm == "theta_diet") {
      if (is.na(sp)) {
        stop(paste0(
          "theta_diet is a predator-indexed DM weight: stratify the linkage ",
          "by species (`by = ~ species`). An unstratified spec would attach ",
          "the prior to a single predator only."), call. = FALSE)
      }
      if (as.integer(data_list$Diet_distribution[sp]) != 1L) {
        stop(sprintf(paste0(
          "theta_diet linkage on predator %s whose Diet_distribution is not ",
          "DirichletMultinomial; the DM weight is not estimated there."), sp),
          call. = FALSE)
      }
    }
  }
  invisible()
}


#' Drop comp priors whose DM weight is fixed in this configuration
#'
#' @description
#' A `comp` linkage row is prior-only, so when the map fixes the DM weight it
#' targets the prior is a constant: it shifts the reported `jnll` without moving
#' an estimate, and makes likelihoods non-comparable across configurations.
#' Such rows are set to `prior_family = "none"`, which the template skips.
#'
#' `.check_comp_linkage_support()` rejects a prior that can never apply to the
#' data at hand (a non-DM `Comp_distribution` / `Diet_distribution`). Whether a
#' weight is estimated *in a given fit* also depends on `msmMode`, `suitMode`,
#' and the fleet setup, and one `compFun` is routinely shared across the
#' single-species and multispecies fits of a stock, so those are reported and
#' ignored rather than rejected.
#'
#' Inertness is read off the finished `map` so it stays in step with
#' [build_map()] and honors a user-supplied `map`. Rows are kept, not dropped:
#' `beta_linkage` is dimensioned by `nrow(linkage_table)`, so dropping them
#' would break `inits` reuse between fits sharing a `compFun`.
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param map the map object from [build_map()] (uses `$mapList`).
#' @param verbose integer; 0 silences the message.
#' @return the linkage table, with inert comp priors neutralized.
#' @keywords internal
#' @noRd
.neutralize_inert_comp_priors <- function(linkage_table, map, verbose = 1) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(linkage_table)
  if (is.null(map) || is.null(map$mapList)) return(linkage_table)

  # comp param -> the map slot holding that DM weight, and how it is indexed.
  slots <- list(theta_comp = "comp_weights",
                theta_caal = "caal_weights",
                theta_diet = "diet_comp_weights")

  inert <- rep(FALSE, nrow(linkage_table))
  for (i in which(linkage_table$process == "comp" &
                  linkage_table$prior_family != "none")) {
    prm  <- linkage_table$param[i]
    slot <- slots[[prm]]
    if (is.null(slot)) next
    m <- map$mapList[[slot]]
    if (is.null(m)) next
    # theta_comp / theta_caal are fleet-indexed, theta_diet species-indexed.
    idx <- if (prm == "theta_diet") linkage_table$species[i] else linkage_table$fleet[i]
    if (is.na(idx) || idx < 1L || idx > length(m)) next
    inert[i] <- is.na(m[[idx]])
  }

  if (any(inert)) {
    if (verbose > 0) {
      message(sprintf(
        paste0("Ignoring %d composition-weighting prior(s) on a DM weight that ",
               "is not estimated in this configuration (%s). A prior on a fixed ",
               "parameter only adds a constant to the objective."),
        sum(inert),
        paste(sprintf("%s[%s]", linkage_table$param[inert],
                      ifelse(linkage_table$param[inert] == "theta_diet",
                             linkage_table$species[inert],
                             linkage_table$fleet[inert])),
              collapse = ", ")))
    }
    linkage_table$prior_family[inert] <- "none"
    linkage_table$prior_p1[inert]     <- NA_real_
    linkage_table$prior_p2[inert]     <- NA_real_
  }
  linkage_table
}
