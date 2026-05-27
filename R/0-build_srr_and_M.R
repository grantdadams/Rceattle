#' Specify the stock-recruit relationship (SRR) for Rceattle
#'
#' @param srr_fun Stock recruit function to be used for hindcast estimation of Rceattle (see @description below). Default = 0
#' @param srr_pred_fun stock recruit function for projection, reference points, and penalties to be used for Rceattle (see below). When \code{srr_fun == 0}, it treats the stock-recruit curve as an additional penalty onto the annualy estimated recruitment from the hindcast (sensu AMAK and Jim Ianelli's pollock model). If \code{srr_fun > 0} then \code{srr_pred_fun = srr_fun} and no additional penalty is included.
#' @param proj_mean_rec Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = SRR(omega, srr_devs)
#' @param srr_hat_styr Integer. The first year used for estimating the recruitment function as an additional penalty. It will add additional penalties sensu AMAK and Jim Ianelli's pollock model when \code{srr_pred_fun > 0} and \code{srr_fun = 0}, starting at \code{styr} + 1. Defaults to $styr + 1$ in $data_list$. Useful if environmental data used to condition stock-recruit relationships is not available until end-year, but projections are desired.
#' @param srr_hat_endyr Integer. The last year used for estimating the recruitment function as an additional penalty. It will add an additional penalties sensu AMAK and Jim Ianelli's pollock model when \code{srr_pred_fun > 0} and \code{srr_fun = 0}. Recruitment Defaults to $endyr$ in $data_list$. Useful if environmental data used to condition stock-recruit relationships is not available for the full time-series, but projections are desired.
#' @param srr_est_mode Switch to determine estimation mode. 0 = fix alpha to prior mean, 1 = freely estimate R0, alpha, and/or beta (default), 2 = use lognormally distributed prior for alpha (Ricker) or steepness (Beverton), 3 = use beta distributed prior for steepness (Beverton) given mean and sd.
#' @param srr_prior mean for normally distributed prior for stock-recruit parameter
#' @param srr_prior_sd Prior standard deviation for stock-recruit parameter
#' @param srr_indices Soft-deprecated. Use the `linkages` argument instead. See `vignette("environmental-linkages-and-priors")`.
#' @param Bmsy_lim Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a likelihood penalty if beta is estimated above this limit. Default `NA` is not used.
#' @param srr_mse_switchyr is used for MSEs to deal with AMAK and Jim Ianelli's estimation where a stock recruit function is estimated as an additional penalty  (srr_fun = 0 and srr_pred_fun > 0). It tells the model in what year to switch to the stock recruit function.
#' @param linkages Optional named list of [linkage_spec()] objects keyed by recruitment parameter name (must be one of `"R0"`, `"alpha"`, `"beta"`). Each spec describes how that parameter depends on environmental covariates and on stratifying factors (species, sex). The offset enters additively (on the log scale) inside the recruitment compute. See `vignette("environmental-linkages-and-priors")` for details.
#'
#' @description
#'
#' **Stock recruitment relationships currently implemented in Rceattle:**
#'
#' - \code{srr_fun = 0} or \code{"mean"}: No stock recruit relationship. Recruitment is a function of R0 and annual deviates (i.e. steepness = 0.99).
#'  \deqn{R_y = exp(R0 + R_{dev,y})}
#'
#' - \code{srr_fun = 2} or \code{"BevertonHolt"}: Beverton-holt stock-recruitment relationship
#'   \deqn{R_y = \frac{\alpha_{srr} * SB_{y-minage}}{1+\beta_{srr} * SB_{y-minage}}}
#'
#' - \code{srr_fun = 4} or \code{"Ricker"}: Ricker stock-recruitment relationship
#'   \deqn{R_y = \alpha_{srr} * SB_{y-minage} * exp(-\beta_{srr} * SB_{y-minage})}
#'
#' When \code{srr_pred_fun > 0} and \code{srr_fun = 0} recruitment in the hindcast is estimated as in \code{srr_fun = 0} \deqn{R_y = exp(R0 + R_{dev,y})}, but an additional stock recruitment relationship defined by \code{srr_pred_fun} is estimated between \code{srr_hat_styr} and \code{srr_hat_endyr} and treated as an additional penalty. The stock recruitment relationship defined by \code{srr_pred_fun} is then used in the projection.
#'
#'
#' @return A \code{list} containing the stock recruitment relationship settings
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
                      srr_indices = NA,
                      Bmsy_lim = NA,
                      linkages = NULL){

  srr_fun      <- .coerce_srr_fun(srr_fun,      "srr_fun")
  srr_pred_fun <- .coerce_srr_fun(srr_pred_fun, "srr_pred_fun")

  # Set pred/RP/penalty to same as SR curve if SR fun > 0
  if(srr_fun > 0){
    srr_pred_fun = srr_fun
  }

  if(!srr_pred_fun %in% c(4,5, "Ricker")){
    Bmsy_lim = -999
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
       srr_indices = srr_indices,
       Bmsy_lim = Bmsy_lim,
       linkages = linkages
  )
}


#' String<->integer mapping for `srr_fun` / `srr_pred_fun` in
#' [build_srr()]
#'
#' Either form is accepted; the canonical integer code is what the
#' TMB template ultimately consumes. Only the structural codes (0,
#' 2, 4) get string aliases. The historical env-driven codes (1, 3,
#' 5) still work with a soft-deprecation warning -- their structural
#' part is identical to 0 / 2 / 4 respectively, and the env effect
#' is now expressed via the `linkages` argument to [build_srr()].
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
  if (length(x) != 1L) {
    stop(sprintf("`%s` must be length 1", what), call. = FALSE)
  }
  if (is.numeric(x)) {
    int <- as.integer(x)
    allowed <- c(unname(.SRR_FUNS), .SRR_DEPRECATED_FUNS)
    if (is.na(int) || !int %in% allowed) {
      stop(sprintf(
        "integer `%s` must be one of: %s (= %s, plus 1/3/5 for legacy env modes)",
        what,
        paste(.SRR_FUNS, collapse = ", "),
        paste(names(.SRR_FUNS), collapse = "/")), call. = FALSE)
    }
    if (int %in% .SRR_DEPRECATED_FUNS) {
      .warn_srr_fun_deprecation(int, what)
    }
    return(int)
  }
  if (is.character(x)) {
    bad <- setdiff(x, names(.SRR_FUNS))
    if (length(bad) > 0L) {
      stop(sprintf(
        "unknown `%s` value(s): %s; allowed: %s",
        what,
        paste(unique(bad), collapse = ", "),
        paste(names(.SRR_FUNS), collapse = ", ")), call. = FALSE)
    }
    return(unname(.SRR_FUNS[x]))
  }
  stop(sprintf("`%s` must be a string or integer; got %s",
               what, class(x)[1]), call. = FALSE)
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
#' `M1` -- with the default log link the offset is added on the log
#' scale to `log_M1` (broadcast across age unless the linkage row
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
#' code is what the TMB template ultimately consumes. The legacy
#' env-driven integer codes 4 and 5 still work with a deprecation
#' warning -- their structural part is identical to 1 and 2
#' respectively, and the env effect is now expressed via the
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
#' env-driven options (4, 5) on M1_model are now better expressed
#' through the `linkages` argument to [build_M1()] but the integer
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
#' For `M1_model` only, the historical env-driven integer codes 4
#' and 5 (controlled by `M1_indices`) are accepted for backwards
#' compatibility but emit a soft-deprecation warning pointing users
#' at the linkage table -- see [build_M1()].
#'
#' @keywords internal
#' @noRd
.coerce_M1_arg <- function(x, map, what) {
  if (length(x) == 0L) {
    stop(sprintf("`%s` must have length >= 1", what), call. = FALSE)
  }
  if (is.numeric(x)) {
    int <- as.integer(x)
    extra <- if (identical(what, "M1_model")) .M1_DEPRECATED_MODELS else integer(0)
    allowed <- c(unname(map), extra)
    if (anyNA(int) || any(!int %in% allowed)) {
      stop(sprintf(
        "integer `%s` must be one of: %s (= %s)",
        what,
        paste(map, collapse = ", "),
        paste(names(map), collapse = "/")), call. = FALSE)
    }
    if (length(extra) > 0L && any(int %in% extra)) {
      .warn_M1_model_deprecation(int)
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
  stop(sprintf(
    "`%s` must be a string or integer; got %s",
    what, class(x)[1]), call. = FALSE)
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


#' Define M1 specifications
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
#' \dontrun{
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


#' Allowed growth-parameter names for `linkages` in [build_growth()]
#'
#' Natural-scale names of the underlying growth-function parameters.
#' Von Bertalanffy uses `K`, `L1`, `Linf`; Richards adds `m`.
#' `sd_L1` / `sd_Linf` are the standard deviations of length-at-age
#' anchored at `L1` and `Linf` (the SD-at-age interpolation endpoints
#' from `growth.hpp`). Only intercept-only specs (`~ 1`) are honored
#' on the SD endpoints -- they thread through `init` / `bounds` /
#' `priors` onto `growth_log_sd` but do not vary by year.
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
#' @param linkages Optional named list of [linkage_spec()] objects
#'   keyed by parameter name (must be one of [GROWTH_LINKAGE_PARAMS]).
#'   The mean-growth keys (`K`, `L1`, `Linf`, `m`)
#'   accept arbitrary one-sided formulas and produce year-varying
#'   offsets applied inside `growth.hpp`. The SD-endpoint keys
#'   (`sd_L1`, `sd_Linf`) only honor intercept-bearing
#'   formulas (typically `~ 1`) -- they thread `init`, `bounds`, and
#'   `priors` onto the underlying `growth_log_sd` parameter, giving
#'   the SDs the same prior/fix/initial-value contract as the mean
#'   parameters. Slope rows on SD specs raise a warning and have no
#'   effect; slope-only formulas (`~ 0 + temp`) error.
#'   Materialization into the global linkage table happens inside
#'   `fit_mod()` once `data_list$env_data` is in scope.
#'
#' @return A list of switches defining the growth model.
#' @export
#'
#' @examples
#' \dontrun{
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
                         minlength = NA,
                         lengthbin = NA,
                         maxlength = NA,
                         linkages = NULL) {
  fun <- .coerce_growth_fun(fun)
  linkages <- .validate_growth_linkages(linkages, fun)
  list(
    fun            = fun,
    linkages       = linkages,
    # Internal integer code consumed by the TMB template until the
    # linkage-driven path replaces it. Vectorized so per-species
    # growth functions (e.g. fun = c("vonBertalanffy", "Richards"))
    # propagate downstream as before.
    growth_model   = unname(.GROWTH_FUN_TO_INT[fun]),
    # VB anchor age (= age at which `l1` is the length). Matches SS3's
    # `Growth_Age_for_L1` ctl input. Scalar or length-nspp; pass NA
    # (default) to inherit max(0.5, minage[sp]) downstream so old
    # configurations stay unchanged at minage >= 1 and minage = 0
    # models pick up an SS3-consistent half-year anchor.
    growth_age_L1  = growth_age_L1,
    # Population length grid (for parametric growth only) -- the ALK and
    # WAA integration happen on this grid. SS3 parity requires a fine
    # grid (~1cm) since the Jensen's-inequality integral over
    # alpha * L^beta converges slowly. Scalars are extended per-species
    # downstream; NA falls back to the data-bin grid (i.e. `lengths` from
    # caal_data), preserving prior behavior for empirical-growth models.
    minlength      = minlength,
    lengthbin      = lengthbin,
    maxlength      = maxlength,
    # Placeholders retained until the random-effects / legacy index
    # paths in 2-build_map.R and src/TMB/growth.hpp are migrated.
    growth_re      = 0L,
    growth_indices = NA
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
.stamp_param <- function(val, param) {
  if (inherits(val, "Rceattle_linkage_spec")) {
    return(.set_linkage_param(val, param))
  }
  lapply(val, .set_linkage_param, param = param)
}
