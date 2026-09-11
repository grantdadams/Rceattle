#' Specify the stock-recruit relationship (SRR) for Rceattle
#'
#' @param srr_fun Stock-recruit function used in the hindcast estimation (see the list below). Default = 0
#' @param srr_pred_fun Stock-recruit function used for projection, reference points, and penalties (see below). When \code{srr_fun == 0}, the stock-recruit curve is added as a penalty on the annually estimated hindcast recruitment (following AMAK and Jim Ianelli's pollock model). If \code{srr_fun > 0}, then \code{srr_pred_fun = srr_fun} and no extra penalty is added.
#' @param proj_mean_rec Recruitment used in the projection: `TRUE`/1 (default) = mean recruitment, the average R over the hindcast; `FALSE`/0 = the stock-recruit relationship given by `srr_pred_fun`. Equilibrium and dynamic reference points follow the curve whenever `srr_pred_fun` is a stock-recruit form, regardless of this switch.
#' @param srr_hat_styr Integer. First year used to estimate the recruitment-penalty function (the AMAK/Ianelli penalty, active when \code{srr_pred_fun > 0} and \code{srr_fun = 0}), starting at \code{styr + 1}. Defaults to \code{styr + 1} in \code{data_list}. Useful when the environmental data conditioning the stock-recruit relationship is not available until the terminal year but projections are still wanted.
#' @param srr_hat_endyr Integer. Last year used to estimate the recruitment-penalty function (the AMAK/Ianelli penalty, active when \code{srr_pred_fun > 0} and \code{srr_fun = 0}). Defaults to \code{endyr} in \code{data_list}. Useful when the environmental data conditioning the stock-recruit relationship does not span the full time series but projections are still wanted.
#' @param srr_est_mode How the curve's built-in prior works, as a code or string: 1 / `"Estimated"` (default, no prior), 2 / `"LognormalPrior"` and 3 / `"BetaPrior"` (a prior on Beverton-Holt steepness, single-species only; mode 2 is a prior on alpha for Ricker), or 0 / `"Fixed"` (alpha held at `srr_prior`) -- for priors on, or fixed values of, `R0`, alpha or beta, use `linkages` instead (see **Priors, fixed values and covariates**).
#' @param srr_prior Prior centre on the natural scale: steepness for a Beverton-Holt curve under modes 2 and 3, alpha for a Ricker curve under mode 2, and the fixed alpha under mode 0.
#' @param srr_prior_sd Prior standard deviation: log scale for the lognormal prior (mode 2), natural scale for the beta prior (mode 3).
#' @param srr_alpha_init,srr_beta_init Optional starting values for alpha and beta, natural scale, one per species, used only when the curve estimates them (see **Starting values**).
#' @param srr_indices Defunct: supplying it is an error, because it has had no effect since 4.4.0. Express an environmental effect through `linkages`; see `vignette("environmental-linkages-and-priors")`.
#' @param Bmsy_lim Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a likelihood penalty if beta is estimated above this limit. Default `NA` is not used.
#' @param srr_mse_switchyr Year at which an MSE switches from the annual recruitment-penalty estimate to the stock-recruit function (the \code{srr_fun = 0}, \code{srr_pred_fun > 0} case).
#' @param linkages Named list of [linkage_spec()] objects keyed by `"R0"`, `"alpha"` or `"beta"`: the recommended way to put a prior on, fix, or add an environmental effect to those parameters (see **Priors, fixed values and covariates**).
#'
#' @description
#' Sets the stock-recruit curve and how recruitment is estimated. Priors, fixed
#' values and environmental effects on \code{R0}, alpha and beta go through
#' \code{linkages}; see **Priors, fixed values and covariates** below.
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
#' **Multispecies models.** Spawning biomass per recruit is undefined when
#' mortality includes predation, so under \code{msmMode > 0} a curve fitted in
#' the hindcast estimates its initial recruitment level (\code{R0}) rather than
#' deriving it; under the fished initModes (3, 4) that level trades off against
#' the initial F. The curve can also enter as the penalty above. A steepness
#' prior is refused and \code{steepness} is reported as 0; priors on alpha or
#' beta go through \code{linkages}.
#'
#' @section Priors, fixed values and covariates:
#' Use \code{linkages} for \code{R0}, alpha and beta. Each entry is a
#' [linkage_spec()], and an intercept-only formula (\code{~ 1}) acts on the
#' parameter itself:
#'
#' - **Prior:** \code{priors = list(`(Intercept)` = prior_lognormal(log(m), s))}
#'   is lognormal with median \code{m} and log-scale SD \code{s};
#'   [prior_normal()] is normal on the natural scale.
#' - **Fixed value:** \code{init = list(`(Intercept)` = v), est_phase = 0}
#'   holds the parameter at \code{v}.
#' - **One species:** add \code{species = 1}; the default applies to every
#'   species.
#' - **Environmental effect:** a covariate formula such as \code{~ temp} adds a
#'   log-scale effect by year; see
#'   \code{vignette("environmental-linkages-and-priors")}.
#'
#' For a Ricker curve the lognormal linkage prior on alpha is identical to
#' \code{srr_est_mode = "LognormalPrior"}. \code{srr_est_mode} and
#' \code{srr_prior} remain for the one prior a linkage cannot express, on
#' Beverton-Holt steepness, which needs spawning biomass per recruit and so
#' exists only in single-species models.
#'
#' @section Starting values:
#' The defaults (\eqn{\alpha = e^3}, \eqn{\beta = 3}) know nothing of the
#' stock's scale. \eqn{\beta} sets the density dependence in
#' \eqn{R = \alpha S / (1 + \beta S)}, so it must be on the order of
#' \eqn{(\alpha - 1/\phi_0) / R_0} -- typically \eqn{10^{-3}} or smaller for a
#' stock measured in tonnes; starting three orders of magnitude away drives
#' predicted recruitment to near zero and the optimizer returns
#' \code{NA/NaN gradient evaluation}. From a steepness \eqn{h} and unfished
#' spawning biomass per recruit \eqn{\phi_0}:
#' \deqn{\alpha = \frac{4h}{\phi_0 (1 - h)}, \qquad
#'       \beta  = \frac{\alpha - 1/\phi_0}{R_0}.}
#' Under predation, where \eqn{\phi_0} is undefined, seed them from a curve
#' fitted to an earlier fit's SSB-recruitment pairs.
#'
#' @return A \code{list} containing the stock recruitment relationship settings
#' @examples
#' # Mean recruitment: no stock-recruit relationship fitted.
#' build_srr(srr_fun = "mean")
#'
#' # Beverton-Holt with a lognormal prior on alpha (median 5, log-scale SD 0.5),
#' # species 1 only.
#' build_srr(srr_fun = "BevertonHolt",
#'           linkages = list(alpha = linkage_spec(~ 1, species = 1,
#'             priors = list(`(Intercept)` = prior_lognormal(log(5), 0.5)))))
#'
#' # Beverton-Holt with alpha fixed at 20.
#' build_srr(srr_fun = "BevertonHolt",
#'           linkages = list(alpha = linkage_spec(~ 1, est_phase = 0,
#'             init = list(`(Intercept)` = 20))))
#'
#' # A temperature effect on alpha (log scale; BTempC is a column of env_data).
#' build_srr(srr_fun = "BevertonHolt",
#'           linkages = list(alpha = linkage_spec(~ BTempC)))
#'
#' # Mean recruitment with the curve as a penalty (Ianelli form).
#' build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt")
#'
#' # A prior on Beverton-Holt steepness, the one prior linkages cannot express
#' # (single-species only).
#' build_srr(srr_fun = "BevertonHolt", srr_est_mode = "LognormalPrior",
#'           srr_prior = 0.8, srr_prior_sd = 0.2)
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

  # Bmsy_lim bounds the Ricker curve only. -999 is the stored "off" value that
  # refits pass back in, so it does not warn.
  if(!srr_pred_fun %in% c(4,5, "Ricker")){
    if (any(!is.na(Bmsy_lim) & Bmsy_lim != -999)) {
      warning("`Bmsy_lim` bounds a Ricker curve only and is ignored for ",
              "`srr_pred_fun = ", srr_pred_fun, "`.", call. = FALSE)
    }
    Bmsy_lim = -999
  }

  # The beta prior is a prior on Beverton-Holt steepness; the template has no
  # Ricker form of it, so a Ricker curve would be estimated with no prior.
  if (isTRUE(srr_est_mode == 3) && srr_pred_fun %in% c(4, 5)) {
    stop("`srr_est_mode = 3` (\"BetaPrior\") is a prior on Beverton-Holt ",
         "steepness and has no Ricker form. For a Ricker curve use ",
         "`srr_est_mode = 2` (\"LognormalPrior\"), a prior on alpha.",
         call. = FALSE)
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

  # `srr_indices` has had no effect since 4.4.0, so it stops rather than fit a model
  # without its covariate. NA or NULL is "not supplied" (.refit_like() passes NULL).
  if (!is.null(srr_indices) && !(length(srr_indices) == 1L && is.na(srr_indices))) {
    .stop_srr_indices_defunct()
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
#' 2, 4) get string aliases. The env-driven codes (1, 3, 5) are errors
#' from 5.32.0: the environmental effect is expressed through the
#' `linkages` argument to [build_srr()].
#'
#' @keywords internal
.SRR_FUNS <- c(
  mean         = 0L,
  BevertonHolt = 2L,
  Ricker       = 4L
)


#' Retired env-driven `srr_fun` / `srr_pred_fun` integer codes: an error in
#' [build_srr()], mapped to 0 / 2 / 4 when a refit reads one off an older fit.
#' @keywords internal
#' @noRd
.SRR_DEPRECATED_FUNS <- c(1L, 3L, 5L)


#' Coerce an `srr_fun` / `srr_pred_fun` value to canonical integer.
#'
#' Accepts either a string from [.SRR_FUNS] (length-1) or a length-1
#' integer 0, 2 or 4. Codes 1, 3 and 5 stop with the linkage that
#' replaces them.
#'
#' @keywords internal
#' @noRd
.coerce_srr_fun <- function(x, what) {
  .coerce_switch_arg(
    x, map = .SRR_FUNS, what = what,
    deprecated = .SRR_DEPRECATED_FUNS,
    warn_fn = function(int) .stop_srr_fun_defunct(int, what),
    length_exact_one = TRUE)
}


# A fit made before 5.32.0 can store code 1, 3 or 5. From 4.4.0 those fitted the
# structural form 0, 2 or 4 with no environmental term, so a refit maps them there.
.srr_fun_structural <- function(x) {
  if (is.null(x)) return(x)
  x   <- as.integer(x)
  old <- !is.na(x) & x %in% .SRR_DEPRECATED_FUNS
  if (any(old)) {
    warning(sprintf(paste0(
      "This fit used srr_fun / srr_pred_fun = %s, whose environmental term has had ",
      "no effect since 4.4.0; refitting it as %s, the model it fitted. Refit with a ",
      "recruitment linkage to include the environmental effect."),
      paste(unique(x[old]), collapse = ", "), paste(unique(x[old] - 1L), collapse = ", ")),
      call. = FALSE)
    x[old] <- x[old] - 1L
  }
  x
}


#' @keywords internal
#' @noRd
.stop_srr_fun_defunct <- function(int, what) {
  form <- switch(as.character(int), "1" = "mean recruitment",
                 "3" = "Beverton-Holt", "Ricker")
  stop(
    sprintf("%s = %d (%s with an environmental effect) is no longer supported: ",
            what, int, form),
    "from 4.4.0 it fitted the model without its environmental term. Express ",
    "the effect as a linkage:\n\n",
    "  build_srr(", what, " = ", int - 1L, ",\n",
    "            linkages = list(", switch(as.character(int), "1" = "R0", "alpha"),
    " = linkage_spec(~ <env_col>)))\n\n",
    "See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.stop_srr_indices_defunct <- function() {
  stop(
    "`srr_indices` is no longer supported: from 4.4.0 it had no effect, so a ",
    "model fitted with it carried no environmental term. Express the effect as ",
    "a linkage:\n\n",
    "  build_srr(srr_fun = ...,\n",
    "            linkages = list(R0 = linkage_spec(~ <env_col>)))\n\n",
    "srr_indices = k referred to env_data column k + 1, counting after Year. ",
    "See vignette('environmental-linkages-and-priors').",
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
