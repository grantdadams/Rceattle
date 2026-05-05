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
#' @param srr_indices vector or single index indicating the columns (excluding "Year") of \code{env_data} to use in a environmentally driven stock recruit curve.
#' @param Bmsy_lim Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a likelihood penalty if beta is estimated above this limit. Default `NA` is not used.
#' @param srr_mse_switchyr is used for MSEs to deal with AMAK and Jim Ianelli's estimation where a stock recruit function is estimated as an additional penalty  (srr_fun = 0 and srr_pred_fun > 0). It tells the model in what year to switch to the stock recruit function.
#'
#' @description
#'
#' **Stock recruitment relationships currently implemented in Rceattle:**
#'
#' - \code{srr_fun = 0} No stock recruit relationship. Recruitment is a function of R0 and annual deviates (i.e. steepness = 0.99).
#'  \deqn{R_y = exp(R0 + R_{dev,y})}
#'
#' - \code{srr_fun = 1} Environmentally driven recruitment without stock recruit relationship
#'  \deqn{R_y = exp(R0 + R_{dev,y} + X * \beta_X)}
#'
#' - \code{srr_fun = 2} Beverton-holt stock-recruitment relationship
#'   \deqn{R_y = \frac{\alpha_{srr} * SB_{y-minage}}{1+\beta_{srr} * SB_{y-minage}}}
#'
#' - \code{srr_fun = 3} Beverton-holt stock-recruitment relationship with environmental covariates impacting larval survival rate and prior is on alpha.
#'   \deqn{R_y = \frac{\alpha_{srr} * e^{X * \beta_X} * SB_{y-minage}}{1+\beta_{srr} * SB_{y-minage}}}
#'
#' - \code{srr_fun = 4} Ricker stock-recruitment relationship
#'   \deqn{R_y = \alpha_{srr} * SB_{y-minage} * exp(-\beta_{srr} * SB_{y-minage})}
#'
#' - \code{srr_fun = 5} Ricker stock-recruitment relationship with environmental covariates impacting larval survival rate and prior is on alpha.
#'   \deqn{R_y = \alpha_{srr} e^{X * \beta_X} * SB_{y-minage} * exp(-\beta_{srr} * SB_{y-minage})}
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
                      Bmsy_lim = NA){

  # Set pred/RP/penalty to same as SR curve if SR fun > 0
  if(srr_fun > 0){
    srr_pred_fun = srr_fun
  }

  if(!srr_pred_fun %in% c(3,4)){
    Bmsy_lim = -999
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
       Bmsy_lim = Bmsy_lim
  )
}



#' Define M1 specifications
#'
#' @param M1_model Vector or scalar specifying M1 fixed effects model (see @description below). 0 = use fixed natural mortality from M1_base in data, 1 = estimate sex- and age-invariant M1, 2 = sex-specific (two-sex model), age-invariant M1, 3 = estimate sex- and age-specific M1, 4 = environmentally driven sex- and age-invariant M1, 5 = environmentally driven age-invariant, but sex-specific M1.
#' @param M1_re Vector or scalar specifying M1 random effects model. See description (default = 0).
#' @param updateM1 If using initial parameters, use M1 fixed effects from data instead (default = FALSE).
#' @param M1_use_prior Vector or scalar specifying if M1 fixed effects come from a lognormal prior
#' @param M2_use_prior Vector or scalar specifying if M1 + M2 come from a lognormal prior in multi-species models (default = FALSE). Lognormal prior for M1 + M2 across species, sexes, ages, and years.
#' @param M_prior Vector or scalar for mean of M prior on natural scale
#' @param M_prior_sd Vector or scalar of SD of lognormal M prior. Used as initial value for random effects variance as well.
#' @param M1_indices vector or single index indicating the columns (excluding Year column) of \code{env_data} to use for environmentally linked M1 when \code{M1_model} is 4 or 5.
#'
#' @description
#'
#' **M1 fixed effects currently implemented in CEATTLE**
#' - \code{M1_model = 0} Fixed based on input \code{M1_base}
#' - \code{M1_model = 1} Single species specific M. Estimates \deqn{M1_{spp}}
#' - \code{M1_model = 2} Sex-specific M. Estimates \deqn{M1_{spp, sex}}
#' - \code{M1_model = 3} Sex- and age-specific M. Estimates \deqn{M1_{spp, sex, age}}
#' - \code{M1_model = 4} Sex-invariant environmental effect. Estimates \deqn{M1_{spp, yr} = M1_{spp} * e^{X * \beta_{X, spp}}}. \code{M1_indices} specifies the environmental indices.
#' - \code{M1_model = 5} Sex- and species specific environmental effect. Estimates \deqn{M1_{spp, sex, yr} = M1_{spp, sex} * e^{X * \beta_{X, spp, sex}}}. \code{M1_indices} specifies the environmental indices.
#' - TODO fit to environmental index
#'
#' **M1 random effects currently implemented in CEATTLE**
#'
#' M1 random effects are applied to each species if \code{M1_model = 1} or each species and sex if \code{M1_model = 2}. Variance and correlation coefficients are species-specific, but sex-invariant.
#'
#' - \code{M1_re = 0}: No random effects (default).
#' - \code{M1_re = 1}: Random effects varies by age, but uncorrelated (IID) and constant over years.
#' - \code{M1_re = 2}: Random effects varies by year, but uncorrelated (IID) and constant over ages.
#' - \code{M1_re = 3}: Random effects varies by year and age, but uncorrelated (IID).
#' - \code{M1_re = 4}: Correlated AR1 random effects varies by age, but constant over years.
#' - \code{M1_re = 5}: Correlated AR1 random effects varies by year, but constant over ages.
#' - \code{M1_re = 6}: Correlated 2D-AR1 random effects varies by year and age.
#'
#' @return A list of switches for defining the M1 model
#' @export
#'
#' Allowed M-parameter names for `linkages` in [build_M1()]
#'
#' Linear-predictor names of the underlying natural-mortality
#' parameters that the linkage system can address. Currently just
#' `log_M` -- the offset is added on the log scale to `ln_M1`
#' (broadcast across age unless the linkage row pins a specific
#' `age_bin`).
#'
#' @keywords internal
M_LINKAGE_PARAMS <- c("log_M")


build_M1 <- function(M1_model = 0,
                     M1_re = 0,
                     updateM1 = FALSE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE,
                     M_prior = 0.40,
                     M_prior_sd = 0.35,
                     M1_indices = NA,
                     linkages = NULL){
  linkages <- .validate_M_linkages(linkages)
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
  if (is.null(linkages)) return(NULL)
  if (!is.list(linkages) || length(linkages) == 0L ||
      is.null(names(linkages)) || any(!nzchar(names(linkages)))) {
    stop("`linkages` must be a named list keyed by M parameter ",
         "(one of: ", paste(M_LINKAGE_PARAMS, collapse = ", "), ")",
         call. = FALSE)
  }
  bad_names <- setdiff(names(linkages), M_LINKAGE_PARAMS)
  if (length(bad_names) > 0) {
    stop("unknown M linkage parameter(s): ",
         paste(bad_names, collapse = ", "),
         "; allowed: ", paste(M_LINKAGE_PARAMS, collapse = ", "),
         call. = FALSE)
  }
  # Each linkages[[param]] may be a single linkage_spec() or a list
  # of them (e.g. species-specific formulas registered together).
  for (nm in names(linkages)) {
    val <- linkages[[nm]]
    if (inherits(val, "Rceattle_linkage_spec")) next
    if (is.list(val) &&
        all(vapply(val, inherits, logical(1),
                   what = "Rceattle_linkage_spec"))) next
    stop("linkages$", nm, " must be a linkage_spec() or a list of ",
         "linkage_spec() objects.", call. = FALSE)
  }
  Map(.stamp_param, linkages, names(linkages))
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
#' Linear-predictor names of the underlying growth-function parameters.
#' Von Bertalanffy uses `log_K`, `log_L1`, `log_Linf`; Richards adds
#' `log_m`. The empirical weight-at-age model admits no linkages.
#'
#' @keywords internal
GROWTH_LINKAGE_PARAMS <- c("log_K", "log_L1", "log_Linf", "log_m")


#' Specify the growth model for Rceattle
#'
#' @param fun Growth function. Either a string ([GROWTH_FUNS]:
#'   `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
#'   equivalent integer code (`0`, `1`, `2`). The canonical string form
#'   is stored on the returned object.
#' @param linkages Optional named list of [linkage_spec()] objects
#'   keyed by parameter name (must be one of [GROWTH_LINKAGE_PARAMS]).
#'   Each spec describes how that growth parameter depends on
#'   environmental covariates and on stratifying factors (species,
#'   sex). The parameter name on each spec is filled in from the list
#'   key. Materialization into the global linkage table happens inside
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
#'     log_K = linkage_spec(
#'       formula = ~ temp,
#'       by      = ~ species + sex,
#'       priors  = list(temp = normal(0, 1))
#'     )
#'   )
#' )
#' }
build_growth <- function(fun = "empirical", linkages = NULL) {
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
  if (is.null(linkages)) return(NULL)
  if (!is.list(linkages) || length(linkages) == 0L ||
      is.null(names(linkages)) || any(!nzchar(names(linkages)))) {
    stop("`linkages` must be a named list keyed by growth parameter ",
         "(one of: ", paste(GROWTH_LINKAGE_PARAMS, collapse = ", "), ")",
         call. = FALSE)
  }
  bad_names <- setdiff(names(linkages), GROWTH_LINKAGE_PARAMS)
  if (length(bad_names) > 0) {
    stop("unknown growth linkage parameter(s): ",
         paste(bad_names, collapse = ", "),
         "; allowed: ", paste(GROWTH_LINKAGE_PARAMS, collapse = ", "),
         call. = FALSE)
  }
  # Each linkages[[param]] may be either a single linkage_spec() or a
  # list of them (e.g. species-specific formulas registered under the
  # same parameter). Normalise both shapes; reject anything else.
  for (nm in names(linkages)) {
    val <- linkages[[nm]]
    if (inherits(val, "Rceattle_linkage_spec")) next
    if (is.list(val) &&
        all(vapply(val, inherits, logical(1),
                   what = "Rceattle_linkage_spec"))) next
    stop("linkages$", nm, " must be a linkage_spec() or a list of ",
         "linkage_spec() objects.", call. = FALSE)
  }
  if (any(fun == "empirical")) {
    warning("at least one species uses fun = 'empirical', which does ",
            "not consume linkages; the supplied specs will be ",
            "retained on the object but ignored at fit time for those ",
            "species.", call. = FALSE)
  }
  if (!all(fun == "Richards") && "log_m" %in% names(linkages)) {
    stop("linkages$log_m is only valid when every species uses ",
         "fun = 'Richards'; von Bertalanffy has no `m` parameter.",
         call. = FALSE)
  }
  # Stamp the parameter name onto each spec (single or list-of-specs)
  # using the linkages-list key, so users don't have to repeat the name.
  Map(.stamp_param, linkages, names(linkages))
}


#' @keywords internal
#' @noRd
.stamp_param <- function(val, param) {
  if (inherits(val, "Rceattle_linkage_spec")) {
    return(.set_linkage_param(val, param))
  }
  lapply(val, .set_linkage_param, param = param)
}
