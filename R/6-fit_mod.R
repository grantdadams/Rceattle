#' Fit the CEATTLE assessment model
#' @description Estimate CEATTLE population parameters by maximum likelihood, and
#'   optionally project the stock and apply a harvest control rule.
#'
#' @param data_list A data list read in via \code{\link{read_data}} or built
#'   directly in R; see \code{vignette("data-without-excel", package = "Rceattle")}.
#' @param inits (Optional) A named list of initial parameter values, as returned by
#'   \code{\link{build_params}} or extracted from a previous fit
#'   (\code{model$estimated_params}). If \code{NULL}, parameters are initialized
#'   from scratch via \code{\link{build_params}}.
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param bounds (Optional) A bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param estimateMode What to fit, given as a string alias or the integer code:
#'   \code{"Estimate"} (0) = fit the hindcast model and the HCR projection
#'   (\code{HCR}); \code{"Hindcast"} (1) = fit the hindcast only (no fitting BRPs/HCR/projection);
#'   \code{"Projection"} (2) = fit the BRPs/HCR/projection only, from the initial
#'   parameters in \code{inits}; \code{"DebugBuild"} (3) = build through
#'   \code{MakeADFun} but not \code{nlminb} -- the returned \code{obj} carries the
#'   real objective and gradient, so \code{obj$fn()} / \code{obj$gr()} are usable
#'   for diagnosing a model before committing to a fit; \code{"DebugOptimize"}
#'   (4) = optimize with all parameters mapped out, so the objective is a
#'   placeholder (\code{dummy^2}), not a likelihood. Defaults to \code{"Estimate"}.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects using the Laplace approximation. The default is FALSE.
#' @param random_q logical. If TRUE, treats annual catchability deviations as random effects using the Laplace approximation. The default is FALSE.
#' @param random_sel logical. If TRUE, treats annual selectivity deviations as random effects using the Laplace approximation. The default is FALSE.
#' @param HCR HCR list object from \code{\link{build_hcr}}
#' @param niter Number of iterations for multispecies model
#' @param recFun The stock recruit-relationship parameterization from \code{\link{build_srr}}.
#' @param M1Fun M1 parameterizations and priors. Use \code{build_M1}.
#' @param growthFun The weight-at-age parameterization from \code{\link{build_growth}}.
#' @param qFun Catchability specification from \code{\link{build_catchability}},
#'   carrying any environmental linkages on q.
#' @param selFun Selectivity specification from \code{\link{build_selectivity}}, carrying any environmental linkages on selectivity parameters.
#' @param compFun Composition-weighting specification from \code{\link{build_composition}}, carrying any priors on the Dirichlet-multinomial weights.
#' @param msmMode The predation-mortality mode, as a string alias or integer code:
#'   \code{"SingleSpecies"} (0, the default, no predation), \code{"MSVPA"}
#'   (1, the Type-II MSVPA predation of Holsman et al. 2015) or
#'   \code{"TypeIIIMSVPA"} (2). Higher integer
#'   codes (Kinzey-Punt, Holling forms) are not yet implemented.
#' @param avgnMode The average abundance-at-age approximation used in the predation-mortality equations. Only mode 0, \eqn{N/Z ( 1 - exp(-Z) )} (the MSVPA form), is currently active; the model always uses it. The alternatives \eqn{N exp(-Z/2)} (1) and \eqn{N} (2) are not currently implemented and setting them has no effect.
#' @param initMode how the population is initialized. 0 = initial age-structure estimated as free parameters; 1 = equilibrium age-structure estimated out from R0 + mortality (M1); 2 = non-equilibrium age-structure estimated out from R0,  mortality (M1), and initial population deviates; 3 = non-equilibrium age-structure estimated out from initial fishing mortality (Finit), R0,  mortality (M1), and initial population deviates; 4 = non-equilibrium age-structure version 2 where initial fishing mortality (Finit) scales R0; 5 = "FishedEquilibrium": F = 0 equilibrium age-structure seeded by the first-year recruitment (\code{exp(rec_pars + rec_dev[year 1])}) decayed by mortality (M1), with initial deviates turned off and no init-dev penalty (the Cole Monnahan / AFSC GOA pollock convention).
#' @param suitMode Switch for suitability derivation for each predator (single value or vector). 0 = empirical based on diet data (Holsman et al. 2015), 1 = length-based gamma suitability, 2 = weight-based gamma suitability, 3 = length-based lognormal suitability, 4 = weight-based lognormal suitability, 5 = length-based normal suitability, 6 = weight-based normal suitability. The length-based modes (1, 3, 5) are not yet implemented and are rejected by `data_check()`; use a weight-based mode (2, 4, 6) or empirical suitability (0).
#' @param suit_styr The first year used to calculate mean suitability. A single integer is applied to every predator, or a vector of length `nspp` sets a distinct start year per predator. Defaults to `styr` in `data_list`. Used when diet data were sampled from a subset of years.
#' @param suit_endyr The last year used to calculate mean suitability. A single integer is applied to every predator, or a vector of length `nspp` sets a distinct end year per predator. Defaults to `endyr` in `data_list`. Used when diet data were sampled from a subset of years.
#' @param fit_control A list returned by [fit_control()] that bundles the
#'   optimizer / sdreport / phasing knobs (`phase`, `bias.correct`, `getsd`,
#'   `getJointPrecision`, `getReportCovariance`, `use_gradient`, `rel_tol`,
#'   `loopnum`, `newtonsteps`, `TMBfilename`, `verbose`, `nlminb_control`).
#'   Defaults to `fit_control()`. See [fit_control()] for the meaning and
#'   defaults of each field.
#' @param config (Optional) An `Rceattle_run_config` from [load_config()] (or
#'   [run_config()]). Its stored `model_config` structure and estimation controls
#'   (`estimateMode`, `random_rec`/`random_q`/`random_sel`, `suit_styr`/
#'   `suit_endyr`, `fit_control`) overlay only the arguments the caller did *not*
#'   pass -- an explicit argument always wins. `NULL` (default) applies no
#'   configuration. Example: `fit_mod(data_list, config = load_config("run.yaml"))`.
#' @param ... Deprecated optimizer / sdreport / phasing arguments
#'   (e.g. `phase`, `getsd`, `bias.correct`, `use_gradient`, `rel_tol`,
#'   `control`, `getJointPrecision`, `getReportCovariance`, `loopnum`,
#'   `newtonsteps`, `verbose`, `TMBfilename`). These are forwarded into
#'   `fit_control` with a deprecation warning; pass them via
#'   [fit_control()] instead.
#'
#' @details
#' CEATTLE is an age-structured population dynamics model that can be fit with or without predation mortality. The default is to exclude predation mortality by setting \code{msmMode} to 0. Predation mortality can be included by setting \code{msmMode} with the following options:
#' \itemize{
#' \item{0. Single species mode}
#' \item{1. Holsman et al. 2015 predation based on multi-species virtual population analysis (MSVPA) based predation formation.}
#' \item{2. MSVPA Holling Type III}
#' }
#'
#' Values 3 through 9 (Kinzey & Punt 2009 functional responses --
#' Holling Type I/II/III, predator interference, predator preemption,
#' Hassell-Varley, Ecosim) are blocked at runtime by \code{data_check()}
#' because the implementations have not been validated against the
#' current parameter set. See \code{src/TMB/predation.hpp}.
#'
#'
#' @return A list of class "Rceattle" including:
#'
#' \itemize{
#'  \item{data_list: List of data inputs}
#'  \item{initial_params: List of starting parameters}
#'  \item{bounds: Parameter bounds used for estimation}
#'  \item{map: List of map used in TMB}
#'  \item{obj: TMB model object}
#'  \item{opt: Optimized model object from `nlminb`}
#'  \item{sdrep: Object of class `sdreport` exported by TMB including the standard errors of estimated parameters}
#'  \item{estimated_params: List of estimated parameters}
#'  \item{quantities: Derived quantities from CEATTLE}
#'  \item{run_time: Model run time}
#'  }
#'
#'
#'
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(
#'   data_list    = BS2017SS,
#'   estimateMode = 0,
#'   msmMode      = 0,
#'   fit_control  = fit_control(phase = FALSE, verbose = 0)
#' )
#' }
#'
#' @export
fit_mod <-
  function(
    data_list = NULL,
    inits = NULL,
    map = NULL,
    bounds = NULL,
    file = NULL,
    estimateMode = 0,
    random_rec = FALSE,
    random_q = FALSE,
    random_sel = FALSE,
    HCR = build_hcr(),
    niter = 3,
    recFun = build_srr(),
    M1Fun = build_M1(),
    growthFun = build_growth(),
    qFun = build_catchability(),
    selFun = build_selectivity(),
    compFun = build_composition(),
    msmMode = 0,
    avgnMode = 0,
    initMode = "NonEquilibrium",
    suitMode = 0,
    suit_styr = NULL,
    suit_endyr = NULL,
    fit_control = NULL,
    config = NULL,
    ...){

    # Whether the caller supplied fit_control -- captured before the default is
    # filled below, since assigning to the argument makes missing() unreliable.
    # (The config= overlay needs this to know if it may fill fit_control.)
    fc_supplied <- !missing(fit_control)

    # Default control bundle. Built inside the body to avoid the
    # "recursive default argument reference" error that occurs when a
    # parameter and the function used in its default share a name.
    if (is.null(fit_control)) fit_control <- Rceattle::fit_control()

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Deprecation catch ----
    # Older versions of fit_mod() accepted the optimizer / sdreport /
    # phasing knobs as individual arguments. They now live exclusively
    # on `fit_control`. Catch the old names from `...`, warn, and
    # forward into the supplied fit_control so existing scripts keep
    # working.
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    .deprecated_ctl_args <- c(
      "phase", "getsd", "bias.correct", "use_gradient", "rel_tol",
      "control", "getJointPrecision", "getReportCovariance",
      "loopnum", "newtonsteps", "verbose", "TMBfilename",
      "projection_uncertainty"
    )
    .extra  <- list(...)
    .legacy <- intersect(names(.extra), .deprecated_ctl_args)
    if (length(.legacy)) {
      warning(
        sprintf(
          paste0(
            "Passing %s directly to fit_mod() is deprecated and will be ",
            "removed in a future release. Bundle these into fit_control() ",
            "instead, e.g. fit_control(%s). Forwarding for now."
          ),
          paste(sQuote(.legacy), collapse = ", "),
          paste(.legacy, "= ...", collapse = ", ")
        ),
        call. = FALSE
      )
      for (.nm in .legacy) {
        .field <- if (.nm == "control") "nlminb_control" else .nm
        fit_control[[.field]] <- .extra[[.nm]]
      }
      .extra <- .extra[setdiff(names(.extra), .deprecated_ctl_args)]
    }
    if (length(.extra)) {
      stop(
        "Unused arguments to fit_mod(): ",
        paste(sQuote(names(.extra)), collapse = ", "),
        call. = FALSE
      )
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Unpack fit_control ----
    # The optimizer / sdreport / phasing knobs all live in fit_control;
    # expand them to local variables so the rest of the function body
    # can refer to them by name without rewriting every reference.
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (!inherits(fit_control, "Rceattle_fit_control")) {
      stop("`fit_control` must be created with fit_control().")
    }
    bias.correct        <- fit_control$bias.correct
    getsd               <- fit_control$getsd
    getJointPrecision   <- fit_control$getJointPrecision
    getReportCovariance <- fit_control$getReportCovariance
    projection_uncertainty <- fit_control$projection_uncertainty
    use_gradient        <- fit_control$use_gradient
    rel_tol             <- fit_control$rel_tol
    loopnum             <- fit_control$loopnum
    newtonsteps         <- fit_control$newtonsteps
    phase               <- fit_control$phase
    TMBfilename         <- fit_control$TMBfilename
    verbose             <- fit_control$verbose
    control             <- fit_control$nlminb_control

    # ---------------------------------------------------------------------
    # Pipeline overview (file prefixes match this execution order):
    #   0-clean_data.R            clean_data() / 0-switches.R switch_check()
    #   1-data_check.R            data_check()      validate inputs
    #   2-build_params.R          build_params()    starting parameter list
    #   3-build_map.R             build_map()       TMB map (fixed vs estimated)
    #   4-build_parameter_bounds.R build_bounds()   lower/upper bounds
    #   5-rearrange_data.R        rearrange_data()  reshape data for TMB
    #   6-fit_mod.R               (this file)       MakeADFun + nlminb + sdreport
    #   6-rename_output.R         rename_output()   label derived quantities
    # HCR map (0-build_hcr.R build_hcr_map()) is applied during projection below.
    # ---------------------------------------------------------------------

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 0 - Start ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    mod_objects <- list()
    start_time  <- Sys.time()

    # Not every estimateMode / HCR combination runs an optimization, and the
    # reporting block below reads `opt$SD`.
    opt <- NULL

    extend_length <- function(x) {
      # A scalar recycles to every species; a length-nspp vector passes through.
      # A multi-element vector of any other length is a misconfiguration (it
      # would otherwise be silently re-recycled, e.g. length 3 with nspp 2 -> 6).
      if (length(x) > 1L && length(x) != data_list$nspp) {
        stop("Expected a scalar or a length-", data_list$nspp,
             " (nspp) vector; got length ", length(x), call. = FALSE)
      }
      if (length(x) == data_list$nspp) x else rep(x, data_list$nspp)
    }

    # Prefer function arg > existing data_list field > fallback year
    resolve_yr <- function(arg, data_val, fallback) {
      if (!is.null(arg)) arg else if (!is.null(data_val)) data_val else fallback
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 1 - Load data and switches ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

    data_list <- Rceattle::clean_data(data_list)

    # Overlay a run configuration (from load_config()) onto arguments the caller
    # did not supply. Its model_config structure is attached to the data list so
    # the model_config overlay below applies it to the omitted structure args;
    # its estimation controls fill the omitted scalar args. An explicitly-passed
    # argument always wins (missing() / fc_supplied is FALSE). config = NULL
    # (default) is a complete no-op, so the fit path is unchanged.
    if (!is.null(config)) {
      if (!inherits(config, "Rceattle_run_config"))
        stop("`config` must be an Rceattle_run_config (from load_config() / run_config()).",
             call. = FALSE)
      if (!is.null(config$model_config)) data_list$model_config <- config$model_config
      if (missing(estimateMode) && !is.null(config$estimateMode)) estimateMode <- config$estimateMode
      if (missing(random_rec)   && !is.null(config$random_rec))   random_rec   <- config$random_rec
      if (missing(random_q)     && !is.null(config$random_q))     random_q     <- config$random_q
      if (missing(random_sel)   && !is.null(config$random_sel))   random_sel   <- config$random_sel
      if (missing(suit_styr)    && !is.null(config$suit_styr))    suit_styr    <- config$suit_styr
      if (missing(suit_endyr)   && !is.null(config$suit_endyr))   suit_endyr   <- config$suit_endyr
      if (!fc_supplied          && !is.null(config$fit_control))  fit_control  <- config$fit_control
    }

    # Overlay a stored model_config (if any) onto arguments the caller did not
    # supply. An explicitly-passed argument always wins (missing() is FALSE); an
    # argument the caller omitted falls back to the stored field. With no
    # model_config slot this is a no-op, so fit_mod() behaves exactly as before.
    # Detected with missing() rather than a sentinel so 0 / FALSE / "" stored
    # values are honoured. NOTE: a call that passes one of these arguments --
    # even at its default -- overrides the slot; omit it to let the config win.
    cfg <- data_list$model_config
    if (!is.null(cfg)) {
      if (missing(msmMode)   && !is.null(cfg$msmMode))   msmMode   <- cfg$msmMode
      if (missing(initMode)  && !is.null(cfg$initMode))  initMode  <- cfg$initMode
      if (missing(avgnMode)  && !is.null(cfg$avgnMode))  avgnMode  <- cfg$avgnMode
      if (missing(suitMode)  && !is.null(cfg$suitMode))  suitMode  <- cfg$suitMode
      if (missing(niter)     && !is.null(cfg$niter))     niter     <- cfg$niter
      if (missing(HCR)       && !is.null(cfg$HCR))       HCR       <- cfg$HCR
      if (missing(recFun)    && !is.null(cfg$recFun))    recFun    <- cfg$recFun
      if (missing(M1Fun)     && !is.null(cfg$M1Fun))     M1Fun     <- cfg$M1Fun
      if (missing(growthFun) && !is.null(cfg$growthFun)) growthFun <- cfg$growthFun
      if (missing(qFun)      && !is.null(cfg$qFun))      qFun      <- cfg$qFun
      if (missing(selFun)    && !is.null(cfg$selFun))    selFun    <- cfg$selFun
      if (missing(compFun)   && !is.null(cfg$compFun))   compFun   <- cfg$compFun
    }

    # Resolve string aliases for the model-mode switches (e.g. "Hindcast" -> 1,
    # "SingleSpecies" -> 0); integers pass through unchanged and an unknown string
    # errors clearly. Done here, after the config overlay, so every downstream
    # integer test on estimateMode / msmMode still applies.
    estimateMode <- .map_switch(estimateMode, estimateMode_map, "estimateMode")
    msmMode      <- .map_switch(msmMode,      msmMode_map,      "msmMode")

    # Add switches from function call
    data_list$random_rec  <- as.numeric(random_rec)
    data_list$estimateMode <- estimateMode
    data_list$niter       <- niter
    data_list$avgnMode    <- avgnMode
    data_list$initMode    <- initMode
    data_list$loopnum     <- loopnum
    data_list$msmMode     <- msmMode
    data_list$suitMode    <- extend_length(suitMode)

    # * Suitability switches ----
    # suit_styr/suit_endyr are per-predator: a scalar is recycled to all species,
    # a length-nspp vector sets a distinct suitability-averaging window per predator.
    data_list$suit_styr  <- extend_length(resolve_yr(suit_styr,  data_list$suit_styr,  data_list$styr))
    data_list$suit_endyr <- extend_length(resolve_yr(suit_endyr, data_list$suit_endyr, data_list$endyr))


    # * Recruitment switches ----
    data_list$srr_fun       <- recFun$srr_fun
    data_list$srr_pred_fun  <- recFun$srr_pred_fun
    data_list$proj_mean_rec <- recFun$proj_mean_rec

    data_list$srr_mse_switchyr <- resolve_yr(recFun$srr_mse_switchyr, data_list$srr_mse_switchyr, data_list$endyr)
    data_list$srr_hat_styr     <- resolve_yr(recFun$srr_hat_styr,     data_list$srr_hat_styr,     data_list$styr + 1)
    data_list$srr_hat_endyr    <- resolve_yr(recFun$srr_hat_endyr,    data_list$srr_hat_endyr,    data_list$endyr)

    data_list$srr_est_mode <- recFun$srr_est_mode
    data_list$srr_prior    <- extend_length(recFun$srr_prior)
    data_list$srr_prior_sd <- extend_length(recFun$srr_prior_sd)
    data_list$srr_linkages <- recFun$linkages
    data_list$Bmsy_lim     <- extend_length(recFun$Bmsy_lim)

    # * M switches ----
    # Fold the deprecated `est_M1` into `M1_model` here -- BEFORE the M1Fun
    # reconciliation below -- so it gets the same treatment as `M1_model`
    # (the switch_check() alias runs later, after M1_model is already resolved
    # from M1Fun, so it cannot help the fit path).
    data_list <- .alias_est_M1(data_list)
    if (!is.null(data_list$M1_model)) {
      if (sum(data_list$M1_model != extend_length(M1Fun$M1_model))) {
        warning("M1_model in data is different than in call `fit_mod`, using switch from 'fit_mod'")
      }
    }

    # FIXME: may want to pull from data here too??
    data_list$M1_model     <- extend_length(M1Fun$M1_model)
    data_list$M1_model     <- ifelse(data_list$nsex == 1 & data_list$M1_model == 2, 1, data_list$M1_model) # sex-specific → sex-invariant for 1-sex model
    data_list$M1_re        <- extend_length(M1Fun$M1_re)
    updateM1               <- M1Fun$updateM1
    data_list$M1_use_prior <- extend_length(M1Fun$M1_use_prior) * (data_list$M1_model > 0) # 0 when M1 is fixed
    data_list$M2_use_prior <- extend_length(M1Fun$M2_use_prior) * (msmMode > 0)            # 0 in single-species mode
    data_list$M_prior      <- extend_length(M1Fun$M_prior)
    data_list$M_prior_sd   <- extend_length(M1Fun$M_prior_sd)
    data_list$M1_indices   <- M1Fun$M1_indices
    data_list$M1_linkages  <- M1Fun$linkages


    # * Growth switches ----
    data_list$growth_fun      <- growthFun$fun
    data_list$growth_linkages <- growthFun$linkages
    data_list$q_linkages      <- qFun$linkages
    data_list$sel_linkages    <- selFun$linkages
    data_list$comp_linkages   <- compFun$linkages
    data_list$growth_model    <- extend_length(growthFun$growth_model)
    data_list$growth_re       <- extend_length(growthFun$growth_re)
    # Plus-group SD-at-age style, resolved like growth_age_L1 below:
    # build_growth() value (if given) > existing data_list$growth_sd_style
    # (so a refit that rebuilds growth via build_growth(fun=) keeps the
    # original SS3/WHAM choice) > WHAM (1) fallback.
    gsd <- extend_length(growthFun$growth_sd_style)
    if (!is.null(data_list$growth_sd_style)) {
      gsd[is.na(gsd)] <- extend_length(data_list$growth_sd_style)[is.na(gsd)]
    }
    gsd[is.na(gsd)] <- 1L   # WHAM
    data_list$growth_sd_style <- gsd
    data_list$growth_indices  <- growthFun$growth_indices
    # VB anchor age per species (= SS3 Growth_Age_for_L1). Resolution
    # order: build_growth() user arg > data_list$growth_age_L1 (e.g. from
    # ss3_to_rceattle converter) > max(0.5, minage[sp]) fallback. The
    # fallback keeps minage >= 1 models backwards-compatible and gives
    # minage = 0 models an SS3-style half-year anchor.
    gal1 <- extend_length(growthFun$growth_age_L1)
    if (!is.null(data_list$growth_age_L1)) {
      gal1[is.na(gal1)] <- extend_length(data_list$growth_age_L1)[is.na(gal1)]
    }
    gal1[is.na(gal1)] <- pmax(0.5, as.numeric(extend_length(data_list$minage)))[is.na(gal1)]
    data_list$growth_age_L1 <- gal1


    # * HCR Switches ----
    data_list$HCR        <- HCR$HCR
    data_list$DynamicHCR <- HCR$DynamicHCR
    if (!HCR$HCR %in% c(2, "ConstantF")) { # Ftarget is also used for fixed F (so may be of length nflts)
      data_list$Ftarget <- extend_length(HCR$Ftarget)
    } else {
      data_list$Ftarget <- HCR$Ftarget
    }
    data_list$Flimit   <- extend_length(HCR$Flimit)
    data_list$Ptarget  <- extend_length(HCR$Ptarget)
    data_list$Plimit   <- extend_length(HCR$Plimit)
    data_list$Alpha    <- extend_length(HCR$Alpha)
    data_list$Pstar    <- extend_length(HCR$Pstar)
    data_list$Sigma    <- extend_length(HCR$Sigma)
    data_list$Fmult    <- extend_length(HCR$Fmult)
    data_list$HCRorder <- extend_length(HCR$HCRorder)
    data_list$QnormHCR <- stats::qnorm(data_list$Pstar, 0, data_list$Sigma)

    # Fill out switches if missing
    data_list <- Rceattle::switch_check(data_list)

    # Check for data error
    data_check(data_list)

    # * Pool process linkages into a global table + design matrix ----
    # No-op when no build_*() supplied a `linkages` list. When
    # linkages are present, this materializes each spec against
    # data_list$env_data and unions design columns by name.
    #
    # Linkages consume env_data POSITIONALLY: row r is applied to model year
    # styr + r - 1 (years beyond env_data get a zero offset). So a `Year`
    # column must be sorted, start at styr, and be contiguous, or a covariate /
    # deviate lands on the wrong year. Validate up front rather than misalign.
    .has_linkage <- any(vapply(
      list(data_list$growth_linkages, data_list$M1_linkages,
           data_list$srr_linkages, data_list$q_linkages, data_list$sel_linkages,
           data_list$comp_linkages),
      function(x) !is.null(x) && length(x) > 0L, logical(1)))
    if (.has_linkage) {
      # Prepend/gap-fill env_data to start at styr (NA for missing years) so a
      # later-starting `observe` covariate aligns; the check then validates the
      # extended table. NA years are skipped in the QAR1 observation (masked).
      data_list$env_data <- .extend_env_data(data_list$env_data, data_list$styr)
      .check_env_data_years(data_list$env_data, data_list$styr)
    }
    .message_auto_fleet_linkages(list(q   = data_list$q_linkages,
                                      sel = data_list$sel_linkages))
    .linkage_pool <- pool_linkages(
      spec_groups = list(growth      = data_list$growth_linkages,
                         M           = data_list$M1_linkages,
                         recruitment = data_list$srr_linkages,
                         q           = data_list$q_linkages,
                         sel         = data_list$sel_linkages,
                         comp        = data_list$comp_linkages),
      env_data    = data_list$env_data,
      # The fleet / species levels carry the model's own names, so a spec built
      # with `linkage_spec(fleet = "Shelikof")` can be resolved to an id.
      strata      = list(
        fleet   = .label_strata(seq_len(nrow(data_list$fleet_control)),
                                data_list$fleet_control$Fleet_name),
        species = .label_strata(seq_len(data_list$nspp), data_list$spnames),
        sex     = if (length(data_list$nsex) > 1L &&
                      length(unique(data_list$nsex)) > 1L) {
          stats::setNames(lapply(seq_len(data_list$nspp),
                                 function(sp) seq_len(data_list$nsex[sp])),
                          as.character(seq_len(data_list$nspp)))
        } else {
          seq_len(max(data_list$nsex))
        },
        age_bin = if (length(data_list$nages) > 1L &&
                      length(unique(data_list$nages)) > 1L) {
          stats::setNames(lapply(seq_len(data_list$nspp),
                                function(sp) seq_len(data_list$nages[sp])),
                         as.character(seq_len(data_list$nspp)))
        } else {
          seq_len(max(data_list$nages))
        }
      )
    )
    data_list$linkage_table <- .linkage_pool$table
    data_list$linkage_X     <- .linkage_pool$X
    # (Fixed-effect covariates with missing years are rejected earlier, in
    # materialize_linkage(), before model.matrix() can silently drop NA rows.)

    # Selectivity linkages are only consumed by the parametric forms wired in
    # the TMB template. Reject a sel linkage on a fleet whose selectivity form
    # is not yet wired, and the non-parametric `coff` param, so the effect is
    # never silently dropped. (Empirical and the RPM random walk cannot carry
    # a covariate offset at all.)
    .check_sel_linkage_support(data_list$linkage_table, data_list$fleet_control)
    .check_q_linkage_support(data_list$linkage_table, data_list$fleet_control)
    .check_comp_linkage_support(data_list$linkage_table, data_list)

    # Random-effect linkage rows (IID `~ (1 | group)`) are now consumed by the
    # TMB template: beta_linkage_re carries the deviations and jnll_comp row 20
    # the N(0, sigma) density. Correlated structures (rw()/ar1()) are still
    # rejected upstream in .materialize_re_design() until their densities land.

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 2: Load/build parameters ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.character(inits) | is.null(inits)) {
      start_par <- suppressWarnings(Rceattle::build_params(data_list = data_list))
    } else {
      start_par <- inits

      # Set F for years with 0 catch to very low number
      zero_catch <- data_list$catch_data |>
        dplyr::filter(.data$Year <= data_list$endyr &
                        .data$Catch == 0) |>
        dplyr::mutate(Year = Year - data_list$styr + 1) |>
        dplyr::select("Fleet_code", "Year") |>
        as.matrix()
      start_par$log_F[zero_catch] <- -999
      rm(zero_catch)

      # Update proj F prop
      start_par$proj_F_prop <- data_list$fleet_control$Proj_F_proportion
    }

    mod_objects$initial_params <- start_par
    if (verbose > 0) { message("Step 1: Parameter build complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 3: Load/build map ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(map)) {
      map <- suppressWarnings(build_map(data_list, start_par,
                                        debug = estimateMode %in% c(2, 4), # turn off hindcast parameters in projection / debug mode
                                        random_rec = random_rec, random_sel = random_sel))
    }
    if (verbose > 0) { message("Step 2: Map build complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 4: Get bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(bounds)) {
      bounds <- Rceattle::build_bounds(param_list = start_par, data_list)
    }
    if (verbose > 0) { message("Step 3: Parameter bounds complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 5: Setup random effects ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Turns on laplace approximation
    random_vars <- c()
    if (random_rec) {
      random_vars <- c(random_vars, "rec_dev")
      # init_dev is a random effect only when it is actually estimated. For
      # initMode = "Equilibrium" build_map() maps ALL of init_dev to NA (the
      # initial age structure is the deterministic equilibrium, init_dev fixed at
      # 0), so adding it to `random` would ask TMB to integrate a fully-mapped
      # parameter -- producing an NA/NaN gradient. Only treat it as random when
      # at least one element is free. build_map() returns a list with $mapList
      # (the raw per-parameter maps, NA = fixed); the init_dev map lives there.
      if (any(!is.na(map$mapList$init_dev))) {
        random_vars <- c(random_vars, "init_dev")
      }
    }
    if (random_q) {
      random_vars <- c(random_vars, "index_q_dev")
    }
    if (random_sel) {
      random_vars <- c(random_vars, "log_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
    }
    if (sum(data_list$M1_re) > 0) {
      random_vars <- c(random_vars, "log_M1_dev")
    }
    # Random-effect linkage deviations enter the Laplace approximation. Only add
    # them when at least one is actually free (mirrors the init_dev guard above):
    # integrating a fully-mapped / length-0 parameter yields an NA/NaN gradient.
    if (any(!is.na(map$mapList$beta_linkage_re))) {
      random_vars <- c(random_vars, "beta_linkage_re")
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 6: Reorganize data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (!is.null(TMBfilename)) {
      TMB::compile(paste0(TMBfilename, ".cpp"),
                   PKG_CXXFLAGS = "",
                   framework = "TMBad",
                   safebounds = FALSE, safeunload = FALSE)
      dyn.load(TMB::dynlib(TMBfilename))
      TMBfilename <- basename(TMBfilename)
    }
    if (is.null(TMBfilename)) {
      TMBfilename <- "ceattle"
    }

    # Composition proportion offset. It lives on `data_list` so that every
    # internal re-fit (projections, retrospective, jitter, run_mse, ...) inherits
    # the same value without threading it through each fit_mod() call. An explicit
    # fit_control(comp_offset=) overrides. Defaults to 1e-5.
    if (!is.null(fit_control$comp_offset)) data_list$comp_offset <- fit_control$comp_offset
    if (is.null(data_list$comp_offset))    data_list$comp_offset <- 1e-5

    # Bias-adjustment flags for lognormal likelihoods (default TRUE/1). Stored on
    # data_list as numeric (0/1) so the DATA_SCALAR is a clean double, so the TMB
    # template can toggle the -sigma^2/2 correction.
    data_list$bias_adjust_obs <- data_list$bias_adjust_proc <- 1
    # TODO(review): a user-supplied NA (e.g. fit_control(bias_adjust_obs = NA))
    # becomes NaN via as.numeric() and propagates a NaN objective. Add a finite
    # check to reject it with a clear error.
    if(!is.null(fit_control$bias_adjust_obs)) data_list$bias_adjust_obs <- as.numeric(fit_control$bias_adjust_obs)
    if(!is.null(fit_control$bias_adjust_proc)) data_list$bias_adjust_proc <- as.numeric(fit_control$bias_adjust_proc)

    # Reorganize data for .cpp file. The OSA observation vector is built on
    # demand by osa_residuals(), so only the fast aggregate metadata is assembled
    # here (build_osa = FALSE, the rearrange_data() default).
    data_list_reorganized <- Rceattle::rearrange_data(data_list)
    data_list_reorganized <- c(list(model = TMBfilename), data_list_reorganized)
    data_list_reorganized$forecast <- rep(0, data_list_reorganized$nspp) # hindcast switch

    # Scrub fields that TMB's dataSanitize cannot recurse into (data
    # frames with character columns, list-of-spec objects, character
    # vectors). They live on data_list for downstream R code; they
    # are re-injected below in TMB-friendly form.
    data_list_reorganized$linkage_table   <- NULL
    data_list_reorganized$growth_linkages <- NULL
    data_list_reorganized$growth_fun      <- NULL
    data_list_reorganized$M1_linkages     <- NULL
    data_list_reorganized$srr_linkages    <- NULL
    data_list_reorganized$q_linkages      <- NULL
    data_list_reorganized$sel_linkages    <- NULL
    data_list_reorganized$comp_linkages   <- NULL

    # OSA residual metadata: obs_ctl maps each obsvec element back to its
    # fleet/species/year/age. It is an R-side data frame (TMB's dataSanitize
    # cannot recurse into it), so stash it for the returned object and scrub it
    # from the TMB data list. The numeric obsvec / *_obsvec_idx / osa_mode
    # fields are TMB-friendly and stay.
    .obs_ctl <- data_list_reorganized$obs_ctl
    data_list_reorganized$obs_ctl <- NULL

    # * Inject linkage-table encoding into the TMB DATA ----
    # Empty when no build_*() supplied a `linkages` list. TMB's
    # DATA_MATRIX can be touchy about 0-dim shapes during ADFun
    # construction, so we substitute a single-cell sentinel when
    # there are no linkages -- the cpp loop short-circuits on
    # n_linkage == 0 without ever indexing into `linkage_X`.
    .linkage_table_now <- data_list$linkage_table %||% new_linkage_table()
    .linkage_X_now     <- if (nrow(.linkage_table_now) > 0L) {
      data_list$linkage_X
    } else {
      matrix(0, nrow = 1L, ncol = 1L)
    }
    .linkage_enc <- encode_linkage_for_tmb(
      table = .linkage_table_now,
      X     = .linkage_X_now
    )
    # The TMB cpp reads these by exact name; do not rename.
    data_list_reorganized$linkage_process      <- .linkage_enc$linkage_process
    data_list_reorganized$linkage_param        <- .linkage_enc$linkage_param
    data_list_reorganized$linkage_species      <- .linkage_enc$linkage_species
    data_list_reorganized$linkage_sex          <- .linkage_enc$linkage_sex
    data_list_reorganized$linkage_age_bin      <- .linkage_enc$linkage_age_bin
    data_list_reorganized$linkage_fleet        <- .linkage_enc$linkage_fleet
    data_list_reorganized$linkage_X_col        <- .linkage_enc$linkage_X_col
    data_list_reorganized$linkage_link         <- .linkage_enc$linkage_link
    data_list_reorganized$linkage_re_index     <- .linkage_enc$linkage_re_index
    data_list_reorganized$linkage_re_sigma      <- .linkage_enc$linkage_re_sigma
    data_list_reorganized$linkage_re_struct     <- .linkage_enc$linkage_re_struct
    data_list_reorganized$linkage_re_rho        <- .linkage_enc$linkage_re_rho
    data_list_reorganized$linkage_re_sigma_prior_family <- .linkage_enc$linkage_re_sigma_prior_family
    data_list_reorganized$linkage_re_sigma_prior_p1     <- .linkage_enc$linkage_re_sigma_prior_p1
    data_list_reorganized$linkage_re_sigma_prior_p2     <- .linkage_enc$linkage_re_sigma_prior_p2
    data_list_reorganized$linkage_re_rho_prior_family   <- .linkage_enc$linkage_re_rho_prior_family
    data_list_reorganized$linkage_re_rho_prior_p1       <- .linkage_enc$linkage_re_rho_prior_p1
    data_list_reorganized$linkage_re_rho_prior_p2       <- .linkage_enc$linkage_re_rho_prior_p2
    data_list_reorganized$linkage_re_obs        <- .linkage_enc$linkage_re_obs
    data_list_reorganized$linkage_re_obs_sd     <- .linkage_enc$linkage_re_obs_sd
    data_list_reorganized$linkage_re_obs_value  <- .linkage_enc$linkage_re_obs_value
    data_list_reorganized$linkage_re_obs_mask   <- .linkage_enc$linkage_re_obs_mask
    data_list_reorganized$linkage_is_intercept <- .linkage_enc$linkage_is_intercept
    data_list_reorganized$linkage_prior_family <- .linkage_enc$linkage_prior_family
    data_list_reorganized$linkage_prior_p1     <- .linkage_enc$linkage_prior_p1
    data_list_reorganized$linkage_prior_p2     <- .linkage_enc$linkage_prior_p2
    data_list_reorganized$linkage_X            <- .linkage_enc$linkage_X

    # Update comp weights, future F (if input), and F_prop from data
    # Age/length composition
    if (!is.null(data_list$fleet_control$Comp_weights)) {
      start_par$comp_weights <- data_list$fleet_control$Comp_weights
    }
    # CAAL
    if (!is.null(data_list$fleet_control$CAAL_weights)) {
      start_par$caal_weights <- data_list$fleet_control$CAAL_weights
    }
    # Diet composition
    if (!is.null(data_list$Diet_comp_weights)) {
      start_par$diet_comp_weights <- data_list$Diet_comp_weights
    }
    # Proportion of projected F to each fleet
    start_par$proj_F_prop <- data_list$fleet_control$Proj_F_proportion
    # Fixed fishing mortality for projections for each species
    if (!is.null(HCR$Ftarget) & HCR$HCR %in% c(2, "ConstantF")) {
      start_par$log_Ftarget <- log(HCR$Ftarget)
    }

    # Update M1 parameter object from data if initial parameter values input
    if (updateM1) {
      m1 <- array(0, dim = c(data_list$nspp,
                             max(data_list$nsex, na.rm = TRUE),
                             max(data_list$nages, na.rm = TRUE)))

      for (i in 1:nrow(data_list$M1_base)) {
        sp  <- as.numeric(as.character(data_list$M1_base$Species[i]))
        sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

        # Handle sex == 0 case for 2-sex species
        sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

        for (j in 1:length(sex_values)) {
          m1[sp, sex_values[j], 1:max(data_list$nages, na.rm = TRUE)] <- as.numeric(data_list$M1_base[i, (1:max(data_list$nages, na.rm = TRUE)) + 2])
        }
      }
      start_par$log_M1 <- log(m1)
    }

    # Update alpha for stock-recruit if fixed/prior and initial parameter values input
    if (data_list$srr_est_mode %in% c(0, 2) & data_list$srr_pred_fun > 3) {
      start_par$rec_pars[, 2] <- log(data_list$srr_prior)
    }

    if (verbose > 0) { message("Step 4: Data rearrange complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 7: Set up parameter bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    L <- c()
    U <- c()
    for (i in 1:length(map$mapFactor)) {
      nm <- names(map$mapFactor)[i]
      if (!nm %in% random_vars) { # no bounds for random effects
        mf   <- unlist(map$mapFactor[[i]])
        keep <- which(!is.na(mf) & !duplicated(mf))
        L <- c(L, unlist(bounds$lower[[nm]])[keep])
        U <- c(U, unlist(bounds$upper[[nm]])[keep])
      }
    }

    # Dimension check
    start_par <- start_par[names(map$mapFactor), drop = FALSE]
    dim_check <- sapply(start_par, function(x) length(unlist(x))) == sapply(map$mapFactor, function(x) length(unlist(x)))
    if (sum(dim_check) != length(dim_check)) {
      stop(paste0("Map and parameter objects are not the same size for: ", names(dim_check)[which(dim_check == FALSE)]))
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 8: Phase hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.logical(phase)) {
      if (phase) {
        phaseList <- set_phases()
      }
    }

    if (!is.logical(phase)) {
      warning("Using input phase. Please set phase = TRUE if using defaults.")
      phaseList <- phase
      phase     <- TRUE
    }


    step <- 5
    if (phase & estimateMode %in% c(0, 1)) {
      if (verbose > 0) { message(paste0("Step ", step, ": Phasing begin")) }; step <- step + 1
      phase_pars <- Rceattle::TMBphase(
        data         = data_list_reorganized,
        parameters   = start_par,
        map          = map$mapFactor,
        random       = random_vars,
        phases       = phaseList,
        model_name   = TMBfilename,
        silent       = verbose != 2,
        use_gradient = use_gradient,
        control      = control
      )

      # Pull the per-phase convergence log (attached by TMBphase) for the
      # phasing diagnostic, then strip it so it doesn't ride along in start_par.
      mod_objects$.conv_phase <- attr(phase_pars, "phase_log")
      attr(phase_pars, "phase_log") <- NULL

      mod_objects$phase_params <- phase_pars
      start_par <- phase_pars

      if (verbose > 0) { message(paste0("Step ", step, ": Phasing complete")) }
      step <- step + 1
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 9: Fit hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # * Build ----
    if (sum(as.numeric(unlist(map$mapFactor)), na.rm = TRUE) == 0) { stop("Map of length 0: all NAs") }
    obj <- TMB::MakeADFun(
      data_list_reorganized,
      parameters = start_par,
      DLL        = TMBfilename,
      map        = map$mapFactor,
      random     = random_vars,
      silent     = verbose != 2
    )

    mod_objects <- c(
      list(
        TMBfilename = TMBfilename,
        bounds      = bounds,
        map         = map
      ),
      mod_objects)

    if (verbose > 0) { message(paste0("Step ", step, ": Hindcast build complete")) }
    step <- step + 1


    # * Optimize hindcast ----
    if (estimateMode %in% c(0, 1, 4)) {
      opt <- suppressMessages(
        .fit_tmb(obj                 = obj,
                 lower               = L,
                 upper               = U,
                 loopnum             = loopnum,
                 newtonsteps         = newtonsteps,
                 getsd               = getsd,
                 control             = control,
                 bias.correct        = bias.correct,
                 getJointPrecision   = getJointPrecision,
                 getReportCovariance = getReportCovariance,
                 quiet               = verbose < 2)
      )

      if (verbose > 0 & estimateMode != 4) {
        message("Step ", step, ": Hindcast optimization complete")
        step <- step + 1
      }
      if (verbose > 0 & estimateMode == 4) {
        message("Step ", step, ": 'dummy' optimization complete")
        step <- step + 1
      }

      # Convergence warnings
      if (estimateMode %in% c(0, 1)) {
        if (is.null(opt$SD) & getsd) {

          message("#################################################")
          message("Model did not converge, check 'identified'")
          message("#################################################")

          identified <- tryCatch(
            { suppressMessages(.check_estimability(obj)) },
            error = function(e) {
              "Some gradients are high, please improve optimization and only then use `Check_Identifiable`"
            }
          )

          # Make into list if gradients were low for diagnostics
          if (!is.character(identified)) {
            identified_param_list <- obj$env$parList(identified$BadParams$Param_check)
            identified_param_list <- rapply(identified_param_list, function(x) ifelse(x == 0, "Not estimated", x), how = "replace")
            identified_param_list <- rapply(identified_param_list, function(x) ifelse(x == 1, "OK",            x), how = "replace")
            identified_param_list <- rapply(identified_param_list, function(x) ifelse(x == 2, "BAD",           x), how = "replace")
            identified$param_list <- identified_param_list
          }
          mod_objects$identified <- identified
        }
      }

      # Capture the hindcast optimizer convergence snapshot now, before any
      # projection re-optimization overwrites `opt` (see R/0-convergence.R).
      if (estimateMode %in% c(0, 1)) {
        mod_objects$.conv_hindcast <- .capture_opt_convergence(
          opt, obj, bounds = bounds, mapFactor = map$mapFactor,
          random_vars = random_vars, getsd = getsd)
      }
    }


    # * Get MLEs ----
    if (estimateMode > 1) { # debugging / projection-only: use initial parameters
      last_par <- start_par
    } else {
      last_par <- obj$env$parList(par=obj$env$last.par.best)
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 10: Run projection ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (estimateMode %in% c(0, 2, 4)) {
      if (data_list$HCR != "NoFishing") { # all HCRs except no F and fixed F

        # * 10.1: Single species mode ----
        if (msmMode == 0) {

          # Turn BRP estimation on within likelihood
          data_list_reorganized$forecast <- rep(1, data_list_reorganized$nspp)

          hcr_map <- Rceattle::build_hcr_map(data_list, map, debug = estimateMode > 3)
          if (sum(!is.na(unlist(hcr_map$mapFactor))) == 0) { stop("HCR map of length 0: all NAs") }

          obj <- TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL        = TMBfilename,
            map        = hcr_map$mapFactor,
            random     = random_vars,
            silent     = verbose != 2
          )

          if (verbose > 0) { message(paste0("Step ", step, ": Projection build complete")) }
          step <- step + 1

          if (data_list$HCR != "ConstantF") { # fixed F does not need estimation
            opt <- suppressMessages(
              .fit_tmb(obj               = obj,
                       loopnum           = loopnum,
                       getsd             = getsd,
                       control           = control,
                       bias.correct      = bias.correct,
                       getJointPrecision = FALSE,
                       quiet             = verbose < 2)
            )
          }
        }


        # * 10.2: Multi-species mode ----
        if (msmMode > 0) {

          for (HCRiter in 1:max(data_list$HCRorder)) {

            hcr_map <- Rceattle::build_hcr_map(data_list, map,
                                               debug         = estimateMode > 3,
                                               all_params_on = FALSE,
                                               HCRiter       = HCRiter)
            if (sum(as.numeric(unlist(hcr_map$mapFactor)), na.rm = TRUE) == 0) { stop("HCR map of length 0: all NAs") }

            # Get SB0: SSB when model is projected forward under no fishing
            data_list_reorganized$forecast <- data_list$HCRorder <= HCRiter # species to include in likelihood
            params_on  <- c(1:data_list$nspp)[which(data_list$HCRorder == HCRiter)] # only update the species being estimated
            quantities <- obj$report(obj$env$last.par.best)
            SB0        <- quantities$ssb[, ncol(quantities$ssb)]
            B0         <- quantities$biomass[, ncol(quantities$biomass)]
            data_list_reorganized$MSSB0[params_on] <- SB0[params_on]
            data_list_reorganized$MSB0[params_on]  <- B0[params_on]

            # Adjust Ftarget inits
            params_off <- c(1:data_list$nspp)[which(data_list$HCRorder > HCRiter)]
            last_par$log_Ftarget[params_on]  <- 0
            last_par$log_Ftarget[params_off] <- -999

            obj <- TMB::MakeADFun(
              data_list_reorganized,
              parameters = last_par,
              DLL        = TMBfilename,
              map        = hcr_map$mapFactor,
              random     = random_vars,
              silent     = verbose != 2
            )

            if (verbose > 0) { message(paste0("Step ", step, " - HCRiter ", HCRiter, ": Projection build complete. Optimizing.")) }
            step <- step + 1

            if (data_list$HCR != "ConstantF") { # fixed F does not need estimation
              opt <- suppressMessages(
                .fit_tmb(obj               = obj,
                         loopnum           = loopnum,
                         getsd             = getsd,
                         bias.correct      = bias.correct,
                         control           = control,
                         getJointPrecision = FALSE,
                         quiet             = verbose < 2)
              )

              # Update F from opt
              last_par$log_Ftarget[params_on] <- opt$par[1:length(params_on)]
            }
          }
        }


        if (verbose > 0) { message(paste0("Step ", step, ": Projection optimization complete")) }
        step <- step + 1


        # * Update MLEs ----
        if (estimateMode > 2) { # debugging: use initial parameters
          last_par <- start_par
        } else {
          last_par <- obj$env$parList(par=obj$env$last.par.best)
        }


        # * Projection uncertainty ----
        # Updates the model with all hindcast and BRP parameters "turned on"
        # to get uncertainty estimates in the projection.
        if (projection_uncertainty) {

          hcr_map_proj <- Rceattle::build_hcr_map(data_list, map, debug = estimateMode > 3, all_params_on = TRUE)
          if (sum(as.numeric(unlist(hcr_map_proj$mapFactor)), na.rm = TRUE) == 0) { stop("HCR projection map of length 0: all NAs") }

          obj <- TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL        = TMBfilename,
            map        = hcr_map_proj$mapFactor,
            random     = random_vars,
            silent     = verbose != 2
          )

          opt$SD <- TMB::sdreport(obj)
        }

      } # End estimable BRP/HCR projections
    } # End projection


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 11: Save output ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    mod_objects$estimated_params <- last_par
    mod_objects$obj              <- obj

    quantities <- obj$report(obj$env$last.par.best)

    # Warning for discontinuous likelihood
    if (estimateMode %in% c(0:2)) {
      if (!(estimateMode == 2 & data_list$HCR == "ConstantF")) { # no optimization of projections with fixed F
        if (!is.null(opt) && !is.null(opt$SD) && random_rec == FALSE) {
          if (abs(opt$objective - quantities$jnll) > rel_tol) {
            message("#################################################")
            message("Convergence warning (8): discontinuous likelihood")
            message("#################################################")
          }
        }
      }
    }

    mod_objects$quantities <- Rceattle::rename_output(data_list = data_list, quantities = quantities)

    mod_objects$data_list <- calc_mcall_ianelli(data_list = data_list, data_list_reorganized = data_list_reorganized, quantities = quantities)
    mod_objects$data_list <- calc_mcall_ianelli_diet(data_list = mod_objects$data_list, quantities = quantities)

    # OSA residual metadata for the aggregate (index/catch) series, mapping
    # obsvec positions to fleet/species/year/age. osa_residuals() rebuilds the
    # full composition / caal / diet metadata on demand, so only the aggregate
    # map is stored here.
    mod_objects$obs_ctl <- .obs_ctl

    mod_objects$run_time <- (Sys.time() - start_time)

    # Record the resolved run configuration so save_config(fit) / run_config(fit)
    # reproduce this fit without re-deriving it. Built from the arguments as
    # actually used (after any config= / model_config overlay).
    mod_objects$run_config <- .rce_run_config(
      mc = model_config(msmMode = msmMode, initMode = initMode, avgnMode = avgnMode,
                        suitMode = suitMode, niter = niter, HCR = HCR, recFun = recFun,
                        M1Fun = M1Fun, growthFun = growthFun, qFun = qFun,
                        selFun = selFun, compFun = compFun),
      estimateMode = estimateMode, random_rec = random_rec, random_q = random_q,
      random_sel = random_sel, suit_styr = suit_styr, suit_endyr = suit_endyr,
      fc = fit_control)

    if (estimateMode < 3) {
      if (!(estimateMode == 2 & data_list$HCR == "ConstantF")) { # no optimization of projections with fixed F
        mod_objects$opt   <- opt
        mod_objects$sdrep <- opt$SD
      }
    }

    class(mod_objects) <- "Rceattle"

    # Convergence diagnostics (R/0-convergence.R): attach the structured object
    # and surface non-OK checks via message(). message() + tryCatch so a
    # non-converged fit is never turned into an error and is always returned.
    if (estimateMode %in% c(0, 1)) {
      mod_objects$convergence <- tryCatch(
        convergence_diagnostics(mod_objects),
        error = function(e) NULL
      )
      tryCatch({
        if (!is.null(mod_objects$convergence) &&
            mod_objects$convergence$status %in% c("WARN", "FAIL")) {
          message("Convergence diagnostics flagged (status: ",
                  mod_objects$convergence$status,
                  "); inspect fit$convergence.")
          for (ch in mod_objects$convergence$checks) {
            if (ch$severity %in% c("WARN", "FAIL")) {
              message(sprintf("  [%s] %s: %s", ch$severity, ch$id, ch$message))
            }
          }
        }
      }, error = function(e) NULL)
    }

    if (!is.null(file)) {
      save(mod_objects, file = paste0(file, ".Rdata"))
    }

    # Free up memory
    if (estimateMode %in% 0:1) {
      TMB::FreeADFun(obj)
    }

    return(mod_objects)
  }
