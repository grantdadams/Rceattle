#' This function runs CEATTLE
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
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
#' @param estimateMode 0 = Fit the hindcast model and projection with HCR specified via \code{HCR}. 1 = Fit the hindcast model only (no projection). 2 = Run the projection only with HCR specified via \code{HCR} given the initial parameters in \code{inits}.  3 = debug mode 1: runs the model through MakeADFun, but not nlminb, 4 = runs the model through MakeADFun and nlminb (will all parameters mapped out).
#' @param projection_uncertainty logical. If TRUE, accounts for hindcast parameter uncertainty in projections when using an HCR. Default is FALSE for speed.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects using the laplace approximation.The default is FALSE.
#' @param random_q logical. If TRUE, treats annual catchability deviations as random effects using the laplace approximation.The default is FALSE.
#' @param random_sel logical. If TRUE, treats annual selectivity deviations as random effects using the laplace approximation.The default is FALSE.
#' @param HCR HCR list object from \code{\link{build_hcr}}
#' @param niter Number of iterations for multispecies model
#' @param recFun The stock recruit-relationship parameterization from \code{\link{build_srr}}.
#' @param M1Fun M1 parameterizations and priors. Use \code{build_M1}.
#' @param growthFun The weight-at-age parameterization from \code{\link{build_growth}}.
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{N/Z ( 1 - exp(-Z) )}, 1 is \eqn{N exp(-Z/2)}, 2 is \eqn{N}.
#' @param initMode how the population is initialized. 0 = initial age-structure estimated as free parameters; 1 = equilibrium age-structure estimated out from R0 + mortality (M1); 2 = non-equilibrium age-structure estimated out from R0,  mortality (M1), and initial population deviates; 3 = non-equilibrium age-structure estimated out from initial fishing mortality (Finit), R0,  mortality (M1), and initial population deviates; 4 = non-equilibrium age-structure version 2 where initial fishing mortality (Finit) scales R0.
#' @param suitMode Switch for suitability derivation for each predator (single value or vector). 0 = empirical based on diet data (Holsman et al. 2015), 1 = length-based gamma suitability, 2 = weight-based gamma suitability, 3 = length-based lognormal suitability, 4 = weight-based lognormal suitability, 5 = length-based normal suitability, 6 = weight-based normal suitability.
#' @param suit_styr Integer. The first year used to calculate mean suitability. Defaults to $styr$ in $data_list$. Used when diet data were sampled from a subset of years.
#' @param suit_endyr Integer. The last year used to calculate mean suitability. Defaults to $endyr$ in $data_list$. Used when diet data were sampled from a subset of years.
#' @param fit_control A list returned by [fit_control()] that bundles the
#'   optimizer / sdreport / phasing knobs (`phase`, `bias.correct`, `getsd`,
#'   `getJointPrecision`, `getReportCovariance`, `use_gradient`, `rel_tol`,
#'   `loopnum`, `newtonsteps`, `TMBfilename`, `verbose`, `nlminb_control`).
#'   Defaults to `fit_control()`. See [fit_control()] for the meaning and
#'   defaults of each field.
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
    projection_uncertainty = FALSE,
    random_rec = FALSE,
    random_q = FALSE,
    random_sel = FALSE,
    HCR = build_hcr(),
    niter = 3,
    recFun = build_srr(),
    M1Fun = build_M1(),
    growthFun = build_growth(),
    msmMode = 0,
    avgnMode = 0,
    initMode = "NonEquilibrium",
    suitMode = 0,
    suit_styr = NULL,
    suit_endyr = NULL,
    fit_control = NULL,
    ...){

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
      "loopnum", "newtonsteps", "verbose", "TMBfilename"
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
    use_gradient        <- fit_control$use_gradient
    rel_tol             <- fit_control$rel_tol
    loopnum             <- fit_control$loopnum
    newtonsteps         <- fit_control$newtonsteps
    phase               <- fit_control$phase
    TMBfilename         <- fit_control$TMBfilename
    verbose             <- fit_control$verbose
    control             <- fit_control$nlminb_control

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 0 - Start ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    mod_objects <- list()
    start_time  <- Sys.time()

    extend_length <- function(x) {
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
    data_list$suit_styr  <- resolve_yr(suit_styr,  data_list$suit_styr,  data_list$styr)
    data_list$suit_endyr <- resolve_yr(suit_endyr, data_list$suit_endyr, data_list$endyr)


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
    data_list$growth_model    <- extend_length(growthFun$growth_model)
    data_list$growth_re       <- extend_length(growthFun$growth_re)
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
    # data_list$env_data and unions design columns by name. The
    # resulting `linkage_table` and `linkage_X` are used downstream
    # once the TMB template knows how to consume them; for now they
    # are computed (and validated) but not yet plumbed into TMB.
    .linkage_pool <- pool_linkages(
      spec_groups = list(growth      = data_list$growth_linkages,
                         M           = data_list$M1_linkages,
                         recruitment = data_list$srr_linkages),
      env_data    = data_list$env_data,
      strata      = list(
        species = seq_len(data_list$nspp),
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
      start_par$proj_F_prop <- data_list$fleet_control$proj_F_prop
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
      if (initMode > 0) {
        random_vars <- c(random_vars, "rec_dev", "init_dev")
      } else {
        random_vars <- c(random_vars, "rec_dev")
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
      TMBfilename <- "ceattle_v01_11"
    }

    # Reorganize data for .cpp file
    data_list_reorganized <- Rceattle::rearrange_dat(data_list)
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
    data_list_reorganized$linkage_X_col        <- .linkage_enc$linkage_X_col
    data_list_reorganized$linkage_link         <- .linkage_enc$linkage_link
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
    start_par$proj_F_prop <- data_list$fleet_control$proj_F_prop
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
      if (!names(map$mapFactor)[i] %in% random_vars) { # no bounds for random effects
        L <- c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
        U <- c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
      }
    }

    # Dimension check
    start_par <- start_par[names(map$mapFactor)]
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
        if (!is.null(opt$SD) & random_rec == FALSE) {
          if (abs(opt$objective - quantities$jnll) > rel_tol) {
            message("#################################################")
            message("Convergence warning (8): discontinuous likelihood")
            message("#################################################")
          }
        }
      }
    }

    mod_objects$quantities <- Rceattle::rename_output(data_list = data_list, quantities = quantities)

    mod_objects$data_list <- Rceattle::calc_mcall_ianelli(data_list = data_list, data_list_reorganized = data_list_reorganized, quantities = quantities)
    mod_objects$data_list <- Rceattle::calc_mcall_ianelli_diet(data_list = mod_objects$data_list, quantities = quantities)

    mod_objects$run_time <- (Sys.time() - start_time)

    if (estimateMode < 3) {
      if (!(estimateMode == 2 & data_list$HCR == "ConstantF")) { # no optimization of projections with fixed F
        mod_objects$opt   <- opt
        mod_objects$sdrep <- opt$SD
      }
    }

    class(mod_objects) <- "Rceattle"

    if (!is.null(file)) {
      save(mod_objects, file = paste0(file, ".Rdata"))
    }

    # Free up memory
    if (estimateMode %in% 0:1) {
      TMB::FreeADFun(obj)
    }

    return(mod_objects)
  }
