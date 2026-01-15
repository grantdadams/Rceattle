#' This functions runs CEATTLE
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#'
#' @param data_list a data_list read in using \code{\link{read_excel}}.
#' @param inits (Optional) Character vector of named initial values from previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters. Can also construct using \code{\link{build_params}}
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param bounds (Optional) A bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param estimateMode 0 = Fit the hindcast model and projection with HCR specified via \code{HCR}. 1 = Fit the hindcast model only (no projection). 2 = Run the projection only with HCR specified via \code{HCR} given the initial parameters in \code{inits}.  3 = debug mode 1: runs the model through MakeADFun, but not nlminb, 4 = runs the model through MakeADFun and nlminb (will all parameters mapped out).
#' @param projection_uncertainty account for hindcast parameter uncertainty in projections when using an HCR? Default is FALSE for speed.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects using the laplace approximation.The default is FALSE.
#' @param random_q logical. If TRUE, treats annual catchability deviations as random effects using the laplace approximation.The default is FALSE.
#' @param random_sel logical. If TRUE, treats annual selectivity deviations as random effects using the laplace approximation.The default is FALSE.
#' @param HCR HCR list object from \code{\link{build_hcr}}
#' @param niter Number of iterations for multispecies model
#' @param recFun The stock recruit-relationship parameterization from \code{\link{build_srr}}.
#' @param M1Fun M1 parameterizations and priors. Use \code{build_M1}.
#' @param growthFun The weight-at-age model from \code{\link{build_growth}}.
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{N/Z ( 1 - exp(-Z) )}, 1 is \eqn{N exp(-Z/2)}, 2 is \eqn{N}.
#' @param initMode how the population is initialized. 0 = initial age-structure estimated as free parameters; 1 = equilibrium age-structure estimated out from R0 + dev-yr1,  mortality (M1); 2 = equilibrium age-structure estimated out from R0,  mortality (M1), and initial population deviates; 3 = non-equilibrium age-structure estimated out from initial fishing mortality (Finit), R0,  mortality (M1), and initial population deviates; 4 = non-equilibrium age-structure version 2 where initial fishing mortality (Finit) scales R0.
#' @param phase TRUE/FALSE If FALSE, will not phase model. If set to \code{"TRUE"}, will use default phasing. Can also accept a list of parameter object names with corresponding phase. See https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R.
#' @param suitMode Switch for suitability derivation for each predator (single value or vector). 0 = empirical based on diet data (Holsman et al. 2015), 1 = length-based gamma suitability, 2 = weight-based gamma suitability, 3 = length-based lognormal suitability, 4 = weight-based lognormal suitability, 5 = length-based normal suitability, 6 = weight-based normal suitability.
#' @param suit_styr Integer. The first year used to calculate mean suitability. Defaults to $styr$ in $data_list$. Used when diet data were sampled from a subset of years.
#' @param suit_endyr Integer. The last year used to calculate mean suitability. Defaults to $endyr$ in $data_list$. Used when diet data were sampled from a subset of years.
#' @param getsd	TRUE/FALSE whether to run standard error calculation (default = TRUE).
#' @param use_gradient use the gradient to phase (default = TRUE).
#' @param rel_tol The relative tolerance for discontinuous likelihood warnings. Set to 1. This evaluates the difference between the TMB object likelihood and the nlminb likelihood.
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param verbose 0 = Silent, 1 = print updates of model fit, 2 = print updates of model fit and TMB estimation progress.
#' @param getJointPrecision return full Hessian of fixed and random effects.
#' @param getReportCovariance return variance covariance of ADREPORT variables
#' @param TMBfilename if a seperate TMB file is to be used for development. Includes location and does not include ".cpp" at the end.
#'
#' @details
#' CEATTLE is an age-structured population dynamics model that can be fit with or without predation mortality. The default is to exclude predation mortality by setting \code{msmMode} to 0. Predation mortality can be included by setting \code{msmMode} with the following options:
#' \itemize{
#' \item{0. Single species mode}
#' \item{1. Holsman et al. 2015 predation based on multi-species virtual population analysis (MSVPA) based predation formation.}
#' \item{2. MSVPA Holling Type III}
#'  \item{3. Kinzey & Punt 2010 Holling Type I (linear) - Deprecated}
#'  \item{4. Kinzey & Punt 2010 Holling Type II - Deprecated}
#'  \item{5. Kinzey & Punt 2010 Holling Type III - Deprecated}
#'  \item{6. Kinzey & Punt 2010 Predator interference - Deprecated}
#'  \item{7. Kinzey & Punt 2010 Predator preemption - Deprecated}
#'  \item{8. Kinzey & Punt 2010 Hassell-Varley - Deprecated}
#'  \item{9. Kinzey & Punt 2010 Ecosim - Deprecated}
#'  }
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
#'  \item{opt: Optimized model object from `nlimb`}
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
#'# Load package and data
#'library(Rceattle)
#'data(BS2017SS) # ?BS2017SS for more information on the data
#'
#'
#'# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
#'ss_run <- fit_mod(data_list = BS2017SS,
#'    inits = NULL, # Initial parameters = 0
#'    file = NULL, # Don't save
#'    estimateMode = 0, # Estimate
#'    random_rec = FALSE, # No random recruitment
#'    msmMode = 0, # Single species mode
#'    avgnMode = 0,
#'    phase = FALSE,
#'    silent = TRUE)
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
    initMode = 2,
    suitMode = 0,
    suit_styr = NULL,
    suit_endyr = NULL,
    phase = FALSE,
    getsd = TRUE,
    bias.correct = FALSE,
    use_gradient = TRUE,
    rel_tol = 1,
    control = list(eval.max = 1e+09,
                   iter.max = 1e+09, trace = 0),
    getJointPrecision = TRUE,
    getReportCovariance = FALSE,
    loopnum = 5,
    verbose = 1,
    newtonsteps = 0,
    catch_hcr = FALSE,
    TMBfilename = NULL){

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Debugging section ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # data_list = GOA2018SS;
    # inits = NULL;
    # map = NULL;
    # bounds = NULL;
    # file = NULL;
    # estimateMode = 0;
    # random_rec = FALSE;
    # random_q = FALSE;
    # random_sel = FALSE;
    # HCR = build_hcr();
    # niter = 3;
    # msmMode = 0;
    # avgnMode = 0;
    # initMode = 2
    # minNByage = 0;
    # suitMode = 0;
    # suit_styr = NULL;
    # suit_endyr = NULL;
    # phase = FALSE;
    # getsd = TRUE;
    # use_gradient = TRUE;
    # rel_tol = 1;
    # control = list(eval.max = 1e+09,
    #                iter.max = 1e+09, trace = 0);
    # getJointPrecision = TRUE;
    # loopnum = 5;
    # verbose = 1;
    # newtonsteps = 0
    # recFun = build_srr()
    # M1Fun = build_M1()
    # growthFun = build_growth()
    # projection_uncertainty = TRUE
    # catch_hcr = FALSE
    # bias.correct = FALSE
    # newtonsteps = 0
    # getReportCovariance = FALSE

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 0 - Start ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    mod_objects <- list() # Objects for saving
    start_time <- Sys.time()

    extend_length <- function(x){
      if(length(x) == data_list$nspp){ return(x)}
      else {return(rep(x, data_list$nspp))}
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 1 - Load data and switches ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

    data_list <- Rceattle::clean_data(data_list)

    # Add switches from function call
    data_list$random_rec <- as.numeric(random_rec)
    data_list$estimateMode <- estimateMode
    data_list$niter <- niter
    data_list$avgnMode <- avgnMode
    data_list$initMode <- initMode
    data_list$loopnum <- loopnum
    data_list$msmMode <- msmMode
    data_list$suitMode <- extend_length(suitMode)

    # * Suitability switches ----
    # -- Start year
    if(is.null(suit_styr) & is.null(data_list$suit_styr)){ # If not provided in data or function, use start year
      data_list$suit_styr <- data_list$styr
    }
    if(!is.null(suit_styr)){ # If provided in function, override data
      data_list$suit_styr <- suit_styr
    }

    # -- End year
    if(is.null(suit_endyr) & is.null(data_list$suit_endyr)){ # If not provided in data or function, use end year
      data_list$suit_endyr <- data_list$endyr
    }
    if(!is.null(suit_endyr)){ # If provided in function, override data
      data_list$suit_endyr <- suit_endyr
    }


    # * Recruitment switches ----
    data_list$srr_fun <- recFun$srr_fun
    data_list$srr_pred_fun <- recFun$srr_pred_fun
    data_list$proj_mean_rec <- recFun$proj_mean_rec
    if(is.null(recFun$srr_mse_switchyr) & is.null(data_list$srr_mse_switchyr)){ # If no meanyear is provided in data or function, use end year
      data_list$srr_mse_switchyr <- data_list$endyr
    }
    if(!is.null(recFun$srr_mse_switchyr)){ # If mean year is provided in function, override data
      data_list$srr_mse_switchyr <- recFun$srr_mse_switchyr
    }

    # -- Start year
    if(is.null(recFun$srr_hat_styr) & is.null(data_list$srr_hat_styr)){ # If not provided in data or function, use start year
      data_list$srr_hat_styr <- data_list$styr + 1
    }
    if(!is.null(recFun$srr_hat_styr)){ # If provided in function, override data
      data_list$srr_hat_styr <- recFun$srr_hat_styr
    }

    # -- End year
    if(is.null(recFun$srr_hat_endyr) & is.null(data_list$srr_hat_endyr)){ # If not provided in data or function, use end year
      data_list$srr_hat_endyr <- data_list$endyr
    }
    if(!is.null(recFun$srr_hat_endyr)){ # If provided in function, override data
      data_list$srr_hat_endyr <- recFun$srr_hat_endyr
    }

    data_list$srr_est_mode <- recFun$srr_est_mode
    data_list$srr_prior <- extend_length(recFun$srr_prior)
    data_list$srr_prior_sd <- extend_length(recFun$srr_prior_sd)
    data_list$srr_indices <- recFun$srr_indices
    data_list$Bmsy_lim <- extend_length(recFun$Bmsy_lim)

    # * M switches ----
    if(!is.null(data_list$M1_model)){
      if(sum(data_list$M1_model != extend_length(M1Fun$M1_model))){
        warning("M1_model in data is different than in call `fit_mod`, using switch from 'fit_mod'")
      }
    }

    # FIXME: may want to pull from data here too??
    data_list$M1_model= extend_length(M1Fun$M1_model)
    data_list$M1_model = ifelse(data_list$nsex == 1 & data_list$M1_model == 2, 1, data_list$M1_model) # Sex specific to sex-invariant if 1-sex model
    data_list$M1_re = extend_length(M1Fun$M1_re)
    updateM1 = M1Fun$updateM1
    data_list$M1_use_prior = extend_length(M1Fun$M1_use_prior) * (data_list$M1_model > 0) # Sets to 0 if M1 is fixed
    data_list$M2_use_prior = extend_length(M1Fun$M2_use_prior) * (msmMode > 0) # Sets to 0 if single-species
    data_list$M_prior = extend_length(M1Fun$M_prior)
    data_list$M_prior_sd = extend_length(M1Fun$M_prior_sd)
    data_list$M1_indices <- M1Fun$M1_indices


    # * Growth switches ----
    data_list$growth_model= extend_length(growthFun$growth_model)
    data_list$growth_re = extend_length(growthFun$growth_re)
    data_list$growth_indices = growthFun$growth_indices


    # * HCR Switches ----
    # - make length of nspp if not
    data_list$HCR = HCR$HCR
    data_list$DynamicHCR = HCR$DynamicHCR
    if(HCR$HCR != 2){ # Ftarget is also used for fixed F (so may be of length nflts)
      data_list$Ftarget = extend_length(HCR$Ftarget)
    } else {
      data_list$Ftarget = HCR$Ftarget
    }
    data_list$Flimit = extend_length(HCR$Flimit)
    data_list$Ptarget = extend_length(HCR$Ptarget)
    data_list$Plimit = extend_length(HCR$Plimit)
    data_list$Alpha = extend_length(HCR$Alpha)
    data_list$Pstar = extend_length(HCR$Pstar)
    data_list$Sigma = extend_length(HCR$Sigma)
    data_list$Fmult = extend_length(HCR$Fmult)
    data_list$HCRorder = extend_length(HCR$HCRorder)
    data_list$QnormHCR = qnorm(data_list$Pstar, 0, data_list$Sigma)

    # if(data_list$HCR == 2 & estimateMode == 2){estimateMode = 4} # If projecting under constant F, run parmeters through obj only

    if(data_list$msmMode > 0 & !data_list$HCR %in% c(0, 1, 2, 3, 6)){
      warning("WARNING:: Only HCRs 1, 2, 3, and 6 work in multi-species mode currently")
    }

    # Fill out switches if missing
    data_list <- Rceattle::switch_check(data_list)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 2: Load/build parameters ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.character(inits) | is.null(inits)) {
      start_par <- suppressWarnings(Rceattle::build_params(data_list = data_list))
    } else{
      start_par <- inits

      # - Set F for years with 0 catch to very low number
      zero_catch <- data_list$catch_data %>%
        dplyr::filter(Year <= data_list$endyr &
                        Catch == 0) %>%
        dplyr::mutate(Year = Year - data_list$styr + 1) %>%
        dplyr::select(Fleet_code, Year) %>%
        as.matrix()
      start_par$ln_F[zero_catch] <- -999
      rm(zero_catch)

      # Update proj F prop
      start_par$proj_F_prop <- data_list$fleet_control$proj_F_prop
    }

    mod_objects$initial_params <- start_par
    if(verbose > 0) {message("Step 1: Parameter build complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 3: Load/build map ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(map)) {
      map <- suppressWarnings(build_map(data_list, start_par, debug = estimateMode %in% c(2, 4), # Turn off hindcast parameters if debugging or projection mode
                                        random_rec = random_rec, random_sel = random_sel))
    } else{
      map <- map
    }
    if(verbose > 0) {message("Step 2: Map build complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 4: Get bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(bounds)) {
      bounds <- Rceattle::build_bounds(param_list = start_par, data_list)
    } else {
      bounds = bounds
    }
    if(verbose > 0) {message("Step 3: Parameter bounds complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 5: Setup random effects ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Turns on laplace approximation
    random_vars <- c()
    if (random_rec) {
      if(initMode > 0){
        random_vars <- c(random_vars , "rec_dev", "init_dev")
      } else{
        random_vars <- c(random_vars , "rec_dev")
      }
    }
    if(random_q){
      random_vars <- c(random_vars , "index_q_dev")
    }
    if(random_sel){
      random_vars <- c(random_vars , "ln_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
    }
    if(sum(data_list$M1_re) > 0){
      random_vars <- c(random_vars, "ln_M1_dev")
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 6: Reorganize data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if(!is.null(TMBfilename)){
      TMB::compile(paste0(TMBfilename,".cpp"),
                   PKG_CXXFLAGS = "",
                   framework = "TMBad",
                   safebounds = FALSE, safeunload = FALSE)
      dyn.load(TMB::dynlib(TMBfilename))
      TMBfilename <- basename(TMBfilename)
    }
    if(is.null(TMBfilename)){
      TMBfilename <- "ceattle_v01_11"
    }

    # - Check for data errors
    Rceattle:::data_check(data_list)

    # - Reorganize data for .cpp file
    data_list_reorganized <- Rceattle::rearrange_dat(data_list)
    data_list_reorganized = c(list(model = TMBfilename), data_list_reorganized)
    data_list_reorganized$forecast <- rep(0, data_list_reorganized$nspp) # Hindcast switch

    # - Update comp weights, future F (if input) and F_prop from data
    # - Age/length composition
    if(!is.null(data_list$fleet_control$Comp_weights)){
      start_par$comp_weights = data_list$fleet_control$Comp_weights
    }
    # - CAAL
    if(!is.null(data_list$fleet_control$CAAL_weights)){
      start_par$caal_weights = data_list$fleet_control$CAAL_weights
    }

    # - Diet composition
    if(!is.null(data_list$Diet_comp_weights)){
      start_par$diet_comp_weights = data_list$Diet_comp_weights
    }

    # - Proportion of projected F to each fleet
    start_par$proj_F_prop = data_list$fleet_control$proj_F_prop

    # - Fixed fishing mortality for projections for each species
    if(!is.null(HCR$Ftarget) & HCR$HCR == 2){
      start_par$ln_Ftarget = log(HCR$Ftarget)
    }

    # - Update M1 parameter object from data if initial parameter values input
    if(updateM1){
      m1 <- array(0, dim = c(data_list$nspp,
                             max(data_list$nsex, na.rm = T),
                             max(data_list$nages, na.rm = T))) # Set up array

      # Initialize from inputs
      for (i in 1:nrow(data_list$M1_base)) {
        sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
        sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

        # Handle sex == 0 case for 2-sex species
        sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

        # Fill in M1 array from fixed values for each sex
        for(j in 1:length(sex_values)){
          m1[sp, sex_values[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = T)) + 2])
        }
      }
      start_par$ln_M1 <- log(m1)
    }

    # - Update alpha for stock-recruit if fixed/prior and initial parameter values input
    if(data_list$srr_est_mode %in% c(0,2) & data_list$srr_pred_fun > 3){
      start_par$rec_pars[,2] <- log(data_list$srr_prior)
    }

    if(verbose > 0) {message("Step 4: Data rearrange complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 7: Set up parameter bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    L <- c()
    U <- c()
    for(i in 1:length(map$mapFactor)){
      if(!names(map$mapFactor)[i] %in% random_vars){ # Dont have bounds for random effects
        L = c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
        U = c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
      }
    }

    # Dimension check
    start_par <- start_par[names(map$mapFactor), drop = F]
    dim_check <- sapply(start_par, function(x) length(unlist(x))) == sapply(map$mapFactor, function(x) length(unlist(x)))
    if(sum(dim_check) != length(dim_check)){
      stop(print(paste0("Map and parameter objects are not the same size for: ", names(dim_check)[which(dim_check == FALSE)])))
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 8: Phase hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Set default phasing
    if(is.logical(phase)){
      if(phase){
        phaseList <- set_phases()
      }
    }

    if(!is.logical(phase)){
      warning("Using input phase. Please set phase = TRUE if using defaults.")
      phaseList = phase
      phase = TRUE
    }


    step = 5
    if(phase & estimateMode %in% c(0,1) ){
      if(verbose > 0) {message(paste0("Step ", step,": Phasing begin"))}; step = step + 1
      phase_pars <- Rceattle::TMBphase(
        data = data_list_reorganized,
        parameters = start_par,
        map = map$mapFactor,
        random = random_vars,
        phases = phaseList,
        model_name = TMBfilename,
        silent = verbose != 2,
        use_gradient = use_gradient,
        control = control
      )

      # Save output
      mod_objects$phase_params <- phase_pars
      start_par <- phase_pars

      if(verbose > 0) {message(paste0("Step ", step,": Phasing complete"))}
      step = step + 1
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 9: Fit hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # * Build ----
    if(sum(as.numeric(unlist(map$mapFactor)), na.rm = TRUE) == 0){stop("Map of length 0: all NAs")}
    obj = TMB::MakeADFun(
      data_list_reorganized,
      parameters = start_par,
      DLL = TMBfilename,
      map = map$mapFactor,
      random = random_vars,
      silent = verbose != 2
    )

    # -- Save objects
    mod_objects <- c(
      list(
        TMBfilename = TMBfilename,
        bounds = bounds,
        map = map
      ),
      mod_objects)

    if(verbose > 0) {message(paste0("Step ",step, ": Hindcast build complete"))}
    step = step + 1


    # * Optimize hindcast ----
    if(estimateMode %in% c(0,1,2,4)){
      opt <- suppressMessages(
        TMBhelper::fit_tmb(obj = obj,
                           fn=obj$fn,
                           gr=obj$gr,
                           startpar=obj$par,
                           lower = L,
                           upper = U,
                           loopnum = loopnum,
                           newtonsteps = newtonsteps,
                           getsd = getsd,
                           control = control,
                           bias.correct = bias.correct,
                           bias.correct.control=list(sd=getsd),
                           getJointPrecision = getJointPrecision,
                           getReportCovariance = getReportCovariance,
                           quiet = verbose < 2)
      )

      if(verbose > 0 & estimateMode != 4) {
        message("Step ",step, ": Hindcast optimization complete")
        step = step + 1
      }
      if(verbose > 0 & estimateMode == 4) {
        message("Step ",step, ": 'dummy' optimization complete")
        step = step + 1
      }

      # -- Convergence warnings
      if(estimateMode %in% c(0,1)){
        if(is.null(opt$SD) & getsd){

          message( "#################################################" )
          message( "Model did not converge, check 'identified'" )
          message( "#################################################" )

          # Bad parameter identification
          identified <- tryCatch({suppressMessages(TMBhelper::check_estimability(obj))
          },
          error = function(e){
            return("Some gradients are high, please improve optimization and only then use `Check_Identifiable`")
          })

          # Make into list if gradients were low for diagnostics
          if(class(identified) != "character"){
            identified_param_list <- obj$env$parList(identified$BadParams$Param_check)
            identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==0,"Not estimated",x), how = "replace")
            identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==1,"OK",x), how = "replace")
            identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==2,"BAD",x), how = "replace")
            identified$param_list <- identified_param_list
          }
          mod_objects$identified <- identified
        }
      }
    }


    # * Get MLEs ----
    if (estimateMode > 1) { # Debugging and projection only: use initial parameters
      last_par <- start_par
    } else{
      # Fixed effects
      if(length(random_vars) == 0){
        last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
      } else { # Random effects
        last_par = try(obj$env$parList())
      }
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 10: Run projection ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if(estimateMode %in% c(0,2,4)){
      if(!data_list$HCR %in% c(0)){ # - All HCRs except no F and fixed F

        # * 10.1: Single species mode ----
        if(msmMode == 0){

          # Turn BRP estimation on within likelihood
          data_list_reorganized$forecast <- rep(1, data_list_reorganized$nspp)

          # -- Update map in obs
          hcr_map <- Rceattle::build_hcr_map(data_list, map, debug = estimateMode > 3)
          if(sum(!is.na(unlist(hcr_map$mapFactor))) == 0){stop("HCR map of length 0: all NAs")}

          obj = TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL = TMBfilename,
            map = hcr_map$mapFactor,
            random = random_vars,
            silent = verbose != 2
          )


          if(verbose > 0) {message(paste0("Step ",step, ": Projection build complete"))}
          step = step + 1

          # -- Optimize
          opt = suppressMessages(
            TMBhelper::fit_tmb(obj = obj,
                               fn=obj$fn,
                               gr=obj$gr,
                               startpar=obj$par,
                               loopnum = loopnum,
                               getsd = getsd,
                               control = control,
                               bias.correct = bias.correct,
                               bias.correct.control=list(sd=getsd),
                               getJointPrecision = FALSE,
                               quiet = verbose < 2,
            )
          )
        }


        # * 10.2: Multi-species mode ----
        if(msmMode > 0){

          # Loop across species orders
          for(HCRiter in 1:max(data_list$HCRorder)){

            # -- Update map in obs
            hcr_map <- Rceattle::build_hcr_map(data_list, map,
                                               debug = estimateMode > 3,
                                               all_params_on = FALSE,
                                               HCRiter = HCRiter)
            if(sum(as.numeric(unlist(hcr_map$mapFactor)), na.rm = TRUE) == 0){stop("HCR map of length 0: all NAs")}

            # -- Get SB0: SSB when model is projected forward under no fishing
            data_list_reorganized$forecast <- data_list$HCRorder <= HCRiter # What species to include in likelihood
            params_on <- c(1:data_list$nspp)[which(data_list$HCRorder == HCRiter)] # Only update the one being estimated (not all)
            quantities <- obj$report(obj$env$last.par.best)
            SB0 <- quantities$ssb[, ncol(quantities$ssb)]
            B0 <- quantities$biomass[, ncol(quantities$biomass)]
            data_list_reorganized$MSSB0[params_on] <- SB0[params_on]
            data_list_reorganized$MSB0[params_on] <- B0[params_on]

            # --- Adjust Ftarget inits
            params_off <- c(1:data_list$nspp)[which(data_list$HCRorder > HCRiter)]
            last_par$ln_Ftarget[params_on] <- 0
            last_par$ln_Ftarget[params_off] <- -999

            # --- Update model object for HCR
            obj = TMB::MakeADFun(
              data_list_reorganized,
              parameters = last_par,
              DLL = TMBfilename,
              map = hcr_map$mapFactor,
              random = random_vars,
              silent = verbose != 2
            )

            if(verbose > 0) {message(paste0("Step ",step," - HCRiter ",HCRiter, ": Projection build complete. Optimizing."))}
            step = step + 1

            # -- Optimize
            if(data_list$HCR != 2){ # Fixed F does not need estimation
              opt = suppressMessages(
                TMBhelper::fit_tmb(obj = obj,
                                   fn=obj$fn,
                                   gr=obj$gr,
                                   startpar=obj$par,
                                   loopnum = loopnum,
                                   getsd = getsd,
                                   bias.correct = bias.correct,
                                   bias.correct.control=list(sd=getsd),
                                   control = control,
                                   getJointPrecision = FALSE,
                                   quiet = verbose < 2,
                )
              )

              # --- Update F from opt
              last_par$ln_Ftarget[params_on] <- opt$par[1:length(params_on)]
            }
          }
        }


        if(verbose > 0) {message(paste0("Step ",step, ": Projection optimization complete"))}
        step = step + 1


        # * Update MLEs ----
        if (estimateMode > 2) { # Debugging, give initial parameters
          last_par <- start_par
        }else{
          if(!random_rec){
            last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
          } else {
            last_par = try(obj$env$parList())
          }
        }


        # * Projection uncertainty ----
        # Updates the model with all hindcast and BRP parameters "turned on" to get out uncertainty estimates in the projection
        if(projection_uncertainty){

          # -- Update both map in to have BRP and hindcast parameters on
          hcr_map_proj <- Rceattle::build_hcr_map(data_list, map, debug = estimateMode > 3, all_params_on = TRUE)
          if(sum(as.numeric(unlist(hcr_map_proj$mapFactor)), na.rm = TRUE) == 0){stop("HCR projection map of length 0: all NAs")}

          # --- Update model object with BRP and hindcast parameters turned for BRP and hindcast
          obj = TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL = TMBfilename,
            map = hcr_map_proj$mapFactor,
            random = random_vars,
            silent = verbose != 2
          )

          # -- Replace sdreport with new sdreport
          opt$SD <- TMB::sdreport(obj)
        }

      } # End estimable BRP/HCR projections
    } # End projection


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 11: Save output ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # - Save estimated parameters
    mod_objects$estimated_params <- last_par
    mod_objects$obj = obj

    # - Get quantities
    quantities <- obj$report(obj$env$last.par.best)

    # -- Warning for discontinuous likelihood
    if(estimateMode %in% c(0:2)){
      if(!is.null(opt$SD) & random_rec == FALSE){
        if(abs(opt$objective - quantities$jnll) > rel_tol){
          message( "#################################################" )
          message( "Convergence warning (8): discontinuous likelihood" )
          message( "#################################################" )
        }
      }
    }

    # -- Save derived quantities
    mod_objects$quantities <- Rceattle::rename_output(data_list = data_list, quantities = quantities)

    # -- Save data w/ mcallister
    mod_objects$data_list <- Rceattle::calc_mcall_ianelli(data_list = data_list, data_list_reorganized = data_list_reorganized, quantities = quantities)
    mod_objects$data_list <- Rceattle::calc_mcall_ianelli_diet(data_list = mod_objects$data_list, quantities = quantities)

    # -- Run time
    mod_objects$run_time = ((Sys.time() - start_time))

    if(estimateMode < 3){
      mod_objects$opt = opt
      mod_objects$sdrep = opt$SD

    }

    class(mod_objects) <- "Rceattle"

    if(!is.null(file)){
      save(mod_objects, file = paste0(file, ".RData"))
    }

    # Free up memory
    if(estimateMode %in% 0:1){
      TMB::FreeADFun(obj) # Free memory if estimated
    }

    return(mod_objects)
  }
