#' This functions runs CEATTLE
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#'
#' @param data_list a data_list created from BSAI CEATTLE dat files \code{\link{build_dat}}, prebuilt data_list \code{\link{BS2017SS}}, or read in using \code{\link{read_excel}}.
#' @param inits (Optional) Character vector of named initial values from previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters. Can also construct using \code{\link{build_params}}
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param bounds (Optional) A bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param estimateMode 0 = Fit the hindcast model and projection with HCR specified via \code{HCR}. 1 = Fit the hindcast model only (no projection). 2 = Run the projection only with HCR specified via \code{HCR} given the initial parameters in \code{inits}.  3 = debug mode 1: runs the model through MakeADFun, but not nlminb, 4 = runs the model through MakeADFun and nlminb (will all parameters mapped out).
#' @param projection_uncertainty account for hindcast parameter uncertainty in projections when using an HCR? Default is FALSE for speed.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE.
#' @param random_q logical. If TRUE, treats annual catchability deviations as random effects.The default is FALSE.
#' @param random_sel logical. If TRUE, treats annual selectivity deviations as random effects.The default is FALSE.
#' @param HCR HCR list object from \code{\link[build_hcr]}
#' @param niter Number of iterations for multispecies model
#' @param recFun The stock recruit-relationship parameterization from \code{\link{build_srr}}.
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{N/Z ( 1 - exp(-Z) )}, 1 is \eqn{N exp(-Z/2)}, 2 is \eqn{N}.
#' @param initMode how the population is initialized. 0 = initial age-structure estimated as free parameters; 1 = equilibrium age-structure estimated out from R0,  mortality (M1), and initial population deviates; 2 = non-equilibrium age-structure estimated out from initial fishing mortality (Finit), R0,  mortality (M1), and initial population deviates.
#' @param minNByage Minimum numbers at age to put in a hard constraint that the number-at-age can not go below.
#' @param phase Optional. List of parameter object names with corresponding phase. See https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R. If NULL, will not phase model. If set to \code{"default"}, will use default phasing.
#' @param suitMode Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = length based gamma suitability, 2 = weight based gamma suitability, 3 = length based lognormal selectivity, 4 = time-varying length based lognormal selectivity.
#' @param meanyr Integer. The last year used to calculate mean suitability and recruitment, starting at \code{styr}. Defaults to $endyr$ in $data_list$. Used for MSE runs where suitability and mean recruitment is held at the value estimated from the years used to condition the OM, but F is estimated for years beyond those used to condition the OM to account for projected catch.
#' @param getsd	TRUE/FALSE whether to run standard error calculation (default = TRUE).
#' @param use_gradient use the gradient to phase (default = TRUE).
#' @param rel_tol The relative tolerance for discontinuous likelihood warnings. Set to 1. This evaluates the difference between the TMB object likelihood and the nlminb likelihood.
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param verbose 0 = Silent, 1 = print updates of model fit, 2 = print updates of model fit and TMB estimation progress.
#' @param M1Fun M1 parameterizations and priors. Use \code{build_M1}.
#' @param getJointPrecision return full Hessian of fixed and random effects.
#' @param TMBfilename TMB filename and location to use if not using the default. Do not include ".cpp" at the end.
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
#'    phase = phaseList,
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
    msmMode = 0,
    avgnMode = 0,
    initMode = 1,
    minNByage = 0,
    suitMode = 0,
    meanyr = NULL,
    phase = NULL,
    getsd = TRUE,
    use_gradient = TRUE,
    rel_tol = 1,
    control = list(eval.max = 1e+09,
                   iter.max = 1e+09, trace = 0),
    getJointPrecision = TRUE,
    loopnum = 5,
    verbose = 1,
    newtonsteps = 0,
    catch_hcr = FALSE,
    TMBfilename = NULL){

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Debugging section ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # data_list = NULL;
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
    # initMode = 1
    # minNByage = 0;
    # suitMode = 0;
    # meanyr = NULL;
    # phase = NULL;
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
    # projection_uncertainty = TRUE

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 0 - Start ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    start_time <- Sys.time()

    extend_length <- function(x){
      if(length(x) == data_list$nspp){ return(x)}
      else {return(rep(x, data_list$nspp))}
    }

    setwd(getwd())


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 1 - Load data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

    data_list <- Rceattle::clean_data(data_list)

    # Switches
    data_list$random_rec <- as.numeric(random_rec)
    data_list$estimateMode <- estimateMode
    data_list$niter <- niter
    data_list$avgnMode <- avgnMode
    if(is.null(data_list$initMode)){
      data_list$initMode <- initMode
    }
    data_list$loopnum <- loopnum
    data_list$msmMode <- msmMode
    data_list$suitMode <- as.numeric(suitMode)
    data_list$minNByage <- as.numeric(minNByage)


    if(is.null(meanyr) & is.null(data_list$meanyr)){ # If no meanyear is provided in data or function, use end year
      data_list$meanyr <- data_list$endyr
    }
    if(!is.null(meanyr)){ # If mean year is provided in function, override data
      data_list$meanyr <- meanyr
    }

    # Recruitment switches
    data_list$srr_fun <- recFun$srr_fun
    data_list$srr_pred_fun <- recFun$srr_pred_fun
    data_list$proj_mean_rec <- recFun$proj_mean_rec
    data_list$srr_est_mode <- recFun$srr_est_mode
    data_list$srr_prior_mean <- extend_length(recFun$srr_prior_mean)
    data_list$srr_prior_sd <- extend_length(recFun$srr_prior_sd)
    data_list$Bmsy_lim <- extend_length(recFun$Bmsy_lim)

    # M1
    if(!is.null(data_list$M1_model)){
      if(sum(data_list$M1_model != extend_length(M1Fun$M1_model))){
        warning("M1_model in data is different than in call `fit_mod`")
      }
    }

    # FIXME: may want to pull from data here too
    data_list$M1_model= extend_length(M1Fun$M1_model)
    updateM1 = M1Fun$updateM1
    data_list$M1_use_prior = extend_length(M1Fun$M1_use_prior) * (data_list$M1_model > 0) # Sets to 0 if M1 is fixed
    data_list$M2_use_prior = extend_length(M1Fun$M2_use_prior) * (msmMode > 0) # Sets to 0 if single-species
    data_list$M1_prior_mean = extend_length(M1Fun$M1_prior_mean)
    data_list$M1_prior_sd = extend_length(M1Fun$M1_prior_sd)


    # HCR Switches (make length of nspp if not)
    data_list$HCR = HCR$HCR
    data_list$DynamicHCR = HCR$DynamicHCR
    if(HCR$HCR != 2){ # FsprTarget is also used for fixed F (so may be of length nflts)
      data_list$FsprTarget = extend_length(HCR$FsprTarget)
    } else {
      data_list$FsprTarget = HCR$FsprTarget
    }
    data_list$FsprLimit = extend_length(HCR$FsprLimit)
    data_list$Ptarget = extend_length(HCR$Ptarget)
    data_list$Plimit = extend_length(HCR$Plimit)
    data_list$Alpha = extend_length(HCR$Alpha)
    data_list$Pstar = extend_length(HCR$Pstar)
    data_list$Sigma = extend_length(HCR$Sigma)
    data_list$Fmult = extend_length(HCR$Fmult)
    data_list$QnormHCR = qnorm(data_list$Pstar, 0, data_list$Sigma)

    if(data_list$HCR == 2 & estimateMode == 2){estimateMode = 4} # If projecting under constant F, run parmeters through obj only
    if(data_list$msmMode > 0 & !data_list$HCR %in% c(0, 2, 3, 6)){
      warning("WARNING:: Only HCRs 2, 3, and 6 work in multi-species mode currently")
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 2: Load/build parameters ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.character(inits) | is.null(inits)) {
      start_par <- suppressWarnings(Rceattle::build_params(data_list = data_list))
    } else{
      start_par <- inits
    }
    if(verbose > 0) {message("Step 1: Parameter build complete")}

    # Set Fdev for years with 0 catch to very low number
    fsh_biom_sub <- data_list$fsh_biom %>%
      filter(Year <= data_list$endyr)
    fsh_ind <- fsh_biom_sub$Fleet_code[which(fsh_biom_sub$Catch == 0)]
    yr_ind <- fsh_biom_sub$Year[which(fsh_biom_sub$Catch == 0)] - data_list$styr + 1
    for(i in 1:length(fsh_ind)){
      start_par$F_dev[fsh_ind[i], yr_ind[i]] <- -999
    }
    rm(fsh_biom_sub)


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 3: Load/build map ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(map)) {
      map <- suppressWarnings(build_map(data_list, start_par, debug = estimateMode == 4, random_rec = random_rec, random_sel = random_sel))
    } else{
      map <- map
    }
    if(verbose > 0) {message("Step 2: Map build complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 4: Get bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(bounds)) {
      bounds <- Rceattle::build_bounds(param_list = start_par, data_list)
    } else {
      bounds = bounds
    }
    if(verbose > 0) {message("Step 3: Param bounds complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 5: Setup random effects ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    random_vars <- c()
    if (random_rec) {
      if(initMode > 0){
        random_vars <- c(random_vars , "rec_dev", "init_dev")
      } else{
        random_vars <- c(random_vars , "rec_dev")
      }
    }
    if(random_q){
      random_vars <- c(random_vars , "ln_srv_q_dev")
    }
    if(random_sel){
      random_vars <- c(random_vars , "ln_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 6: Reorganize data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if(!is.null(TMBfilename)){
      TMB::compile(paste0(TMBfilename,".cpp"))
      dyn.load(dynlib(TMBfilename))
      TMBfilename <- basename(TMBfilename)

    }
    if(is.null(TMBfilename)){
      TMBfilename <- "ceattle_v01_10"
    }

    Rceattle:::data_check(data_list)

    data_list_reorganized <- Rceattle::rearrange_dat(data_list)
    data_list_reorganized = c(list(model = TMBfilename), data_list_reorganized)
    if(msmMode > 0 & !catch_hcr){
      data_list_reorganized$HCR = 0 # Estimate model with F = 0 for the projection if multispecies
    }
    data_list_reorganized$forecast <- FALSE # Don't include BRPs in likelihood of hindcast

    # - Update comp weights, future F (if input) and F_prop from data
    if(!is.null(data_list$fleet_control$Comp_weights)){
      start_par$comp_weights = data_list$fleet_control$Comp_weights
    }
    start_par$proj_F_prop = data_list$fleet_control$proj_F_prop

    nyrs_proj <- data_list$projyr - data_list$styr + 1
    if(!is.null(HCR$FsprTarget) & HCR$HCR == 2){
      start_par$ln_Ftarget = log(HCR$FsprTarget) # Fixed fishing mortality for projections for each species
    }

    # - Update M1 parameter object from data if initial parameter values input
    if(updateM1){
      m1 <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T))) # Set up array

      # Initialize from inputs
      for (i in 1:nrow(data_list$M1_base)) {
        sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
        sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

        # Fill in M1 array from fixed values for each sex
        if(sex == 0){ sex = c(1, 2)} # If sex = combined/both males and females, fill in both dimensions
        for(j in 1:length(sex)){
          m1[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = T)) + 2])
        }
      }
      start_par$ln_M1 <- log(m1)
    }

    # - Update alpha for stock-recruit if fixed/prior and initial parameter values input
    if(data_list$srr_est_mode %in% c(0,2)){
      start_par$rec_pars[,2] <- log(data_list$srr_prior_mean)
    }

    if(verbose > 0) {message("Step 4: Data rearranged complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 7: Set up parameter bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    L <- c()
    U <- c()
    for(i in 1:length(map$mapFactor)){
      if(names(map$mapFactor)[i] %!in% random_vars){ # Dont have bounds for random effects
        L = c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
        U = c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
      }
    }

    # Dimension check
    dim_check <- sapply(start_par, unlist(length)) == sapply(map$mapFactor, unlist(length))
    if(sum(dim_check) != length(dim_check)){
      stop(print(paste0("Map and parameter objects are not the same size for: ", names(dim_check)[which(dim_check == FALSE)])))
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 8: Phase hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Set default phasing
    if(!is.null(phase)){
      if(class(phase) == "character"){
        if(tolower(phase) == "default"){
          phase = list(
            dummy = 1,
            ln_pop_scalar = 4, # Scalar for input numbers-at-age
            rec_pars = 1, # Stock-recruit parameters or log(mean rec) if no stock-recruit relationship
            ln_rec_sigma = 2, # Variance for annual recruitment deviats
            rec_dev = 2, # Annual recruitment deviats
            init_dev = 2, # Age specific initial age-structure deviates or parameters
            ln_sex_ratio_sigma = 3, # Variance of sex ratio (usually fixed)
            ln_M1 = 4, #  Estimated natural or residual mortality
            ln_mean_F = 1, # Mean fleet-specific fishing mortality
            ln_Flimit = 3, # Estimated F limit
            ln_Ftarget = 3, # Estimated F target
            ln_Finit = 3, # Estimated fishing mortality for non-equilibrium initial age-structure
            proj_F_prop = 1, # Fixed fleet-specific proportion of Flimit and Ftarget apportioned within each species
            F_dev = 1, # Annual fleet specific fishing mortality deviates
            ln_srv_q = 3, # Survey catchability
            ln_srv_q_dev = 5, # Annual survey catchability deviates (if time-varying)
            ln_sigma_srv_q = 4, # Prior SD for survey catchability deviates
            ln_sigma_time_varying_srv_q = 4, # SD for annual survey catchability deviates (if time-varying)
            sel_coff = 3, # Non-parametric selectivity coefficients
            sel_coff_dev = 4, # Annual deviates for non-parametric selectivity coefficients
            ln_sel_slp = 3, # Slope parameters for logistic forms of selectivity
            sel_inf = 3, # Asymptote parameters for logistic forms of selectivity
            ln_sel_slp_dev = 5, # Annual deviates for slope parameters for logistic forms of selectivity (if time-varying)
            sel_inf_dev = 5, # Annual deviates for asymptote parameters for logistic forms of selectivity (if time-varying)
            ln_sigma_sel = 4, # SD for annual selectivity deviates (if time-varying)
            sel_curve_pen = 4, # Penalty for non-parametric selectivity
            ln_sigma_srv_index = 2, # Log SD for survey lognormal index likelihood (usually input)
            ln_sigma_fsh_catch = 2, # Log SD for lognormal catch likelihood (usually input)
            comp_weights = 4, # Weights for multinomial comp likelihood
            logH_1 = 6,  # Functional form parameter (not used in MSVPA functional form)
            logH_1a = 6, # Functional form parameter (not used in MSVPA functional form)
            logH_1b = 6, # Functional form parameter (not used in MSVPA functional form)
            logH_2 = 6, # Functional form parameter (not used in MSVPA functional form)
            logH_3 = 6, # Functional form parameter (not used in MSVPA functional form)
            H_4 = 6, # Functional form parameter (not used in MSVPA functional form)
            log_gam_a = 5, # Suitability parameter (not used in MSVPA style)
            log_gam_b = 5, # Suitability parameter (not used in MSVPA style)
            log_phi = 5 # Suitability parameter (not used in MSVPA style)
          )


          # debugphase = list(
          #   dummy = 1,
          #   ln_pop_scalar = 5, # Scalar for input numbers-at-age
          #   rec_pars = 1, # Stock-recruit parameters or log(mean rec) if no stock-recruit relationship
          #   ln_rec_sigma = 4, # Variance for annual recruitment deviats
          #   rec_dev = 2, # Annual recruitment deviats
          #   init_dev = 3, # Age specific initial age-structure deviates or parameters
          #   ln_sex_ratio_sigma = 3, # Variance of sex ratio (usually fixed)
          #   ln_M1 = 4, #  Estimated natural or residual mortality
          #   ln_mean_F = 6, # Mean fleet-specific fishing mortality
          #   ln_Flimit = 15, # Estimated F limit
          #   ln_Ftarget = 15, # Estimated F target
          #   ln_Finit = 7, # Estimated fishing mortality for non-equilibrium initial age-structure
          #   proj_F_prop = 14, # Fixed fleet-specific proportion of Flimit and Ftarget apportioned within each species
          #   F_dev = 7, # Annual fleet specific fishing mortality deviates
          #   ln_srv_q = 10, # Survey catchability
          #   ln_srv_q_dev = 11, # Annual survey catchability deviates (if time-varying)
          #   ln_sigma_srv_q = 15, # Prior SD for survey catchability deviates
          #   ln_sigma_time_varying_srv_q = 15, # SD for annual survey catchability deviates (if time-varying)
          #   sel_coff = 8, # Non-parametric selectivity coefficients
          #   sel_coff_dev = 11, # Annual deviates for non-parametric selectivity coefficients
          #   ln_sel_slp = 9, # Slope parameters for logistic forms of selectivity
          #   sel_inf = 9, # Asymptote parameters for logistic forms of selectivity
          #   ln_sel_slp_dev = 11, # Annual deviates for slope parameters for logistic forms of selectivity (if time-varying)
          #   sel_inf_dev = 11, # Annual deviates for asymptote parameters for logistic forms of selectivity (if time-varying)
          #   ln_sigma_sel = 12, # SD for annual selectivity deviates (if time-varying)
          #   sel_curve_pen = 13, # Penalty for non-parametric selectivity
          #   ln_sigma_srv_index = 14, # Log SD for survey lognormal index likelihood (usually input)
          #   ln_sigma_fsh_catch = 14, # Log SD for lognormal catch likelihood (usually input)
          #   comp_weights = 15, # Weights for multinomial comp likelihood
          #   logH_1 = 15,  # Functional form parameter (not used in MSVPA functional form)
          #   logH_1a = 15, # Functional form parameter (not used in MSVPA functional form)
          #   logH_1b = 15, # Functional form parameter (not used in MSVPA functional form)
          #   logH_2 = 15, # Functional form parameter (not used in MSVPA functional form)
          #   logH_3 = 15, # Functional form parameter (not used in MSVPA functional form)
          #   H_4 = 15, # Functional form parameter (not used in MSVPA functional form)
          #   log_gam_a = 15, # Suitability parameter (not used in MSVPA style)
          #   log_gam_b = 15, # Suitability parameter (not used in MSVPA style)
          #   log_phi = 15 # Suitability parameter (not used in MSVPA style)
          # )
        }
      }

      if(class(phase) == "character"){
        if(tolower(phase) != "default"){
          warning("phase misspecified: please set to 'default' or list with the same order as parameters.")
        }
      }
    }

    step = 5
    if(!is.null(phase) & estimateMode %in% c(0,1) ){
      if(verbose > 0) {message(paste0("Step ", step,": Phasing begin"))}
      phase_pars <- Rceattle::TMBphase(
        data = data_list_reorganized,
        parameters = start_par,
        map = map$mapFactor,
        random = random_vars,
        phases = phase,
        model_name = TMBfilename,
        silent = verbose != 2,
        use_gradient = use_gradient,
        control = control
      )

      start_par <- phase_pars

      if(verbose > 0) {message(paste0("Step ", step,": Phasing complete - getting final estimates"))}
      step = step + 1
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 9: Fit hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if(estimateMode != 2){ # dont build if projection and estimating HCR parameters
      if(sum(as.numeric(unlist(map$mapFactor)), na.rm = TRUE) == 0){stop("Map of length 0: all NAs")}
      obj = TMB::MakeADFun(
        data_list_reorganized,
        parameters = start_par,
        DLL = TMBfilename,
        map = map$mapFactor,
        random = random_vars,
        silent = verbose != 2
      )
    }

    # -- Save objects
    mod_objects <-
      list(
        TMBfilename = TMBfilename,
        initial_params = start_par,
        bounds = bounds,
        map = map
      )

    if(verbose > 0) {message(paste0("Step ",step, ": final build complete. Optimizing."))}
    step = step + 1


    # -- Optimize hindcast
    if(estimateMode %in% c(0,1,4)){
      opt = Rceattle::fit_tmb(obj = obj,
                              fn=obj$fn,
                              gr=obj$gr,
                              startpar=obj$par,
                              lower = L,
                              upper = U,
                              loopnum = loopnum,
                              getsd = getsd,
                              control = control,
                              getJointPrecision = getJointPrecision,
                              quiet = verbose < 2,
      )
      if(verbose > 0) {message("Step ",step, ": Final optimization complete")
        step = step + 1
      }

      # -- Convergence warnings
      if(estimateMode %in% c(0,1)){
        # Bad parameter identification
        if(is.null(opt$SD) & getsd){
          identified <- suppressMessages(TMBhelper::check_estimability(obj))

          # Make into list
          identified_param_list <- obj$env$parList(identified$BadParams$Param_check)
          identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==0,"Not estimated",x), how = "replace")
          identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==1,"OK",x), how = "replace")
          identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==2,"BAD",x), how = "replace")
          identified$param_list <- identified_param_list
          mod_objects$identified <- identified
        }
      }
    }

    # -- Get MLEs
    if (estimateMode > 1) { # Debugging and projection only: use initial parameters
      last_par <- start_par
    } else{
      if(!random_rec){
        last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
      } else {
        last_par = try(obj$env$parList())
      }
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 10: Run HCR projections ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if(estimateMode %in% c(0,2,4)){
      if(data_list$HCR > 2){
        data_list_reorganized$forecast <- TRUE # Turn BRP estimation on within likelihood

        # - Single species mode
        if(msmMode == 0){

          # -- Update map in obs
          hcr_map <- build_hcr_map(data_list, map, debug = estimateMode > 3)
          if(sum(!is.na(unlist(hcr_map$mapFactor))) == 0){stop("HCR map of length 0: all NAs")}

          obj = TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL = TMBfilename,
            map = hcr_map$mapFactor,
            random = random_vars,
            silent = verbose != 2
          )

          # -- Optimize
          opt = Rceattle::fit_tmb(obj = obj,
                                  fn=obj$fn,
                                  gr=obj$gr,
                                  startpar=obj$par,
                                  loopnum = loopnum,
                                  getsd = getsd,
                                  control = control,
                                  getJointPrecision = FALSE,
                                  quiet = verbose < 2,
          )
        }


        # - Multi-species mode
        if(msmMode > 0){

          # -- Update map in obs
          hcr_map <- build_hcr_map(data_list, map, debug = estimateMode > 3)
          if(sum(as.numeric(unlist(hcr_map$mapFactor)), na.rm = TRUE) == 0){stop("HCR map of length 0: all NAs")}

          # -- Get quantities
          if(estimateMode == 2){ # Build obj if we havent done so already
            obj = TMB::MakeADFun(
              data_list_reorganized,
              parameters = last_par,
              DLL = TMBfilename,
              map = hcr_map$mapFactor,
              random = random_vars,
              silent = verbose != 2
            )
          }
          quantities <- obj$report(obj$env$last.par.best)

          # -- Get SB0: SSB when model is projected forward under no fishing
          SB0 <- quantities$biomassSSB[, ncol(quantities$biomassSSB)]
          B0 <- quantities$biomass[, ncol(quantities$biomass)]
          data_list_reorganized$MSSB0 <- SB0
          data_list_reorganized$MSB0 <- B0

          # -- Set HCR back to original
          data_list_reorganized$HCR <- data_list$HCR

          # --- Update model object for HCR
          obj = TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL = TMBfilename,
            map = hcr_map$mapFactor,
            random = random_vars,
            silent = verbose != 2
          )

          # -- Optimize
          opt = Rceattle::fit_tmb(obj = obj,
                                  fn=obj$fn,
                                  gr=obj$gr,
                                  startpar=obj$par,
                                  loopnum = loopnum,
                                  getsd = getsd,
                                  control = control,
                                  getJointPrecision = FALSE,
                                  quiet = verbose < 2,
          )
        }

        if(verbose > 0) {message("Step ",step, ": Projections complete")}

        # -- Update MLEs
        if (estimateMode > 2) { # Debugging, give initial parameters
          last_par <- start_par
        }else{
          if(!random_rec){
            last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
          } else {
            last_par = try(obj$env$parList())
          }
        }

        # Updates the model with all hindcast and BRP parameters "turned on" to get out uncertainty estimates in the projection
        if(projection_uncertainty){

          # -- Update both map in to have BRP and hindcast parameters on
          hcr_map_proj <- build_hcr_map_projection(data_list, map, debug = estimateMode > 3)
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
    # STEP 11: Save output ----
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

    # -- Rename jnll
    colnames(quantities$jnll_comp) <- paste0("Sp/Srv/Fsh_", 1:ncol(quantities$jnll_comp))
    rownames(quantities$jnll_comp) <- c(
      "Survey biomass",
      "Total catch",
      "Age/length composition data",
      "Sex ratio",
      "Non-parametric selectivity",
      "Selectivity deviates",
      "Selectivity normalization",
      "Catchability prior",
      "Catchability deviates",
      "Stock-recruit prior",
      "Recruitment deviates",
      "Initial abundance deviates",
      "Fishing mortality deviates",
      "SPR Calculation",
      "Zero n-at-age penalty",
      "M prior",
      "Ration",
      "Ration penalties",
      "Stomach content proportion by weight"
    )


    colnames(quantities$biomassSSB) <- data_list$styr:data_list$projyr
    colnames(quantities$R) <- data_list$styr:data_list$projyr

    rownames(quantities$biomassSSB) <- data_list$spnames
    rownames(quantities$R) <- data_list$spnames

    # -- Save derived quantities
    mod_objects$quantities <- quantities


    # - Calculate Mcallister-Iannelli coefficients
    # Effective sample size for the length data for year y

    eff_n_mcallister <- rowSums(quantities$comp_hat * (1 - quantities$comp_hat), na.rm = TRUE)/rowSums((data_list_reorganized$comp_obs - quantities$comp_hat)^2, na.rm = TRUE) # sum_length (p_hat * (1 - p_hat))/ sum_length ((p - p_hat) ^ 2)


    # Loop fleets and take harmonic mean
    data_list$fleet_control$Est_weights_mcallister <- NA
    for(flt in unique(data_list$comp_data$Fleet_code)){
      comp_sub <- which(data_list$comp_data$Fleet_code == flt & data_list$comp_data$Year > 0)
      data_list$fleet_control$Est_weights_mcallister[which(data_list$fleet_control$Fleet_code == flt)] <- ((1/length(comp_sub))*sum((eff_n_mcallister[comp_sub]/data_list$comp_data$Sample_size[comp_sub])^-1))^-1
    }

    # -- Save data w/ mcallister
    mod_objects$data_list <- data_list


    # - Save objects
    mod_objects$run_time = ((Sys.time() - start_time))

    if(estimateMode < 3){
      mod_objects$opt = opt
      mod_objects$sdrep = opt$SD

    }

    class(mod_objects) <- "Rceattle"

    if(!is.null(file)){
      save(mod_objects, file = paste0(file, ".RData"))
    }

    # suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
    return(mod_objects)

    # Free up memory
    TMB::FreeADFun(obj)
  }



#' Function to clean data for Rceattle runs
#'
#' @param data_list
#'
#' @export
#'
clean_data <- function(data_list){

  # - Remove years of data previous to start year
  data_list$UobsWtAge <- as.data.frame(data_list$UobsWtAge)
  data_list$UobsAge <- as.data.frame(data_list$UobsAge)
  data_list$wt <- data_list$wt[which(data_list$wt$Year == 0 | data_list$wt$Year >= data_list$styr),]
  data_list$UobsAge <- data_list$UobsAge[which(data_list$UobsAge$Year == 0 | data_list$UobsAge$Year >= data_list$styr),]
  data_list$UobsWtAge <- data_list$UobsWtAge[which(data_list$UobsWtAge$Year == 0 | data_list$UobsWtAge$Year >= data_list$styr),]
  data_list$srv_biom <- data_list$srv_biom[which(abs(data_list$srv_biom$Year) >= data_list$styr),]
  data_list$fsh_biom <- data_list$fsh_biom[which(abs(data_list$fsh_biom$Year) >= data_list$styr),]
  data_list$comp_data <- data_list$comp_data[which(abs(data_list$comp_data$Year) >= data_list$styr),]
  data_list$emp_sel <- data_list$emp_sel[which(data_list$emp_sel$Year == 0 | data_list$emp_sel$Year >= data_list$styr),]
  data_list$NByageFixed <- data_list$NByageFixed[which(data_list$NByageFixed$Year == 0 | data_list$NByageFixed$Year >= data_list$styr),]
  data_list$Pyrs <- data_list$Pyrs[which(data_list$Pyrs$Year == 0 | data_list$Pyrs$Year >= data_list$styr),]

  # - Add temp multi-species SB0
  if(is.null(data_list$MSSB0)){
    data_list$MSSB0 <- rep(999, data_list$nspp)
    data_list$MSB0 <- rep(999, data_list$nspp)
  }

  # - Remove years of data after to proj year
  data_list$wt <- data_list$wt[which(data_list$wt$Year <= data_list$projyr),]
  data_list$UobsAge <- data_list$UobsAge[which(data_list$UobsAge$Year <= data_list$projyr),]
  data_list$UobsWtAge <- data_list$UobsWtAge[which(data_list$UobsWtAge$Year <= data_list$projyr),]
  data_list$srv_biom <- data_list$srv_biom[which(abs(data_list$srv_biom$Year) <= data_list$projyr),]
  data_list$fsh_biom <- data_list$fsh_biom[which(abs(data_list$fsh_biom$Year) <= data_list$projyr),]
  data_list$comp_data <- data_list$comp_data[which(abs(data_list$comp_data$Year) <= data_list$projyr),]
  data_list$emp_sel <- data_list$emp_sel[which(data_list$emp_sel$Year <= data_list$projyr),]
  data_list$NByageFixed <- data_list$NByageFixed[which(data_list$NByageFixed$Year <= data_list$projyr),]
  data_list$Pyrs <- data_list$Pyrs[which(data_list$Pyrs$Year <= data_list$projyr),]


  # - Extend catch data to proj year for projections
  if(data_list$projyr > data_list$endyr){
    # yrs_proj <- (data_list$endyr + 1):data_list$projyr
    # proj_fsh_biom <- data_list$fsh_biom %>%
    #   group_by(Fleet_code) %>%
    #   slice(rep(n(),  length(yrs_proj))) %>%
    #   mutate(Year = yrs_proj, Catch = NA)
    # data_list$fsh_biom <- rbind(data_list$fsh_biom, proj_fsh_biom)

    for(flt in (unique(data_list$fsh_biom$Fleet_code))){
      fsh_biom_sub <- data_list$fsh_biom[which(data_list$fsh_biom$Fleet_code == flt),]
      yrs_proj <- (data_list$endyr + 1):data_list$projyr
      yrs_proj <- yrs_proj[which(yrs_proj %!in% fsh_biom_sub$Year)]
      nyrs_proj <- length(yrs_proj)
      proj_fsh_biom <- data.frame(Fleet_name = rep(fsh_biom_sub$Fleet_name[1], nyrs_proj),
                                  Fleet_code = rep(flt, nyrs_proj),
                                  Species = rep(fsh_biom_sub$Species[1], nyrs_proj),
                                  Year = yrs_proj,
                                  Month = rep(fsh_biom_sub$Month[length(fsh_biom_sub$Month)], nyrs_proj),
                                  Selectivity_block = rep(fsh_biom_sub$Selectivity_block[length(fsh_biom_sub$Selectivity_block)], nyrs_proj),
                                  Catch = rep(NA, nyrs_proj),
                                  Log_sd = rep(fsh_biom_sub$Log_sd[length(fsh_biom_sub$Log_sd)], nyrs_proj))
      data_list$fsh_biom <- rbind(data_list$fsh_biom, proj_fsh_biom)
    }
  }
  data_list$fsh_biom <- data_list$fsh_biom[
    with(data_list$fsh_biom, order(Fleet_code, Year)),]

  return(data_list)
}

#' Not in function
#'
#' @param x
#' @param y
#'
#' @export
#'
'%!in%' <- function(x,y){!('%in%'(x,y))}

