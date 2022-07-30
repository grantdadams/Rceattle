#' This functions runs CEATTLE
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#' @param data_list a data_list created from BSAI CEATTLE dat files \code{\link{build_dat}}, prebuilt data_list \code{\link{BS2017SS}}, or read in using \code{\link{read_excel}}.
#' @param inits (Optional) Character vector of named initial values from previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters. Can also construct using \code{\link{build_params}}
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param bounds (Optional) A bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param estimateMode 0 = Fit the hindcast model and projection with HCR specified via \code{HCR}. 1 = Fit the hindcast model only (no projection). 2 = Run the projection only with HCR specified via \code{HCR} given the initial parameters in \code{inits}.  3 = debug mode 1: runs the model through MakeADFun, but not nlminb, 4 = runs the model through MakeADFun and nlminb (will all parameters mapped out).
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE.
#' @param random_q logical. If TRUE, treats annual catchability deviations as random effects.The default is FALSE.
#' @param random_rec logical. If TRUE, treats annual selectivity deviations as random effects.The default is FALSE.
#' @param HCR HCR list object from \code{\link[build_hcr]}
#' @param niter Number of iterations for multispecies model
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{N/Z ( 1 - exp(-Z) )}, 1 is \eqn{N exp(-Z/2)}, 2 is \eqn{N}.
#' @param minNByage Minimum numbers at age to put in a hard constraint that the number-at-age can not go below.
#' @param phase Optional. List of parameter object names with corresponding phase. See https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R. If NULL, will not phase model. If set to \code{"default"}, will use default phasing.
#' @param suitMode Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = length based gamma suitability, 2 = weight based gamma suitability, 3 = length based lognormal selectivity, 4 = time-varying length based lognormal selectivity.
#' @param meanyr Integer. The last year used to calculate mean suitability and recruitment, starting at \code{styr}. Defaults to $endyr$ in $data_list$. Used for MSE runs where suitability and mean recruitment is held at the value estimated from the years used to condition the OM, but F is estimated for years beyond those used to condition the OM to account for projected catch.
#' @param updateM1 TRUE/FALSE whether to update M1 from data, if inits are used (default = TRUE). Useful for phasing in multi-species models from single-species models, but the data have updated residual mortality (M1) for the multi-species model.
#' @param getsd	TRUE/FALSE whether to run standard error calculation (default = TRUE).
#' @param use_gradient use the gradient to phase (default = TRUE).
#' @param rel_tol The relative tolerance for discontinuous likelihood warnings. Set to 1. This evaluates the difference between the TMB object likelihood and the nlminb likelihood.
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param verbose 0 = Silent, 1 = print updates of model fit, 2 = print updates of model fit and TMB estimation progress.
#' @param getJointPrecision return full Hessian of fixed and random effects.
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
#'# Set up phases, also the default
# phase = list(
#'   dummy = 1,
#'   ln_pop_scalar = 4,
#'   ln_mean_rec = 1,
#'   ln_rec_sigma = 2,
#'   rec_dev = 2,
#'   init_dev = 2,
#'   ln_mean_F = 1,
#'   ln_Flimit = 3,
#'   ln_Ftarget = 3,
#'   proj_F_prop = 1,
#'   F_dev = 1,
#'   ln_srv_q = 3,
#'   srv_q_pow = 4,
#'   ln_srv_q_dev = 5,
#'   ln_sigma_srv_q = 4,
#'   ln_sigma_time_varying_srv_q = 4,
#'   sel_coff = 3,
#'   sel_coff_dev = 4,
#'   sel_curve_pen = 4,
#'   ln_sex_ratio_sigma = 3,
#'   ln_sel_slp = 3,
#'   ln_M1 = 4,
#'   sel_inf = 3,
#'   ln_sel_slp_dev = 5,
#'   sel_inf_dev = 5,
#'   ln_sigma_sel = 4,
#'   ln_sigma_srv_index = 2,
#'   ln_sigma_fsh_catch = 2,
#'   log_gam_a = 4,
#'   log_gam_b = 4,
#'   log_phi = 4,
#'   comp_weights = 4
#')
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
    random_rec = FALSE,
    random_q = FALSE,
    random_sel = FALSE,
    HCR = build_hcr(),
    niter = 3,
    msmMode = 0,
    avgnMode = 0,
    minNByage = 0,
    suitMode = 0,
    meanyr = NULL,
    updateM1 = TRUE,
    phase = NULL,
    getsd = TRUE,
    use_gradient = TRUE,
    rel_tol = 1,
    control = list(eval.max = 1e+09,
                   iter.max = 1e+09, trace = 0),
    getJointPrecision = TRUE,
    loopnum = 5,
    verbose = 1,
    newtonsteps = 0){


    # ### For debugging
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
    # updateM1 = TRUE;
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


    start_time <- Sys.time()

    setwd(getwd())

    #--------------------------------------------------
    # 1. DATA and MODEL PREP
    #--------------------------------------------------

    # STEP 1 - LOAD DATA
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

    data_list <- Rceattle::clean_data(data_list)

    # Switches
    data_list$random_rec <- as.numeric(random_rec)
    data_list$estimateMode <- estimateMode
    data_list$niter <- niter
    data_list$avgnMode <- avgnMode
    data_list$msmMode <- msmMode
    data_list$suitMode <- as.numeric(suitMode)
    data_list$minNByage <- as.numeric(minNByage)
    if(is.null(meanyr) & is.null(data_list$meanyr)){ # If no meanyear is provided in data or function, use end year
      data_list$meanyr <- data_list$endyr
    }
    if(!is.null(meanyr)){ # If mean year is provided in function, override data
      data_list$meanyr <- meanyr
    }


    # HCR Switches (make length of nspp if not)
    data_list$HCR = HCR$HCR
    data_list$DynamicHCR = HCR$DynamicHCR
    data_list$FsprTarget = HCR$FsprTarget[1] # NOTE: if fixed F, inputs it again below
    data_list$FsprLimit = ifelse(length(HCR$FsprLimit) == data_list$nspp, HCR$FsprLimit, rep(HCR$FsprLimit, data_list$nspp))
    data_list$Ptarget = ifelse(length(HCR$Ptarget) == data_list$nspp, HCR$Ptarget, rep(HCR$Ptarget, data_list$nspp))
    data_list$Plimit = ifelse(length(HCR$Plimit) == data_list$nspp, HCR$Plimit, rep(HCR$Plimit, data_list$nspp))
    data_list$Alpha = ifelse(length(HCR$Alpha) == data_list$nspp, HCR$Alpha, rep(HCR$Alpha, data_list$nspp))
    data_list$Pstar = ifelse(length(HCR$Pstar) == data_list$nspp, HCR$Pstar, rep(HCR$Pstar, data_list$nspp))
    data_list$Sigma = ifelse(length(HCR$Sigma) == data_list$nspp, HCR$Sigma, rep(HCR$Sigma, data_list$nspp))
    data_list$QnormHCR = ifelse(HCR$HCR == 6, qnorm(HCR$Pstar, 0, HCR$Sigma), rep(0, data_list$nspp)) # Pstar HCR

    if(data_list$HCR == 2 & estimateMode == 2){estimateMode = 4} # If projecting under constant F, run parmeters through obj only

    # STEP 1 - LOAD PARAMETERS
    if (is.character(inits) | is.null(inits)) {
      start_par <- suppressWarnings(Rceattle::build_params(
        data_list = data_list,
        inits = inits
      ))
    } else{
      # inits$proj_F <- data_list$fleet_control$proj_F
      start_par <- inits
    }
    if(verbose > 0) {message("Step 1: Parameter build complete")}



    # STEP 2 - BUILD MAP
    if (is.null(map)) {
      map <-
        suppressWarnings(build_map(data_list, start_par, debug = estimateMode > 3, random_rec = random_rec))
    } else{
      map <- map
    }
    if(verbose > 0) {message("Step 2: Map build complete")}


    # STEP 3 - Get bounds
    if (is.null(bounds)) {
      bounds <- Rceattle::build_bounds(param_list = start_par, data_list)
    } else {
      bounds = bounds
    }
    if(verbose > 0) {message("Step 3: Param bounds complete")}


    # STEP 4 - Setup random effects
    random_vars <- c()
    if (random_rec) {
      random_vars <- c(random_vars , "rec_dev", "init_dev")
    }
    if(random_q){
      random_vars <- c(random_vars , "ln_srv_q_dev")
    }
    if(random_sel){
      random_vars <- c(random_vars , "ln_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
    }



    # Set default phasing
    if(!is.null(phase)){
      if(class(phase) == "character"){
        if(tolower(phase) == "default"){
          phase = list(
            dummy = 1,
            ln_pop_scalar = 4,
            ln_mean_rec = 1,
            ln_rec_sigma = 2,
            rec_dev = 2,
            init_dev = 2,
            ln_sex_ratio_sigma = 3,
            ln_M1 = 4,
            ln_mean_F = 1,
            ln_Flimit = 3,
            ln_Ftarget = 3,
            proj_F_prop = 1,
            F_dev = 1,
            ln_srv_q = 3,
            # srv_q_pow = 4,
            ln_srv_q_dev = 5,
            ln_sigma_srv_q = 4,
            ln_sigma_time_varying_srv_q = 4,
            sel_coff = 3,
            sel_coff_dev = 4,
            ln_sel_slp = 3,
            sel_inf = 3,
            ln_sel_slp_dev = 5,
            sel_inf_dev = 5,
            ln_sigma_sel = 4,
            sel_curve_pen = 4,
            ln_sigma_srv_index = 2,
            ln_sigma_fsh_catch = 2,
            comp_weights = 4,
            logH_1 = 6,
            logH_1a = 6,
            logH_1b = 6,
            logH_2 = 6,
            logH_3 = 6,
            H_4 = 6,
            log_gam_a = 5,
            log_gam_b = 5,
            log_phi = 5
          )
        }
      }

      if(class(phase) == "character"){
        if(tolower(phase) != "default"){
          warning("phase misspecified: please set to 'default' or list with the same order as parameters.")
        }
      }
    }


    # STEP 5 - Compile CEATTLE is providing cpp file
    # - Get cpp file if not provided
    TMBfilename <- "ceattle_v01_09"


    # STEP 6 - Reorganize data and build model object
    Rceattle:::data_check(data_list)
    data_list_reorganized <- Rceattle::rearrange_dat(data_list)
    data_list_reorganized = c(list(model = "ceattle_v01_09"),data_list_reorganized)
    if(msmMode > 0 & data_list$HCR == 3){
      data_list_reorganized$HCR = 0 # Estimate model with F = 0 for the projection if multispecies
    }

    # - Update comp weights, future F (if input) and F_prop from data
    if(!is.null(data_list$fleet_control$Comp_weights)){
      start_par$comp_weights = data_list$fleet_control$Comp_weights
    }
    start_par$proj_F_prop = data_list$fleet_control$proj_F_prop

    nyrs_proj <- data_list$projyr - data_list$styr + 1
    if(!is.null(HCR$FsprTarget) & HCR$HCR == 2){
      start_par$ln_Ftarget = matrix(log(HCR$FsprTarget), nrow = data_list$nspp, ncol = nyrs_proj) # Fixed fishing mortality for projections for each species
    }

    # - Update M1 for inits
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

    if(verbose > 0) {message("Step 4: Data rearranged complete")}

    # STEP 7 - Set up parameter bounds
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


    # STEP 8 - Phase hindcast
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


    # STEP 9 - Fit final hindcast model
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


    # STEP 10 - Run HCR projections
    if(estimateMode %in% c(0,2,4)){
      if(data_list$HCR > 2){

        # - Single species mode
        if(msmMode == 0){

          # -- Update map in obs
          hcr_map <- build_hcr_map(data_list, map, debug = estimateMode > 3)
          if(sum(as.numeric(unlist(hcr_map$mapFactor)), na.rm = TRUE) == 0){stop("HCR map of length 0: all NAs")}

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
        # obj$report()$DynamicSPRtarget/obj$report()$DynamicSPR0
        # obj$report()$DynamicSPRlimit/obj$report()$DynamicSPR0
        #
        # obj$report()$SPRtarget/obj$report()$SPR0
        # obj$report()$SPRlimit/obj$report()$SPR0
        #
        # obj$report()$SPR0
        # obj$report()$SB0

        if(verbose > 0) {message("Step ",step, ": Projections complete")}

        # -- Update MLEs
        if (estimateMode > 2) { # Debugging, give initial parameters
          last_par <- start_par
        }
        else{
          if(!random_rec){
            last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
          } else {
            last_par = try(obj$env$parList())
          }
        }
      }
    }

    # - Save estimated parameters
    mod_objects$estimated_params <- last_par
    mod_objects$obj = obj

    # - Get quantities
    quantities <- obj$report(obj$env$last.par.best)

    # -- Warning for discontinuous likelihood
    if(estimateMode %in% c(0:2)){
      if(!is.null(opt$SD) & random_rec == FALSE){
        if(abs(opt$objective - quantities$jnll) > rel_tol){
          message( "#########################" )
          message( "Convergence warning (8)" )
          message( "#########################" )
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
      "NA",
      "Selectivity normalization",
      "Catchability prior",
      "Catchability deviates",
      "Recruitment deviates",
      "Initial abundance deviates",
      "Fishing mortality deviates",
      "SPR Calculation",
      "Zero n-at-age penalty",
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
  data_list$MSSB0 <- rep(999, data_list$nspp)

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

