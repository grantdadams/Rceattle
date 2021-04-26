#' This functions runs CEATTLE
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#'
#' @param TMBfilename (Optional) A version of the cpp CEATTLE \code{cpp_directory}. If NULL, uses the deafult and built .cpp file
#' @param cpp_directory (Optional) The directory where the cpp file is found
#' @param data_list a data_list created from BSAI CEATTLE dat files \code{\link{build_dat}}, prebuilt data_list \code{\link{BS2017SS}}, or read in using \code{\link{read_excel}}.
#' @param inits (Optional) Character vector of named initial values from previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters. Can also construct using \code{\link{build_params}}
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param bounds (Optional) A bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE.
#' @param niter Number of iterations for multispecies model
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{N/Z ( 1 - exp(-Z) )}, 1 is \eqn{N exp(-Z/2)}, 2 is \eqn{N}.
#' @param minNByage Minimum numbers at age to put in a hard constraint that the number-at-age can not go below.
#' @param phase Optional. List of parameter object names with corresponding phase. See https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R. If NULL, will not phase model. If set to \code{"default"}, will use default phasing.
#' @param silent logical. IF TRUE, includes TMB estimation progress
#' @param suitMode Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = length based gamma selectivity from Kinzey and Punt (2009), 2 = time-varying length based gamma selectivity from Kinzey and Punt (2009), 3 = time-varying weight based gamma selectivity from Kinzey and Punt (2009), 4 = length based lognormal selectivity, 5 = time-varying length based lognormal selectivity, 6 = time-varying weight based lognormal selectivity,
#' @param getsd	Boolean whether to run standard error calculation
#' @param use_gradient use the gradient to phase. Default = TRUE
#' @param rel_tol The relative tolerance for discontinuous likelihood warnings. Set to 1. This evaluates the difference between the TMB object likelihood and the nlminb likelihood.
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param getHessian return Hessian for usage in later code
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
#'  \item{opt: Optimized model object from `nlimb``}
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
#'
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
#'   ln_FSPR = 3,
#'   proj_F_prop = 1,
#'   F_dev = 1,
#'   ln_srv_q = 3,
#'   srv_q_pow = 4,
#'   ln_srv_q_dev = 5,
#'   ln_srv_q_dev_re = 4,
#'   ln_sigma_srv_q = 4,
#'   ln_sigma_time_varying_srv_q = 4,
#'   sel_coff = 3,
#'   sel_curve_pen = 4,
#'   ln_sex_ratio_sigma = 3,
#'   ln_sel_slp = 3,
#'   ln_M1 = 4,
#'   sel_inf = 3,
#'   ln_sel_slp_dev = 5,
#'   sel_inf_dev = 5,
#'   ln_sel_slp_dev_re = 4,
#'   sel_inf_dev_re = 4,
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
#'    debug = 0, # Estimate
#'    random_rec = FALSE, # No random recruitment
#'    msmMode = 0, # Single species mode
#'    avgnMode = 0,
#'    phase = phaseList,
#'    silent = TRUE)
#'
#' @export

fit_mod <-
  function(TMBfilename = "ceattle_v01_07",
           cpp_directory = NULL,
           data_list = NULL,
           inits = NULL,
           map = NULL,
           bounds = NULL,
           file = NULL,
           debug = FALSE,
           random_rec = TRUE,
           niter = 3,
           msmMode = 0,
           avgnMode = 0,
           minNByage = 0,
           suitMode = 0,
           phase = NULL,
           silent = FALSE,
           recompile = FALSE,
           getsd = TRUE,
           use_gradient = TRUE,
           rel_tol = 1,
           control = list(eval.max = 1e+09,
                          iter.max = 1e+09, trace = 0),
           getHessian = TRUE,
           loopnum = 5,
           newtonsteps = 0) {
    start_time <- Sys.time()

    setwd(getwd())
    '%!in%' <- function(x,y)!('%in%'(x,y))

    #--------------------------------------------------
    # 1. DATA and MODEL PREP
    #--------------------------------------------------

    # STEP 1 - LOAD DATA
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

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

    # Switches
    data_list$random_rec <- as.numeric(random_rec)
    data_list$debug <- debug
    data_list$niter <- niter
    data_list$avgnMode <- avgnMode
    data_list$msmMode <- msmMode
    data_list$suitMode <- as.numeric(suitMode)
    data_list$minNByage <- as.numeric(minNByage)


    # STEP 1 - LOAD PARAMETERS
    if (is.character(inits) | is.null(inits)) {
      params <- suppressWarnings(Rceattle::build_params(
        data_list = data_list,
        inits = inits
      ))
    } else{
      # inits$proj_F <- data_list$fleet_control$proj_F
      params <- inits
    }
    start_par <- params
    message("Step 1: Parameter build complete")



    # STEP 2 - BUILD MAP
    if (is.null(map)) {
      map <-
        suppressWarnings(build_map(data_list, params, debug = FALSE, random_rec = random_rec))
    } else{
      map <- map
    }
    message("Step 2: Map build complete")


    # STEP 3 - Get bounds
    if (is.null(bounds)) {
      bounds <- Rceattle::build_bounds(param_list = params, data_list)
    } else {
      bounds = bounds
    }
    message("Step 3: Param bounds complete")


    # STEP 4 - Setup random effects
    random_vars <- c("ln_srv_q_dev_re", "ln_sel_slp_dev_re", "sel_inf_dev_re")
    if (random_rec == TRUE) {
      random_vars <- c(random_vars , "rec_dev", "init_dev")
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
            ln_mean_F = 1,
            ln_FSPR = 3,
            proj_F_prop = 1,
            F_dev = 1,
            ln_srv_q = 3,
            srv_q_pow = 4,
            ln_srv_q_dev = 5,
            ln_srv_q_dev_re = 4,
            ln_sigma_srv_q = 4,
            ln_sigma_time_varying_srv_q = 4,
            sel_coff = 3,
            sel_curve_pen = 4,
            ln_sex_ratio_sigma = 3,
            ln_sel_slp = 3,
            ln_M1 = 4,
            sel_inf = 3,
            ln_sel_slp_dev = 5,
            sel_inf_dev = 5,
            ln_sel_slp_dev_re = 4,
            sel_inf_dev_re = 4,
            ln_sigma_sel = 4,
            ln_sigma_srv_index = 2,
            ln_sigma_fsh_catch = 2,
            log_gam_a = 4,
            log_gam_b = 4,
            log_phi = 4,
            comp_weights = 4
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
    if(is.null(TMBfilename) | is.null(cpp_directory)){
      cpp_directory <- system.file("executables",package="Rceattle")
      TMBfilename <- "ceattle_v01_07"
    } else{
      cpp_directory <- cpp_directory
      TMBfilename <- TMBfilename
    }
    cpp_file <- paste0(cpp_directory, "/", TMBfilename)

    # - Remove compiled files if not compatible with system
    version_files <-
      list.files(path = cpp_directory, pattern = TMBfilename)
    if (Sys.info()[1] == "Windows" &
        paste0(TMBfilename, ".so") %in% version_files) {
      suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE))
      suppressWarnings(file.remove(paste0(cpp_file, ".so")))
      suppressWarnings(file.remove(paste0(cpp_file, ".o")))
    }
    if (Sys.info()[1] != "Windows" &
        paste0(TMBfilename, ".dll") %in% version_files) {
      suppressMessages(suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE)))
      suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
      suppressWarnings(file.remove(paste0(cpp_file, ".o")))
    }
    if (Sys.info()[1] != "Windows" &
        paste0(TMBfilename, ".so") %!in% version_files) {
      suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE))
      suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
      suppressWarnings(file.remove(paste0(cpp_file, ".o")))
    }
    if(recompile){
      suppressMessages(suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE)))
      suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
      suppressWarnings(file.remove(paste0(cpp_file, ".so")))
      suppressWarnings(file.remove(paste0(cpp_file, ".o")))
    }

    old_wd <- getwd()
    setwd(cpp_directory)
    TMB::compile(paste0(TMBfilename, ".cpp"), CPPFLAGS="-Wno-ignored-attributes")
    dyn.load(TMB::dynlib(paste0(TMBfilename)), silent = TRUE)
    setwd(old_wd)


    message("Step 4: Compile CEATTLE complete")


    # STEP 6 - Reorganize data and build model object
    Rceattle:::data_check(data_list)
    data_list_reorganized <- Rceattle::rearrange_dat(data_list)

    # - Update comp weights from data
    if(!is.null(data_list$fleet_control$Comp_weights)){
      params$comp_weights = data_list$fleet_control$Comp_weights
    }


    # STEP 7 - Set up parameter bounds
    L <- c()
    U <- c()
    for(i in 1:length(map[[1]])){
      if(names(map[[1]])[i] %!in% random_vars){ # Dont have bounds for random effects
        L = c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map[[1]][[i]])) & !duplicated(unlist(map[[1]][[i]])))])
        U = c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map[[1]][[i]])) & !duplicated(unlist(map[[1]][[i]])))])
      }
    }


    # STEP 8 - Fit model object
    step = 5
    # If phased
    if(!is.null(phase) & debug == FALSE ){
      message(paste0("Step ", step,": Phasing begin"))
      phase_pars <- Rceattle::TMBphase(
        data = data_list_reorganized,
        parameters = params,
        map = map[[1]],
        random = random_vars,
        phases = phase,
        model_name = TMBfilename,
        silent = silent,
        use_gradient = use_gradient,
        control = control
      )

      start_par <- phase_pars

      message(paste0("Step ", step,": Phasing complete - getting final estimates"))
      step = step + 1
    }


    # STEP 9 - Fit final model
    obj = TMB::MakeADFun(
      data_list_reorganized,
      parameters = start_par,
      DLL = TMBfilename,
      map = map[[1]],
      random = random_vars,
      silent = silent
    )

    message(paste0("Step ",step, ": final build complete"))
    step = step + 1

    # Optimize
    if(debug == FALSE){
      opt = Rceattle::fit_tmb(obj = obj,
                              fn=obj$fn,
                              gr=obj$gr,
                              startpar=obj$par,
                              lower = L,
                              upper = U,
                              loopnum = loopnum,
                              getsd = getsd,
                              control = control,
                              getHessian = getHessian,
                              quiet = silent,
      )
    }

    message("Step ",step, ": Final optimization complete")

    # Get quantities
    quantities <- obj$report(obj$env$last.par.best)


    # Rename jnll
    colnames(quantities$jnll_comp) <- paste0("Sp/Srv/Fsh_", 1:ncol(quantities$jnll_comp))
    rownames(quantities$jnll_comp) <- c(
      "Survey biomass",
      "Total catch",
      "Age/length composition data",
      "Sex ratio",
      "Non-parametric selectivity",
      "Selectivity random walk deviates",
      "Selectivity random effect deviates",
      "Selectivity normalization",
      "Catchability random walk deviates",
      "Catchability random effect deviates",
      "Recruitment deviates",
      "Initial abundance deviates",
      "Fishing mortality deviates",
      "SPR Calculation",
      "Zero n-at-age penalty",
      "Ration",
      "Ration penalties",
      "Stomach content weight",
      "Stomach content numbers"
    )


    colnames(quantities$biomassSSB) <- data_list$styr:data_list$projyr
    colnames(quantities$R) <- data_list$styr:data_list$projyr

    rownames(quantities$biomassSSB) <- data_list$spnames
    rownames(quantities$R) <- data_list$spnames

    # Calculate Mcallister-Iannelli coefficients
    # Effective sample size for the length data for year y

    eff_n_mcallister <- rowSums(quantities$comp_hat * (1 - quantities$comp_hat), na.rm = TRUE)/rowSums((data_list_reorganized$comp_obs - quantities$comp_hat)^2, na.rm = TRUE) # sum_length (p_hat * (1 - p_hat))/ sum_length ((p - p_hat) ^ 2)


    # Loop fleets and take harmonic mean
    data_list$fleet_control$Est_weights_mcallister <- NA
    for(flt in unique(data_list$comp_data$Fleet_code)){
      comp_sub <- which(data_list$comp_data$Fleet_code == flt & data_list$comp_data$Year > 0)
      data_list$fleet_control$Est_weights_mcallister[which(data_list$fleet_control$Fleet_code == flt)] <- ((1/length(comp_sub))*sum((eff_n_mcallister[comp_sub]/data_list$comp_data$Sample_size[comp_sub])^-1))^-1
    }


    if (debug) {
      last_par <- params
    }
    else if (random_rec == F) {
      last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))
    }
    else{
      last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))
    }

    run_time = ((Sys.time() - start_time))


    # Return objects
    mod_objects <-
      list(
        TMBfilename = TMBfilename,
        cpp_directory = cpp_directory,
        data_list = data_list,
        initial_params = params,
        bounds = bounds,
        map = map,
        obj = obj,
        estimated_params = last_par,
        quantities = quantities,
        run_time = run_time
      )

    if(debug == FALSE){
      mod_objects$opt = opt
      mod_objects$sdrep = opt$SD

    }

    if(debug == FALSE){
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


    if(debug == FALSE){
      if(!is.null(opt$SD) & random_rec == FALSE){
        # Warning for discontinuous likelihood
        if(abs(opt$objective - quantities$jnll) > rel_tol){
          message( "#########################" )
          message( "Convergence warning (8)" )
          message( "#########################" )
        }
      }
    }




    class(mod_objects) <- "Rceattle"

    if(!is.null(file)){
      save(mod_objects, file = paste0(file, ".RData"))
    }

    # suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
    return(mod_objects)
  }


