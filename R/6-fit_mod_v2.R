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
#' @details
#' CEATTLE is an age-structured population dynamics model that can be fit with or without predation mortality. The default is to exclude predation mortality by setting \code{msmMode} to 0. Predation mortality can be included by setting \code{msmMode} with the following options:
#' \itemize{
#' \item{0. Single species mode}
#' \item{1. Holsman et al. 2015 predation based on multi-species virtual population analysis (MSVPA) based predation formation.}
#' \item{2. MSVPA Holling Type III}
#'  \item{3. Kinzey & Punt 2010 Holling Type I (linear)}
#'  \item{4. Kinzey & Punt 2010 Holling Type II}
#'  \item{5. Kinzey & Punt 2010 Holling Type III}
#'  \item{6. Kinzey & Punt 2010 Predator interference}
#'  \item{7. Kinzey & Punt 2010 Predator preemption}
#'  \item{8. Kinzey & Punt 2010 Hassell-Varley}
#'  \item{9. Kinzey & Punt 2010 Ecosim}
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
#' `quantities` from the returned `Rceattle` object includes the following:
#' \itemize{
#'  \item 1. Population components
#'  \item mn_rec: Mean recruitment; dim = [1, nspp]
#'  \item{Zed: Total mortality at age; dim = [nspp, nages, nyrs] }
#'  \item{NByage: Numbers at age; dim = [nspp, nages, nyrs] }
#'  \item{AvgN: Average numbers-at-age; dim = [nspp, nages, nyrs] }
#'  \item{S: Survival at age; dim = [nspp, nages, nyrs] }
#'  \item{biomassByage: Estimated biomass-at-age (kg); dim = [nspp, nages, nyrs] }
#'  \item{biomassSSBByage: Spawning biomass at age (kg); dim = [nspp, nages, nyrs] }
#'  \item{biomass: Estimated biomass (kg); dim = [nspp, nyrs] }
#'  \item{biomassSSB: Estimated spawning stock biomass (kg); dim = [nspp, nyrs] }
#'  \item{pmature: Estimated recruitment (n); dim = [nspp, nyrs] }
#'  \item{r_sigma: Standard deviation of recruitment variation}
#'  \item{R: Estimated recruitment (n); dim = [nspp, nyrs] }
#'  \item{M1: Base natural mortality; dim = [nspp, nages] }
#'  \item{M2: Predation mortality at age; dim = [nyrs, nages, nspp] }
#'  \item{M: Total natural mortality at age; dim = [nyrs, nages, nspp] }
#'  \item{}
#'  \item{2. Survey components}
#'  \item{srv_age_obs:}
#'  \item{srv_bio_hat: Estimated BT survey biomass (kg); dim = [nspp, nyrs] }
#'  \item{srv_hat: Estimated BT survey total abundance (n); dim = [nspp, nyrs] }
#'  \item{srv_age_hat: Estimated BT age comp; dim = [nspp, nages, nyrs] }
#'  \item{eit_hat: Estimated EIT survey biomass (kg); dim = [nyrs] }
#'  \item{eit_age_hat: Estimated EIT catch-at-age ; dim = [nyrs, srv_age_bins(0)] }
#'  \item{eit_age_comp_hat : Estimated EIT age comp ; dim = [nyrs, srv_age_bins(0)] }
#'  \item{obs_eit_age:}
#'  \item{eit_age_comp: Eit age comp; dim = [n_eit, srv_age_bins(0)] }
#'  \item{avgsel_srv: Average survey selectivity; dim = [1, nspp] }
#'  \item{srv_sel: Estimated survey selectivity at age; dim = [nspp, nyrs] }
#'  \item{}
#'  \item{3. Fishery components}
#'  \item{F: Estimated fishing mortality; dim = [nspp, nages, nyrs] }
#'  \item{F_dev:}
#'  \item{fsh_sel: Log estimated fishing selectivity; dim = [nyrs, srv_age_bins(0)] }
#'  \item{avgsel_fsh: Average fishery selectivity }
#'  \item{tc_biom_hat: Estimated total yield (kg); dim = [nspp, nyrs] }
#'  \item{catch_hat: Estimated catch-at-age (n); dim = [nspp, nages, nyrs] }
#'  \item{tc_hat: Estimated total catch (n); dim = [nspp, nyrs] }
#'  \item{fsh_age_hat: Estimated fishery age comp; dim = [nspp, nages, nyrs] }
#'  \item{fsh_age_obs: Observed fishery age comp; dim = [nyrs_fsh_comp, fsh_age_bins, nspp] }
#'  \item{}
#'  \item{3. Likelihood components}
#'  \item{jnll_comp: Matrix of negative log-likelihood components (See below) }
#'  \item{offset: Offsets for multinomial likelihood }
#'  \item{}
#'  \item{4. Ration components}
#'  \item{ConsumAge: Pre-allocated indiviudal consumption in grams per predator-age; dim = [nyrs, nages, nspp] }
#'  \item{Consum_livingAge: Pre-allocated indiviudal consumption in grams per predator-age; dim = [nyrs, nages, nspp] }
#'  \item{S2Age: pre-allocate mean stomach weight as a function of sp_age }
#'  \item{LbyAge: Length by age from LW regression }
#'  \item{mnWt_obs: Mean observed weight at age (across years); dim = [nspp, nages] }
#'  \item{fT:  Pre-allocation of temperature function of consumption; dim = [nspp, nTyrs]}
#'  \item{TempC: Bottom temperature; dim = [1, nTyrs] }
#'  \item{ration2Age: Annual ration at age (kg/yr); dim = [nyrs, nages, nspp] }
#'  \item{}
#'  \item{5. Suitability components}
#'  \item{suma_suit: Sum of suitabilities; dim = [nyrs, nages, nspp] }
#'  \item{suit_main: Suitability/gamma selectivity of predator age u on prey age a; dim = [nspp, nspp, nages, nages] }
#'  \item{suit_other: Suitability not accounted for by the included prey; dim = [nspp, nages] }
#'  \item{stom_div_bio2: // Stomach proportion over biomass; U/ (W * N) ; dim = [nspp, nspp, nages, nages, nyrs] }
#'  \item{stomKir: Stomach proportion U; dim = [nspp, nspp, nages, nages, nyrs] }
#'  \item{avail_food: Available food to predator; dim = [nyrs, nages, nspp] }
#'  \item{of_stomKir: Other food stomach content; dim = [nyrs, nages, nspp] }
#'  \item{B_eaten: Biomass of prey eaten via predation; dim = [nyrs, nages, nspp] }
#'  \item{}
#'  \item{6. Kinzey predation functions}
#'  \item{H_1: Functional response parameters from Kinzey & Punt (2009) }
#'  \item{H_1a: Functional response parameters from Kinzey & Punt (2009) }
#'  \item{H_1b: Functional response parameters from Kinzey & Punt (2009) }
#'  \item{H_2: Functional response parameters from Kinzey & Punt (2009) }
#'  \item{H_3: Functional response parameters from Kinzey & Punt (2009) }
#'  \item{gam_a: Predator gamma selectivity parameters }
#'  \item{gam_b: Predator gamma selectivity parameters }
#'  \item{N_pred_yrs: Effective numbers of predators for each age of prey }
#'  \item{N_prey_yrs: Effective numbers of prey for each age of prey }
#'  \item{N_pred_eq: Effective numbers of predators for each age of prey in equilibrium (styr_pred) }
#'  \item{N_prey_eq: Effective numbers of prey for each age of predator in equilibrium (styr_pred) }
#'  \item{pred_resp: Predator functional response }
#'  \item{Pred_r: Pred_ratio values }
#'  \item{Prey_r: Prey_ratio values }
#'  \item{Vmort_ua: Predation mortality on prey age a by single predator age u }
#'  \item{eaten_la: Number of prey of age a eaten by predator length l }
#'  \item{eaten_ua: Number of prey of age a eaten by predator age u }
#'  \item{Q_mass_l: Mass of each prey sp consumed by predator at length // FIXME: make into 4D array }
#'  \item{Q_mass_u: Mass of each prey sp consumed by predator at age // FIXME: make into 4D array }
#'  \item{Q_other_u: Mass of other prey consumed by predator at age }
#'  \item{Q_hat: Fraction for each prey type of total mass eaten by predator length }
#'  \item{T_hat: Fraction of prey of length m in predator of length l }
#'  \item{omega_hat: Estimated daily ration by predator age each year }
#'  \item{omega_hat_ave: Estimated daily ration by predator age averaged over years }
#' }
#'
#'
#' @examples
#'
#'# Load package and data
#'library(Rceattle)
#'data(BS2017SS) # ?BS2017SS for more information on the data
#'
#'# Set up phases, also the default
#'phaseList = list(
#'dummy = 1,
#'ln_pop_scalar = 11,
#'ln_mean_rec = 1,
#'ln_rec_sigma = 2,
#'rec_dev = 2,
#'init_dev = 2,
#'ln_mean_F = 1,
#'ln_FSPR = 1,
#'F_dev = 1,
#'log_srv_q = 3,
#'ln_srv_q_dev = 4,
#'ln_srv_q_dev_re = 4,
#'ln_sigma_srv_q = 4,
#'sel_coff = 3,
#'ln_sel_slp = 3,
#'sel_inf = 3,
#'ln_sel_slp_dev = 4,
#'sel_inf_dev = 4,
#'ln_sel_slp_dev_re = 4,
#'sel_inf_dev_re = 4,
#'ln_sigma_sel = 4,
#'ln_sigma_srv_index = 2,
#'ln_sigma_fsh_catch = 2,
#'logH_1 = 6,
#'logH_1a = 6,
#'logH_1b = 6,
#'logH_2 = 6,
#'logH_3 = 6,
#'H_4 = 6,
#'log_gam_a = 5,
#'log_gam_b = 5,
#'log_phi = 5
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
  function(TMBfilename = "ceattle_v01_06",
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
           rel_tol = 1) {
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


    # Remove years of data previous to start year
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


    # Extend catch data to proj year for projections
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


    # Get cpp file if not provided
    if(is.null(TMBfilename) | is.null(cpp_directory)){
      cpp_directory <- system.file("executables",package="Rceattle")
      TMBfilename <- "ceattle_v01_06"
    } else{
      cpp_directory <- cpp_directory
      TMBfilename <- TMBfilename
    }


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
    message("Step 1: Parameter build complete")



    # STEP 2 - BUILD MAP
    if (is.null(map)) {
      map <-
        suppressWarnings(build_map(data_list, params, debug = debug, random_rec = random_rec))
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
          sel_inf = 3,
          ln_sel_slp_dev = 5,
          sel_inf_dev = 5,
          ln_sel_slp_dev_re = 4,
          sel_inf_dev_re = 4,
          ln_sigma_sel = 4,
          ln_sigma_srv_index = 2,
          ln_sigma_fsh_catch = 2,
          logH_1 = 4,
          logH_1a = 4,
          logH_1b = 4,
          logH_2 = 4,
          logH_3 = 4,
          H_4 = 4,
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


    # STEP 5 - Compile CEATTLE is providing cpp file
    cpp_file <- paste0(cpp_directory, "/", TMBfilename)

    # Remove compiled files if not compatible with system
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
    if(!is.null(phase) & debug == FALSE & debug == FALSE){
      message(paste0("Step ", step,": Phasing begin"))
      phase_pars <- Rceattle::TMBphase(
        data = data_list_reorganized,
        parameters = params,
        map = map[[1]],
        random = random_vars,
        phases = phase,
        cpp_directory = cpp_directory,
        model_name = TMBfilename,
        silent = silent,
        use_gradient = use_gradient
      )

      start_par <- phase_pars

      message(paste0("Step ", step,": Phasing complete - getting final estimates"))
      step = step + 1
    }

    # Not phased
    if(is.null(phase) | debug == TRUE){
      start_par <- params
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
                              loopnum = 5,
                              getsd = getsd,
                              control = list(eval.max = 1e+09,
                                             iter.max = 1e+09, trace = 0)
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

    # Calculate MaCallister-Iannelli coefficients
    # Effective sample size for the length data for year y

    eff_n_macallister <- rowSums(quantities$comp_hat * (1 - quantities$comp_hat), na.rm = TRUE)/rowSums((data_list_reorganized$comp_obs - quantities$comp_hat)^2, na.rm = TRUE) # sum_length (p_hat * (1 - p_hat))/ sum_length ((p - p_hat) ^ 2)


    # Loop fleets and take harmonic mean
    weights_macallister <- rep(NA, length(unique(data_list$comp_data$Fleet_code)))
    data_list$fleet_control$Est_weights_macallister <- NA
    for(flt in unique(data_list$comp_data$Fleet_code)){
      comp_sub <- which(data_list$comp_data$Fleet_code == flt & data_list$comp_data$Year > 0)
      data_list$fleet_control$Est_weights_macallister[which(data_list$fleet_control$Fleet_code == flt)] <- ((1/length(comp_sub))*sum((eff_n_macallister[comp_sub]/data_list$comp_data$Sample_size[comp_sub])^-1))^-1
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
        identified <- suppressMessages(TMBhelper::Check_Identifiable(obj))

        # Make into list
        identified_param_list <- obj$env$parList(as.numeric(identified$BadParams$Param_check))
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


