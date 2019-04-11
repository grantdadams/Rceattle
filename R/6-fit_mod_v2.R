#' This functions runs CEATTLE
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#'
#' @param TMBfilename (Optional) A version of the cpp CEATTLE \code{cpp_directory}. If NULL, uses the deafult and built .cpp file
#' @param cpp_directory (Optional) The directory where the cpp file is found
#' @param data_list a data_list created from BSAI CEATTLE dat files \code{\link{build_dat}}, prebuilt data_list \code{\link{BS2017SS}}, or read in using \code{\link{read_excel}}.
#' @param inits (Optional) Character vector of named initial values from ADMB or list of previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters. Can also consturct using \code{\link{build_params}}
#' @param map (Optional) A prebuilt map object from \code{\link{build_map}}.
#' @param bounds (Optional) A prebuild bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE.
#' @param niter Number of iterations for multispecies model
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{\frac{N}{Z} \left( 1 - exp^{-Z} \right)}, 1 is \eqn{N e^{-Z/2}}, 2 is \eqn{N}.
#' @param silent logical. IF TRUE, includes TMB estimation progress
#' @param suitMode Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = length based gamma selectivity from Kinzey and Punt (2009), 2 = time-varing length based gamma selectivity from Kinzey and Punt (2009), 3 = time-varying weight based gamma selectivity from Kinzey and Punt (2009), 4 = length based lognormal selectivity, 5 = time-varing length based lognormal selectivity, 6 = time-varying weight based lognormal selectivity,
#' @details
#' CEATTLE is an age-structured population dynamics model that can be fit with or without predation mortality. The default is to exclude predation mortality by setting \code{msmMode} to 0. Predation mortality can be included by setting \code{msmMode} with the following options:
#' \itemize{
#' \item{0. Single species mode}
#' \item{1. Holsman et al. 2015 predation based on multi-species virtual population analysis (MSVPA) based predation formation.}
#'  \item{2. Kinzey & Punt 2010 Holling Type I (linear)}
#'  \item{3. Kinzey & Punt 2010 Holling Type II}
#'  \item{4. Kinzey & Punt 2010 Holling Type III}
#'  \item{5. Kinzey & Punt 2010 Predator interference}
#'  \item{6. Kinzey & Punt 2010 Predator preemption}
#'  \item{7. Kinzey & Punt 2010 Hassell-Varley}
#'  \item{8. Kinzey & Punt 2010 Ecosim}
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
#' `jnll_comp` includes:
#' \itemize{
#'  \item{-- Data components}
#' \item{Slot 0 -- Survey biomass}
#' \item{Slot 1 -- Survey age/length composition}
#' \item{Slot 2 -- Fishery selectivity}
#' \item{Slot 3 -- Fishery selectivity normalization}
#' \item{Slot 4 -- Total catch -- Fishery observer data}
#' \item{Slot 5 -- Fishery age/length composition -- Fishery observer data}
#' \item{Slot 6 -- Survey selectivity}
#' \item{Slot 7 -- Survey selectivity normalization}
#' \item{-- Priors}
#' \item{Slot 8 -- Tau -- Annual recruitment deviation}
#' \item{Slot 9 -- init_dev -- Initial abundance-at-age}
#' \item{Slot 10 -- Epsilon -- Annual fishing mortality deviation}
#' \item{-- M2 likelihood components}
#' \item{Slot 13 -- Ration likelihood}
#' \item{Slot 14 -- Ration penalties}
#' \item{Slot 15 -- Diet weight likelihood}
#' \item{Slot 16 -- Stomach content of prey length ln in predator length a likelihood}
#' }
#'
#' @examples
#'library(Rceattle)
#'data(BS2017SS) # ?BS2017SS for more information on the data
#'
#'# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
#'ss_run <- fit_mod(data_list = BS2017SS,
#'    inits = NULL, # Initial parameters = 0
#'    file = NULL, # Don't save
#'    debug = 0, # Estimate
#'    random_rec = FALSE, # No random recruitment
#'    msmMode = 0, # Single species mode
#'    avgnMode = 0,
#'    silent = TRUE)
#'
#' @export

fit_mod <-
  function(TMBfilename = "ceattle_v01_04",
           cpp_directory = NULL,
           data_list = NULL,
           inits = NULL,
           map = NULL,
           bounds = NULL,
           file = NULL,
           debug = T,
           random_rec = FALSE,
           niter = 3,
           msmMode = 0,
           avgnMode = 0,
           suitMode = 0,
           silent = FALSE,
           recompile = FALSE) {
    start_time <- Sys.time()

    setwd(getwd())

    #--------------------------------------------------
    # 1. DATA and MODEL PREP
    #--------------------------------------------------
    # # Check if require packages are installed and install if not
    # if ("TMB" %in% rownames(installed.packages()) == FALSE) {
    #  install.packages("TMB")
    # }
    # if ("TMBhelper" %in% rownames(installed.packages()) == FALSE) {
    #  install.packages("TMBhelper")
    # }
    # library(TMB)
    # library(TMBhelper)


    # STEP 1 - LOAD DATA
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }


    # Switches
    data_list$random_rec <- as.numeric(random_rec)
    data_list$debug <- debug
    data_list$niter <- niter
    data_list$avgnMode <- avgnMode
    data_list$msmMode <- msmMode
    data_list$suitMode <- as.numeric(suitMode)


    # Get cpp file if not provided
    if(is.null(TMBfilename) | is.null(cpp_directory)){
      cpp_directory <- system.file("executables",package="Rceattle")
      TMBfilename <- "ceattle_v01_04"
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
      inits$proj_F <- data_list$fsh_control$proj_F
      params <- inits
    }
    message("Step 1: Parameter build complete")



    # STEP 2 - BUILD MAP
    if (is.null(map)) {
      map <-
        suppressWarnings(Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec))
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
    random_vars <- c()
    if (random_rec == TRUE) {
      random_vars <- c("rec_dev")
    }

    '%!in%' <- function(x,y)!('%in%'(x,y))


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
      suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
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
    TMB::compile(paste0(TMBfilename, ".cpp"))
    dyn.load(TMB::dynlib(paste0(TMBfilename)), silent = TRUE)
    setwd(old_wd)


    message("Step 4: Compile CEATTLE complete")


    # STEP 6 - Reorganize data
    data_list2 <- rearrange_dat(data_list)

    # STEP 7 - Build and fit model object
    obj = TMB::MakeADFun(
      data_list2,
      parameters = params,
      DLL = TMBfilename,
      map = map[[1]],
      random = random_vars,
      silent = silent

    )
    message(paste0("Step 5: Build object complete"))
    # opt <- nlminb(obj$par, obj$fn, obj$gr)
    # methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
    # opt_list <- list()
    # for(i in 1:length(methods)){
    #  opt_list[i] = optimx(obj$par, function(x) as.numeric(obj$fn(x)), obj$gr, control = list(maxit = 10000), method = methods[i])
    # }


    # Remove inactive parameters from bounds and vectorize
    L = unlist(bounds$lower)[which(!is.na(unlist(map[[1]])))]
    U = unlist(bounds$upper)[which(!is.na(unlist(map[[1]])))]

    # Remove random effects from bounds
    if (random_rec == TRUE) {
      L <- L[-grep(random_vars, names(L))]
      U <- U[-grep(random_vars, names(U))]
    }

    # Optimize
    opt = Rceattle::Optimize(obj = obj,
                             fn=obj$fn,
                             gr=obj$gr,
                             startpar=obj$par,
                             lower = L,
                             upper = U,
                             loopnum = 8,
                             control = list(eval.max = 1e+08,
                                            iter.max = 1e+08, trace = 0)
    )

    message("Step 6: Optimization complete")

    # Get quantities
    quantities <- obj$report(obj$env$last.par.best)

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
        data_list = data_list,
        initial_params = params,
        bounds = bounds,
        map = map,
        obj = obj,
        opt = opt,
        sdrep = opt$SD,
        estimated_params = last_par,
        quantities = quantities,
        run_time = run_time
      )

    if(debug == 0){
      # Check identifiability
      identified <- suppressMessages(TMBhelper::Check_Identifiable(obj))
#
#       # Make into list
#       identified_param_list <- obj$env$parList(as.numeric(identified$BadParams$Param_check))
#       identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==0,"Not estimated",x), how = "replace")
#       identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==1,"OK",x), how = "replace")
#       identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==2,"BAD",x), how = "replace")
#
#       identified$param_list <- identified_param_list

      mod_objects$identified <- identified
    }

    class(mod_objects) <- "Rceattle"

    if(!is.null(file)){
      save(mod_objects, file = paste0(file, ".RData"))
    }

    # suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
    return(mod_objects)
  }
