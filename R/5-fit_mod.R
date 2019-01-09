#' This functions runs CEATTLE
#' @description  This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#'
#' @param TMBfilename The version of the cpp CEATTLE file found in the src folder
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param inits Character vector of named initial values from ADMB or list of previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters.
#' @param file_name Filename where files will be saved. If NULL, no file is saved.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE.
#' @param niter Number of iterations for multispecies model
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param avgnMode The average abundance-at-age approximation to be used for predation mortality equations. 0 (default) is the \eqn{\frac{N}{Z} \left( 1 - exp^{-Z} \right)}, 1 is \eqn{N e^{-Z/2}}, 2 is \eqn{N}.
#' @param silent logical.  IF TRUE, includes TMB estimation progress
#' @param est_diet logical. If FALSE, does not include diet in the likelihood.The default is FALSE.
#' @param suitMode logical. If FALSE, does not estimate suitability parameters. If TRUE, estimates gamma selectivity parameters. The default is FALSE.
#'
#' @details
#' CEATTLE is an age-structured population dynamics model that can be fit with or without predation mortality. The default is to exclude predation mortality by setting \code{msmMode} to 0. Predation mortality can be included by setting \code{msmMode} with the following options:
#' \describe{
#' \item{0. Single species mode}
#' \item{1. Holsman et al. 2015 predation based on multi-species virtual population analysis (MSVPA) based predation formation.}
#'   \item{2. Kinzey & Punt 2010 Holling Type I (linear)}
#'   \item{3. Kinzey & Punt 2010 Holling Type II}
#'   \item{4. Kinzey & Punt 2010 Holling Type III}
#'   \item{5. Kinzey & Punt 2010 Predator interference}
#'   \item{6. Kinzey & Punt 2010 Predator preemption}
#'   \item{7. Kinzey & Punt 2010 Hassell-Varley}
#'   \item{8. Kinzey & Punt 2010 Ecosim}
#'   }

Rceattle <-
  function(TMBfilename = "CEATTLE_BSAI_MS_v01_02",
           data_list = NULL,
           inits = NULL,
           file_name = NULL,
           debug = T,
           random_rec = FALSE,
           niter = 3,
           msmMode = 0,
           avgnMode = 0,
           est_diet = FALSE,
           suitMode = FALSE,
           silent = FALSE) {
    start_time <- Sys.time()

    setwd(getwd())

    #--------------------------------------------------
    # 1. DATA and MODEL PREP
    #--------------------------------------------------
    # Check if require packages are installed and install if not
    if ("TMB" %in% rownames(installed.packages()) == FALSE) {
      install.packages("TMB")
    }
    if ("TMBhelper" %in% rownames(installed.packages()) == FALSE) {
      install.packages("TMBhelper")
    }
    library(TMB)
    library(TMBhelper)

    # Load data
    # source("R/2-build_params.R")
    # source("R/3-build_map.R")


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
    data_list$suitMode <- suitMode
    data_list$est_diet <- est_diet


    # STEP 1 - LOAD PARAMETERS
    if (is.character(inits) | is.null(inits)) {
      params <- Rceattle::build_params(
        data_list = data_list,
        nselages = 8,
        inits = inits,
        TMBfilename = TMBfilename
      )
    } else{
      params <- inits
    }
    print("Step 1: Parameter build complete")



    # STEP 2 - BUILD MAP
    map  <-
      Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec)
    print("Step 2: Map build complete")


    # STEP 3 - Get bounds
    bounds <- Rceattle::build_bounds(param_list = params)



    # STEP 4 - Setup random effects
    random_vars <- c()
    if (random_rec == TRUE) {
      random_vars <- c("rec_dev")
    }


    # STEP 5 - Compile CEATTLE
    version <- TMBfilename
    cpp_directory <- "inst"
    cpp_file <- paste0(cpp_directory, "/", version)

    # Remove compiled files if not compatible with system
    version_files <-
      list.files(path = cpp_directory, pattern = version)
    if (Sys.info()[1] == "Windows" &
        paste0(version, ".so") %in% version_files) {
      try(dyn.unload(TMB::dynlib(paste0(cpp_file))))
      file.remove(paste0(cpp_file, ".so"))
      file.remove(paste0(cpp_file, ".o"))
    }
    if (Sys.info()[1] != "Windows" &
        paste0(version, ".dll") %in% version_files) {
      try(dyn.unload(TMB::dynlib(paste0(cpp_file))))
      file.remove(paste0(cpp_file, ".dll"))
      file.remove(paste0(cpp_file, ".o"))
    }

    TMB::compile(paste0(cpp_file, ".cpp"))
    dyn.load(TMB::dynlib(paste0(cpp_file)))
    print("Step 4: Compile CEATTLE complete")



    # STEP 6 - Build and fit model object
    obj = TMB::MakeADFun(
      data_list,
      parameters = params,
      DLL = version,
      map = map,
      random = random_vars,
      silent = silent

    )
    print(paste0("Step 5: Optimizing model"), hessian = TRUE)
    # opt <- nlminb(obj$par, obj$fn, obj$gr)
    # methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
    # opt_list <- list()
    # for(i in 1:length(methods)){
    #   opt_list[i] = optimx(obj$par, function(x) as.numeric(obj$fn(x)), obj$gr, control = list(maxit = 10000), method = methods[i])
    # }
    opt = tryCatch(
      TMBhelper::Optimize(obj),
      lower = bounds$lower,
      upper = bounds$upper,
      error = function(e)
        NULL,
      loopnum = 3
    )

    # Get quantities
    sdrep = TMB::sdreport(obj)
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

    run_time = paste((Sys.time() - start_time) / 60)

    # Return objects
    mod_objects <-
      list(
        data_list = data_list,
        initial_params = params,
        estimated_params = last_par,
        map = map,
        sdrep = sdrep,
        obj = obj,
        opt = opt,
        quantities = quantities,
        bounds = bounds,
        run_time = run_time
      )

    if(!is.null(file_name)){
      save(mod_objects, file = paste0(file_name, ".RData"))
    }

    # dyn.unload(TMB::dynlib(paste0(cpp_file)))
    return(mod_objects)
  }
