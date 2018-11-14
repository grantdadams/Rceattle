# This functions runs CEATTLE

#' This function estimates population parameters of CEATTLE using maximum likelihood in TMB. The currently uses data from \code{.dat} files used in the ADMB version of CEATTLE.
#'
#' @param ctlFilename The ADMB control file used for CEATTLE
#' @param TMBfilename The version of the cpp CEATTLE file found in the src folder
#' @param dat_dir The directory where dat files are stored
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#' @param data_list (Optional) a data_list from a previous model run.
#' @param inits Character vector of named initial values from ADMB or list of previous parameter estimates from single species model.
#' @param plot_trajectory Boolian of whether to include plotting functions
#' @param random_rec Boolian of whether to treat recruitment deviations as random effects.
#' @param niter Number of iterations for multispecies model
#' @param file_name Filename where files will be saved
#'
#' @return
#' @export
#'
#' @examples

Rceattle <- function(data_list = NULL, ctlFilename = NULL, TMBfilename = NULL, dat_dir = NULL, inits = NULL, debug = T, random_rec = FALSE, niter = 3, file_name = NULL){

  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}
  if("TMBhelper" %in% rownames(installed.packages()) == FALSE) {install.packages("TMBhelper")}
  library(TMB)
  library(TMBhelper)

  # Load data
  source("R/1-build_dat.R")
  source("R/2-build_params.R")
  source("R/3-build_map.R")



  # STEP 1 - LOAD DATA
  if(is.null(data_list)){
    data_list <- build_dat(ctlFilename = ctlFilename,
                           TMBfilename = TMBfilename,
                           dat_dir = dat_dir,
                           nspp = 3,
                           debug = debug,
                           random_rec = random_rec)
    print("Step 1: Data build complete")
  }
  data_list$random_rec <- as.numeric(random_rec)
  data_list$debug <- debug
  data_list$niter <- niter



  # STEP 2 - LOAD PARAMETERS
  if(is.character(inits) | is.null(inits)){
    params <- build_params(data_list = data_list,
                           nselages = 8,
                           incl_prev = ifelse(is.null(inits), FALSE, TRUE),
                           Rdata_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/CEATTLE_results.Rdata"),
                           param_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/", inits),
                           TMBfilename = TMBfilename)
  } else{
    params <- inits
  }
  print("Step 2: Parameter build complete")



  # STEP 3 - BUILD MAP
  map  <- build_map(data_list, params, debug = debug, random_rec = random_rec)
  print("Step 3: Map build complete")



  # STEP 4 - Setup random effects
  random_vars <- c()
  if(random_rec == TRUE){
    random_vars <- c("rec_dev")
  }



  # STEP 5 - Compile CEATTLE
  version <- TMBfilename
  cpp_directory <- "inst"
  cpp_file <- paste0(cpp_directory, "/", version)

  # Remove compiled files if not compatible with system
  version_files <- list.files(path = cpp_directory, pattern = version)
  if(Sys.info()[1] == "Windows" & paste0(version,".so") %in% version_files){
    try(  dyn.unload(TMB::dynlib(paste0(cpp_file))))
    file.remove(paste0(cpp_file,".so"))
    file.remove(paste0(cpp_file,".o"))
  }
  if(Sys.info()[1] != "Windows" & paste0(version,".dll") %in% version_files){
    try(  dyn.unload(TMB::dynlib(paste0(cpp_file))))
    file.remove(paste0(cpp_file,".dll"))
    file.remove(paste0(cpp_file,".o"))
  }

  TMB::compile(paste0(cpp_file, ".cpp"))
  dyn.load(TMB::dynlib(paste0(cpp_file)))
  print("Step 4: Compile CEATTLE complete")



  # STEP 6 - Build object

  obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map, random = random_vars, silent = FALSE)
  print(paste0("Optimizing model"))
  # opt <- nlminb(obj$par, obj$fn, obj$gr)
  # methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
  # opt_list <- list()
  # for(i in 1:length(methods)){
  #   opt_list[i] = optimx(obj$par, function(x) as.numeric(obj$fn(x)), obj$gr, control = list(maxit = 10000), method = methods[i])
  # }
  opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)



  # STEP 7 -  Refit - if not debugging

  if(debug == FALSE){
    iter <- 2
    for(i in 1:iter){
      last_par = obj$env$parList(obj$env$last.par.best)
      last_par$dummy = 0
      print(paste0("Re-running model ", i))
      obj = TMB::MakeADFun(data_list, parameters = last_par,  DLL = version, map = map, random = random_vars, silent = TRUE)
      opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
    }
  }
  # Unload model
  rep = TMB::sdreport(obj)
  dyn.unload(TMB::dynlib(paste0(cpp_file)))


  # Return objects
  mod_objects <- list(data_list = data_list, params = params, map = map, rep = rep, obj = obj, opt = opt)
  save(mod_objects, file = paste0(file_name, ".RData"))

  return(mod_objects)
}
