# This functions runs CEATTLE

#' This function estimates population parameters of CEATTLE using maximum likelihood in TMB. The currently uses data from \code{.dat} files used in the ADMB version of CEATTLE.
#'
#' @param ctlFilename The ADMB control file used for CEATTLE
#' @param TMBfilename The version of the cpp CEATTLE file found in the src folder
#' @param dat_dir The directory where dat files are stored
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#' @param data_list (Optional) a data_list from a previous model run.
#'
#' @return
#' @export
#'
#' @examples
Rceattle <- function( data_list = NULL, ctlFilename, TMBfilename, dat_dir, debug = T){
  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}

  # Load data
  source("R/1-build_dat.R")
  source("R/2-build_params.R")
  source("R/3-build_map.R")

  # STEP 1 - LOAD DATA
  if(is.null(data_list)){
    data_list <- build_dat(ctlFilename = ctlFilename, TMBfilename = TMBfilename, dat_dir = dat_dir, debug = debug)
    print("Step 1: Data build complete")
  }
  data_list$debug <- debug

  # STEP 2 - LOAD PARAMETERS
  params <- build_params(data_list, nselages = 8, incl_prev = TRUE, Rdata_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/CEATTLE_results.Rdata"),  std_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/ceattle_est.std"), TMBfilename = TMBfilename)
  print("Step 2: Parameter build complete")

  # STEP 3 - BUILD MAP
  map  <- build_map(data_list, params, debug = debug)
  print("Step 3: Map build complete")


  # Compile CEATTLE
  version <- TMBfilename
  setwd("inst")
  TMB::compile(paste0(version, ".cpp"))
  dyn.load(TMB::dynlib(paste0(version)))
  setwd("..")
  print("Step 4: Compile CEATTLE complete")

  # Build object
  obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map)
  opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
  rep = obj$report()

  # Refit - if not debugging
  if(debug == FALSE){
    for(i in 1:10){
      last_par = obj$env$parList(opt$par)
      print("Re-running model ", i)
      obj = TMB::MakeADFun(data_list, parameters = last_par,  DLL = version, map = map, silent = TRUE)
      opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
    }
  }
  rep = obj$report()


  # Return objects
  mod_objects <- list(data_list = data_list, params = params, map = map, rep = rep, obj = obj, opt = opt)
  return(mod_objects)
}
