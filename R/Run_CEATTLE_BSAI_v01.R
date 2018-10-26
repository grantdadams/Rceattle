# This functions runs CEATTLE

#' This function estimates population parameters of CEATTLE using maximum likelihood in TMB. The currently uses data from \code{.dat} files used in the ADMB version of CEATTLE.
#'
#' @param ctlFilename The ADMB control file used for CEATTLE
#' @param TMBfilename The version of the cpp CEATTLE file found in the src folder
#' @param dat_dir The directory where dat files are stored
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#' @param data_list (Optional) a data_list from a previous model run.
#' @param inits Boolian of wether to have the model start optimization at ADMBs MLEs \code{TRUE} or all 0s \code{FALSE}
#'
#' @return
#' @export
#'
#' @examples
Rceattle <- function( data_list = NULL, ctlFilename, TMBfilename, dat_dir, debug = T, inits = NULL){
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
  source("R/2-build_params_par_file.R")
  source("R/3-build_map.R")

  # STEP 1 - LOAD DATA
  if(is.null(data_list)){
<<<<<<< HEAD
    data_list <- build_dat(ctlFilename = ctlFilename,
                           TMBfilename = TMBfilename,
                           dat_dir = dat_dir,
                           nspp = 3,
                           debug = debug)
=======
    data_list <- build_dat(ctlFilename = ctlFilename, TMBfilename = TMBfilename, dat_dir = dat_dir, debug = debug)
>>>>>>> parent of 17f0802... Added par file into parameter creation
    print("Step 1: Data build complete")
  }
  data_list$debug <- debug

  # STEP 2 - LOAD PARAMETERS
<<<<<<< HEAD
  params <- build_params(data_list = data_list,
                         nselages = 8,
                         incl_prev = ifelse(is.null(inits), FALSE, TRUE),
                         Rdata_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/CEATTLE_results.Rdata"),
                         param_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/", inits),
                         TMBfilename = TMBfilename)
=======
  params <- build_params(data_list, nselages = 8, incl_prev = inits, Rdata_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/CEATTLE_results.Rdata"),  std_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/ceattle.par"), TMBfilename = TMBfilename)
>>>>>>> parent of 17f0802... Added par file into parameter creation
  print("Step 2: Parameter build complete")

  # STEP 3 - BUILD MAP
  map  <- build_map(data_list, params, debug = debug)
  print("Step 3: Map build complete")


  # STEP 4 - Compile CEATTLE
  version <- TMBfilename
  cpp_directory <- "inst"
  cpp_file <- paste0(cpp_directory, "/", version)

  # Remove compiled files if not compatible with system
  version_files <- list.files(path = cpp_directory, pattern = version)
  if(Sys.info()[1] == "Windows" & paste0(version,".so") %in% version_files){
    file.remove(paste0(cpp_file,".so"))
    file.remove(paste0(cpp_file,".o"))
  }
  if(Sys.info()[1] != "Windows" & paste0(version,".dll") %in% version_files){
    file.remove(paste0(cpp_file,".dll"))
    file.remove(paste0(cpp_file,".o"))
  }

  TMB::compile(paste0(cpp_file, ".cpp"))
  dyn.load(TMB::dynlib(paste0(cpp_file)))
  print("Step 4: Compile CEATTLE complete")

  # STEP 5 - Build object
  obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map)
  opt = TMBhelper::Optimize( obj ) ; #tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
  rep = obj$report()

  # Refit - if not debugging
  if(debug == TRUE){ iter = 1}
  if(debug == FALSE){ iter = 3}
  for(i in 1:iter){
    last_par = obj$env$parList(opt$par)
    last_par$dummy = 0
    print("Re-running model ", i)
    obj = TMB::MakeADFun(data_list, parameters = last_par,  DLL = version, map = map, silent = TRUE)
    opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
  }

  obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map)
  opt = TMBhelper::Optimize( obj ) ; #tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
  rep = obj$report()


  # Return objects
  mod_objects <- list(data_list = data_list, params = params, map = map, rep = rep, obj = obj, opt = opt)
  return(mod_objects)
}
