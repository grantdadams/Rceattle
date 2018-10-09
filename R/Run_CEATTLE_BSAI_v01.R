# This functions runs CEATTLE

#' This function estimates population parameters of CEATTLE using maximum likelihood in TMB. The currently uses data from \code{.dat} files used in the ADMB version of CEATTLE.
#'
#' @param ctlFilename The ADMB control file used for CEATTLE
#' @param TMBfilename The version of the cpp CEATTLE file found in the src folder
#' @param dat_dir The directory where dat files are stored
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values.
#'
#' @return
#' @export
#'
#' @examples
Rceattle <- function( ctlFilename, TMBfilename, dat_dir, debug = T){
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
  data_list <- build_dat(ctlFilename = ctlFilename, TMBfilename = TMBfilename, dat_dir = dat_dir, debug = debug)
  print("Step 1: Data build complete")

  # STEP 2 - LOAD PARAMETERS
  params <- build_params(data_list, nselages = 8, incl_prev = TRUE, Rdata_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/CEATTLE_results.Rdata"),  std_file = paste0(strsplit(dat_dir, "/dat")[[1]][1], "/ceattle_est.std"), TMBfilename = TMBfilename)
   print("Step 2: Parameter build complete")

  # STEP 3 - BUILD MAP
  map  <- build_map(data_list, params, debug = debug)
  print("Step 3: Map build complete")


  # Compile CEATTLE
  library(TMBhelper)
  library(TMB)
  version <- TMBfilename
  setwd("src")
  compile(paste0(version, ".cpp"))
  dyn.load(dynlib(paste0(version)))
  setwd("..")
  print("Step 4: Compile CEATTLE complete")

  # Build object
  obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map)
  opt = Optimize( obj )

  # Refit
  for(i in 1:10){
    last_par = obj$env$parList(opt$par)
    obj = TMB::MakeADFun(data_list, parameters = last_par,  DLL = version, map = map)
    opt = Optimize( obj )
  }
  rep = obj$report()
  #Opt$opt$diagnostics

  mod_objects <- list(data_list = data_list, params = params, map = map, rep = rep, obj = obj, opt = opt)
  return(mod_objects)
}

mod_objects <- Rceattle( ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", debug = TRUE)
rep <- mod_objects$rep
data_list <- mod_objects$data_list
params <- mod_objects$params
opt <- mod_objects$opt
rep <- mod_objects$obj$report()
rep$jnll
opt$objective
