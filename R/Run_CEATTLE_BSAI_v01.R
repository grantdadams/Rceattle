# This functions runs CEATTLE

Rceattle <- function( ctlFilename, TMBfilename, dat_dir ){
  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}

  # Load data
  source("R/1-build_dat.R")
  source("R/2-build_params.R")
  source("R/3-build_map.R")

  data_list <- build_dat(ctlFilename = ctlFilename, TMBfilename = TMBfilename, dat_dir = dat_dir)
  print("Step 1: Data build complete")

  params <- build_params(data_list, nselages = 8, incl_prev = F, Rdata_file = "data/CEATTLE_results.Rdata",  std_file = "data/ceattle_est.std", TMBfilename = TMBfilename )
  print("Step 2: Parameter build complete")

  map  <- build_map(data_list, params)
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
  last_par = obj$env$parList(opt$par)
  obj = TMB::MakeADFun(data_list, parameters = last_par,  DLL = version, map = map)
  opt = Optimize( obj )
  rep = obj$report()
  #Opt$opt$diagnostics

  mod_objects <- list(data_list = data_list, params = params, map = map, rep = rep, obj = obj, opt = opt)
  return(mod_objects)
}

mod_objects <- Rceattle( ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_v01_1_0", dat_dir =  "data/dat files/" )
rep <- mod_objects$rep
data_list <- mod_objects$data_list
params <- mod_objects$params
opt <- mod_objects$opt
rep <- mod_objects$obj$report()
rep$jnll
opt$objective
