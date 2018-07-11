# This functions runs CEATTLE

Rceattle <- function( ctlFilename = "asmnt2017_0", TMBfilename = "CEATTLE_BSAI_v01", dat_dir =  "data/dat files/" ){
  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}

  # Load data
  source("R/1-build_dat.R")
  source("R/2-build_params.R")
  source("R/3-build_map.R")
  data_list <- build_dat(ctlFilename, TMBfilename, dat_dir)
  params <- build_params(data_list, nselages = 8, incl_prev = T, Rdata_file = "data/CEATTLE_results.Rdata",  std_file = "data/ceattle_est.std")
  map  <- build_map(data_list, params)


  # Compile CEATTLE
  library(TMBhelper)
  library(TMB)
  setwd("src")
  version <- "CEATTLE_BSAI_v01"
  compile(paste0(version, ".cpp"))
  dyn.load(dynlib(version))
  setwd("../")

  # Build object
  Obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map)
  Rep = Obj$report()
  #Opt = Optimize( Obj )
  #Opt$opt$diagnostics

  return(Rep)

}

rep <- Rceattle( )
