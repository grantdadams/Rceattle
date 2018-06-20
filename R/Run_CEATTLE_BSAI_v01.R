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
  data_list <- build_dat(ctlFilename, TMBfilename, dat_dir)
  params <- build_params(data_list, nselages = 8)


  # Compile CEATTLE
  library(TMB)
  setwd("src")
  version <- "CEATTLE_BSAI_v01"
  compile(paste0(version, ".cpp"))
  dyn.load(dynlib(version))
  setwd("../")

  # Build object
  Obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version)


}
