# This functions runs CEATTLE

Rceattle <- function( ctlFilename ){
  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}
  # Load required packages


  source("R/build_dat.R")
  data_list <- build_dat(ctlFilename = "asmnt2017_0", TMBfilename = "CEATTLE_BSAI_v01", dat_dir = "data/dat files/")



  # Compile CEATTLE
  library(TMBdebug)
  setwd("tests")
  version <- "data_load"
  compile(paste0(version, ".cpp"))
  dyn.load(dynlib(version))
  setwd("../")

  # Build object
  Obj = TMBdebug::MakeADFun(data_list, parameters = list(),  DLL = version, type = "Fun")


}
