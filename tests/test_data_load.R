# This functions tests the loading of data into TMB.

data_load <- function( ctlFilename = "asmnt2017_0", TMBfilename = "CEATTLE_BSAI_v01", dat_dir =  "data/dat files/"){
  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}
  # Load data
  source("R/1-build_dat.R")
  data_list <- build_dat(ctlFilename, TMBfilename, dat_dir)


  #--------------------------------------------------
  # 2. Build CPP file data entry
  #--------------------------------------------------
  cpp_fn<-file(paste("src/", TMBfilename,".cpp",sep=""))
  cpp_file <- readLines(cpp_fn)
  nrow <- grep('PARAMETER SECTION', cpp_file) # Last line of data files
  cpp_file <- cpp_file[1:(nrow - 2)] # Retain data section
  cpp_file <- c(cpp_file, "return 0;", "}", "")
  writeLines(cpp_file, con = file("tests/data_load.cpp"))


  #--------------------------------------------------
  # 3. Run TMB with data only
  #--------------------------------------------------
  library(TMBdebug)
  library(TMB)
  setwd("tests")
  version <- "data_load"
  compile(paste0(version, ".cpp"))
  dyn.load(dynlib(version))
  setwd("../")

  # Build object
  Obj = MakeADFun(data_list, parameters = list(),  DLL = version, type = "Fun")
}
