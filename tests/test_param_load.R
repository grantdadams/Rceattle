# This functions tests the loading of parameters into TMB.

param_load <- function( ctlFilename = "asmnt2017_0", TMBfilename = "CEATTLE_BSAI_v01", dat_dir =  "data/dat files/"){
  version <- "param_load"
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


  #--------------------------------------------------
  # 2. Build CPP file data entry
  #--------------------------------------------------
  cpp_fn<-file(paste("src/", TMBfilename,".cpp",sep=""))
  cpp_file <- readLines(cpp_fn)
  nrow <- grep('6.2 EIT Components', cpp_file) # Last line of data files
  cpp_file <- cpp_file[1:(nrow - 1)] # Retain data section
  cpp_file <- c(cpp_file, "return 0;", "}", "")
  writeLines(cpp_file, con = file(paste0("tests/", version, ".cpp")))


  #--------------------------------------------------
  # 3. Run TMB with data only
  #--------------------------------------------------
  library(TMBdebug)
  library(TMB)
  setwd("tests")
  compile(paste0(version, ".cpp"))
  dyn.load(dynlib(version))
  setwd("../")

  # Build object
  Obj = TMBdebug::MakeADFun(data_list, parameters = params,  DLL = version)
  Rep <- Obj$report()
  Rep$srv_sel
}
