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

  params <- build_params(data_list, nselages = 8, incl_prev = T, Rdata_file = "data/CEATTLE_results.Rdata",  std_file = "data/ceattle_est.std", TMBfilename = TMBfilename )
  print("Step 2: Parameter build complete")

  map  <- build_map(data_list, params)
  print("Step 3: Map build complete")

  # Compile CEATTLE
  library(TMBhelper)
  library(TMB)
  version <- "CEATTLE_BSAI_v01"
  compile(paste0("src/", version, ".cpp"))
  dyn.load(dynlib(paste0("src/", version)))
  print("Step 4: Compile CEATTLE complete")

  # Build object
  Obj = TMB::MakeADFun(data_list, parameters = params,  DLL = version, map = map)
  rep = Obj$report()
  #Opt = Optimize( Obj )
  #Opt$opt$diagnostics

  return(rep)

}

rep <- Rceattle( ctlFilename = "asmnt2017_0", TMBfilename = "CEATTLE_BSAI_v01", dat_dir =  "data/dat files/" )
