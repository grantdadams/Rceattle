
library(TMB)
library(TMBhelper)
library(TMBdebug)

# Load data
# source("R/2-build_params.R")
# source("R/3-build_map.R")
data_list_ss$nages
TMBfilename = "ceattle_v01_02"
data_list = data_list_ss
inits = NULL # Initial parameters = 0
file_name = NULL # Don't save
debug = 1 # Estimate
random_rec = FALSE # No random recruitment
msmMode = 3 # Single species mode
avgnMode = 0
silent = FALSE
niter = 10


# Switches
data_list$random_rec <- as.numeric(random_rec)
data_list$debug <- debug
data_list$niter <- niter
data_list$avgnMode <- avgnMode
data_list$msmMode <- msmMode

# STEP 1 - LOAD PARAMETERS
if (is.character(inits) | is.null(inits)) {
  params <- Rceattle::build_params(
    data_list = data_list,
    nselages = 8,
    inits = inits,
    TMBfilename = TMBfilename
  )
} else{
  params <- inits
}
print("Step 1: Parameter build complete")



# STEP 2 - BUILD MAP
map  <-
  Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec)
print("Step 2: Map build complete")


# STEP 3 - Get bounds
bounds <- Rceattle::build_bounds(param_list = params)



# STEP 4 - Setup random effects
random_vars <- c()
if (random_rec == TRUE) {
  random_vars <- c("rec_dev")
}


# STEP 5 - Compile CEATTLE
version <- TMBfilename
cpp_directory <- "inst"
cpp_file <- paste0(cpp_directory, "/", version)


dyn.unload(TMB::dynlib(paste0(cpp_file)))
TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)))
print("Step 4: Compile CEATTLE complete")


obj = TMBdebug::MakeADFun(
  data_list,
  parameters = params,
  DLL = version,
  map = map,
  random = random_vars,
  silent = silent
)
