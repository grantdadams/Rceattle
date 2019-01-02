library(Rceattle)

data_list_ms <- build_dat(ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01_02", dat_dir = "data/BSAI/BS_MS_Files/dat files/", nspp = 3)

TMBfilename = "CEATTLE_BSAI_MS_v01_02";
data_list = data_list_ms;
inits = NULL; # Initial parameters = 0
file_name = NULL; # Don't save
debug = 0; # Estimate
random_rec = FALSE; # No random recruitment
niter = 10; # Number of iterations around predation/pop dy functions
msmMode = 2; # Multi-species holsman mode
avgnMode = 0


setwd(getwd())

#--------------------------------------------------
# 1. DATA and MODEL PREP
#--------------------------------------------------
library(TMB)
library(TMBhelper)


# Switches
data_list$random_rec <- as.numeric(random_rec)
data_list$debug <- debug
data_list$niter <- niter
data_list$avgnMode <- avgnMode
data_list$msmMode <- msmMode

# STEP 2 - LOAD PARAMETERS
  params <- Rceattle::build_params(
    data_list = data_list,
    nselages = 8,
    inits = inits,
    TMBfilename = TMBfilename
  )





# STEP 3 - BUILD MAP
map  <-
  Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec)



# STEP 4 - Setup random effects
random_vars <- c()
if (random_rec == TRUE) {
  random_vars <- c("rec_dev")
}


# STEP 5 - Compile CEATTLE
version <- TMBfilename
cpp_directory <- "inst"
cpp_file <- paste0(cpp_directory, "/", version)

# Remove compiled files if not compatible with system
version_files <-
  list.files(path = cpp_directory, pattern = version)
if (Sys.info()[1] == "Windows" &
    paste0(version, ".so") %in% version_files) {
  try(dyn.unload(TMB::dynlib(paste0(cpp_file))))
  file.remove(paste0(cpp_file, ".so"))
  file.remove(paste0(cpp_file, ".o"))
}
if (Sys.info()[1] != "Windows" &
    paste0(version, ".dll") %in% version_files) {
  try(dyn.unload(TMB::dynlib(paste0(cpp_file))))
  file.remove(paste0(cpp_file, ".dll"))
  file.remove(paste0(cpp_file, ".o"))
}

TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)))
print("Step 4: Compile CEATTLE complete")



# STEP 6 - Build and fit model object
obj = TMB::MakeADFun(
  data_list,
  parameters = params,
  DLL = version,
  map = map,
  random = random_vars,
  silent = FALSE
)
