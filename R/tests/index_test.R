
library(TMB)
library(TMBhelper)
library(TMBdebug)
library(Rceattle)

# Load data
# source("R/2-build_params.R")
# source("R/3-build_map.R")
data_list_ss$nages
cpp_directory <- "inst/executables"
TMBfilename <- "ceattle_v01_02"
data("BS2017SS")
data_list = BS2017SS
inits = NULL # Initial parameters = 0
file_name = NULL # Don't save
debug = 0 # Estimate
random_rec = FALSE # No random recruitment
msmMode = 0 # Single species mode
avgnMode = 0
silent = FALSE
niter = 10
est_diet = FALSE
suitMode = FALSE
map <- NULL
bounds <- NULL

setwd(getwd())

#--------------------------------------------------
# 1. DATA and MODEL PREP
#--------------------------------------------------
# # Check if require packages are installed and install if not
# if ("TMB" %in% rownames(installed.packages()) == FALSE) {
#  install.packages("TMB")
# }
# if ("TMBhelper" %in% rownames(installed.packages()) == FALSE) {
#  install.packages("TMBhelper")
# }
# library(TMB)
# library(TMBhelper)


# STEP 1 - LOAD DATA
if (is.null(data_list)) {
  stop("Missing data_list object")
}


# Switches
data_list$random_rec <- as.numeric(random_rec)
data_list$debug <- debug
data_list$niter <- niter
data_list$avgnMode <- avgnMode
data_list$msmMode <- msmMode
data_list$suitMode <- as.numeric(suitMode)


# Get cpp file if not provided
if(is.null(TMBfilename) | is.null(cpp_directory)){
  cpp_directory <- system.file("executables",package="Rceattle")
  TMBfilename <- "ceattle_v01_02"
} else{
  cpp_directory <- cpp_directory
  TMBfilename <- TMBfilename
}


# STEP 1 - LOAD PARAMETERS
if (is.character(inits) | is.null(inits)) {
  params <- suppressWarnings(Rceattle::build_params(
    data_list = data_list,
    inits = inits,
    TMBfilename = TMBfilename,
    cpp_directory = cpp_directory
  ))
} else{
  params <- inits
}
print("Step 1: Parameter build complete")



# STEP 2 - BUILD MAP
if (is.null(map)) {
  map <-
    suppressWarnings(Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec))
} else{
  map <- map
}
print("Step 2: Map build complete")


# STEP 3 - Get bounds
if (is.null(bounds)) {
  bounds <- Rceattle::build_bounds(param_list = params, data_list)
} else {
  bounds = bounds
}
print("Step 3: Param bounds complete")


# STEP 4 - Setup random effects
random_vars <- c()
if (random_rec == TRUE) {
  random_vars <- c("rec_dev")
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# STEP 5 - Compile CEATTLE is providing cpp file
cpp_file <- paste0(cpp_directory, "/", TMBfilename)

# Remove compiled files if not compatible with system
version_files <-
  list.files(path = cpp_directory, pattern = TMBfilename)
if (Sys.info()[1] == "Windows" &
    paste0(TMBfilename, ".so") %in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
  suppressWarnings(file.remove(paste0(cpp_file, ".so")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}
if (Sys.info()[1] != "Windows" &
    paste0(TMBfilename, ".dll") %in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
  suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}
if (Sys.info()[1] != "Windows" &
    paste0(TMBfilename, ".so") %!in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
  suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}
if(recompile){
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
  suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
  suppressWarnings(file.remove(paste0(cpp_file, ".so")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}

dyn.unload(TMB::dynlib(paste0(cpp_file)), silent = TRUE)
TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)), silent = TRUE)
print("Step 4: Compile CEATTLE complete")



# STEP 6 - Build and fit model object
obj = TMB::MakeADFun(
  data_list,
  parameters = params,
  DLL = TMBfilename,
  map = map,
  random = random_vars,
  silent = silent

)



# Remove inactive parameters from bounds and vectorize
L = unlist(bounds$lower)[which(!is.na(unlist(map)))]
U = unlist(bounds$upper)[which(!is.na(unlist(map)))]

# Remove random effects from bounds
if (random_rec == TRUE) {
  L <- L[-grep(random_vars, names(L))]
  U <- U[-grep(random_vars, names(U))]
}

# Optimize
opt = Rceattle::Optimize(obj = obj,
                         fn=obj$fn,
                         gr=obj$gr,
                         startpar=obj$par,
                         lower = L,
                         upper = U,
                         loopnum = 5
)
Check_Identifiable(obj)
