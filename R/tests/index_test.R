
library(TMB)
library(TMBhelper)
library(TMBdebug)

# Load data
# source("R/2-build_params.R")
# source("R/3-build_map.R")
data_list_ss$nages
TMBfilename = "ceattle_v01_02"
data("BS2017SS")
data_list = BS2017SS
inits = NULL # Initial parameters = 0
file_name = NULL # Don't save
debug = 1 # Estimate
random_rec = FALSE # No random recruitment
msmMode = 1 # Single species mode
avgnMode = 0
silent = FALSE
niter = 10
est_diet = FALSE
suitMode = FALSE

#--------------------------------------------------
# 1. DATA and MODEL PREP
#--------------------------------------------------
# # Check if require packages are installed and install if not
# if ("TMB" %in% rownames(installed.packages()) == FALSE) {
#   install.packages("TMB")
# }
# if ("TMBhelper" %in% rownames(installed.packages()) == FALSE) {
#   install.packages("TMBhelper")
# }
# library(TMB)
# library(TMBhelper)
#
# # Load data
# # source("R/2-build_params.R")
# # source("R/3-build_map.R")


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
data_list$suitMode <- suitMode
data_list$est_diet <- est_diet


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

# Remove compiled files if not compatible with system
version_files <-
  list.files(path = cpp_directory, pattern = version)
if (Sys.info()[1] == "Windows" &
    paste0(version, ".so") %in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
  file.remove(paste0(cpp_file, ".so"))
  file.remove(paste0(cpp_file, ".o"))
}
if (Sys.info()[1] != "Windows" &
    paste0(version, ".dll") %in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
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
  DLL = "Rceattle",
  map = map,
  random = random_vars,
  silent = silent

)


print(paste0("Step 5: Optimizing model"), hessian = TRUE)
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
# opt_list <- list()
# for(i in 1:length(methods)){
#   opt_list[i] = optimx(obj$par, function(x) as.numeric(obj$fn(x)), obj$gr, control = list(maxit = 10000), method = methods[i])
# }
opt = tryCatch(
  TMBhelper::Optimize(obj),
  lower = bounds$lower,
  upper = bounds$upper,
  error = function(e)
    NULL,
  loopnum = 3
)

# Get quantities
sdrep = TMB::sdreport(obj)
quantities <- obj$report(obj$env$last.par.best)

if (debug) {
  last_par <- params
}
else if (random_rec == F) {
  last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))
}
else{
  last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))
}

run_time = ((Sys.time() - start_time))

# Return objects
mod_objects <-
  list(
    data_list = data_list,
    initial_params = params,
    estimated_params = last_par,
    map = map,
    sdrep = sdrep,
    obj = obj,
    opt = opt,
    quantities = quantities,
    bounds = bounds,
    run_time = run_time
  )
