library(Rceattle)
library(TMB)
library(TMBhelper)
library(TMBdebug)

# Load data
# source("R/2-build_params.R")
# source("R/3-build_map.R")
TMBfilename = "ceattle_v01_04"
cpp_directory <- "inst/executables"
data("BS2017SS")
BS2017SS$endyr <- 2017
data_list = BS2017SS
inits = NULL # Initial parameters = 0
map = NULL
bounds = NULL
file_name = NULL # Don't save
debug = 0 # Estimate
random_rec = FALSE # No random recruitment
msmMode = 1 # Type 1 species mode
avgnMode = 0
silent = TRUE
niter = 5
minNByage = 0
est_diet = FALSE
suitMode = FALSE
recompile = FALSE



# #--------------------------------------------------
# # 1. DATA and MODEL PREP
# #--------------------------------------------------
inits <- ms_run$estimated_params

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
data_list$minNByage <- as.numeric(minNByage)


# Get cpp file if not provided
if(is.null(TMBfilename) | is.null(cpp_directory)){
  cpp_directory <- system.file("executables",package="Rceattle")
  TMBfilename <- "ceattle_v01_04"
} else{
  cpp_directory <- cpp_directory
  TMBfilename <- TMBfilename
}


# STEP 1 - LOAD PARAMETERS
if (is.character(inits) | is.null(inits)) {
  params <- suppressWarnings(Rceattle::build_params(
    data_list = data_list,
    inits = inits
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


dyn.unload(TMB::dynlib(paste0(cpp_file)))
TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)), silent = TRUE)
print("Step 4: Compile CEATTLE complete")


# STEP 6 - Reorganize data
data_list2 <- rearrange_dat(data_list)
data_list2$msmMode = 3

# STEP 7 - Build and fit model object
obj = TMB::MakeADFun(
  data_list2,
  parameters = params,
  DLL = TMBfilename,
  random = random_vars,
  silent = FALSE,
  map = map[[1]]
)

obj$gr()
obj$fn()

opt <- nlminb(objective = obj$fn,
              start=obj$par)

print(paste0("Step 5: Build object complete"), hessian = TRUE)


# Remove inactive parameters from bounds and vectorize
L = unlist(bounds$lower)[which(!is.na(unlist(map[[1]])))]
U = unlist(bounds$upper)[which(!is.na(unlist(map[[1]])))]

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

# Get quantities
quantities <- obj$report(obj$env$last.par.best)
round(quantities$jnll_comp,3)

Check_Identifiable(obj)

round(mod_objects$quantities$jnll_comp, 3)


# fit <- tmbstan(obj, chains=1)
