library(Rceattle)
library(TMB)
library(TMBhelper)
library(TMBdebug)

# Load data
# source("R/2-build_params.R")
# source("R/3-build_map.R")
TMBfilename = "ceattle_v01_05"
cpp_directory <- "inst/executables"
data("BS2017SS")
BS2017SS$endyr <- 2014
data_list = BS2017SS
inits = NULL # Initial parameters = 0
map = NULL
bounds = NULL
file_name = NULL # Don't save
debug = 1 # Estimate
random_rec = FALSE # No random recruitment
msmMode = 2 # Type 1 species mode
avgnMode = 0
silent = TRUE
niter = 5
est_diet = FALSE
suitMode = 1
recompile = FALSE



# #--------------------------------------------------
# # 1. DATA and MODEL PREP
# #--------------------------------------------------

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
data_list$minNByage <- 0



# STEP 1 - LOAD PARAMETERS
params <- suppressWarnings(Rceattle::build_params(
  data_list = data_list,
  inits = NULL
))


# STEP 2 - BUILD MAP
map <-
  suppressWarnings(Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec))

print("Step 2: Map build complete")


# STEP 3 - Get bounds
bounds <- Rceattle::build_bounds(param_list = params, data_list)


# STEP 4 - Setup random effects
random_vars <- c("ln_srv_q_dev_re", "srv_sel_slp_dev_re", "srv_sel_inf_dev_re", "fsh_sel_slp_dev_re", "fsh_sel_inf_dev_re")
if (random_rec == TRUE) {
  random_vars <- c(random_vars , "rec_dev")
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

# STEP 7 - Build and fit model object
obj = TMB::MakeADFun(
  data_list2,
  parameters = params,
  DLL = TMBfilename,
  random = random_vars,
  silent = silent,
  map = map[[1]]
)
print(paste0("Step 5: Build object complete"), hessian = TRUE)
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
# opt_list <- list()
# for(i in 1:length(methods)){
#  opt_list[i] = optimx(obj$par, function(x) as.numeric(obj$fn(x)), obj$gr, control = list(maxit = 10000), method = methods[i])
# }


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
