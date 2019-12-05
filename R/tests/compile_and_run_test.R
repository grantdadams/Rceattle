library(Rceattle)
library(TMB)
library(TMBhelper)
library(TMBdebug)

# Load data
# source("R/2-build_params.R")
# source("R/3-build_map.R")
TMBfilename = "ceattle_v01_06"
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
msmMode = 4 # Type 1 species mode
avgnMode = 0
silent = TRUE
niter = 5
est_diet = FALSE
suitMode = FALSE
recompile = FALSE


data_list$UobsWtAge <- as.data.frame(data_list$UobsWtAge)
data_list$UobsAge <- as.data.frame(data_list$UobsAge)

# Step 0 - remove years prior to start year
data_list$wt <- data_list$wt[which(data_list$wt$Year == 0 | data_list$wt$Year >= data_list$styr),]
data_list$UobsAge <- data_list$UobsAge[which(data_list$UobsAge$Year == 0 | data_list$UobsAge$Year >= data_list$styr),]
data_list$UobsWtAge <- data_list$UobsWtAge[which(data_list$UobsWtAge$Year == 0 | data_list$UobsWtAge$Year >= data_list$styr),]
data_list$srv_biom <- data_list$srv_biom[which(data_list$srv_biom$Year >= data_list$styr),]
data_list$fsh_biom <- data_list$fsh_biom[which(data_list$fsh_biom$Year >= data_list$styr),]
data_list$comp_data <- data_list$comp_data[which(data_list$comp_data$Year >= data_list$styr),]
data_list$emp_sel <- data_list$emp_sel[which(data_list$emp_sel$Year == 0 | data_list$emp_sel$Year >= data_list$styr),]
data_list$NByageFixed <- data_list$NByageFixed[which(data_list$NByageFixed$Year == 0 | data_list$NByageFixed$Year >= data_list$styr),]
data_list$Pyrs <- data_list$Pyrs[which(data_list$Pyrs$Year == 0 | data_list$Pyrs$Year >= data_list$styr),]


# #--------------------------------------------------
# # 1. DATA and MODEL PREP
# #--------------------------------------------------


# Switches
data_list$random_rec <- as.numeric(random_rec)
data_list$debug <- debug
data_list$niter <- niter
data_list$avgnMode <- avgnMode
data_list$msmMode <- msmMode
data_list$suitMode <- as.numeric(suitMode)
data_list$minNByage <- 0


# Get cpp file if not provided
if(is.null(TMBfilename) | is.null(cpp_directory)){
  cpp_directory <- system.file("executables",package="Rceattle")
  TMBfilename <- "ceattle_v01_04"
} else{
  cpp_directory <- cpp_directory
  TMBfilename <- TMBfilename
}


# STEP 1 - LOAD PARAMETERS

params <- suppressWarnings(Rceattle::build_params(
  data_list = data_list,
  inits = inits
))
print("Step 1: Parameter build complete")



# STEP 2 - BUILD MAP
map <-
  suppressWarnings(Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec))

print("Step 2: Map build complete")


# STEP 3 - Get bounds
bounds <- Rceattle::build_bounds(param_list = params, data_list)

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
