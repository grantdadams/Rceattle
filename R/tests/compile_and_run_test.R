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
data_list = BS2017SS
inits = NULL # Initial parameters = 0
map = NULL
bounds = NULL
file_name = NULL # Don't save
debug = 0 # Estimate
random_rec = FALSE # No random recruitment
msmMode = 0 # Single species mode
avgnMode = 0
silent = FALSE
niter = 5
est_diet = FALSE
suitMode = FALSE
recompile = FALSE

# #--------------------------------------------------
# # 1. DATA and MODEL PREP
# #--------------------------------------------------
load("C:/Users/Grant Adams/Documents/GitHub/RceattleRuns/2019 Think Tank/Models/ss_no_re.RData")
initial_params <- mod_objects$estimated_params

initial_params$srv_sel_coff
initial_params$srv_sel_inf <- cbind(initial_params$srv_sel_inf, matrix(0, ncol = 1, nrow = 2))
initial_params$srv_sel_slp <- cbind(initial_params$srv_sel_slp, matrix(0, ncol = 1, nrow = 2))
initial_params$log_srv_q <- c(initial_params$log_srv_q, initial_params$log_eit_q)
inits <- initial_params
inits$log_eit_q <- NULL
inits$phi <- NULL

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
  TMBfilename <- "ceattle_v01_04"
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
  params <- suppressWarnings(Rceattle::build_params(
    data_list = data_list,
    inits = inits,
    TMBfilename = TMBfilename,
    cpp_directory = cpp_directory
  ))
  inits$srv_sel_coff <- params$srv_sel_coff
  inits$fsh_sel_inf <- params$fsh_sel_inf
  inits$fsh_sel_slp <- params$fsh_sel_slp
  inits$log_phi <- params$log_phi
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


Check_Identifiable(obj)

print("Step 6: Optimization complete")

# Get quantities
quantities <- obj$report(obj$env$last.par.best)

if (debug) {
  last_par <- params
} else if (random_rec == F) {
  last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))
} else{
  last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))
}

run_time = ((Sys.time() - start_time))

# # Return objects
# mod_objects <-
#   list(
#     data_list = data_list,
#     initial_params = params,
#     bounds = bounds,
#     map = map,
#     obj = obj,
#     opt = opt,
#     sdrep = opt$SD,
#     estimated_params = last_par,
#     quantities = quantities,
#     run_time = run_time
#   )
#
# class(mod_objects) <- "Rceattle"
#
# if(!is.null(file_name)){
#   save(mod_objects, file = paste0(file_name, ".RData"))
# }
#
# quantities$fsh_sel[,,1]
# mod_objects$quantities$fsh_sel
#
# quantities$F[,,1]
# mod_objects$quantities$F[,,1]
#
# quantities$M[,,1]
# mod_objects$quantities$M[,,1]
#
# quantities$Zed[,,1]
# mod_objects$quantities$Zed[,,1]
#
# quantities$S[,,1]
# mod_objects$quantities$S[,,1]
#
# quantities$R
# mod_objects$quantities$R
#
# quantities$NByage[,,6]
# mod_objects$quantities$NByage[,,6]
#
# quantities$biomassSSB
# mod_objects$quantities$biomassSSB
#
#
# quantities$srv_sel[,,1]
# mod_objects$quantities$srv_sel
#
#
# quantities$srv_bio_hat
# c(t(mod_objects$quantities$srv_bio_hat))
#
#
# quantities$srv_hat[37:46]
# c(t(mod_objects$quantities$srv_hat))[40:49]
#
# quantities$srv_comp_hat[1:10, 1:12]
# t(mod_objects$quantities$srv_age_hat[1,1:12,1:10])
#
# quantities$srv_comp_hat[37:46, 1:12]
# t(mod_objects$quantities$srv_age_hat[2,1:12,1:10])
#
# quantities$srv_comp_hat[73:82, 1:12]
# t(mod_objects$quantities$srv_age_hat[3,1:12,1:10])
#
# quantities$fsh_bio_hat
# mod_objects$quantities$tc_biom_hat
#
# quantities$fsh_hat
# mod_objects$quantities$tc_hat
#
# quantities$fsh_comp_hat[1:10, 1:12]
# t(mod_objects$quantities$fsh_age_hat[1,1:12,1:10])
#
# quantities$fsh_comp_hat[39:48, 1:12]
# t(mod_objects$quantities$fsh_age_hat[2,1:12,1:10])
#
# quantities$fsh_comp_hat[78:87, 1:12]
# t(mod_objects$quantities$fsh_age_hat[3,1:12,1:10])
#
# round(quantities$jnll_comp, 4)
# round(mod_objects$quantities$jnll_comp,4)
#
#
# # round(quantities$offset, 4)
# # round(mod_objects$quantities$offset_fsh,4)
# # round(mod_objects$quantities$offset_eit,4)
# # round(mod_objects$quantities$offset_srv,4)
