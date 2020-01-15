library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data



TMBfilename = "ceattle_v01_06";
cpp_directory = NULL;
data_list = NULL;
inits = NULL;
map = NULL;
bounds = NULL;
file = NULL;
debug = T;
random_rec = FALSE;
niter = 3;
msmMode = 0;
avgnMode = 0;
minNByage = 0;
suitMode = 0;
phase = NULL;
silent = FALSE;
recompile = FALSE


data_list = BS2017SS
TMBfilename = "ceattle_v01_06";
cpp_directory = "inst/executables";
inits = NULL; # Initial parameters = 0
file = NULL; # Don't save
debug = 0; # Estimate
random_rec = FALSE; # No random recruitment
msmMode = 0; # Single species mode
phase = "default";
silent = TRUE;
recompile = TRUE


setwd(getwd())

#--------------------------------------------------
# 1. DATA and MODEL PREP
#--------------------------------------------------

# STEP 1 - LOAD DATA
if (is.null(data_list)) {
  stop("Missing data_list object")
}


# Remove years of data previous to start year
data_list$UobsWtAge <- as.data.frame(data_list$UobsWtAge)
data_list$UobsAge <- as.data.frame(data_list$UobsAge)

data_list$wt <- data_list$wt[which(data_list$wt$Year == 0 | data_list$wt$Year >= data_list$styr),]
data_list$UobsAge <- data_list$UobsAge[which(data_list$UobsAge$Year == 0 | data_list$UobsAge$Year >= data_list$styr),]
data_list$UobsWtAge <- data_list$UobsWtAge[which(data_list$UobsWtAge$Year == 0 | data_list$UobsWtAge$Year >= data_list$styr),]
data_list$srv_biom <- data_list$srv_biom[which(data_list$srv_biom$Year >= data_list$styr),]
data_list$fsh_biom <- data_list$fsh_biom[which(data_list$fsh_biom$Year >= data_list$styr),]
data_list$comp_data <- data_list$comp_data[which(data_list$comp_data$Year >= data_list$styr),]
data_list$emp_sel <- data_list$emp_sel[which(data_list$emp_sel$Year == 0 | data_list$emp_sel$Year >= data_list$styr),]
data_list$NByageFixed <- data_list$NByageFixed[which(data_list$NByageFixed$Year == 0 | data_list$NByageFixed$Year >= data_list$styr),]
data_list$Pyrs <- data_list$Pyrs[which(data_list$Pyrs$Year == 0 | data_list$Pyrs$Year >= data_list$styr),]


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
  TMBfilename <- "ceattle_v01_06"
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
  inits$proj_F <- data_list$fleet_control$proj_F
  params <- inits
}
message("Step 1: Parameter build complete")



# STEP 2 - BUILD MAP
if (is.null(map)) {
  map <-
    suppressWarnings(Rceattle::build_map(data_list, params, debug = debug, random_rec = random_rec))
} else{
  map <- map
}
message("Step 2: Map build complete")


# STEP 3 - Get bounds
if (is.null(bounds)) {
  bounds <- Rceattle::build_bounds(param_list = params, data_list)
} else {
  bounds = bounds
}
message("Step 3: Param bounds complete")


# STEP 4 - Setup random effects
random_vars <- c("ln_srv_q_dev_re", "sel_slp_dev_re", "sel_inf_dev_re")
if (random_rec == TRUE) {
  random_vars <- c(random_vars , "rec_dev")
}

'%!in%' <- function(x,y)!('%in%'(x,y))



# Set default phasing
if(class(phase) == "character"){
  if(tolower(phase) == "default"){
    phase = list(
      dummy = 1,
      ln_pop_scalar = 11,
      ln_mn_rec = 1,
      ln_rec_sigma = 2,
      rec_dev = 2,
      init_dev = 2,
      ln_mean_F = 1,
      proj_F = 1,
      proj_F_prop = 1,
      F_dev = 1,
      log_srv_q = 5,
      srv_q_pow = 6,
      ln_srv_q_dev = 7,
      ln_srv_q_dev_re = 7,
      ln_sigma_srv_q = 8,
      sel_coff = 3,
      sel_slp = 3,
      sel_inf = 3,
      sel_slp_dev = 4,
      sel_inf_dev = 4,
      sel_slp_dev_re = 4,
      sel_inf_dev_re = 4,
      ln_sigma_sel = 4,
      ln_sigma_srv_index = 2,
      ln_sigma_fsh_catch = 2,
      logH_1 = 10,
      logH_1a = 10,
      logH_1b = 10,
      logH_2 = 10,
      logH_3 = 10,
      H_4 = 10,
      log_gam_a = 9,
      log_gam_b = 9,
      log_phi = 9
    )
  }
}

if(class(phase) == "character"){
  if(tolower(phase) != "default"){
    warning("phase misspecified: please set to 'default' or list with the same order as parameters.")
  }
}


# STEP 5 - Compile CEATTLE is providing cpp file
cpp_file <- paste0(cpp_directory, "/", TMBfilename)

# Remove compiled files if not compatible with system
version_files <-
  list.files(path = cpp_directory, pattern = TMBfilename)
if (Sys.info()[1] == "Windows" &
    paste0(TMBfilename, ".so") %in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE))
  suppressWarnings(file.remove(paste0(cpp_file, ".so")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}
if (Sys.info()[1] != "Windows" &
    paste0(TMBfilename, ".dll") %in% version_files) {
  suppressMessages(suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE)))
  suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}
if (Sys.info()[1] != "Windows" &
    paste0(TMBfilename, ".so") %!in% version_files) {
  suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE))
  suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}
if(recompile){
  suppressMessages(suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file))), silent = TRUE)))
  suppressWarnings(file.remove(paste0(cpp_file, ".dll")))
  suppressWarnings(file.remove(paste0(cpp_file, ".so")))
  suppressWarnings(file.remove(paste0(cpp_file, ".o")))
}

old_wd <- getwd()
setwd(cpp_directory)
TMB::compile(paste0(TMBfilename, ".cpp"))
dyn.load(TMB::dynlib(paste0(TMBfilename)), silent = TRUE)
setwd(old_wd)


message("Step 4: Compile CEATTLE complete")


# STEP 6 - Reorganize data and build model object
Rceattle:::data_check(data_list)
data_list_reorganized <- Rceattle::rearrange_dat(data_list)


# STEP 7 - Set up parameter bounds
L <- c()
U <- c()
for(i in 1:length(map[[1]])){
  L = c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map[[1]][[i]])) & !duplicated(unlist(map[[1]][[i]])))])
  U = c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map[[1]][[i]])) & !duplicated(unlist(map[[1]][[i]])))])
}

# Remove random effects from bounds
for(i in 1:length(random_vars)){
  remove_bounds <- grep(random_vars[i], names(L))
  if(length(remove_bounds) > 0){
    L <- L[-remove_bounds]
    U <- U[-remove_bounds]
  }
}


# STEP 8 - Fit model object
step = 5
# If phased
if(!is.null(phase)){
  phase_pars <- Rceattle::TMBphase(
    data = data_list_reorganized,
    parameters = params,
    map = map[[1]],
    random = random_vars,
    phases = phase,
    cpp_directory = cpp_directory,
    model_name = TMBfilename,
    silent = silent
  )

  start_par <- phase_pars

  message(paste0("Step ", step,": Phasing complete - getting final estimates"))
  step = step + 1
}

# Not phased
if(is.null(phase)){
  start_par <- params
}


# STEP 9 - Fit final model
obj = TMB::MakeADFun(
  data_list_reorganized,
  parameters = start_par,
  DLL = TMBfilename,
  map = map[[1]],
  random = random_vars,
  silent = silent
)

message(paste0("Step ",step, ": final build complete"))
step = step + 1

# Optimize
opt = Rceattle::fit_tmb(obj = obj,
                        fn=obj$fn,
                        gr=obj$gr,
                        startpar=obj$par,
                        lower = L,
                        upper = U,
                        loopnum = 5,
                        control = list(eval.max = 1e+09,
                                       iter.max = 1e+09, trace = 0)
)

est_params <- suppressWarnings(obj$env$parList(obj$env$last.par.best))


message("Step ",step, ": Final optimization complete")


# STEP 10 - Get Fspr
f_spr <- data.frame(matrix(NA, nrow = 3, ncol = data_list$nspp))
f_spr[1,] <- data_list$spnames
f_spr[2,] <- exp(est_params$ln_mn_rec)




# STEP 11 - Get quantities
quantities <- obj$report(obj$env$last.par.best)

# Rename jnll
colnames(quantities$jnll_comp) <- paste0("Sp/Srv/Fsh_", 1:ncol(quantities$jnll_comp))
rownames(quantities$jnll_comp) <- c(
  "Survey biomass",
  "Total catch",
  "Age/length composition data",
  "Non-parametric selectivity",
  "Selectivity random walk deviates",
  "Selectivity random effect deviates",
  "Selectivity normalization",
  "Catchability random walk deviates",
  "Catchability random effect deviates",
  "Recruitment deviates",
  "Initial abundance deviates",
  "Fishing mortality deviates",
  "Empty",
  "Ration",
  "Ration penalties",
  "Stomach content weight",
  "Stomach content numbers"
)


colnames(quantities$biomassSSB) <- data_list$styr:data_list$projyr
colnames(quantities$R) <- data_list$styr:data_list$projyr

rownames(quantities$biomassSSB) <- data_list$spnames
rownames(quantities$R) <- data_list$spnames


if (debug) {
  last_par <- params
}
else{
  last_par = suppressWarnings(est_params)
}

run_time = ((Sys.time() - start_time))


# Return objects
mod_objects <-
  list(
    data_list = data_list,
    initial_params = params,
    bounds = bounds,
    map = map,
    obj = obj,
    opt = opt,
    sdrep = opt$SD,
    estimated_params = last_par,
    quantities = quantities,
    run_time = run_time
  )

if(is.null(opt$SD)){
  identified <- suppressMessages(TMBhelper::Check_Identifiable(obj))

  # Make into list
  identified_param_list <- obj$env$parList(as.numeric(identified$BadParams$Param_check))
  identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==0,"Not estimated",x), how = "replace")
  identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==1,"OK",x), how = "replace")
  identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==2,"BAD",x), how = "replace")

  identified$param_list <- identified_param_list

  mod_objects$identified <- identified
}

class(mod_objects) <- "Rceattle"

if(!is.null(file)){
  save(mod_objects, file = paste0(file, ".RData"))
}

# suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
