library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")


################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
mydata$fleet_control$proj_F_prop <-rep(1,7)
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)





### For debugging
data_list = NULL;
inits = NULL;
map = NULL;
bounds = NULL;
file = NULL;
estimateMode = 0;
random_rec = FALSE;
random_q = FALSE;
random_sel = FALSE;
HCR = build_hcr();
niter = 3;
msmMode = 0;
avgnMode = 0;
minNByage = 0;
suitMode = 0;
suityr = NULL;
phase = NULL;
getsd = TRUE;
use_gradient = TRUE;
rel_tol = 1;
control = list(eval.max = 1e+09,
               iter.max = 1e+09, trace = 0);
getJointPrecision = TRUE;
loopnum = 5;
verbose = 1;
newtonsteps = 0



#### Dynamic
data_list = mydata
inits = ss_run$estimated_params # Initial parameters from ss_run
estimateMode = 2 # Run projection only
HCR = build_hcr(HCR = 5, # Tier3 HCR
                DynamicHCR = TRUE, # Use dynamic reference points
                FsprTarget = 0.4, # F40%
                FsprLimit = 0.35, # F35%
                Plimit = 0.2, # No fishing when SB<SB20
                Alpha = 0.05)
msmMode = 0 # Single species mode
verbose = 1
# ### For debugging
# data_list = NULL;
# inits = NULL;
# map = NULL;
# bounds = NULL;
# file = NULL;
# estimateMode = 0;
# random_rec = FALSE;
# random_q = FALSE;
# random_sel = FALSE;
# HCR = build_hcr();
# niter = 3;
# msmMode = 0;
# avgnMode = 0;
# minNByage = 0;
# suitMode = 0;
# suityr = NULL;
# phase = NULL;
# getsd = TRUE;
# use_gradient = TRUE;
# rel_tol = 1;
# control = list(eval.max = 1e+09,
#                iter.max = 1e+09, trace = 0);
# getJointPrecision = TRUE;
# loopnum = 5;
# verbose = 1;
# newtonsteps = 0


start_time <- Sys.time()

setwd(getwd())
'%!in%' <- function(x,y)!('%in%'(x,y))

#--------------------------------------------------
# 1. DATA and MODEL PREP
#--------------------------------------------------

# STEP 1 - LOAD DATA
if (is.null(data_list)) {
  stop("Missing data_list object")
}

# - Remove years of data previous to start year
data_list$UobsWtAge <- as.data.frame(data_list$UobsWtAge)
data_list$UobsAge <- as.data.frame(data_list$UobsAge)
data_list$wt <- data_list$wt[which(data_list$wt$Year == 0 | data_list$wt$Year >= data_list$styr),]
data_list$UobsAge <- data_list$UobsAge[which(data_list$UobsAge$Year == 0 | data_list$UobsAge$Year >= data_list$styr),]
data_list$UobsWtAge <- data_list$UobsWtAge[which(data_list$UobsWtAge$Year == 0 | data_list$UobsWtAge$Year >= data_list$styr),]
data_list$srv_biom <- data_list$srv_biom[which(abs(data_list$srv_biom$Year) >= data_list$styr),]
data_list$fsh_biom <- data_list$fsh_biom[which(abs(data_list$fsh_biom$Year) >= data_list$styr),]
data_list$comp_data <- data_list$comp_data[which(abs(data_list$comp_data$Year) >= data_list$styr),]
data_list$emp_sel <- data_list$emp_sel[which(data_list$emp_sel$Year == 0 | data_list$emp_sel$Year >= data_list$styr),]
data_list$NByageFixed <- data_list$NByageFixed[which(data_list$NByageFixed$Year == 0 | data_list$NByageFixed$Year >= data_list$styr),]
data_list$Pyrs <- data_list$Pyrs[which(data_list$Pyrs$Year == 0 | data_list$Pyrs$Year >= data_list$styr),]


# - Remove years of data after to proj year
data_list$wt <- data_list$wt[which(data_list$wt$Year <= data_list$projyr),]
data_list$UobsAge <- data_list$UobsAge[which(data_list$UobsAge$Year <= data_list$projyr),]
data_list$UobsWtAge <- data_list$UobsWtAge[which(data_list$UobsWtAge$Year <= data_list$projyr),]
data_list$srv_biom <- data_list$srv_biom[which(abs(data_list$srv_biom$Year) <= data_list$projyr),]
data_list$fsh_biom <- data_list$fsh_biom[which(abs(data_list$fsh_biom$Year) <= data_list$projyr),]
data_list$comp_data <- data_list$comp_data[which(abs(data_list$comp_data$Year) <= data_list$projyr),]
data_list$emp_sel <- data_list$emp_sel[which(data_list$emp_sel$Year <= data_list$projyr),]
data_list$NByageFixed <- data_list$NByageFixed[which(data_list$NByageFixed$Year <= data_list$projyr),]
data_list$Pyrs <- data_list$Pyrs[which(data_list$Pyrs$Year <= data_list$projyr),]


# - Extend catch data to proj year for projections
if(data_list$projyr > data_list$endyr){
  # yrs_proj <- (data_list$endyr + 1):data_list$projyr
  # proj_fsh_biom <- data_list$fsh_biom %>%
  #   group_by(Fleet_code) %>%
  #   slice(rep(n(),  length(yrs_proj))) %>%
  #   mutate(Year = yrs_proj, Catch = NA)
  # data_list$fsh_biom <- rbind(data_list$fsh_biom, proj_fsh_biom)

  for(flt in (unique(data_list$fsh_biom$Fleet_code))){
    fsh_biom_sub <- data_list$fsh_biom[which(data_list$fsh_biom$Fleet_code == flt),]
    yrs_proj <- (data_list$endyr + 1):data_list$projyr
    yrs_proj <- yrs_proj[which(yrs_proj %!in% fsh_biom_sub$Year)]
    nyrs_proj <- length(yrs_proj)
    proj_fsh_biom <- data.frame(Fleet_name = rep(fsh_biom_sub$Fleet_name[1], nyrs_proj),
                                Fleet_code = rep(flt, nyrs_proj),
                                Species = rep(fsh_biom_sub$Species[1], nyrs_proj),
                                Year = yrs_proj,
                                Month = rep(fsh_biom_sub$Month[length(fsh_biom_sub$Month)], nyrs_proj),
                                Selectivity_block = rep(fsh_biom_sub$Selectivity_block[length(fsh_biom_sub$Selectivity_block)], nyrs_proj),
                                Catch = rep(NA, nyrs_proj),
                                Log_sd = rep(fsh_biom_sub$Log_sd[length(fsh_biom_sub$Log_sd)], nyrs_proj))
    data_list$fsh_biom <- rbind(data_list$fsh_biom, proj_fsh_biom)
  }
}
data_list$fsh_biom <- data_list$fsh_biom[
  with(data_list$fsh_biom, order(Fleet_code, Year)),]

# Switches
data_list$random_rec <- as.numeric(random_rec)
data_list$debug <- estimateMode
data_list$niter <- niter
data_list$avgnMode <- avgnMode
data_list$msmMode <- msmMode
data_list$suitMode <- as.numeric(suitMode)
data_list$minNByage <- as.numeric(minNByage)
if(is.null(data_list$suityr)){
  data_list$suityr <- data_list$endyr
}

# HCR Switches
data_list$HCR = HCR$HCR
data_list$DynamicHCR = HCR$DynamicHCR
data_list$FsprTarget = HCR$FsprTarget
data_list$FsprLimit = HCR$FsprLimit
data_list$Ptarget = HCR$Ptarget
data_list$Plimit = HCR$Plimit
data_list$Alpha = HCR$Alpha
data_list$QnormHCR = ifelse(HCR$HCR == 6, qnorm(HCR$Pstar, 0, HCR$Sigma), 0) # Pstar HCR

# STEP 1 - LOAD PARAMETERS
if (is.character(inits) | is.null(inits)) {
  start_par <- suppressWarnings(Rceattle::build_params(
    data_list = data_list,
    inits = inits
  ))
} else{
  # inits$proj_F <- data_list$fleet_control$proj_F
  start_par <- inits
}
if(verbose > 0) {message("Step 1: Parameter build complete")}



# STEP 2 - BUILD MAP
if (is.null(map)) {
  map <-
    suppressWarnings(build_map(data_list, start_par, debug = estimateMode > 3, random_rec = random_rec))
} else{
  map <- map
}
if(verbose > 0) {message("Step 2: Map build complete")}


# STEP 3 - Get bounds
if (is.null(bounds)) {
  bounds <- Rceattle::build_bounds(param_list = start_par, data_list)
} else {
  bounds = bounds
}
if(verbose > 0) {message("Step 3: Param bounds complete")}


# STEP 4 - Setup random effects
random_vars <- c()
if (random_rec) {
  random_vars <- c(random_vars , "rec_dev", "init_dev")
}
if(random_q){
  random_vars <- c(random_vars , "ln_srv_q_dev")
}
if(random_sel){
  random_vars <- c(random_vars , "ln_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
}



# Set default phasing
if(!is.null(phase)){
  if(class(phase) == "character"){
    if(tolower(phase) == "default"){
      phase = list(
        dummy = 1,
        ln_pop_scalar = 4,
        ln_mean_rec = 1,
        ln_rec_sigma = 2,
        rec_dev = 2,
        init_dev = 2,
        ln_sex_ratio_sigma = 3,
        ln_M1 = 4,
        ln_mean_F = 1,
        ln_Flimit = 3,
        ln_Ftarget = 3,
        proj_F_prop = 1,
        F_dev = 1,
        ln_srv_q = 3,
        # srv_q_pow = 4,
        ln_srv_q_dev = 5,
        ln_sigma_srv_q = 4,
        ln_sigma_time_varying_srv_q = 4,
        sel_coff = 3,
        sel_coff_dev = 4,
        ln_sel_slp = 3,
        sel_inf = 3,
        ln_sel_slp_dev = 5,
        sel_inf_dev = 5,
        ln_sigma_sel = 4,
        sel_curve_pen = 4,
        ln_sigma_srv_index = 2,
        ln_sigma_fsh_catch = 2,
        comp_weights = 4,
        logH_1 = 6,
        logH_1a = 6,
        logH_1b = 6,
        logH_2 = 6,
        logH_3 = 6,
        H_4 = 6,
        log_gam_a = 5,
        log_gam_b = 5,
        log_phi = 5
      )
    }
  }

  if(class(phase) == "character"){
    if(tolower(phase) != "default"){
      warning("phase misspecified: please set to 'default' or list with the same order as parameters.")
    }
  }
}


# STEP 5 - Compile CEATTLE is providing cpp file
# - Get cpp file if not provided
TMBfilename <- "ceattle_v01_09"




# STEP 6 - Reorganize data and build model object
Rceattle:::data_check(data_list)
data_list_reorganized <- Rceattle::rearrange_dat(data_list)
data_list_reorganized = c(list(model = "ceattle_v01_09"),data_list_reorganized)
data_list_reorganized$HCR = 0 # Estimate model with F = 0 for the projection

# - Update comp weights and F_prop from data
if(!is.null(data_list$fleet_control$Comp_weights)){
  start_par$comp_weights = data_list$fleet_control$Comp_weights
}
start_par$proj_F_prop = data_list$fleet_control$proj_F_prop

if(verbose > 0) {message("Step 4: Data rearranged complete")}

# STEP 7 - Set up parameter bounds
L <- c()
U <- c()
for(i in 1:length(map$mapFactor)){
  if(names(map$mapFactor)[i] %!in% random_vars){ # Dont have bounds for random effects
    L = c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
    U = c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
  }
}


# STEP 8 - Fit hindcast
step = 5
# If phase: phase hindcast
if(!is.null(phase) & estimateMode %in% c(0,1) ){
  message(paste0("Step ", step,": Phasing begin"))
  phase_pars <- Rceattle::TMBphase(
    data = data_list_reorganized,
    parameters = start_par,
    map = map$mapFactor,
    random = random_vars,
    phases = phase,
    model_name = TMBfilename,
    silent = verbose != 2,
    use_gradient = use_gradient,
    control = control
  )

  start_par <- phase_pars

  if(verbose > 0) {message(paste0("Step ", step,": Phasing complete - getting final estimates"))}
  step = step + 1
}


# STEP 9 - Fit final hindcast model
if(estimateMode != 2){ # dont build if projection
  obj = TMB::MakeADFun(
    data_list_reorganized,
    parameters = start_par,
    DLL = TMBfilename,
    map = map$mapFactor,
    random = random_vars,
    silent = verbose != 2
  )
}

# -- Save objects
mod_objects <-
  list(
    TMBfilename = TMBfilename,
    initial_params = start_par,
    bounds = bounds,
    map = map
  )

if(verbose > 0) {message(paste0("Step ",step, ": final build complete. Optimizing."))}
step = step + 1


# -- Optimize hindcast
if(estimateMode %in% c(0,1,4)){
  opt = Rceattle::fit_tmb(obj = obj,
                          fn=obj$fn,
                          gr=obj$gr,
                          startpar=obj$par,
                          lower = L,
                          upper = U,
                          loopnum = loopnum,
                          getsd = getsd,
                          control = control,
                          getJointPrecision = getJointPrecision,
                          quiet = verbose == 2,
  )
  if(verbose > 0) {message("Step ",step, ": Final optimization complete")
    step = step + 1
  }

  # -- Convergence warnings
  if(estimateMode %in% c(0,1)){
    # Bad parameter identification
    if(is.null(opt$SD) & getsd){
      identified <- suppressMessages(TMBhelper::check_estimability(obj))

      # Make into list
      identified_param_list <- obj$env$parList(identified$BadParams$Param_check)
      identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==0,"Not estimated",x), how = "replace")
      identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==1,"OK",x), how = "replace")
      identified_param_list <- rapply(identified_param_list,function(x) ifelse(x==2,"BAD",x), how = "replace")
      identified$param_list <- identified_param_list
      mod_objects$identified <- identified
    }
  }
}

# -- Get MLEs
if (estimateMode > 1) { # Debugging and projection only: use initial parameters
  last_par <- start_par
} else{
  if(!random_rec){
    last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
  } else {
    last_par = try(obj$env$parList())
  }
}


# STEP 10 - Run HCR projections
if(estimateMode %in% c(0,2,4)){
  if(data_list$HCR > 2){
    data_list_reorganized$HCR = data_list$HCR # Set HCR back to original


    # START HERE
    # -- Update map in obs
    hcr_map <- build_hcr_map(data_list, map, debug = estimateMode > 3)

    obj = TMB::MakeADFun(
      data_list_reorganized,
      parameters = last_par,
      DLL = TMBfilename,
      map = hcr_map$mapFactor,
      random = random_vars,
      silent = verbose != 2
    )

    # -- Optimize
    opt = Rceattle::fit_tmb(obj = obj,
                            fn=obj$fn,
                            gr=obj$gr,
                            startpar=obj$par,
                            loopnum = loopnum,
                            getsd = getsd,
                            control = control,
                            getJointPrecision = FALSE,
                            quiet = verbose == 2,
    )

    # obj$report()$DynamicSPRtarget/obj$report()$DynamicSPR0
    # obj$report()$DynamicSPRlimit/obj$report()$DynamicSPR0
    #
    # obj$report()$SPRtarget/obj$report()$SPR0
    # obj$report()$SPRlimit/obj$report()$SPR0
    #
    # obj$report()$SPR0
    # obj$report()$SB0

    if(verbose > 0) {message("Step ",step, ": Projections complete")}

    # -- Update MLEs
    if (estimateMode > 2) { # Debugging, give initial parameters
      last_par <- start_par
    }
    else{
      if(!random_rec){
        last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
      } else {
        last_par = try(obj$env$parList())
      }
    }
  }
}


obj$report()$DynamicSPRtarget/obj$report()$DynamicSPR0
obj$report()$DynamicSPRlimit/obj$report()$DynamicSPR0

obj$report()$SPRtarget/obj$report()$SPR0
obj$report()$SPRlimit/obj$report()$SPR0

obj$report()$SPR0
obj$report()$SB0


obj$report(opt$opt$par)$dynamicsprlimitcheck

data_list_reorganized$DynamicHCR
round(obj$report(opt$par)$jnll_comp[12:15,1:3])
round(obj$report(opt$par-1)$jnll_comp[12:15,1:3])
round(obj$report(opt$par+1)$jnll_comp[12:15,1:3])
