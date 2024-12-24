library(Rceattle)
data("BS2017SS")
BS2017SS$fleet_control$Comp_loglike <- -1


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Debugging section ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
data_list = BS2017SS;
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
initMode = 2
minNByage = 0;
suitMode = 0;
suit_styr = NULL;
suit_endyr = NULL;
phase = FALSE;
getsd = TRUE;
use_gradient = TRUE;
rel_tol = 1;
control = list(eval.max = 1e+09,
               iter.max = 1e+09, trace = 0);
getJointPrecision = TRUE;
loopnum = 5;
verbose = 1;
newtonsteps = 0
recFun = build_srr()
M1Fun = build_M1()
projection_uncertainty = TRUE
catch_hcr = FALSE

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 0 - Start ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
start_time <- Sys.time()

compiler::enableJIT(0)
extend_length <- function(x){
  if(length(x) == data_list$nspp){ return(x)}
  else {return(rep(x, data_list$nspp))}
}

setwd(getwd())


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 1 - Load data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if (is.null(data_list)) {
  stop("Missing data_list object")
}

data_list <- Rceattle::clean_data(data_list)

# Add switches from function call
data_list$random_rec <- as.numeric(random_rec)
data_list$estimateMode <- estimateMode
data_list$niter <- niter
data_list$avgnMode <- avgnMode
data_list$initMode <- initMode
data_list$loopnum <- loopnum
data_list$msmMode <- msmMode
data_list$suitMode <- as.numeric(suitMode)

# - Suitability
# -- Start year
if(is.null(suit_styr) & is.null(data_list$suit_styr)){ # If not provided in data or function, use start year
  data_list$suit_styr <- data_list$styr
}
if(!is.null(suit_styr)){ # If provided in function, override data
  data_list$suit_styr <- suit_styr
}

# -- End year
if(is.null(suit_endyr) & is.null(data_list$suit_endyr)){ # If not provided in data or function, use end year
  data_list$suit_endyr <- data_list$endyr
}
if(!is.null(suit_endyr)){ # If provided in function, override data
  data_list$suit_endyr <- suit_endyr
}


# * Recruitment switches ----
data_list$srr_fun <- recFun$srr_fun
data_list$srr_pred_fun <- recFun$srr_pred_fun
data_list$proj_mean_rec <- recFun$proj_mean_rec
if(is.null(recFun$srr_meanyr) & is.null(data_list$srr_meanyr)){ # If no meanyear is provided in data or function, use end year
  data_list$srr_meanyr <- data_list$endyr
}
if(!is.null(recFun$srr_meanyr)){ # If mean year is provided in function, override data
  data_list$srr_meanyr <- recFun$srr_meanyr
}

# -- Start year
if(is.null(recFun$srr_hat_styr) & is.null(data_list$srr_hat_styr)){ # If not provided in data or function, use start year
  data_list$srr_hat_styr <- data_list$styr
}
if(!is.null(recFun$srr_hat_styr)){ # If provided in function, override data
  data_list$srr_hat_styr <- recFun$srr_hat_styr
}

# -- End year
if(is.null(recFun$srr_hat_endyr) & is.null(data_list$srr_hat_endyr)){ # If not provided in data or function, use end year
  data_list$srr_hat_endyr <- data_list$endyr
}
if(!is.null(recFun$srr_hat_endyr)){ # If provided in function, override data
  data_list$srr_hat_endyr <- recFun$srr_hat_endyr
}

data_list$srr_est_mode <- recFun$srr_est_mode
data_list$srr_prior <- extend_length(recFun$srr_prior)
data_list$srr_prior_sd <- extend_length(recFun$srr_prior_sd)
data_list$srr_env_indices <- recFun$srr_env_indices
data_list$Bmsy_lim <- extend_length(recFun$Bmsy_lim)

# * M1 switches ----
if(!is.null(data_list$M1_model)){
  if(sum(data_list$M1_model != extend_length(M1Fun$M1_model))){
    warning("M1_model in data is different than in call `fit_mod`")
  }
}

# FIXME: may want to pull from data here too
data_list$M1_model= extend_length(M1Fun$M1_model)
updateM1 = M1Fun$updateM1
data_list$M1_use_prior = extend_length(M1Fun$M1_use_prior) * (data_list$M1_model > 0) # Sets to 0 if M1 is fixed
data_list$M2_use_prior = extend_length(M1Fun$M2_use_prior) * (msmMode > 0) # Sets to 0 if single-species
data_list$M_prior = extend_length(M1Fun$M_prior)
data_list$M_prior_sd = extend_length(M1Fun$M_prior_sd)


# - HCR Switches (make length of nspp if not)
data_list$HCR = HCR$HCR
data_list$DynamicHCR = HCR$DynamicHCR
if(HCR$HCR != 2){ # FsprTarget is also used for fixed F (so may be of length nflts)
  data_list$FsprTarget = extend_length(HCR$FsprTarget)
} else {
  data_list$FsprTarget = HCR$FsprTarget
}
data_list$FsprLimit = extend_length(HCR$FsprLimit)
data_list$Ptarget = extend_length(HCR$Ptarget)
data_list$Plimit = extend_length(HCR$Plimit)
data_list$Alpha = extend_length(HCR$Alpha)
data_list$Pstar = extend_length(HCR$Pstar)
data_list$Sigma = extend_length(HCR$Sigma)
data_list$Fmult = extend_length(HCR$Fmult)
data_list$HCRorder = extend_length(HCR$HCRorder)
data_list$QnormHCR = qnorm(data_list$Pstar, 0, data_list$Sigma)

if(data_list$HCR == 2 & estimateMode == 2){estimateMode = 4} # If projecting under constant F, run parmeters through obj only
if(data_list$msmMode > 0 & !data_list$HCR %in% c(0, 1, 2, 3, 6)){
  warning("WARNING:: Only HCRs 1, 2, 3, and 6 work in multi-species mode currently")
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 2: Load/build parameters ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if (is.character(inits) | is.null(inits)) {
  start_par <- suppressWarnings(Rceattle::build_params(data_list = data_list))
} else{
  start_par <- inits

  # - Adjust srr parameters
  if(ncol(start_par$beta_rec_pars) != length(data_list$srr_env_indices)){
    start_par$beta_rec_pars <- matrix(0, nrow = data_list$nspp, ncol = length(data_list$srr_env_indices))
  }
}
if(verbose > 0) {message("Step 1: Parameter build complete")}

# Set Fdev for years with 0 catch to very low number
catch_data_sub <- data_list$catch_data %>%
  filter(Year <= data_list$endyr)
fsh_ind <- catch_data_sub$Fleet_code[which(catch_data_sub$Catch == 0)]
yr_ind <- catch_data_sub$Year[which(catch_data_sub$Catch == 0)] - data_list$styr + 1
for(i in 1:length(fsh_ind)){
  start_par$F_dev[fsh_ind[i], yr_ind[i]] <- -999
}
rm(catch_data_sub)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 3: Load/build map ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if (is.null(map)) {
  map <- suppressWarnings(build_map(data_list, start_par, debug = estimateMode == 4, random_rec = random_rec, random_sel = random_sel))
} else{
  map <- map
}
if(verbose > 0) {message("Step 2: Map build complete")}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 4: Get bounds ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if (is.null(bounds)) {
  bounds <- Rceattle::build_bounds(param_list = start_par, data_list)
} else {
  bounds = bounds
}
if(verbose > 0) {message("Step 3: Param bounds complete")}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 5: Setup random effects ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# FIXME: this should be controled by fleet_control
random_vars <- c()
if (random_rec) {
  if(initMode > 0){
    random_vars <- c(random_vars , "rec_dev", "init_dev")
  } else{
    random_vars <- c(random_vars , "rec_dev")
  }
}
if(random_q){
  random_vars <- c(random_vars , "index_q_dev")
}
if(random_sel){
  random_vars <- c(random_vars , "ln_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 6: Reorganize data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
Rceattle:::data_check(data_list)

data_list_reorganized <- Rceattle::rearrange_dat(data_list)
data_list_reorganized$forecast <- rep(0, data_list_reorganized$nspp) # Don't include BRPs in likelihood of hindcast

# - Update comp weights, future F (if input) and F_prop from data
if(!is.null(data_list$fleet_control$Comp_weights)){
  start_par$comp_weights = data_list$fleet_control$Comp_weights
}
start_par$proj_F_prop = data_list$fleet_control$proj_F_prop

nyrs_proj <- data_list$projyr - data_list$styr + 1
if(!is.null(HCR$FsprTarget) & HCR$HCR == 2){
  start_par$ln_Ftarget = log(HCR$FsprTarget) # Fixed fishing mortality for projections for each species
}

# - Update M1 parameter object from data if initial parameter values input
if(updateM1){
  m1 <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T))) # Set up array

  # Initialize from inputs
  for (i in 1:nrow(data_list$M1_base)) {
    sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
    sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

    # Fill in M1 array from fixed values for each sex
    if(sex == 0){ sex = c(1, 2)} # If sex = combined/both males and females, fill in both dimensions
    for(j in 1:length(sex)){
      m1[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = T)) + 2])
    }
  }
  start_par$ln_M1 <- log(m1)
}

# - Update alpha for stock-recruit if fixed/prior and initial parameter values input
if(data_list$srr_est_mode %in% c(0,2) & data_list$srr_pred_fun > 3){
  start_par$rec_pars[,2] <- log(data_list$srr_prior)
}

if(verbose > 0) {message("Step 4: Data rearranged complete")}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 7: Set up parameter bounds ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
L <- c()
U <- c()
for(i in 1:length(map$mapFactor)){
  if(names(map$mapFactor)[i] %!in% random_vars){ # Dont have bounds for random effects
    L = c(L, unlist(bounds$lower[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
    U = c(U, unlist(bounds$upper[[i]])[which(!is.na(unlist(map$mapFactor[[i]])) & !duplicated(unlist(map$mapFactor[[i]])))])
  }
}

# Dimension check
dim_check <- sapply(start_par, unlist(length)) == sapply(map$mapFactor, unlist(length))
if(sum(dim_check) != length(dim_check)){
  stop(print(paste0("Map and parameter objects are not the same size for: ", names(dim_check)[which(dim_check == FALSE)])))
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STEP 9: Fit hindcast ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
rtmb_ceattle(start_par, data_list_reorganized)
cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
obj <- RTMB::MakeADFun(cmb(rtmb_ceattle, data_list_reorganized),
                       parameters = start_par,
                       map = map$mapFactor,
                       random = random_vars,
                       silent = verbose != 2
)
