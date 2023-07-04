library(Rceattle)
library(dplyr)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060

data("GOA2018SS") # Note: the only difference is the residual mortality (M1_base) is lower
GOA2018SS$projyr <- 2060

BS2017SS$fleet_control$proj_F_prop <-rep(1,7)
GOA2018SS$fleet_control$proj_F_prop <- rep(1, length(GOA2018SS$fleet_control$proj_F_prop))


DATA_SET <-BS2017SS
DATA_SET <-GOA2018SS

# For EBS MS:
alpha = exp(c(4.121, 2.119, 1.553))
# For GOA MS:
alpha = exp(c(3.143, 1.975, 1.44))

#--------------------------------------------------------------------
# SINGLE-SPECIES FIX M
#--------------------------------------------------------------------
ss_run <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  phase = "default",
  initMode = 2,
  verbose = 1)

ss_run_ricker1 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  M1Fun = build_M1(M1_model = 0),
  phase = "default",
  verbose = 1,
  initMode = 0)

ss_run_ricker2 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  M1Fun = build_M1(M1_model = 0),
  phase = "default",
  verbose = 1,
  initMode = 1)

ss_run_ricker3 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  M1Fun = build_M1(M1_model = 0),
  phase = "default",
  verbose = 1,
  initMode = 2)



#--------------------------------------------------------------------
# SINGLE-SPECIES ESTIMATE M
#--------------------------------------------------------------------
ss_run_M <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL,#ss_run$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  M1Fun = build_M1(M1_model = c(1,2,1)), # Estimate M
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  phase = "default",
  initMode = 2,
  verbose = 1)


ss_run_M_ricker1 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 0)

ss_run_M_ricker2 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 1)

ss_run_M_ricker3 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 2)



#--------------------------------------------------------------------
# MULTI-SPECIES
#--------------------------------------------------------------------
ms_run <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = ss_run_M$estimated_params, # Initial parameters from single species ests
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  M1Fun = build_M1(M1_model = c(1,2,1)),
  niter = 3, # 3 iterations around population and predation dynamics
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # MSVPA based
  suitMode = 0, # empirical suitability
  initMode = 2,
  phase = "default",
  verbose = 1)


# Ricker with prior and initial age structure as free parameters
ms_run_ricker1 <- Rceattle::fit_mod(
  data_list = DATA_SET,#;
  inits = ss_run_M$estimated_params,#; # Initial parameters = 0
  file = NULL,#; # Don't save
  estimateMode = 1,#; # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),#;
  random_rec = FALSE,#; # No random recruitment
  msmMode = 1,#; # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),#;
  phase = default,#;
  verbose = 1,#;
  initMode = 0
)


# Ricker with prior and initial age structure unfished equilibrium
ms_run_ricker2 <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = ms_run$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 1)


# Ricker with prior and initial age structure fished equilibrium
ms_run_ricker3 <- Rceattle::fit_mod(
  data_list = DATA_SET,#;
  inits = ms_run$estimated_params,#; # Initial parameters = 0
  file = NULL,#; # Don't save
  estimateMode = 1,#; # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),#;
  random_rec = FALSE,#; # No random recruitment
  msmMode = 1,#; # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),#;
  phase = "default",#;
  verbose = 1,#;
  initMode = 2)


# Ricker with prior and initial age structure as fished equilibrium
ms_run_ricker3_prior <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = ss_run_M$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 2,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 2)


# Ricker with fixed alpha and initial age structure as fished equilibrium
ms_run_ricker3_fixed <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = ss_run_M$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 0,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 2)


# Ricker via Ianelli penalty and initial age structure as fished equilibrium
ms_run_ricker3_Ianelli <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = ms_run$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 0,
                     srr_pred_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1, # Freely estimate (prior doesnt gives the bomb)
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 2)


# Ricker via Ianelli penalty and initial age structure as fished equilibrium (prior on alpha)
ms_run_ricker3_Ianelli_prior <- Rceattle::fit_mod(
  data_list = DATA_SET,
  inits = ms_run$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  recFun = build_srr(srr_fun = 0,
                     srr_pred_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 2, # Freely estimate (prior doesnt gives the bomb)
                     srr_prior_mean = alpha,
                     srr_prior_sd = 1),
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Single species mode
  M1Fun = build_M1(M1_model = c(1,2,1)),
  phase = "default",
  verbose = 1,
  initMode = 2)


# Plot it
pars <- exp(ms_run_ricker3_Ianelli_prior$estimated_params$rec_pars)
ssb <- (ms_run_ricker3_Ianelli_prior$quantities$biomassSSB)
r <- (ms_run_ricker3_Ianelli_prior$quantities$R)

par(mfrow = c(3,1), mar = c(3,2,0.5,0.5))
for(i in 1:3){
  sp <- i

  # no sr
  ssb <- ms_run$quantities$biomassSSB[,1:41]
  r <- ms_run$quantities$R[,2:42]
  mod <- lm(log(r[sp,]/ssb[sp,]) ~ c(ssb[sp,]))


  plot(x = ssb[sp,], y = r[sp,], xlab = "SSB", ylab = "R", pch = 16, col = 2, xlim = c(0, max(ssb[sp,])), ylim = c(0, max(r[sp,])))
  curve(exp(mod$coefficients[1]) * x * exp(mod$coefficients[2] * x), from = 0, add=TRUE, to = max(ssb[sp,]), col = 2, lwd = 2)

# sr penalty
  ssb <- ms_run_ricker3_Ianelli_prior$quantities$biomassSSB[,1:41]
  r <- ms_run_ricker3_Ianelli_prior$quantities$R[,2:42]
  points(x = ssb[sp,], y = r[sp,], pch = 16, col = 1)
  curve(pars[sp,2] * x * exp(-pars[sp,3]/1000000 * x), from = 0, to = max(ssb[sp,]), add = TRUE, lwd = 2)
  legend("topleft", legend = c("GOA Pollock", "GOA ATF", "GOA Cod")[sp], bty = "n")
}
