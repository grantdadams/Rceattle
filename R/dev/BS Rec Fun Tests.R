library(Rceattle)
library(dplyr)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060

data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
BS2017MS$projyr <- 2060

BS2017SS$fleet_control$proj_F_prop <-rep(1,7)
BS2017MS$fleet_control$proj_F_prop <- rep(1, 7)


ss_run_ricker1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 1, # Estimate hindcast only
                                   recFun = build_srr(srr_fun = 3,
                                                      proj_mean_rec = FALSE,
                                                      srr_use_prior = TRUE,
                                                      srr_prior_mean = 4,
                                                      srr_prior_sd = 4),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   phase = "default",
                                   verbose = 1,
                                   initMode = 0)

ss_run_ricker2 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 1, # Estimate hindcast only
                                   recFun = build_srr(srr_fun = 3,
                                                      proj_mean_rec = FALSE,
                                                      srr_use_prior = TRUE,
                                                      srr_prior_mean = 4,
                                                      srr_prior_sd = 2),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   phase = "default",
                                   verbose = 1,
                                   initMode = 1)

ss_run_ricker3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 1, # Estimate hindcast only
                                   recFun = build_srr(srr_fun = 3,
                                                      proj_mean_rec = FALSE,
                                                      srr_use_prior = TRUE,
                                                      srr_prior_mean = 4,
                                                      srr_prior_sd = 2),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   phase = "default",
                                   verbose = 1,
                                   initMode = 2)
