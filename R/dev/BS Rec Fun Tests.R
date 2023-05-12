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


ss_run_M_ricker1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 1, # Estimate hindcast only
                                   recFun = build_srr(srr_fun = 3,
                                                      proj_mean_rec = FALSE,
                                                      srr_use_prior = TRUE,
                                                      srr_prior_mean = 4,
                                                      srr_prior_sd = 1),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   M1Fun = build_M1(M1_model = 1),
                                   phase = "default",
                                   verbose = 1,
                                   initMode = 0)

ss_run_M_ricker2 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 1, # Estimate hindcast only
                                   recFun = build_srr(srr_fun = 3,
                                                      proj_mean_rec = FALSE,
                                                      srr_use_prior = TRUE,
                                                      srr_prior_mean = 4,
                                                      srr_prior_sd = 1),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   M1Fun = build_M1(M1_model = 1),
                                   phase = "default",
                                   verbose = 1,
                                   initMode = 1)

ss_run_M_ricker3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 0, # Estimate hindcast only
                                   recFun = build_srr(srr_fun = 3,
                                                      proj_mean_rec = FALSE,
                                                      srr_use_prior = TRUE,
                                                      srr_prior_mean = 4,
                                                      srr_prior_sd = 1),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   M1Fun = build_M1(M1_model = 1),
                                   phase = "default",
                                   verbose = 1,
                                   initMode = 2)




ss_run_M_ricker1g <- Rceattle::fit_mod(data_list = GOA2018SS,
                                    inits = NULL, # Initial parameters = 0
                                    file = NULL, # Don't save
                                    estimateMode = 1, # Estimate hindcast only
                                    recFun = build_srr(srr_fun = 3,
                                                       proj_mean_rec = FALSE,
                                                       srr_use_prior = TRUE,
                                                       srr_prior_mean = 4,
                                                       srr_prior_sd = 1),
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 0, # Single species mode
                                    M1Fun = build_M1(M1_model = 1),
                                    phase = "default",
                                    verbose = 1,
                                    initMode = 0)

ss_run_M_ricker2g <- Rceattle::fit_mod(data_list = GOA2018SS,
                                    inits = NULL, # Initial parameters = 0
                                    file = NULL, # Don't save
                                    estimateMode = 1, # Estimate hindcast only
                                    recFun = build_srr(srr_fun = 3,
                                                       proj_mean_rec = FALSE,
                                                       srr_use_prior = TRUE,
                                                       srr_prior_mean = 4,
                                                       srr_prior_sd = 1),
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 0, # Single species mode
                                    M1Fun = build_M1(M1_model = 1),
                                    phase = "default",
                                    verbose = 1,
                                    initMode = 1)

ss_run_M_ricker3g <- Rceattle::fit_mod(data_list = GOA2018SS,
                                    inits = NULL, # Initial parameters = 0
                                    file = NULL, # Don't save
                                    estimateMode = 0, # Estimate hindcast only
                                    recFun = build_srr(srr_fun = 3,
                                                       proj_mean_rec = FALSE,
                                                       srr_use_prior = TRUE,
                                                       srr_prior_mean = 4,
                                                       srr_prior_sd = 1),
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 0, # Single species mode
                                    M1Fun = build_M1(M1_model = 1),
                                    phase = "default",
                                    verbose = 1,
                                    initMode = 2)



