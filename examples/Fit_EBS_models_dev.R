# Grant Adams, Kirstin Holsman, Andre Punt - April 2019
# Code to run Bering Sea CEATTLE model in TMB
# Citation:
# Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., and Moffitt, E. A. 2015. A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134: 360–378.

library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data("BS2017SS") # Single-species data. ?BS2017SS for more information on the data
data("BS2017MS") # Multi-species data. Note: the only difference is the residual mortality (M1_base) is lower

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")


################################################
# Estimation
################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
mydata$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)
load("~/GitHub/ss.RData")
plot_biomass(list(ss_run, mod_objects))




# Single-species, but estimate M
load("~/GitHub/ssm.RData")
inits <- mod_objects$estimated_params
inits$rec_pars <- matrix(9, nrow = mod_objects$data_list$nspp, ncol = 2)
inits$rec_pars[,1] <- inits$ln_mean_rec
inits$ln_mean_rec <- NULL

ss_run_M <- Rceattle::fit_mod(data_list = mydata,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Estimate
                              M1Fun = build_M1(M1_model = 1,
                                               M1_use_prior = FALSE,
                                               M2_use_prior = FALSE),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = "default",
                              initMode = 1,
                              verbose = 1)

plot_biomass(list(ss_run_M, mod_objects))
plot_recruitment(list(ss_run_M, mod_objects))
plot_ssb(list(ss_run_M, mod_objects))
round(ss_run_M$quantities$jnll_comp,3)[1:3,] - round(mod_objects$quantities$jnll_comp,3)[1:3,]
round(ss_run_M$quantities$NByage[1,1,1:12,1:5]-mod_objects$quantities$NByage[1,1,1:12,1:5],3)
round(ss_run_M$quantities$NByage[1,1,1:12,1]-mod_objects$quantities$NByage[1,1,1:12,1],3)

ss_run_M$quantities$NByage[1,,12,1]
ss_run_M$quantities$R0[1] * exp(- sum(ss_run_M$quantities$M1[1,1,1:11]) + ss_run_M$estimated_params$init_dev[1,11])/(1-exp(-ss_run_M$quantities$M1[1,1,12]))

# - Multi-species
# For the a multispecies model we from the single species parameters.
BS2017MS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            M1Fun = build_M1(M1_model = 0,
                                             updateM1 = TRUE,
                                             M1_use_prior = FALSE,
                                             M2_use_prior = FALSE),
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

load("~/GitHub/ms.RData")
plot_biomass(list(ms_run, mod_objects))


################################################
# Plotting
################################################
# We can plot all runs
mod_list <- list(ss_run, ss_run_M, ms_run)
mod_names <- c("Single-species", "Single-species estimate M", "Multi-species")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

# Plot mortality and predation
plot_b_eaten(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey by each predator
plot_mort(Rceattle = ms_run, type = 3) # Mortality-at-age time series

# Run diagnostics
plot_selectivity(Rceattle = ms_run)
plot_comp(ms_run) # Fitted survey composition data
plot_index(ms_run) # Fitted indices of abundance
plot_catch(ms_run) # Fitted catch series
