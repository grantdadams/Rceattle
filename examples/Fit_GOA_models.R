# Code to run Gulf of Alaska CEATTLE model in TMB
# Citation:
# Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn, M. W., Ianelli, J. N., Spies, I., ... & Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fisheries Research, 251, 106303.

library(Rceattle)
library(dplyr)

################################################
# Data
################################################
# Example
# To run the 2018 single species assessment for the Gulf of Alaska, a data file must first be loaded:
data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data


################################################
# Estimation
################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
# GOA2018SS$fleet_control$proj_F_prop <- rep(1, nrow(GOA2018SS$fleet_control))
ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


# Single-species, but estimate M
ss_run_M <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              M1Fun = build_M1(M1_model = c(1,2,1)), # Estimate M
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

plot_biomass(ss_run_M)

# - Multi-species
# For the a multispecies model we from the single species parameters.
ms_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                            inits = ss_run_M$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            M1Fun = build_M1(M1_model = c(1,2,1)),
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

plot_biomass(ms_run)


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
plot_mortality(Rceattle = ms_run, type = 3) # Mortality-at-age time series

# Run diagnostics
plot_selectivity(Rceattle = ms_run)
plot_comp(ms_run) # Fitted survey composition data
plot_index(ms_run) # Fitted indices of abundance
plot_catch(ms_run) # Fitted catch series
