#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Setup ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Rceattle)
library(dplyr)
library(ggplot2)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Estimation ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            file = NULL,
                            phase = TRUE, # Phase the model
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1)


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
# Or we can use the previous data set
BS2017MS$projyr <- 2060

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species MLEs
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 5, # 5 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1,  # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Simulation ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Survey and composition data can be simulated from the estimated quantities using `sim_mod`:
# - simulates all historical data, assuming historical variance and distributions
ss_sim_data <- sim_mod(ss_run, simulate = TRUE)

ss_sim_run <- Rceattle::fit_mod(
  data_list = ss_sim_data,
  inits = NULL, # Initial parameters = 0
  phase = TRUE,
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  verbose = 1)

# Simulate data for the multispecies model
ms_sim_data <- sim_mod(ms_run, simulate = TRUE)

ms_sim_run <- Rceattle::fit_mod(
  data_list = ms_sim_data,
  inits = ms_run$estimated_params, # Initial parameters at previous MLEs (quicker estimation)
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 1,
  verbose = 1)

# We can plot the simulated runs
mod_list <- list(ss_run, ss_sim_run, ms_run, ms_sim_run)
mod_names <- c("SS", "SS sim", "MS", "MS sim")

plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names)
