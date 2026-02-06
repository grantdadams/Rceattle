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
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Time-series plots ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# We can plot both runs as well:
mod_list <- list(ss_run, ms_run)
mod_names <- c("Single-species", "Multi-species")

# Plot trajectories
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_b_eaten(Rceattle = mod_list, model_names = mod_names)
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_m_at_age(Rceattle = mod_list, model_names = mod_names, age = 2)
plot_m2_at_age_prop(Rceattle = mod_list, model_names = mod_names) # Predation M by each species
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names)
plot_depletion(Rceattle = mod_list, model_names = mod_names)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Diagnostics plots ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# - Most functions only take a single model
plot_comp(ms_run) # Fitted composition data
plot_index(ms_run) # Fitted indices of abundance
plot_logindex(ms_run) # Fitted log indices of abundance
plot_indexresidual(ms_run) # Residuals of log indices of abundance
plot_catch(ms_run) # Fitted catch series


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Retrospective ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# we can also do retrospective peels and calculate Mohn's Rho on the CEATTLE
# NOTE: this is using mean historical F for projections as we changed it above
ss_run_retro <- Rceattle::retrospective(Rceattle = ss_run, peels = 5)

# Look at Mohns rho
ss_run_retro$mohns

# Plot retrospectives
plot_biomass(ss_run_retro$Rceattle_list)

# See how forecast changes
plot_biomass(ss_run_retro$Rceattle_list, incl_proj = TRUE)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Jitter ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
jitters <- Rceattle::jitter(Rceattle = ss_run, njitter = 10, phase = TRUE)
hist(log(jitters$nll - min(jitters$nll)))


# Plot the jitters
plot_biomass(jitters$Rceattle_list)
length(jitters$Rceattle_list) # Some models did not converge, and are not returned
