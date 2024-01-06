################################################
# Setup
################################################
library(Rceattle)
library(dplyr)
library(ggplot2)


################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2030

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")


################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            file = NULL,
                            phase = "default", # Phase the model
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1)
# Type ?fit_mod for more details

# The you can plot the model results using using
plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run, add_ci = TRUE)

# Note: fitting the model outputs the composition weights using the harmonic mean MacCallister-Ianelli method
# - setting the comp weights to the Mcallister weights and refitting the model will change the MLE based on the new weighting
ss_run$data_list$fleet_control$Comp_weights <- ss_run$data_list$fleet_control$Est_weights_mcallister
ss_run_weighted <- Rceattle::fit_mod(data_list = ss_run$data_list,
                                     inits = ss_run$estimated_params, # Initial parameters from previous MLEs
                                     file = NULL, # Don't save
                                     estimateMode = 0, # Estimate
                                     random_rec = FALSE, # No random recruitment
                                     msmMode = 0, # Single species mode
                                     verbose = 1)

# The you can plot the model results using using
plot_biomass(Rceattle =  list(ss_run, ss_run_weighted))
plot_recruitment(Rceattle =  list(ss_run, ss_run_weighted))


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
# Or we can use the previous data set
BS2017MS$projyr <- 2030

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species MLEs
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 5, # 5 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

# We can plot both runs as well:
mod_list <- list(ss_run, ms_run)
mod_names <- c("Single-species", "Multi-species")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_selectivity(Rceattle = mod_list, model_names = mod_names)
plot_mort(Rceattle = mod_list, model_names = mod_names, age = 2)

# Run diagnostics
plot_fsh_comp(ms_run) # Fitted fishery composition data
plot_srv_comp(ms_run) # Fitted survey composition data
plot_index(ms_run) # Fitted indices of abundance
plot_catch(ms_run) # Fitted catch series


# Save results
write_results(Rceattle = ms_run, file = "yourresults.xlsx")


################################################
# Projection
################################################

# PROJECTION 1: CHANGING F
# Rceattle automatically projects the population forward when estimating under no fishing
BS2017MS$projyr # Year the population is projected forward

# Re-run, without estimating
ms_run_proj <- Rceattle::fit_mod(data_list = BS2017MS,
                                 inits = ms_run$estimated_params, # Initial parameters from single species ests
                                 file = NULL, # Don't save
                                 estimateMode = 0, # Run in projection mode
                                 HCR = build_hcr(HCR = 2,
                                                 FsprTarget = c(0.2342936, 0.513, 0.0774777)), # Set projection F mean historical F
                                 niter = 5, # 5 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 suitMode = 0, # empirical suitability
                                 verbose = 1)




# PROJECTION 1: CHANGING FUTURE RECRUITMENT
# Change recruitment deviations - NOTE: the workflow for this may change in the future

# Get years from projection
nyrs <- BS2017MS$endyr - BS2017MS$styr + 1
nyrs_proj <- BS2017MS$projyr - BS2017MS$styr + 1
yrs_proj <- (nyrs + 1):nyrs_proj

# Replace future rec_devs with numbers
ms_run$estimated_params$rec_dev[,yrs_proj] <- replace(
  ms_run$estimated_params$rec_dev[,yrs_proj],
  values = rnorm( length(ms_run$estimated_params$rec_dev[,yrs_proj]),
                  mean = 0,
                  sd = 0.707) # Assumed value from penalized likelihood
)

ms_run_proj2 <- Rceattle::fit_mod(data_list = BS2017MS,
                                  inits = ms_run$estimated_params, # Initial parameters from single species ests
                                  file = NULL, # Don't save
                                  estimateMode = 3, # Do not estimate. Not changing parameters right now
                                  niter = 5, # 5 iterations around population and predation dynamics
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 1, # MSVPA based
                                  suitMode = 0, # empirical suitability
                                  verbose = 1)


# plot
mod_list <- list(ss_run, ms_run, ms_run_proj, ms_run_proj2)
mod_names <- c("SS", "MS no F", "MS mean historical F", "MS mean historical F w random rec")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)



################################################
# Retrospective
################################################
# we can also do retrospective peels and calculate Mohn's Rho on the CEATTLE
# NOTE: this is using mean historical F for projections as we changed it above
ms_run_proj_retro <- retrospective(Rceattle = ms_run_proj, peels = 5)

# Look at Mohns rho
ms_run_proj_retro$mohns

# Plot retrospectives
plot_biomass(ms_run_proj_retro$Rceattle_list)

# See how forecast changes
plot_biomass(ms_run_proj_retro$Rceattle_list, incl_proj = TRUE)


################################################
# Simulation
################################################
# Survey and composition data can be simulated from the estimated quantities using `sim_mod`:
ss_sim <- sim_mod(ss_run, simulate = TRUE)

ss_sim_run <- Rceattle::fit_mod(
  data_list = ss_sim,
  inits = NULL, # Initial parameters = 0
  phase = "default",
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  verbose = 1)

# Simulate data for the multispecies model
ms_sim <- sim_mod(ms_run, simulate = TRUE)

ms_sim_run <- Rceattle::fit_mod(
  data_list = ms_sim,
  inits = ms_run$estimated_params, # Initial parameters at previous MLEs
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Holsman MS mode
  verbose = 1)

# We can plot the simulated runs
mod_list <- list(ss_run, ss_sim_run, ms_run, ms_sim_run)
mod_names <- c("SS", "SS sim", "MS", "MS sim")

plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names)


################################################
# Model variants
################################################

# For recruitment, the model can estimate recruitment deviates as random effects
ss_re <- Rceattle::fit_mod(
  data_list = mydata,
  inits = ss_re$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = TRUE, # Turn of recruitment deviations as random effects
  msmMode = 0, # Single species mode
  verbose = 1,
  phase = NULL)

# Diet estimation
ms_gamma <- Rceattle::fit_mod(
  data_list = mydata,
  inits = ss_run$estimated_params, # Initial parameters from single species ests
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  niter = 5, # 10 iterations around population and predation dynamics
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # MSVPA based
  suitMode = 1, # Have a gamma function based on predator-prey length ratio for suitability. Includes diet proportion by weight in likelihood as multinomial
  verbose = 1)
# Can try different functions. Look at ?fit_mod suitmode

