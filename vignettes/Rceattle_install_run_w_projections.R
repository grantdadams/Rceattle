# Grant Adams, Kirstin Holsman, Andre Punt - April 2019
# Function to run CEATTLE model in TMB
# Citation:
# Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., and Moffitt, E. A. 2015. A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134: 360–378.

################################################
# Installation
################################################
# Install devtools if you don't already have it
install.packages("devtools")
# Install TMB and rtools
# Instructions can be found here for non-pc: https://github.com/kaskr/adcomp/wiki/Download
install.packages('TMB', type = 'source')
# Try "TMB::runExample(all = TRUE)" to see if TMB works

# Install Rceattle
devtools::install_github("grantdadams/Rceattle", ref = "grant/time_varying_q_sel_2sex2", auth_token = "4925b42ac46f1e0aefd671e9dc0c1cf1b3157017")


################################################
# Setup
################################################
library(Rceattle)


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
                             inits = NULL, # Initial parameters = 0
                             file = NULL, # Don't save
                             debug = 0, # Estimate
                             random_rec = FALSE, # No random recruitment
                             msmMode = 0, # Single species mode
                             silent = TRUE,
                             recompile = FALSE)
# Type ?fit_mod for more details

# The you can plot the model results using using
plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run, add_ci = TRUE)
plot_mortality(Rceattle = ms_run, maxage = 21)

# Note: fitting the model updates the composition weights using the harmonic mean MacCallister-Ianelli method, refitting the model will change the MLE based on the new weighting
ss_run_weighted <- Rceattle::fit_mod(data_list = ss_run$data_list,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            silent = TRUE,
                            recompile = FALSE)

# The you can plot the model results using using
plot_biomass(Rceattle =  list(ss_run, ss_run_weighted))
plot_recruitment(Rceattle =  list(ss_run, ss_run_weighted))


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
# Or we can use the previous data set
BS2017MS$projyr <- 2030

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            niter = 5, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            silent = TRUE)

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
# Rceattle automatically projects the population forward when estimating
# To check to see what the F rates are for each fishery we can check:
BS2017MS$fsh_control$proj_F
BS2017MS$projyr # Year the population is projected forward

# We can then change the F - For example, mean historical F
BS2017MS$fsh_control$proj_F <- c(0.2342936, 0.513, 0.0774777)

# Re-run, without estimating
ms_run_proj <- Rceattle::fit_mod(data_list = BS2017MS,
                                 inits = ms_run$estimated_params, # Initial parameters from single species ests
                                 file = NULL, # Don't save
                                 debug = TRUE, # Do not estimate. Not changing parameters right now
                                 niter = 10, # 10 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 suitMode = 0, # empirical suitability
                                 silent = TRUE)




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
                                  debug = TRUE, # Do not estimate. Not changing parameters right now
                                  niter = 10, # 10 iterations around population and predation dynamics
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 1, # MSVPA based
                                  suitMode = 0, # empirical suitability
                                  silent = TRUE)


# plot
mod_list <- list(ss_run, ms_run, ms_run_proj, ms_run_proj2)
mod_names <- c("SS", "MS no F", "MS mean historical F", "MS mean historical F w random rec")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, incl_proj = TRUE)



################################################
# Retrospective
################################################
# we can also do retrospective peels and calculate Mohn's Rho on the CEATTLE
# NOTE: this is using mean historical F for projections as we changed it above
ms_run_proj_retro <- retrospective(Rceattle = ms_run_proj, peels = 10)

# Look at Mohns rho
ms_run_proj_retro$mohns

# Plot retrospectives
plot_biomass(ms_run_proj_retro$Rceattle_list)

# Seet how forecast changes
plot_biomass(ms_run_proj_retro$Rceattle_list, incl_proj = TRUE, mohns = ms_run_proj_retro$mohns)


################################################
# Simulation
################################################
# Data can be simulated from the estimated quantities using `sim_mod`:
ss_sim <- sim_mod(ss_run)

ss_sim_run <- Rceattle::fit_mod(
  data_list = ss_sim,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  debug = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  silent = TRUE)

# Simulate the multispecies model
ms_sim <- sim_mod(ms_run)

ms_sim_run <- Rceattle::fit_mod(
  data_list = ms_sim,
  inits = ss_run$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  debug = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # Holsman MS mode
  silent = TRUE)

# We can plot the simulated runs
mod_list <- list(ss_sim_run, ms_sim_run)
mod_names <- c("SS sim", "MS sim")

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
  debug = 0, # Estimate
  random_rec = TRUE, # Turn of recruitment deviations as random effects
  msmMode = 0, # Single species mode
  silent = TRUE,
  phase = NULL)

# Diet estimation
ms_gamma <- Rceattle::fit_mod(
  data_list = mydata,
  inits = ss_run$estimated_params, # Initial parameters from single species ests
  file = NULL, # Don't save
  debug = 0, # Estimate
  niter = 10, # 10 iterations around population and predation dynamics
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # MSVPA based
  suitMode = 1, # Have a gamma function with time-independent length ratio for suitability. Includes diet proportion by weight in likelihood as multinomial
  silent = TRUE)
# Can try different functions. Look at ?fit_mod suitmode


ms_gamma2 <- Rceattle::fit_mod(
  data_list = mydata,
  inits = ms_run$estimated_params, # Initial parameters from single species ests
  file = NULL, # Don't save
  debug = 1, # Estimate
  niter = 10, # 10 iterations around population and predation dynamics
  random_rec = FALSE, # No random recruitment
  msmMode = 4, # MSVPA based
  suitMode = 0, # Have a gamma function with time-independent length ratio for suitability. Includes diet proportion by weight in likelihood as multinomial
  silent = FALSE)

ms_gamma2$quantities$jnll_comp
sum(is.nan(ms_gamma2$quantities$T_hat))

################################################
# Included files
################################################
# If we want to extract the cpp file
cpp_directory <- system.file("executables",package="Rceattle")
TMBfilename <- "ceattle_v01_04"
cpp_file <- paste0(cpp_directory, "/", TMBfilename, ".cpp")
cpp_file <- file(cpp_file)
cpp_file <- readLines(cpp_file)

