# Install Rceattle
devtools::install_github("grantdadams/Rceattle", auth_token = "4925b42ac46f1e0aefd671e9dc0c1cf1b3157017")

# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
library(Rceattle)
data(BS2017SS) # ?BS2017SS for more information on the data 

# Write data to excel
Rceattle::write_excel(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_excel( file = "BS2017SS.xlsx")

# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = mydata,
                   inits = NULL, # Initial parameters = 0
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   random_rec = FALSE, # No random recruitment
                   msmMode = 0, # Single species mode
                   silent = TRUE)
# Type ?fit_mod for more details

# The you can plot the model results using using
plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run)


# For the a multispecies model starting from the single species parameters, the following can be specified:
data("BS2017MS")
# Or we can use the previous data set

ms_run <- Rceattle::fit_mod(data_list = mydata,
                   inits = ss_run$estimated_params, # Initial parameters from single species ests
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   niter = 10, # 10 iterations around population and predation dynamics
                   random_rec = FALSE, # No random recruitment
                   msmMode = 1, # MSVPA based
                   suitMode = 0, # empirical suitability
                   silent = TRUE)

# We can plot both runs as well:
mod_list <- list(ss_run, ms_run)
mod_names <- c("SS", "MS")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names)

plot_selectivity(Rceattle = mod_list, model_names = mod_names)
plot_mort(Rceattle = mod_list, model_names = mod_names, age = 2)


################################################
# Simulation
################################################
# Data can be simulated from the estimated quantities using `sim_mod`:
ss_sim <- sim_mod(ss_run)

ss_sim_run <- Rceattle::fit_mod(
  data_list = ss_sim,
  inits = NULL, # Initial parameters = 0
  file_name = NULL, # Don't save
  debug = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  silent = TRUE)

# Simulate the multispecies model
ms_sim <- sim_mod(ms_run)

ms_sim_run <- Rceattle::fit_mod(
  data_list = ms_sim,
  inits = NULL, # Initial parameters = 0
  file_name = NULL, # Don't save
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
  data_list = ss_sim,
  inits = NULL, # Initial parameters = 0
  file_name = NULL, # Don't save
  debug = 0, # Estimate
  random_rec = TRUE, # Random recruitment
  msmMode = 0, # Single species mode
  silent = TRUE)

# Diet estimation
ms_gamma <- Rceattle::fit_mod(
  data_list = mydata2,
  inits = ss_run$estimated_params, # Initial parameters from single species ests
  file_name = NULL, # Don't save
  debug = 0, # Estimate
  niter = 10, # 10 iterations around population and predation dynamics
  random_rec = FALSE, # No random recruitment
  msmMode = 1, # MSVPA based
  suitMode = 1, # Have a gamma function with time-independent length ratio for suitability. Includes diet in likelihood as multinomial
  silent = TRUE)
# Can try different functions. Look at ?fit_mod suitmode




# If we want to extract the cpp file
cpp_directory <- system.file("executables",package="Rceattle")
TMBfilename <- "ceattle_v01_04"
cpp_file <- paste0(cpp_directory, "/", TMBfilename, ".cpp")
cpp_file <- file(cpp_file)
cpp_file <- readLines(cpp_file)
