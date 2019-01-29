
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
library(Rceattle)
data(BS2017SS) # ?BS2017SS for more information on the data

# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle(TMBfilename = "ceattle_v01_02",
                   cpp_directory = "inst/executables",
                   data_list = BS2017SS,
                   inits = NULL, # Initial parameters = 0
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   random_rec = FALSE, # No random recruitment
                   msmMode = 0, # Single species mode
                   avgnMode = 0,
                   silent = TRUE)
# The you can plot the model results using using
plot_selectivity(ss_run)
  plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run)


# For the 2017 multispecies model starting from the single species parameters, the following can be specified:
data(BS2017MS) # ?BS2017MS for more information on the data
library(Rceattle)
ms_run <- Rceattle(TMBfilename = "ceattle_v01_02",
                    cpp_directory = "inst/executables",
                    data_list = BS2017MS,
                    inits = ss_run$estimated_params, # Initial parameters from single species ests
                    file_name = NULL, # Don't save
                    debug = 0, # Estimate
                    niter = 10,
                    random_rec = FALSE, # No random recruitment
                    msmMode = 1, # Holsman et al empirical suitability
                    avgnMode = 0,
                    silent = TRUE)

plot_biomass(list(ss_run, ms_run, mod_objects), model_names = c("SS", "MS", "MS_2"),  file_name = "tests/example")
plot_recruitment(list(ss_run, ms_run, mod_objects), model_names = c("SS", "MS", "MS_2"),  file_name = "tests/example")

ms_run2 <- Rceattle(TMBfilename = "ceattle_v01_02",
                   cpp_directory = "inst/executables",
                   data_list = BS2017MS,
                   inits = ss_run$estimated_params, # Initial parameters from single species ests
                   file_name = "tests/example", # Don't save
                   debug = 0, # Estimate
                   niter = 10,
                   random_rec = FALSE, # No random recruitment
                   msmMode = 1, # Holsman et al empirical suitability
                   avgnMode = 0,
                   silent = TRUE)

plot_biomass(list(ss_run, ms_run, ms_run2), model_names = c("SS", "MS", "MS_diet"),  file_name = "tests/example")

# Estimate diet
ms_run_diet <- Rceattle(TMBfilename = "ceattle_v01_02",
                        cpp_directory = "inst/executables",
                        data_list = BS2017MS,
                        inits = ms_run$estimated_params, # Initial parameters from single species ests
                        file_name = "tests/example", # Don't save
                        debug = 0, # Estimate
                        niter = 10,
                        random_rec = FALSE, # No random recruitment
                        msmMode = 1, # Holsman et al
                        avgnMode = 0,
                        suitMode = 2, # Weight based lognormal suit
                        silent = FALSE,
                        est_diet = TRUE)

plot_biomass(list(ss_run, ms_run, ms_run2), model_names = c("SS", "MS", "MS_diet"),  file_name = "tests/example")
plot_recruitment(list(ss_run, ms_run, ms_run_diet), model_names = c("SS", "MS", "MS_diet"),  file_name = "tests/example")
