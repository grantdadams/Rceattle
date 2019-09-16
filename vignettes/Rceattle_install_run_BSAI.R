# Grant Adams, Kirstin Holsman, Andre Punt - April 2019
# Function to run CEATTLE model in TMB
# Citation:
# Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., and Moffitt, E. A. 2015. A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134: 360–378.

# Install devtools if you don't already have it
install.packages("devtools")
# Install TMB and rtools https://cran.r-project.org/bin/windows/Rtools/
# Instructions can be found here for non-pc: https://github.com/kaskr/adcomp/wiki/Download
install.packages('TMB', type = 'source')
# Try "TMB::runExample(all = TRUE)" to see if TMB works

# Install Rceattle
devtools::install_github("grantdadams/Rceattle", auth_token = "4925b42ac46f1e0aefd671e9dc0c1cf1b3157017")



library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data

# Write data to excel
# Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")
#
# # Change the data how you want in excel
# # Read the data back in
# mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")
phaseList = list(
      dummy = 1,
      ln_mn_rec = 1,
      ln_rec_sigma = 2,
      rec_dev = 2,
      init_dev = 2,
      ln_mean_F = 1,
      proj_F = 1,
      F_dev = 1,
      log_srv_q = 3,
      ln_srv_q_dev = 4,
      ln_srv_q_dev_re = 4,
      ln_sigma_srv_q = 4,
      sel_coff = 3,
      sel_slp = 3,
      sel_inf = 3,
      sel_slp_dev = 4,
      sel_inf_dev = 4,
      sel_slp_dev_re = 4,
      sel_inf_dev_re = 4,
      ln_sigma_sel = 4,
      ln_sigma_srv_index = 2,
      ln_sigma_fsh_catch = 2,
      logH_1 = 6,
      logH_1a = 6,
      logH_1b = 6,
      logH_2 = 6,
      logH_3 = 6,
      H_4 = 6,
      log_gam_a = 5,
      log_gam_b = 5,
      log_phi = 5
  )



################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            TMBfilename = "ceattle_v01_06",
                            cpp_directory = "inst/executables",
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = phaseList,
                            silent = TRUE,
                            recompile = FALSE)

# The you can plot the model results using using
plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run)


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            TMBfilename = "ceattle_v01_06",
                            cpp_directory = "inst/executables",
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            niter = 10, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            silent = TRUE,
                            recompile = FALSE)


# We can plot both runs as well:
mod_list <- list(ss_run, ms_run)
mod_names <- c("SS", "MS")

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
