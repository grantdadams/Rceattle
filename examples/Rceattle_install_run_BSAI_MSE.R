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
devtools::install_github("grantdadams/Rceattle")



library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")



################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
mydata$fleet_control$proj_F_prop <-rep(1,7) # All species have one fishing fleet
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = FALSE, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 2)

# The you can plot the model results using using
plot_biomass(Rceattle =  ss_run), incl_proj = F)
plot_recruitment(Rceattle =  ss_run, add_ci = TRUE)
plot_catch(Rceattle =  ss_run, incl_proj = FALSE)


# Estimate M
mydata_M <- mydata
mydata_M$est_M1 <- c(1,1,1)
ss_run_M <- Rceattle::fit_mod(data_list = mydata_M,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              debug = FALSE, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = "default",
                              verbose = 1)



# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
BS2017MS$est_M1 <- c(1,1,1)
BS2017MS$projyr <- 2060
ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run_M$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            niter = 3, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

# We can plot both runs as well:
mod_list <- list(ss_run, ms_run)
mod_names <- c("SS", "MS")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_selectivity(Rceattle = mod_list, model_names = mod_names)
plot_mort(Rceattle = mod_list, model_names = mod_names, age = 2)

# Run diagnostics
plot_comp(ms_run) # Fitted fishery composition data
plot_index(ms_run) # Fitted indices of abundance
plot_catch(ms_run) # Fitted catch series

# Run MSE
# - MS-OM: SS-EM
mse1 <- mse_run(operating_model = ms_run, estimation_model = ss_run, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# Run MSE
# - SS-OM: SS-EM
mse2 <- mse_run(operating_model = ss_run_M, estimation_model = ss_run, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))
