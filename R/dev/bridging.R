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

load("~/GitHub/BS_ss_run.RData")

inits <- mod_objects$estimated_params
inits$rec_pars <- matrix(9, nrow = mod_objects$data_list$nspp, ncol = 2)
inits$rec_pars[,1] <- inits$ln_mean_rec

inits$ln_mean_rec <- NULL

################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            file = NULL,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            phase = "default")
# Type ?fit_mod for more details

plot_biomass(list(ss_run, mod_objects))
plot_recruitment(list(ss_run, mod_objects))

ss_run$estimated_params$init_dev


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
# Or we can use the previous data set
BS2017MS$projyr <- 2030

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 5, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1,
                            initMode = 1)
