
# Load
library(Rceattle)

# Run the 2017 single species assessment for the Bering Sea
data(BS2017SS)
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


# Update composition weights using McAllister-Ianelli method
# - I estimate the weights after model optimization in the fit_mod function
BS2017SS_weighted <- BS2017SS
BS2017SS_weighted$fleet_control$Comp_weights <- ss_run$data_list$fleet_control$Est_weights_mcallister

ss_run_reweighted <- Rceattle::fit_mod(data_list = BS2017SS_weighted,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)
