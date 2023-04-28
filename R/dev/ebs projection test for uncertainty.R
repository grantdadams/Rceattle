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
mydata$fleet_control$proj_F_prop <-rep(1,7)
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)

# Estimate single-species and estimate M
mydata_M <- mydata
mydata_M$est_M1 <- c(1,1,1)
ss_run_M <- Rceattle::fit_mod(data_list = mydata_M,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = "default",
                              verbose = 1)

# -- PFMC Category 1
ss_run_M_Cat1 <- Rceattle::fit_mod(data_list = mydata_M,#;
                                   inits = ss_run_M$estimated_params,#; # Initial parameters from ss_run
                                   estimateMode = 0,#; # Run projection only
                                   updateM1 = FALSE,#;
                                   HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                   FsprLimit = 0.45, # F45%
                                                   Ptarget = 0.4, # Target is 40% B0
                                                   Plimit = 0.1, # No fishing when SB<SB10
                                                   Pstar = 0.45,
                                                   Sigma = 0.5),#;
                                   msmMode = 0,#; # Single species mode
                                   verbose = 1,#;
                                   projection_uncertainty = TRUE)

plot_biomass(list(ss_run_M, ss_run_M_Cat1), add_ci = TRUE, incl_proj = TRUE)
