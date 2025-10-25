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
                            file = "ss_run_dev", # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)

# Estimate single-species and estimate M
ss_run_M <- Rceattle::fit_mod(data_list = mydata,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = "ss_run_M_dev", # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              M1Fun = build_M1(M1_model = 1),
                              msmMode = 0, # Single species mode
                              phase = "default",
                              verbose = 1)

# -- PFMC Category 1
ss_run_M_Cat1 <- Rceattle::fit_mod(data_list = mydata,#;
                                   inits = ss_run_M$estimated_params,#; # Initial parameters from ss_run
                                   estimateMode = 0,#; # Run projection only
                                   file = "ss_run_M_Cat1_dev.Rdata",
                                   HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                   FsprLimit = 0.45, # F45%
                                                   Ptarget = 0.4, # Target is 40% B0
                                                   Plimit = 0.1, # No fishing when SB<SB10
                                                   Pstar = 0.45,
                                                   Sigma = 0.5),#;
                                   M1Fun = build_M1(M1_model = 1),
                                   msmMode = 0,#; # Single species mode
                                   verbose = 1,#;
                                   projection_uncertainty = TRUE)


# Add prior to M
ss_run_M_Cat1_prior <- Rceattle::fit_mod(data_list = mydata,#;
                                   inits = ss_run_M$estimated_params,#; # Initial parameters from ss_run
                                   estimateMode = 0,#; # Run projection only
                                   HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                   FsprLimit = 0.45, # F45%
                                                   Ptarget = 0.4, # Target is 40% B0
                                                   Plimit = 0.1, # No fishing when SB<SB10
                                                   Pstar = 0.45,
                                                   Sigma = 0.5),#;
                                   M1Fun = build_M1(M1_model = 1,
                                                    M1_use_prior = TRUE,
                                                    M2_use_prior = FALSE,
                                                    M1_prior_mean = 0.20,
                                                    M1_prior_sd = 0.05),
                                   msmMode = 0,#; # Single species mode
                                   verbose = 1,#;
                                   projection_uncertainty = TRUE)


# Plot
plot_biomass(list(ss_run_M, ss_run_M_Cat1, ss_run_M_Cat1_prior), add_ci = TRUE, incl_proj = TRUE, model_names = c("No F", "PFMC", "PFMC w/ prior"))

# Look at M1
ss_run$quantities$M1[,1,1] # fixed
ss_run_M_Cat1$quantities$M1[,1,1] # no prior
ss_run_M_Cat1_prior$quantities$M1[,1,1] # w/ prior

