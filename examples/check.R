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

################################################
# Alternative harvest control rules
################################################
# -- NPFMC Tier 3
ss_run_Tier3 <- Rceattle::fit_mod(data_list = mydata,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 2, # Run projection only
                                  HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                  FsprTarget = 0.4, # F40%
                                                  FsprLimit = 0.35, # F35%
                                                  Plimit = 0.2, # No fishing when SB<SB20
                                                  Alpha = 0.05),
                                  msmMode = 0, # Single species mode
                                  verbose = 1)



ss_run_dynamicTier3 <- Rceattle::fit_mod(data_list = mydata,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 2, # Run projection only
                                         HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         FsprTarget = 0.4, # F40%
                                                         FsprLimit = 0.35, # F35%
                                                         Plimit = 0.2, # No fishing when SB<SB20
                                                         Alpha = 0.05),
                                         msmMode = 0, # Single species mode
                                         verbose = 1)
ss_run_dynamicTier3$quantities$DynamicSPRtarget/ss_run_dynamicTier3$quantities$DynamicSPR0

data_list = mydata
inits = ss_run$estimated_params # Initial parameters from ss_run
estimateMode = 2 # Run projection only
HCR = build_hcr(HCR = 5, # Tier3 HCR
                DynamicHCR = TRUE, # Use dynamic reference points
                FsprTarget = 0.4, # F40%
                FsprLimit = 0.35, # F35%
                Plimit = 0.2, # No fishing when SB<SB20
                Alpha = 0.05)
msmMode = 0 # Single species mode
verbose = 1
