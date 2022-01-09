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



# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
BS2017MS$est_M1 <- c(1,1,1) # Estimate residual M
BS2017MS$projyr <- 2060
ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run_M$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            niter = 3, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

# We can plot both runs as well:
mod_list <- list(ss_run, ss_run_M, ms_run)
mod_names <- c("SS", "SS-M", "MS")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)


################################################
# Alternative harvest control rules
################################################
# -- Constant F as a percentage of SB0
ss_run_fb0 <- Rceattle::fit_mod(data_list = mydata,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 2, # Run projection only
                                  HCR = build_hcr(HCR = 3, # Tier3 HCR
                                                  DynamicHCR = FALSE, # Use dynamic reference points
                                                  FsprTarget = 0.4), # F that achieves 40% SB0
                                  msmMode = 0, # Single species mode
                                  verbose = 1)


ss_run_dynamicfb0 <- Rceattle::fit_mod(data_list = mydata,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 2, # Run projection only
                                         HCR = build_hcr(HCR = 3, # Tier3 HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         FsprTarget = 0.4), # F that achieves 40% SB0
                                         msmMode = 0, # Single species mode
                                         verbose = 1)


# -- Constant Fspr
ss_run_Fspr <- Rceattle::fit_mod(data_list = mydata,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 2, # Run projection only
                                  HCR = build_hcr(HCR = 4, # Tier3 HCR
                                                  FsprTarget = 0.4 # F40%
                                                  ),
                                  msmMode = 0, # Single species mode
                                  verbose = 1)


ss_run_dynamicFspr <- Rceattle::fit_mod(data_list = mydata,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 2, # Run projection only
                                         HCR = build_hcr(HCR = 4, # Tier3 HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         FsprTarget = 0.4 # F40%
                                                         ),
                                         msmMode = 0, # Single species mode
                                         verbose = 1)


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

# -- PFMC Category 1
ss_run_Cat1 <- Rceattle::fit_mod(data_list = mydata,
                                 inits = ss_run$estimated_params, # Initial parameters from ss_run
                                 estimateMode = 2, # Run projection only
                                 HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                 FsprLimit = 0.45, # F45%
                                                 Ptarget = 0.4, # Target is 40% B0
                                                 Plimit = 0.1, # No fishing when SB<SB10
                                                 Pstar = 0.45,
                                                 Sigma = 0.5),
                                 msmMode = 0, # Single species mode
                                 verbose = 1)

ss_run_dynamicCat1 <- Rceattle::fit_mod(data_list = mydata,
                                        inits = ss_run$estimated_params, # Initial parameters from ss_run
                                        estimateMode = 2, # Run projection only
                                        HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                        DynamicHCR = TRUE, # Use dynamic reference points
                                                        FsprLimit = 0.45, # F45%
                                                        Ptarget = 0.4, # Target is 40% SB0
                                                        Plimit = 0.1, # No fishing when SB<SB10
                                                        Pstar = 0.45,
                                                        Sigma = 0.5),
                                        msmMode = 0, # Single species mode
                                        verbose = 1)

# -- SESSF Tier 1
ss_run_Tier1 <- Rceattle::fit_mod(data_list = mydata,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 2, # Run projection only
                                  HCR = build_hcr(HCR = 7, # Tier 1 HCR
                                                  FsprTarget = 0.48, # F40%
                                                  FsprLimit = 0.20, # F20%
                                                  Ptarget = 0.35, # Target is 35% SSB0
                                                  Plimit = 0.20, # No fishing when B<B20
                                  ),
                                  msmMode = 0, # Single species mode
                                  verbose = 1)


ss_run_dynamicTier1 <- Rceattle::fit_mod(data_list = mydata,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 2, # Run projection only
                                         HCR = build_hcr(HCR = 7, # Tier 1 HCR
                                                         DynamicHCR = TRUE,
                                                         FsprTarget = 0.48, # F40%
                                                         FsprLimit = 0.20, # F20%
                                                         Ptarget = 0.35, # Target is 35% SSB0
                                                         Plimit = 0.20, # No fishing when B<B20
                                         ),
                                         msmMode = 0, # Single species mode
                                         verbose = 1)

# -- Plot
mod_list <- list(ss_run, ss_run_fb0, ss_run_Fspr, ss_run_Tier3, ss_run_Cat1, ss_run_Tier1)
model_names <- c("F=0","F 40% B0", "Fspr 40%", "NPFMC Tier 3", "PFMC Cat 1", "SESSF Tier 1")
plot_biomass(mod_list, model_names = model_names, incl_proj = TRUE)
plot_ssb(mod_list, model_names = model_names, incl_proj = TRUE)


dynamic_mod_list <- list(ss_run, ss_run_dynamicfb0, ss_run_dynamicFspr, ss_run_dynamicTier3, ss_run_dynamicCat1, ss_run_dynamicTier1)
dynamic_model_names <- c("F=0","F 40% B0", "Fspr 40%", "NPFMC Tier 3", "PFMC Cat 1", "SESSF Tier 1")
plot_biomass(dynamic_mod_list, model_names = dynamic_model_names, incl_proj = TRUE)
plot_ssb(dynamic_mod_list, model_names = dynamic_model_names, incl_proj = TRUE)

################################################
# Management strategy evaluation
################################################
# -- No F
# - MS-OM: SS-EM No F
mse1 <- mse_run(om = ms_run, em = ss_run, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM No F
mse2 <- mse_run(om = ss_run_M, em = ss_run, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))


# -- NPFMC Tier 3 HCRs
# - MS-OM: SS-EM Tier 3 HCR
mse3 <- mse_run(om = ms_run, em = ss_run_Tier3, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM Tier 3 HCR
mse4 <- mse_run(om = ss_run_M, em = ss_run_Tier3, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - MS-OM: SS-EM dynamic Tier 3 HCR
mse5 <- mse_run(om = ms_run, em = ss_run_dynamicTier3, nsim = 1, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM dynamic Tier 3 HCR
mse6 <- mse_run(om = ss_run_M, em = ss_run_dynamicTier3, nsim = 1, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))


# -- PFMC Category 1 HCRs
# - MS-OM: SS-EM Tier 3 HCR
mse7 <- mse_run(om = ms_run, em = ss_run_Cat1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM Tier 3 HCR
mse8 <- mse_run(om = ss_run_M, em = ss_run_Cat1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - MS-OM: SS-EM dynamic Tier 3 HCR
mse9 <- mse_run(om = ms_run, em = ss_run_dynamicCat1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM dynamic Tier 3 HCR
mse10 <- mse_run(om = ss_run_M, em = ss_run_dynamicCat1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))


# -- SESSF Tier 1 HCRs
# - MS-OM: SS-EM Tier 3 HCR
mse11 <- mse_run(om = ms_run, em = ss_run_Tier1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM Tier 3 HCR
mse12 <- mse_run(om = ss_run_M, em = ss_run_Tier1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - MS-OM: SS-EM dynamic Tier 3 HCR
mse13 <- mse_run(om = ms_run, em = ss_run_dynamicTier1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))

# - SS-OM: SS-EM dynamic Tier 3 HCR
mse14 <- mse_run(om = ss_run_M, em = ss_run_dynamicTier1, nsim = 50, assessment_period = 2, sampling_period = 2, simulate = TRUE, cap = c(1500000))
