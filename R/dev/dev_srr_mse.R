library(Rceattle)
library(dplyr)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060

data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
BS2017MS$projyr <- 2060

BS2017SS$fleet_control$proj_F_prop <-rep(1,7)
BS2017MS$fleet_control$proj_F_prop <- rep(1, 7)

#
# BS2017SS$fsh_biom <- BS2017SS$fsh_biom %>%
#   mutate(Catch = ifelse(Year > 2000 & grep( "Cod", BS2017SS$fsh_biom$Fleet_name), Catch * 4,
#                         Catch)
#   )


################################################
# Estimate OMs
################################################
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1,
                            initMode = 2)

# - Single-species fixed M
alpha = exp(c(4.121, 2.119, 1.553))
ss_run_ricker <- Rceattle::fit_mod(
  data_list = BS2017SS,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  M1Fun = build_M1(M1_model = 0,
                   M1_use_prior = FALSE,
                   M2_use_prior = FALSE),
  recFun = build_srr(srr_fun = 0,
                     srr_pred_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 0.2),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  phase = "default",
  verbose = 1,
  initMode = 2)



# -- NPFMC Tier 3
ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = ss_run$estimated_params,
                                  estimateMode = 0, # Run projection only
                                  HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                  FsprTarget = 0.4, # F40%
                                                  FsprLimit = 0.35, # F35%
                                                  Plimit = c(0.2, 0.2, 0), # No fishing when SB<SB20
                                                  Alpha = 0.05),
                                  msmMode = 0, # Single species mode
                                  verbose = 1,
                                  initMode = 2)


ss_run_dynamicTier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 0, # Run projection only
                                         HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         FsprTarget = 0.4, # F40%
                                                         FsprLimit = 0.35, # F35%
                                                         Plimit = c(0.2, 0.2, 0), # No fishing when SB<SB20
                                                         Alpha = 0.05),
                                         msmMode = 0, # Single species mode
                                         verbose = 1,
                                         initMode = 2)


plot_depletionSSB(list(ss_run, ss_run_Tier3, ss_run_dynamicTier3), incl_proj = TRUE, model_names = 1:3)
plot_ssb(list(ss_run, ss_run_Tier3, ss_run_dynamicTier3), incl_proj = TRUE)

check <- ss_run_dynamicTier3
check$quantities$biomassSSB <- check$quantities$DynamicSB0
plot_ssb(list(ss_run, ss_run_Tier3, check, ss_run_dynamicTier3), incl_proj = TRUE)

ss_run_dynamicTier3$quantities$DynamicNByage0[1,1,1,1:10]
ss_run_dynamicTier3$quantities$R[1,1:10]





# Update OM to have save HCR as EM to calculate performance metrics
om <- ss_run_ricker
em <- ss_run_dynamicTier3

om <- Rceattle::fit_mod(
  data_list = om$data_list,
  inits = om$estimated_params,
  map =  NULL,
  bounds = NULL,
  file = NULL,
  estimateMode = 2, # Run projection only
  HCR = build_hcr(HCR = em$data_list$HCR, # Tier3 HCR
                  DynamicHCR = em$data_list$DynamicHCR,
                  FsprTarget = em$data_list$FsprTarget,
                  FsprLimit = em$data_list$FsprLimit,
                  Ptarget = em$data_list$Ptarget,
                  Plimit = em$data_list$Plimit,
                  Alpha = em$data_list$Alpha,
                  Pstar = em$data_list$Pstar,
                  Sigma = em$data_list$Sigma
  ),
  recFun = build_srr(srr_fun = om$data_list$srr_fun,
                     srr_pred_fun  = om$data_list$srr_pred_fun ,
                     proj_mean_rec  = om$data_list$proj_mean_rec ,
                     srr_est_mode  = om$data_list$srr_est_mode ,
                     srr_prior_mean  = om$data_list$srr_prior_mean,
                     srr_prior_sd   = om$data_list$srr_prior_sd,
                     Bmsy_lim = om$data_list$Bmsy_lim),
  M1Fun =     build_M1(M1_model= om$data_list$M1_model,
                       updateM1 = FALSE,
                       M1_use_prior = om$data_list$M1_use_prior,
                       M2_use_prior = om$data_list$M2_use_prior,
                       M1_prior_mean = om$data_list$M1_prior_mean,
                       M1_prior_sd = om$data_list$M1_prior_sd),
  random_rec = om$data_list$random_rec,
  niter = om$data_list$niter,
  msmMode = om$data_list$msmMode,
  avgnMode = om$data_list$avgnMode,
  minNByage = om$data_list$minNByage,
  suitMode = om$data_list$suitMode,
  initMode = om$data_list$initMode,
  phase = NULL,
  loopnum = 3,
  getsd = FALSE,
  verbose = 0)


print(paste0("Running OM ",om, " and EM ", em))

# Run MSE
mse <- mse_run_parallel(om = om,
                        em = em,
                        nsim = 1,
                        start_sim = 1,
                        seed = 666,
                        regenerate_seed = 666,
                        assessment_period = 1,
                        simulate_data = TRUE,
                        sample_rec = TRUE,
                        cap = NULL,
                        dir = "check",
                        file = NULL,
                        regenerate_past = TRUE)
EMs_from_OM_Sim_1 <- readRDS("check/EMs_from_OM_Sim_1.rds")
mod_list <- list(ss_run_dynamicTier3, ss_run_ricker, EMs_from_OM_Sim_1$OM, EMs_from_OM_Sim_1$EM$`OM_Sim_1. EM_yr_2060`)
plot_ssb(mod_list, incl_proj = TRUE)
plot_stock_recruit(mod_list)
plot_depletionSSB(mod_list, incl_proj = TRUE)
