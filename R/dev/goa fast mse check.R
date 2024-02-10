## Load packages and data ----
pacman::p_load(Rceattle, readxl, dplyr, tidyr, writexl)
combined_data <- read_data(file = "~/GitHub/Rceattle_MSE/Data/GOA_23_1_1_data_1977_2023_edited.xlsx")
combined_data$projyr <- 2100
alpha = exp(c(3.143, 1.975, 1.44))


## Estimate OMs ----
# OM 1) Single-spp fix M ----
# * Density-independent recruitment ----
# - Climate naive
ss_mod <- Rceattle::fit_mod(data_list = combined_data,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            phase = "default",
                            initMode = 1)

# - Adjust Fprop for Cod
avg_F <- (exp(ss_mod$estimated_params$ln_mean_F+ss_mod$estimated_params$F_dev)) # Average F from last 2 years
avg_F <- rowMeans(avg_F[,(ncol(avg_F)-2) : ncol(avg_F)])
f_ratio <- avg_F[14:16]
f_ratio <- f_ratio/sum(f_ratio)

# Adjust future F proportion to each fleet
ss_mod$data_list$fleet_control$proj_F_prop <- c(rep(0, 7), 1,0,0,1, 0,0, f_ratio, 0, 0)
combined_data$fleet_control$proj_F_prop <- c(rep(0, 7), 1,0,0,1, 0,0, f_ratio, 0, 0)
ss_mod$estimated_params$proj_F_prop <- ss_mod$data_list$fleet_control$proj_F_prop


## Estimate EMs ----
# EM 1) Single-spp fix M ----
# - Tier-3
ss_run_Tier3 <- Rceattle::fit_mod(data_list = combined_data,
                                  inits = ss_mod$estimated_params, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 0, # Single species mode
                                  verbose = 1,
                                  phase = NULL,
                                  initMode = 1,
                                  HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                  DynamicHCR = FALSE, # Dont dynamic reference points
                                                  FsprTarget = 0.4, # F40%
                                                  FsprLimit = 0.35, # F35%
                                                  Plimit = c(0.2, 0, 0.2), # No fishing when SB<SB20
                                                  Alpha = 0.05))

plot_biomass(list(ss_mod, ss_run_Tier3), incl_proj = T)


ms_run_fb40iter <- Rceattle::fit_mod(data_list = combined_data,
                                     inits = ss_mod$estimated_params, # Initial parameters = 0
                                     file = NULL, # Don't save
                                     estimateMode = 0, # Estimate
                                     random_rec = FALSE, # No random recruitment
                                     msmMode = 1, # Multi species mode
                                     verbose = 2,
                                     niter = 3,
                                     suit_meanyr = 2018,
                                     phase = "default",
                                     M1Fun = build_M1(M1_model = c(1,2,1),
                                                      M1_use_prior = FALSE,
                                                      M2_use_prior = FALSE),
                                     HCR = build_hcr(HCR = 3, # Constant F HCR
                                                     DynamicHCR = FALSE, # Use dynamic reference points
                                                     FsprTarget = 0.4,
                                                     HCRorder = c(1,2,1)), # F that achieves 40% SB0
                                     initMode = 1)


## MSE ----
