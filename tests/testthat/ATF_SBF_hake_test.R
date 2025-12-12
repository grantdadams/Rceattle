library(Rceattle)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)

#access toke: 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

SBF_ATF_hakedata <- Rceattle::read_data(file = "241025_SBF_ATF_Hake.xlsx")

ss_run <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run$quantities$jnll #2217.028
#save(ss_run, file = "results/models/ATF_SBF/ATF_SBF_ss_run.Rdata")

ms_run <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            M1Fun = build_M1(M1_model = 0,  #do not estimate mortality!
                                             updateM1 = FALSE,
                                             M1_use_prior = FALSE,
                                             M2_use_prior = FALSE),
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            initMode = 2, # Fished start with init devs
                            verbose = 1)

ms_run$quantities$jnll #2229.18
ms_run$quantities$jnll_comp
#save(ms_run, file = "results/models/ATF_SBF/ATF_SBF_ms_run.Rdata")

plot_diet_comp(ms_run) #diet estimate hake canni is good, in this model we do not have ATF predation, hence = 0
plot_biomass(ms_run)
# Correct way to access the data sent to TMB
ms_run$obj$env$data$n_stomach_obs
head(ms_run$obj$env$data$stomach_id)

plot_biomass(ms_run)

#####################################################################
# FIX SUITABILITY AND SUM across prey ages (NEW PART)
# Create initial parameter list:
test_data <- SBF_ATF_hakedata

#test_data$diet_data <- test_data$diet_data %>%
#  filter(!(Pred == 2 & Pred_sex == 2)) %>% # Removes the duplicate male ATF data
#  mutate(Pred_sex = ifelse(Pred == 2, 0, Pred_sex)) # Relabels the remaining ATF data as combined sex

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
# Set weight ratio parameters
#inits$log_gam_a = c(0, 3.006)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
#inits$log_gam_b = c(0, 1.887)

#ORIGINAL (TRUE) VALUES
inits$log_gam_a = c(0, 3.7, 3.1)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, 1.8, 1.120)

#ORIGINAL (TRUE) VALUES
#inits$log_gam_a = c(0, 3.3)  # Median log weight ratio for ATF, 0 for other species (pred/prey)
#inits$log_gam_b = c(0, 1.8)

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on ATF
inits$log_phi[2,2] <- -999 # Set ATF do not feed on ATF
inits$log_phi[1,3] <- -999 #Fixing so hake do not prey on SBF
inits$log_phi[3,3] <- -999 # Set SBF do not feed on SBF
inits$log_phi[2,3] <- -999 # Set ATF do not feed on SBF
inits$log_phi[3,2] <- -999 # Set SBF do not feed on ATF

# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on atf
map$mapList$log_phi[2,2] <- NA #so we dont estimate atf on atf
map$mapList$log_phi[1,3] <- NA #so we dont estimate on hake on SBF
map$mapList$log_phi[3,3] <- NA #so we dont estimate on SBF on SBF
map$mapList$log_phi[2,3] <- NA #so we dont estimate on atf on sbf
map$mapList$log_phi[3,2] <- NA #so we dont estimate sbf on atf

map$mapFactor$log_phi <- factor(map$mapList$log_phi)

run_ms_LN_3 <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 1,
                                                  updateM1 = TRUE,
                                                  M1_use_prior = FALSE,
                                                  M2_use_prior = FALSE),
                                 file = NULL, # Don't save
                                 estimateMode = 0, # estimate
                                 niter = 3, # 3 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 loopnum = 5,
                                 phase = TRUE,
                                 suitMode = c(0, 4, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_3$quantities$jnll #2762.643
run_ms_LN_3$quantities$jnll_comp
plot_diet_comp(run_ms_LN_3)
plot_diet_comp2(run_ms_LN_3)

save(run_ms_LN_3, file = "results/models/ATF_SBF/ATF_SBF_ms_LN_estM.Rdata")

mod_list <- list(ss_run, ms_run, run_ms_LN_3, run_ms_LN_4, run_ms_LN_4_prior)
mod_names <- c("ss_run", "Cannib", "run_ms_LN_estM", "comb", "comb_P")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)

########
## Sex-combined diet
test_data$diet_data <- test_data$diet_data %>%
  filter(!(Pred %in% c(2, 3) & Pred_sex == 2)) %>%  # Removes the duplicate male ATF data
  mutate(Pred_sex = ifelse(Pred %in% c(2, 3), 0, Pred_sex))  # Relabels remaining ATF data as combined sex


run_ms_LN_4 <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 1,
                                                  updateM1 = TRUE,
                                                  M1_use_prior = FALSE,
                                                  M2_use_prior = FALSE),
                                 file = NULL, # Don't save
                                 estimateMode = 0, # estimate
                                 niter = 3, # 3 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 loopnum = 5,
                                 phase = TRUE,
                                 suitMode = c(0, 4, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_4$quantities$jnll #2762.643
run_ms_LN_4$quantities$jnll_comp
plot_diet_comp(run_ms_LN_4)
plot_diet_comp2(run_ms_LN_4)
save(run_ms_LN_4, file = "results/models/ATF_SBF/ATF_SBF_ms_LN_estM_SexComb.Rdata")

run_ms_LN_4_prior <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 1,
                                                  updateM1 = TRUE,
                                                  M1_use_prior = 0.22,
                                                  M2_use_prior = 0.31),
                                 file = NULL, # Don't save
                                 estimateMode = 0, # estimate
                                 niter = 3, # 3 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 loopnum = 5,
                                 phase = TRUE,
                                 suitMode = c(0, 4, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_4_prior$quantities$jnll #2463.555
run_ms_LN_4$quantities$jnll_comp
plot_diet_comp(run_ms_LN_4)
plot_diet_comp2(run_ms_LN_4_prior)

### Biological reference points
run_ms_LN_5 <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 1,
                                                  updateM1 = TRUE,
                                                  M1_use_prior = 0.22,
                                                  M2_use_prior = 0.31),
                                 file = NULL, # Don't save
                                 estimateMode = 0, # estimate
                                 niter = 3, # 3 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 loopnum = 5,
                                 phase = TRUE,
                                 suitMode = c(0, 4, 4), # empirical + LN suitability
                                 initMode = 2,
                                 HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                 DynamicHCR = FALSE,
                                                 FsprLimit = 0.4, # F40%
                                                 Ptarget = 0.4, # Target is 40% B0
                                                 Plimit = 0.1, # No fishing when SB<SB10
                                                 Pstar = 0.5,
                                                 Sigma = 0.5),
                                 verbose = 1)
