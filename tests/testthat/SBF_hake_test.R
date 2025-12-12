library(Rceattle)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)

SBF_hakedata <- Rceattle::read_data(file = "070725_Sablefish_Hake_model_Neffhigh.xlsx")
SBF_hakedata$endyr<- 2014

## SINGLE SPSS MODEL
ss_run <- Rceattle::fit_mod(data_list = SBF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run2$quantities$jnll #2109.864
#save(ss_run, file = "results/models/SBF/SBF_ss_run.Rdata")

## MULISPECIS MODEL - CANNIBALISM
ms_run <- Rceattle::fit_mod(data_list = SBF_hakedata,
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

ms_run$quantities$jnll #2122.016
ms_run$quantities$jnll #2122.016
#save(ms_run, file = "results/models/SBF/SBF_ms_run.Rdata")

plot_diet_comp(ms_run) #diet estimate hake canni is good (SBF on hake = 0)

# Correct way to access the data sent to TMB
ms_run$obj$env$data$n_stomach_obs
head(ms_run$obj$env$data$stomach_id)

#####################################################################
# FIX SUITABILITY AND SUM across prey ages (NEW PART)
# Create initial parameter list:
test_data <- SBF_hakedata

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
# Set weight ratio parameters
#inits$log_gam_a = c(0, 2.5)  # Mean log weight ratio for SBF, 0 for other species (pred/prey)
#inits$log_gam_b = c(0, 1.120)  # Standard deviation of log weight ratio for SBF, 0 for other species

##ORIGINAL (TRUE) VALUES
inits$log_gam_a = c(0, 3.0)  # Mean log weight ratio for SBF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, 1.120)

##TEST values
#inits$log_gam_a = c(0, 2.94)  # Median log weight ratio for SBF, 0 for other species (pred/prey)
#inits$log_gam_b = c(0, 1.120)

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on SBF
inits$log_phi[2,2] <- -999 # Set SBF do not feed on SBF
#inits$log_phi[2,1] <- 4

# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on SBF
map$mapList$log_phi[2,2] <- NA #so we dont estimate SBF on SBF

map$mapFactor$log_phi <- factor(map$mapList$log_phi)

#In this model we are estimating estimateMode = 0
run_ms_LN_1 <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 0,
                                                  updateM1 = FALSE,
                                                  M1_use_prior = FALSE,
                                                  M2_use_prior = FALSE),
                                 file = NULL, # Don't save
                                 estimateMode = 0, # estimate
                                 niter = 3, # 3 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 loopnum = 5,
                                 phase = TRUE,
                                 suitMode = c(0, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_1$quantities$jnll #2345.017
run_ms_LN_1$quantities$jnll_comp #Stomach content data 46.33596 / 210.558429

plot_diet_comp2(run_ms_LN_1)
plot_b_eaten_prop(run_ms_LN_1)

#In this model we are updating M1
run_ms_LN_2 <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 0,
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
                                 suitMode = c(0, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_2$quantities$jnll #2407.375
run_ms_LN_2$quantities$jnll_comp #Stomach content data 46.33596 /210.558429
run_ms_LN_2$quantities$diet_hat
run_ms_LN_2$quantities$diet_prop_hat
run_ms_LN_2$quantities$jnll
run_ms_LN_1$estimated_params$log_phi
run_ms_LN_3$quantities$vulnerability
run_ms_LN_2$estimated_params$diet_comp_weights
run_ms_LN_1$data_list$Diet_weights_mcallister

plot_diet_comp(run_ms_LN_2)
plot_b_eaten_prop(run_ms_LN_2)

#In this model we are estimating M1
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
                                 suitMode = c(0, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_3$quantities$jnll # 2387.54
save(run_ms_LN_3, file = "results/models/SBF/SBF_run_ms_LN_estM.Rdata")
run_ms_LN_3$quantities$jnll_comp #Stomach content data 46.33821 /235.935094
plot_diet_comp2(run_ms_LN_3)

##NOTE: IN THIS MODEL, DIET IS NOT SEX SPECIFIC. WE ASSUME THAT MALE AND FEMALE HAS THE SAME PROP OF HAKE IN DIET.
##BUT FEMALE SBF IS BIGGER THAN MALE SBF, SO THERE IS LIKELY SOME DIFFERENCES THERE WE CAN NOT PICK.

##Solution?
#### AVERAGE ACORSS SEX ########
test_data$diet_data <- test_data$diet_data %>%
  filter(!(Pred == 2 & Pred_sex == 2)) %>% # Removes the duplicate male SBF data
  mutate(Pred_sex = ifelse(Pred == 2, 0, Pred_sex)) # Relabels the remaining SBF data as combined sex


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
                                 suitMode = c(0, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_4$quantities$jnll #2263.066
run_ms_LN_4$quantities$jnll_comp #Stomach content data 46.33612 /98.938432
save(run_ms_LN_4, file = "results/models/SBF/SBF_run_ms_LN_estM_SexComb.Rdata")

plot_diet_comp2(run_ms_LN_4)

## BIOLOGICAL REFERENCE POINTS
run_ms_LN_BRP <- Rceattle::fit_mod(data_list = test_data,
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
                                 suitMode = c(0, 4), # empirical + LN suitability
                                 initMode = 2,
                                 HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                 DynamicHCR = FALSE,
                                                 Ftarget = 0.4,
                                                 Flimit  = 0.4, # F40%
                                                 Ptarget = 0.4, # Target is 40% B0
                                                 Plimit = 0.1, # No fishing when SB<SB10
                                                 Pstar = 0.5,
                                                 Sigma = 0.5),
                                 verbose = 1)

run_ms_LN_BRP$quantities$jnll #2263.161
plot_diet_comp2(run_ms_LN_BRP)
run_ms_LN_BRP$quantities$B0 #Unfished total biomass
run_ms_LN_BRP$quantities$SB0 #Unfished spawning biomass
run_ms_LN_BRP$quantities$R0 #Unfished recruitment
run_ms_LN_BRP$quantities$Ftarget #Target fishing mortality rate


###plots
##plot single species models
load("results/models/SBF/SBF_run_ms_LN_estM.Rdata")
load("results/models/SBF/SBF_ss_run.Rdata")
load("results/models/SBF/SBF_ms_run.Rdata")

mod_list <- list(ss_run, ms_run, run_ms_LN_3, run_ms_LN_1)
mod_names <- c("ss_run", "Cannib", "run_ms_LN_estM", "run_ms_LN_noEst")

models <- list(ss_run, ms_run, run_ms_LN_3)
model_names <-  c("ss_run", "Cannib", "run_ms_LN_estM")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)

# Index and catch fits can show multiple models
plot_index(mod_list, model_names= mod_names)
plot_logindex(mod_list, model_names= mod_names)
plot_catch(mod_list, model_names= mod_names)

# Composition plots can be only show one model at a time
plot_comp(run_ms_LN_3) # - Produces many plots!
plot_comp(run_ms_LN_4)

# Plot ration
plot_ration(run_ms_LN_3)
plot_ration(run_ms_LN_4)

plot_selectivity(run_ms_LN_3)
plot_selectivity(run_ms_LN_4)

# Plot mortality-at-age
library(ggplot2)
plot_mortality(run_ms_LN_3)
plot_mortality(run_ms_LN_4)

# Look at jnll composition
# - pull “dev-name-change” branch to get unweighted likelihood components!
run_ms_LN_3$quantities$jnll_comp
sum(run_ms_LN_3$quantities$jnll_comp[-20,]) # - Sum across everything except diet likelihood
sum(run_ms_LN_3$quantities$unweighted_jnll_comp[-20,]) # - Excludes weights

sum(ms_run$quantities$jnll_comp[-20,]) # - Sum across everything except diet likelihood
sum(ms_run$quantities$unweighted_jnll_comp[-20,]) # - Excludes weights

# Get residual mortality (time variant M1)
run_ms_LN_3$quantities$M1_at_age[1,1,,]

## other exploration
# Predator: 1=Hake_F, 2=SBF_F, 3=Hake_M, 4=SBF_M
# Prey:     1=Hake_F, 2=SBF_F, 3=Hake_M, 4=SBF_M
sum(run_ms_LN_3$quantities$suitability[2,1,1:5 ,1:5 ,1])
run_ms_LN_3$quantities$suit_other[2,1,1:5 ,1]
run_ms_LN_3$quantities$avail_food[2,1,1:5 ,1]
run_ms_LN_3$quantities$vulnerability
run_ms_LN_3$quantities$vulnerability_other
run_ms_LN_3$data_list$other_food
run_ms_LN_3$quantities$diet_hat

run_ms_LN_3$quantities$consumption_at_age[1,1,1:5 ,1]

plot_mortality(run_ms_LN_3, type = 3)

