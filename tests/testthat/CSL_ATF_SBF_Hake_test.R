library(Rceattle)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)

#hake = 1
#CSL = 2
#SBF = 3
#ATF= 4


CSL_SBF_ATF_hakedata <- Rceattle::read_data(file = "241025_CSL_SBF_ATF_Hake.xlsx")
CSL_SBF_ATF_hakedata$endyr

ss_run <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run$quantities$jnll #1958.761
#save(ss_run, file = "results/models/ATF_SBF/ATF_SBF_ss_run.Rdata")

ms_run <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
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

ms_run$quantities$jnll #1972.422
ms_run$quantities$jnll_comp
#save(ms_run, file = "results/models/ATF_SBF/ATF_SBF_ms_run.Rdata")

plot_diet_comp2(ms_run) #diet estimate hake canni is good, in this model we do not have ATF predation, hence = 0


# Correct way to access the data sent to TMB
ms_run$obj$env$data$n_stomach_obs
head(ms_run$obj$env$data$stomach_id)

plot_biomass(ms_run)
plot_biomass(ss_run)

mod_list <- list(ss_run, ms_run)
mod_names <- c("ss_run", "ms_run")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)

#####################################################################
# FIX SUITABILITY AND SUM across prey ages (NEW PART)
# Create initial parameter list:

test_data <- CSL_SBF_ATF_hakedata

#test_data$diet_data <- test_data$diet_data %>%
#  filter(!(Pred == 2 & Pred_sex == 2)) %>% # Removes the duplicate male ATF data
#  mutate(Pred_sex = ifelse(Pred == 2, 0, Pred_sex)) # Relabels the remaining ATF data as combined sex

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here

#ORIGINAL (TRUE) VALUES
inits$log_gam_a = c(0, 0, 3.1, 3.7)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, 0, 1.120, 1.8)

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on ATF
inits$log_phi[2,2] <- -999 # Set ATF do not feed on ATF
inits$log_phi[1,3] <- -999 #Fixing so hake do not prey on SBF
inits$log_phi[3,3] <- -999 # Set SBF do not feed on SBF
inits$log_phi[2,3] <- -999 # Set ATF do not feed on SBF
inits$log_phi[3,2] <- -999 # Set SBF do not feed on ATF
inits$log_phi[2,4] <- -999
inits$log_phi[3,4] <- -999
inits$log_phi[3,4] <- -999
inits$log_phi[4,2] <- -999
inits$log_phi[4,3] <- -999
inits$log_phi[4,4] <- -999
inits$log_phi[1,4] <- -999

# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on atf
map$mapList$log_phi[1,3] <- NA #so we dont estimate on hake on SBF
map$mapList$log_phi[1,4] <- NA

map$mapList$log_phi[2,1] <- NA #so we dont estimate atf on atf
map$mapList$log_phi[2,2] <- NA #so we dont estimate atf on atf
map$mapList$log_phi[2,3] <- NA #so we dont estimate on SBF on SBF
map$mapList$log_phi[2,4] <- NA #so we dont estimate on atf on sbf

map$mapList$log_phi[3,2] <- NA #so we dont estimate sbf on atf
map$mapList$log_phi[3,3] <- NA
map$mapList$log_phi[3,4] <- NA

map$mapList$log_phi[4,2] <- NA
map$mapList$log_phi[4,3] <- NA
map$mapList$log_phi[4,4] <- NA


map$mapFactor$log_phi <- factor(map$mapList$log_phi)

run_ms_CSL_debugg3 <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                       inits = inits, # Initial parameters from single species ests
                                       map = map,
                                       M1Fun = build_M1(M1_model = 0,
                                                        updateM1 = FALSE,
                                                        M1_use_prior = FALSE,
                                                        M2_use_prior = FALSE),
                                       file = NULL, # Don't save
                                       estimateMode = 3, # estimate
                                       niter = 3, # 3 iterations around population and predation dynamics
                                       random_rec = FALSE, # No random recruitment
                                       msmMode = 1, # MSVPA based
                                       loopnum = 5,
                                       phase = TRUE,
                                       suitMode = c(0, 0, 4, 4), # empirical + LN suitability
                                       initMode = 2,
                                       verbose = 1)

run_ms_CSL_debugg3$quantities$jnll_comp #NA in CSL Stomach content data

run_ms_CSL_debugg3$estimated_params$d

# Check both sexes for suitability
# 4. Check suitability
suitability <- run_ms_CSL_debugg3$quantities$suitability
suitability[4, 1, 1:10, 1:5, 1]  # Female CSL eating Pollock
suitability[8, 1, 1:10, 1:5, 1]  # Male CSL eating Pollock

# Sum suitability across all prey for BOTH sexes
apply(suitability[4, , , , 1], c(1), sum)  # Female by predator age
apply(suitability[8, , , , 1], c(1), sum)  # Male by predator age

# Check diet_prop_hat for both sexes
diet_prop_hat <- run_ms_CSL_debugg3$quantities$diet_prop_hat
diet_prop_hat[4, 1, 1:10, 1:5, 1]  # Female
diet_prop_hat[8, 1, 1:10, 1:5, 1]  # Male

# Sum diet_prop_hat by predator age for both sexes
apply(diet_prop_hat[4, , , , 1], c(2), sum)  # Female
apply(diet_prop_hat[8, , , , 1], c(2), sum)  # Male

# Check avail_food for both sexes
avail_food <- run_ms_CSL_debugg3$quantities$avail_food
avail_food[4, 1, 1:20, 1]  # Female CSL
avail_food[4, 2, 1:20, 1]  # Male CSL

# Check diet data sex structure
# Filter for CSL
csl_diet <- diet_data_with_pred %>%
  filter(Pred == 4) %>%
  select(Pred, Prey, Pred_sex, Pred_age, Prey_age, Year,
         Stomach_proportion_by_weight, Est_diet)


print(csl_diet)
table(csl_diet_check$Pred_sex)

# If Pred_sex shows 1 and 2, then you have sex-specific diet data
# If Pred_sex shows 0, then diet data is sex-aggregated but model is 2-sex

# Check avgN_at_age for CSL
avgN_at_age <- run_ms_CSL_debugg3$quantities$avgN_at_age
avgN_at_age[4, 1, 1:24, 1]  # Female CSL
avgN_at_age[4, 2, 1:24, 1]  # Male CSL

# Check ration for CSL
run_ms_CSL_debugg3$quantities$ration[4, 1, 1:24, 1]  # Female
run_ms_CSL_debugg3$quantities$ration[4, 2, 1:24, 1]  # Male

# Get CSL diet rows
csl_diet_rows <- which(run_ms_CSL_debugg3$data$diet_data$Pred %in% c(4, 8))

# Look at the structure of predicted vs observed
diet_check <- data.frame(
  Pred = run_ms_CSL_debugg3$data$diet_data$Pred[csl_diet_rows],
  Pred_sex = run_ms_CSL_debugg3$data$diet_data$Pred_sex[csl_diet_rows],
  Pred_age = run_ms_CSL_debugg3$data$diet_data$Pred_age[csl_diet_rows],
  Prey_age = run_ms_CSL_debugg3$data$diet_data$Prey_age[csl_diet_rows],
  Observed = run_ms_CSL_debugg3$data$diet_data$Stomach_proportion_by_weight[csl_diet_rows],
  Predicted = diet_hat[csl_diet_rows, 2]
)

# Check for NAs in predicted
sum(is.na(diet_check$Predicted))
sum(is.infinite(diet_check$Predicted))

# Look at summary by predator age
diet_summary <- diet_check %>%
  group_by(Pred, Pred_sex, Pred_age) %>%
  summarise(
    n = n(),
    sum_obs = sum(Observed),
    sum_pred = sum(Predicted),
    any_na = any(is.na(Predicted)),
    .groups = 'drop'
  )

print(diet_summary, n = 50)

# Find problematic combinations
diet_summary %>%
  filter(sum_pred < 0.01 | sum_pred > 0.99 | any_na)

run_ms_CSL_debugg3$

run_ms_CSL_debugg4 <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                        inits = inits, # Initial parameters from single species ests
                                        map = map,
                                        M1Fun = build_M1(M1_model = 0,
                                                         updateM1 = FALSE,
                                                         M1_use_prior = FALSE,
                                                         M2_use_prior = FALSE),
                                        file = NULL, # Don't save
                                        estimateMode = 4, # estimate
                                        niter = 1, # 3 iterations around population and predation dynamics
                                        random_rec = FALSE, # No random recruitment
                                        msmMode = 1, # MSVPA based
                                        loopnum = 5,
                                        phase = TRUE,
                                        suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                        initMode = 2,
                                        verbose = 1)

plot_diet_comp(run_ms_CSL_debugg3) #no estimate ATF & SBF diet
plot_diet_comp(run_ms_CSL_debugg4) #no estimate ATF & SBF diet

run_ms_CSL_debugg3$quantities$jnll_comp # CSL = NA Stomach content data
run_ms_CSL_debugg4$quantities$jnll_comp # CSL = NA Stomach content data

##BUT this doesn't converge..
run_ms_CSL_Mnoest <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                       inits = inits, # Initial parameters from single species ests
                                       map = map,
                                       M1Fun = build_M1(M1_model = 0,
                                                        updateM1 = FALSE,
                                                        M1_use_prior = FALSE,
                                                        M2_use_prior = FALSE),
                                       file = NULL, # Don't save
                                       estimateMode = 0, # estimate
                                       niter = 1, # 3 iterations around population and predation dynamics
                                       random_rec = FALSE, # No random recruitment
                                       msmMode = 1, # MSVPA based
                                       loopnum = 5,
                                       phase = TRUE,
                                       suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                       initMode = 2,
                                       verbose = 1)

run_ms_CSL_Mest <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
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
                                 suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)


run_ms_CSL_Mest$quantities$jnll #
run_ms_CSL_Mnoest$quantities$jnll #

run_ms_LN_3$quantities$jnll_comp
plot_diet_comp(run_ms_CSL_Mnoest)
plot_diet_comp2(run_ms_CSL_Mnoest)

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


devtools::document()
