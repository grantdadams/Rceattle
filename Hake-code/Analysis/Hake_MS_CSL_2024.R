library(Rceattle)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)


#############################################################################################################
## LOAD DATA
#############################################################################################################

#CSL_SBF_ATF_hakedata <- Rceattle::read_data(file = "hake_yr24_SBF_ATF_CSL_Final.xlsx")
CSL_SBF_ATF_hakedata <- Rceattle::read_data(file = "hake_yr24_SBF_ATF_CSL_JUNE26.xlsx")

#############################################################################################################
## RUN MODELS - single spss and multispss MSVPA (hake cannibalims only)
#############################################################################################################

##Single spss
ss_run <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run$quantities$jnll #2469.184
ss_run$quantities$jnll_comp
save(ss_run, file = "results/models_all_pred_2024_prior2020/ss_run.Rdata")

ss_run$data_list$fleet_control


## Run multispecies model estimating M1
ms_run <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            M1Fun = build_M1(M1_model = 1,  #estimate mortality!
                                             updateM1 = FALSE,
                                             M1_use_prior = FALSE,
                                             M2_use_prior = FALSE),
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability'
                            suit_styr = c(1980, 1980, 1980, 1980),  # hake, ATF, sablefish
                            suit_endyr = c(2019, 2019, 2019, 2019), # hake, ATF, sablefish
                            initMode = 2, # Fished start with init devs
                            verbose = 1)


ms_run$quantities$jnll # 2549.503
ms_run$quantities$jnll_comp
ms_run$quantities$M1
save(ms_run, file = "results/models_all_pred_2024_prior2020/ms_run.Rdata")

plot_diet_comp(ms_run) #diet estimate hake canni is good, in this model we do not have ATF predation, hence = 0
plot_diet_comp2(ms_run)
plot_biomass(ms_run)
plot_ssb(ms_run)


#############################################################################################################
## LOGNORMAL SUITABILITY MODEL + DM LIKELIHOOD DISTRIBUTION
#############################################################################################################

# Create initial parameter list:
test_data <- CSL_SBF_ATF_hakedata

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
#ORIGINAL (TRUE) VALUES
inits$log_gam_a = c(0, log(3.7), log(3.1), 0)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, log(1.83), log(1.120), 0)

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on ATF
inits$log_phi[1,3] <- -999 #Fixing so hake do not prey on SBF
inits$log_phi[1,4] <- -999 #Fixing so hake do not prey on CSL

inits$log_phi[2,2] <- -999 # Set ATF do not feed on ATF
inits$log_phi[2,3] <- -999 # Set ATF do not feed on SBF
inits$log_phi[2,4] <- -999 # Set ATF do not feed on CSL

inits$log_phi[3,3] <- -999 # Set SBF do not feed on SBF
inits$log_phi[3,2] <- -999 # Set SBF do not feed on ATF
inits$log_phi[3,4] <- -999 # Set SBF do not feed on CSL

inits$log_phi[4,4] <- -999 # Set CSL do not feed on CSL
inits$log_phi[4,2] <- -999 # Set CSL do not feed on ATF
inits$log_phi[4,3] <- -999 # Set CSL do not feed on SBF


# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on ATF
map$mapList$log_phi[1,3] <- NA #so we dont estimate on hake on SBF
map$mapList$log_phi[1,4] <- NA #so we dont estimate on hake on CSL

map$mapList$log_phi[2,2] <- NA #so we dont estimate ATF on ATF
map$mapList$log_phi[2,3] <- NA #so we dont estimate on SBF on SBF
map$mapList$log_phi[2,4] <- NA #so we dont estimate on ATF on CSL

map$mapList$log_phi[3,2] <- NA #so we dont estimate SBF on ATF
map$mapList$log_phi[3,3] <- NA #so we dont estimate SBF on SBF
map$mapList$log_phi[3,4] <- NA #so we dont estimate SBF on CSL

map$mapList$log_phi[4,1] <- NA #so we dont estimate CSL on hake
map$mapList$log_phi[4,2] <- NA #so we dont estimate CSL on ATF
map$mapList$log_phi[4,3] <- NA #so we dont estimate CSL on SBF
map$mapList$log_phi[4,4] <- NA #so we dont estimate CSL on CSL

map$mapFactor$log_phi <- factor(map$mapList$log_phi)

# RUN MODEL
run_ms_CSL_Mest_prior0<- Rceattle::fit_mod(data_list = test_data,
                                          inits = inits, # Initial parameters from single species ests
                                          map = map,
                                          M1Fun = build_M1(M1_model = 1,
                                                           #updateM1 = TRUE,
                                                           M1_use_prior = TRUE,
                                                           M_prior = 0.2,
                                                           M_prior_sd = 0.1),
                                          file = NULL, # Don't save
                                          estimateMode = 0, # estimate
                                          niter = 3, # 3 iterations around population and predation dynamics
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 1, # MSVPA based
                                          loopnum = 5,
                                          phase = TRUE,
                                          suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                          suit_styr = c(1980, 1980, 1980, 1980),  # hake, ATF, sablefish
                                          suit_endyr = c(2019, 2019, 2019, 2019), # hake, ATF, sablefish
                                          initMode = 2,
                                          verbose = 1)

run_ms_CSL_Mest_prior_1<- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_prior0$data_list,
                                            inits = run_ms_CSL_Mest_prior0$estimated_params, # Initial parameters from single species ests
                                            map = run_ms_CSL_Mest_prior0$map,
                                            M1Fun = build_M1(M1_model = 1,
                                                             #updateM1 = TRUE,
                                                             M1_use_prior = TRUE,
                                                             M_prior = 0.2,
                                                             M_prior_sd = 0.1),
                                            file = NULL, # Don't save
                                            estimateMode = 0, # estimate
                                            niter = 3, # 3 iterations around population and predation dynamics
                                            random_rec = FALSE, # No random recruitment
                                            msmMode = 1, # MSVPA based
                                            loopnum = 5,
                                            phase = TRUE,
                                            suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                            suit_styr = c(1980, 1980, 1980, 1980),  # hake, ATF, sablefish
                                            suit_endyr = c(2019, 2019, 2019, 2019), # hake, ATF, sablefish
                                            initMode = 2,
                                            verbose = 1)



run_ms_CSL_Mest_prior<- Rceattle::fit_mod(data_list = test_data,
                                          inits = inits, # Initial parameters from single species ests
                                          map = map,
                                          M1Fun = build_M1(M1_model = 1,
                                                           #updateM1 = TRUE,
                                                           M1_use_prior = TRUE,
                                                           M_prior = 0.2,
                                                           M_prior_sd = 0.1),
                                          file = NULL, # Don't save
                                          estimateMode = 0, # estimate
                                          niter = 3, # 3 iterations around population and predation dynamics
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 1, # MSVPA based
                                          loopnum = 5,
                                          phase = TRUE,
                                          suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                          suit_styr = c(1980, 2013, 2005, 1980),  # hake, ATF, sablefish
                                          suit_endyr = c(2019, 2018, 2008, 2019), # hake, ATF, sablefish
                                          initMode = 2,
                                          verbose = 1)

run_ms_CSL_Mest_prior$quantities$jnll #2394.82/ 13142.93
run_ms_CSL_Mest_prior$quantities$jnll_comp
save(run_ms_CSL_Mest_prior, file = "results/models_all_pred_2024_prior2020/ms_LN_run.Rdata")

### PRIOR DIFFERENCES #######
prior2020<- run_ms_CSL_Mest_prior
load(file = "results/models_all_pred_2024/ms_LN_run_refit_weigth.Rdata")
priorHammelCope<- run_ms_CSL_Mest_prior

mod_list <- list(ss_run, prior2020, priorHammelCope)
mod_names <- c("ss_run", "prior2020", "prior_HammelCope")

mod_list <- list(ss_run, ms_run, run_ms_CSL_Mest_prior)
mod_names <- c("ss_run", "ms_run", "run_ms_CSL_Mest_prior")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names, species =1)
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species =1)
plot_ssb(Rceattle = mod_list, model_names = mod_names, species =1)

#NOTE: We see that when we use the prior on M1 the outputs looks closer to the Stock assessment model

# Compare M1 estimates between models
run_ms_CSL_Mest_prior$quantities$M1
run_ms_CSL_Mest_prior$quantities$vulnerability


plot_diet_comp(run_ms_CSL_Mest_prior) #diet estimate hake canni is good, in this model we do not have ATF predation, hence = 0
plot_diet_comp2(run_ms_CSL_Mest_prior)


## Here we see that ln_M1 is the parameter with higher gradient
gr <- run_ms_CSL_Mest_prior$obj$gr() # Get gradients (unnamed vector)
pars <- run_ms_CSL_Mest_prior$obj$par # Get parameters and names (same order as gradient)
grcheck2 <- data.frame(names = names(pars), par = pars, gr = abs(gr[1,]))

###########################################################################################
## Lets try to estabilize the model and reduce gradient for M1 by re-fitting the model again.
############################################################################################
run_ms_CSL_Mest_test<- run_ms_CSL_Mest_prior_1


run_ms_CSL_Mest_test <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_prior$data_list,
                                          inits = run_ms_CSL_Mest_prior$estimated_params, # Initial parameters from single species ests
                                          map = run_ms_CSL_Mest_prior$map,
                                          M1Fun = build_M1(M1_model = 1,
                                                           #updateM1 = TRUE,
                                                           M1_use_prior = TRUE,
                                                           M_prior = 0.2,
                                                           M_prior_sd = 0.1),
                                          file = NULL, # Don't save
                                          estimateMode = 0, # estimate
                                          niter = 3, # 3 iterations around population and predation dynamics
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 1, # MSVPA based
                                          loopnum = 5,
                                          phase = TRUE,
                                          suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                          suit_styr = c(1980, 1980, 1980, 1980),  # hake, ATF, sablefish
                                          suit_endyr = c(2019, 2019, 2019, 2019), # hake, ATF, sablefish
                                          initMode = 2,
                                          verbose = 1)

run_ms_CSL_Mest_test$quantities$jnll #13142.93
run_ms_CSL_Mest_test$quantities$jnll_comp
save(run_ms_CSL_Mest_test, file = "results/models_all_pred_2024_prior2020/ms_LN_run_refit.Rdata")

gr <- run_ms_CSL_Mest_test$obj$gjnll_compgr <- run_ms_CSL_Mest_test$obj$gr() # Get gradients (unnamed vector)
pars <- run_ms_CSL_Mest_test$obj$par # Get parameters and names (same order as gradient)
grcheck <- data.frame(names = names(pars), par = pars, gr = abs(gr[1,]))

plot_diet_comp(run_ms_CSL_Mest_test)
plot_diet_comp2(run_ms_CSL_Mest_test)
plot_biomass(run_ms_CSL_Mest_test)

run_ms_CSL_Mest_test$data_list$Diet_weights_mcallister


###########################################################################################
## test model weigthing
###########################################################################################
map<- run_ms_CSL_Mest_test$map
map$mapList$diet_comp_weights <- c(NA, 2, 3, NA) # spp 1 fixed, spp 2&3 estimated
map$mapFactor$diet_comp_weights <- factor(map$mapList$diet_comp_weights)

run_ms_CSL_Mest_test_W0 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test$data_list,
                                            inits = run_ms_CSL_Mest_test$estimated_params, # Initial parameters from single species ests
                                            map = map,
                                            M1Fun = build_M1(M1_model = 1,
                                                             #updateM1 = TRUE,
                                                             M1_use_prior = TRUE,
                                                             M_prior = 0.2,
                                                             M_prior_sd = 0.1),
                                            file = NULL, # Don't save
                                            estimateMode = 0, # estimate
                                            niter = 3, # 3 iterations around population and predation dynamics
                                            random_rec = FALSE, # No random recruitment
                                            msmMode = 1, # MSVPA based
                                            loopnum = 5,
                                            phase = TRUE,
                                            suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                            suit_styr = c(1980, 2013, 2005, 1980),  # hake, ATF, sablefish
                                            suit_endyr = c(2019, 2018, 2008, 2019), # hake, ATF, sablefish
                                            initMode = 2,
                                            verbose = 1)
run_ms_CSL_Mest_test_W_new<- run_ms_CSL_Mest_test_W

run_ms_CSL_Mest_test_W$quantities$jnll #13140.69
run_ms_CSL_Mest_test_W$quantities$jnll_comp
save(run_ms_CSL_Mest_test_W, file = "results/models_all_pred_2024_prior2020/ms_LN_run_refit_weigth.Rdata")

load("results/models_all_pred_2024_prior2020/ms_LN_run_refit_weigth.Rdata")

plot_diet_comp(run_ms_CSL_Mest_test_W)
plot_diet_comp2(run_ms_CSL_Mest_test_W)
run_ms_CSL_Mest_test_W$estimated_params$diet_comp_weights # 1.000000 1.242114 2.763260 1.000000
run_ms_CSL_Mest_test_W$data_list$Diet_weights_mcallister #1.37801148 0.08223869 0.31101459 2.07549537
run_ms_CSL_Mest_test_W$sdrep
run_ms_CSL_Mest_test_W$quantities$diet_ESS  #0.00000 46.04381 65.87318

mod_list <- list(ss_run, run_ms_CSL_Mest_test_W)
mod_names <- c("ss model", "ms model")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)

plot_mortality(run_ms_CSL_Mest_test_W0, type = 3, width = 8,
               height = 8), "figures/Manuscript2/M1_M2_plot_final")

#######################################################################
## DIAGNOSTICS
#######################################################################

fit_out <- function(run) {
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient

  fit <- cbind(objective, jnll, K, AIC, gradient)
  return(fit)
}

fit_out(run_ms_CSL_Mest_prior)
fit_out(run_ms_CSL_Mest_test_W)


# PLOT POPULATION DYNAMICS
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, incl_proj = TRUE)
plot_biomass(run_ms_CSL_Mest_test_W, species =2:4) #Now biomass looks alike

plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1, width = 5,
             height = 4,
             lwd = 2,              # slightly thinner looks cleaner at 300 dpi
             mod_cex = 1,        # slightly larger text
             alpha = 0.3, file = "figures/Manuscript1/B_plot")

plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1, width = 5,
                 height = 4,
                 lwd = 2,              # slightly thinner looks cleaner at 300 dpi
                 mod_cex = 1,        # slightly larger text
                 alpha = 0.3, file = "figures/Manuscript1/Rec_plot")

plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1, width = 5,
         height = 4,
         lwd = 2,              # slightly thinner looks cleaner at 300 dpi
         mod_cex = 1,        # slightly larger text
         alpha = 0.3, file = "figures/Manuscript1/SBB_plot")

plot_b_eaten(run_ms_CSL_Mest_test_W, add_ci = TRUE, width = 7,
             height = 7,
             lwd = 2,              # slightly thinner looks cleaner at 300 dpi
             mod_cex = 1,        # slightly larger text
             alpha = 0.3, file = "figures/Manuscript1/Biomass_consumed_plot")

plot_b_eaten_prop(run_ms_CSL_Mest_test_W, add_ci = TRUE, width = 7,
                  height = 7,
                  lwd = 2,              # slightly thinner looks cleaner at 300 dpi
                  file = "figures/Manuscript1/Biomass_consumed_prop_plot")


plot_m_at_age(Rceattle = mod_list, model_names = mod_names, age = 2, species =1, "figures/Manuscript1/M_Age2_plot")
plot_m2_at_age_prop(run_ms_CSL_Mest_test_W, species =1) # Predation M by each species
plot_mortality(run_ms_CSL_Mest_test_W, type = 3, width = 7,
               height = 7, "figures/Manuscript1/M1_M2_plot")

plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species =1)
plot_depletion(Rceattle = mod_list, model_names = mod_names, species =1)

## DIAGNOSTICS
plot_diet_comp2(run_ms_CSL_Mest_test_W, file = "figures/Manuscript1/diagnostics/diet_fit")

plot_catch(Rceattle = mod_list, model_names = mod_names, width = 6,
           height = 5,
           lwd = 2,              # slightly thinner looks cleaner at 300 dpi
           alpha = 0.3, file = "figures/Manuscript1/diagnostics/catch_plot")

plot_index(Rceattle = mod_list, model_names = mod_names, width = 6,
           height = 5, file = "figures/Manuscript1/diagnostics/survey_index_plot")

plot_selectivity(run_ms_CSL_Mest_test_W, width = 6,
                 height = 5, file = "figures/Manuscript1/diagnostics/selectivity_plot")

plot_comp(run_ms_CSL_Mest_test_W, file = "figures/Manuscript1/diagnostics/comps_plot")

plot_comp(run_ms_CSL_Mest_test_W) # Fitted composition data
plot_index(run_ms_CSL_Mest_test_W) # Fitted indices of abundance

plot_logindex(Rceattle = mod_list, model_names = mod_names, width = 6,
              height = 5, file = "figures/Manuscript1/diagnostics/survey_index_plot")  # Fitted log indices of abundance

plot_indexresidual(Rceattle = mod_list, model_names = mod_names, width = 6,
                   height = 5, file = "figures/Manuscript1/diagnostics/survey_index_residual_plot") # Residuals of log indices of abundance

plot_catch(run_ms_CSL_Mest_test_W) # Fitted catch series

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Retrospective ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# we can also do retrospective peels and calculate Mohn's Rho on the CEATTLE
# NOTE: this is using mean historical F for projections as we changed it above
run_ms_CSL_Mest_test_W_retro <- Rceattle::retrospective(Rceattle = run_ms_CSL_Mest_test_W, peels = 5)

# Look at Mohns rho
run_ms_CSL_Mest_test_W_retro$mohns

# Plot retrospectives
plot_biomass(run_ms_CSL_Mest_test_W_retro$Rceattle_list, species = 1, legend.pos = "topright", width = 6,
             height = 5, file = "figures/Manuscript1/diagnostics/retros_plot")

# See how forecast changes
plot_biomass(run_ms_CSL_Mest_test_W_retro$Rceattle_list, incl_proj = TRUE, species = 1, legend.pos = "topright", width = 6,
             height = 5, file = "figures/Manuscript1/diagnostics/retros_forecast_plot")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Jitter ---- Not converging
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
jitters <- Rceattle::jitter(Rceattle = run_ms_CSL_Mest_test_W, njitter = 10, phase = TRUE)
hist(log(jitters$nll - min(jitters$nll)))

# Plot the jitters
plot_biomass(jitters$Rceattle_list)
length(jitters$Rceattle_list) # Some models did not converge, and are not returned


