library(Rceattle)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)


#############################################################################################################
## LOAD DATA
#############################################################################################################

SBF_ATF_hakedata <- Rceattle::read_data(file = "300426_SBF_ATF_Hake_Final.xlsx")

## DEFINE DM for compositional data
SBF_ATF_hakedata$fleet_control$Comp_loglike <- c(0,0) # default: DM for age comps

# Manually add diet_ll_type -- length must equal nspp
SBF_ATF_hakedata$diet_ll_type <- rep(0L, SBF_ATF_hakedata$nspp)  # default: DM for diet comps


run_ms_CSL_Mest_test_W$data_list$diet_ll_type
run_ms_CSL_Mest_test_W$data_list$diet_ll_type

#############################################################################################################
## RUN MODELS - single spss and multispss MSVPA (hake cannibalims only)
#############################################################################################################

##Single spss
ss_run <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


ss_run$quantities$jnll #2164.514
ss_run$quantities$jnll_comp
ss_run$quantities$M1 #Fixed at 0.21
save(ss_run, file = "results/Models_July11/Comps_multinomial/ss_run.Rdata")

ss_run_M <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              M1Fun = build_M1(M1_model = 1,
                                               #updateM1 = TRUE,
                                               M1_use_prior = TRUE,
                                               M_prior = 0.2,
                                               M_prior_sd = 0.1),
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

ss_run_M$quantities$jnll #2214.202
ss_run_M$quantities$jnll_comp
ss_run_M$quantities$M1 #estimated: 0.2276400
save(ss_run_M, file = "results/Models_July11/Comps_multinomial/ss_run_M.Rdata")

ss_run_M_noPrior <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              M1Fun = build_M1(M1_model = 1,
                                               #updateM1 = TRUE,
                                               M1_use_prior = FALSE),
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

ss_run_M_noPrior$quantities$jnll #2162.233
ss_run_M_noPrior$quantities$jnll_comp
ss_run_M_noPrior$quantities$M1 #estimated: 0.268

inits_warm <- ss_run_M$estimated_params
ss_M1_re <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
  inits = inits_warm,
  map = NULL,
  estimateMode = 0,
  random_rec = FALSE,        # keep consistent with your other runs
  M1Fun = build_M1(
    M1_model = 1,             # species-specific M1
    M1_re = 2,                # AR1 random effect varying by year, constant across ages
    M1_use_prior = TRUE,      # keep the same informative prior as your other estimated-M1 runs
    M_prior = 0.2,
    M_prior_sd = 0.1
  ),
  msmMode = 0,                # single-species mode
  phase = TRUE,               # recommended for random effects models, improves convergence
  verbose = 1,
  loopnum = 5
)

ss_M1_re$quantities$jnll

mod_list <- list(ss_run, run_ms_CSL_Mest_test_W, ss_run_M, ss_run_M_noPrior)
mod_names <- c("ss_model_Fix", "MS_model_Est", "ss_model_Est", "ss_run_M_noPrior")

# Plot biomass trajectory
plot_ssb(Rceattle = mod_list, model_names = mod_names, species = 1) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species = 1)

plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species = 1)

#############################################################################################
#### MULTISPECIES MODELS ###################################################################
## Run multispecies model estimating M1
ms_run <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata,
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
                            suitMode = 0, # empirical suitability
                            suit_styr = c(1980, 1980, 1980),  # hake, ATF, sablefish
                            suit_endyr = c(2019, 2019, 2019), # hake, ATF, sablefish
                            initMode = 2, # Fished start with init devs
                            verbose = 1)


ms_run$quantities$jnll #2215.4
ms_run$quantities$M1 #0.2976292
save(ms_run, file = "results/models_final/ms_run.Rdata")

#############################################################################################################
## LOGNORMAL SUITABILITY MODEL + DM LIKELIHOOD DISTRIBUTION
#############################################################################################################

# Create initial parameter list:
test_data <- SBF_ATF_hakedata

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
#ORIGINAL (TRUE) VALUES
inits$log_gam_a = c(0, log(3.7), log(3.1))  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, log(1.83), log(1.120))

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

# RUN MODEL
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
                                            suitMode = c(0, 4, 4), # empirical + LN suitability
                                            suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                            suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                            initMode = 2,
                                            verbose = 1)

run_ms_CSL_Mest_prior$quantities$jnll_comp #2305.217
run_ms_CSL_Mest_prior$quantities$M1 #0.2661498
run_ms_CSL_Mest_prior$sdrep
run_ms_CSL_Mest_prior$estimated_params$log_phi
run_ms_CSL_Mest_prior$quantities$vulnerability #0.8450781 and 0.7672554
plot_diet_comp2(run_ms_CSL_Mest_prior)

## Here we see that ln_M1 is the parameter with higher gradient
gr <- run_ms_CSL_Mest_prior$obj$gr() # Get gradients (unnamed vector)
pars <- run_ms_CSL_Mest_prior$obj$par # Get parameters and names (same order as gradient)
grcheck2 <- data.frame(names = names(pars), par = pars, gr = abs(gr[1,]))

###########################################################################################
## Lets try to estabilize the model and reduce gradient for M1 by re-fitting the model again.
############################################################################################
## Diet_comps_Multinomial
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
                                          suitMode = c(0, 4, 4), # empirical + LN suitability
                                          suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                          suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                          initMode = 2,
                                          verbose = 1)

run_ms_CSL_Mest_test$quantities$jnll #2305.217
run_ms_CSL_Mest_test$quantities$M1 #0.2661498
run_ms_CSL_Mest_test$quantities$jnll_comp
run_ms_CSL_Mest_test$sdrep
plot_diet_comp2(run_ms_CSL_Mest_test)

save(run_ms_CSL_Mest_test, file = "results/models_final/ms_LN_run_refit.Rdata")

gr <- run_ms_CSL_Mest_test$obj$gr() # Get gradients (unnamed vector)
pars <- run_ms_CSL_Mest_test$obj$par # Get parameters and names (same order as gradient)
grcheck <- data.frame(names = names(pars), par = pars, gr = abs(gr[1,]))

save(run_ms_CSL_Mest_test, file = "results/Models_July11/Comps_multinomial/run_ms_CSL_Mest_test.Rdata")
load("results/Models_July11/Comps_multinomial/run_ms_CSL_Mest_test.Rdata")

###########################################################################################
## test model weigthing (Diet comps Dirichlet multinomial)
###########################################################################################
#load("results/models_final/ms_LN_run_refit.Rdata")
map<- run_ms_CSL_Mest_test$map
map$mapList$diet_comp_weights <- c(NA, 2, 3) # spp 1 fixed, spp 2&3 estimated
map$mapFactor$diet_comp_weights <- factor(map$mapList$diet_comp_weights)

run_ms_CSL_Mest_test_W <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test$data_list,
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
                                            suitMode = c(0, 4, 4), # empirical + LN suitability
                                            suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                            suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                            initMode = 2,
                                            verbose = 1)


run_ms_CSL_Mest_test_W$quantities$jnll #2293.756
run_ms_CSL_Mest_test_W$quantities$jnll_comp
run_ms_CSL_Mest_test_W$quantities$M1 # 0.2642167
run_ms_CSL_Mest_test_W$sdrep
run_ms_CSL_Mest_test_W$estimated_params$log_phi
run_ms_CSL_Mest_test_W$quantities$vulnerability #0.8450781 and 0.7672554
run_ms_CSL_Mest_test_W$quantities$unweighted_jnll_comp
save(run_ms_CSL_Mest_test_W, file = "results/models_final/ms_LN_run_refit_weigth_test.Rdata")


mod_list <- list(ss_run, run_ms_CSL_Mest_test_W, ss_run_M, ss_run_M_noPrior)
mod_names <- c("ss_model_Fix", "MS_model_Est", "ss_model_Est", "ss_run_M_noPrior")

# Plot biomass trajectory
plot_ssb(Rceattle = mod_list, model_names = mod_names, species = 1) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species = 1)

plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species = 1)

###########################################################################################
## test diet multinomial instead of DM
###########################################################################################
Data<- run_ms_CSL_Mest_test$data_list
Data$Diet_comp_weights<- c(20, 20, 20)

inits_large_theta <- run_ms_CSL_Mest_test$estimated_params
print(inits_large_theta$diet_comp_weights)  # should show whatever it currently is

# Overwrite
inits_large_theta$diet_comp_weights <- c(20, 20, 20)
# Confirm the overwrite actually took
print(inits_large_theta$diet_comp_weights)

run_ms_diet_multinomial_equiv <- Rceattle::fit_mod(data_list = Data,
                                            inits = inits_large_theta, # Initial parameters from single species ests
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
                                            suitMode = c(0, 4, 4), # empirical + LN suitability
                                            suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                            suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                            initMode = 2,
                                            verbose = 1)

run_ms_diet_multinomial_equiv$quantities$jnll
run_ms_diet_multinomial_equiv$quantities$jnll_comp
run_ms_diet_multinomial_equiv$quantities$M1
run_ms_diet_multinomial_equiv$sdrep
run_ms_diet_multinomial_equiv$quantities$unweighted_jnll_comp
run_ms_diet_multinomial_equiv$initial_params$diet_comp_weights

mod_list <- list(run_ms_CSL_Mest_test_W, run_ms_CSL_Mest_test, run_ms_diet_multinomial_equiv)
mod_names <- c("model_DM", "model_1", "model_M20")

# Plot biomass trajectory
plot_ssb(Rceattle = mod_list, model_names = mod_names, species = 1) #Now biomass looks alike
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species = 1)

################################################################################################
load("results/models_final/ms_LN_run_refit_weigth_test.Rdata")
##############################################################################################
run_ms_CSL_Mest_0<- Rceattle::fit_mod(data_list = test_data,
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
                                          suitMode = c(0, 4, 4), # empirical + LN suitability
                                          suit_styr = c(1980, 1980, 1980),  # hake, ATF, sablefish
                                          suit_endyr = c(2019, 2019, 2019), # hake, ATF, sablefish
                                          initMode = 2,
                                          verbose = 1)

run_ms_CSL_Mest_1<- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_0$data_list,
                                      inits = run_ms_CSL_Mest_0$estimated_params, # Initial parameters from single species ests
                                      map = run_ms_CSL_Mest_0$map,
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
                                      suitMode = c(0, 4, 4), # empirical + LN suitability
                                      suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                      suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                      initMode = 2,
                                      verbose = 1)

run_ms_CSL_Mest_1$quantities$jnll #2301.866
run_ms_CSL_Mest_1$quantities$jnll_comp
run_ms_CSL_Mest_1$quantities$M1 # 0.2617886
run_ms_CSL_Mest_1$sdrep

map<- run_ms_CSL_Mest_1$map
map$mapList$diet_comp_weights <- c(NA, 2, 3) # spp 1 fixed, spp 2&3 estimated
map$mapFactor$diet_comp_weights <- factor(map$mapList$diet_comp_weights)

run_ms_CSL_Mest_test_W0 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_1$data_list,
                                            inits = run_ms_CSL_Mest_1$estimated_params, # Initial parameters from single species ests
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
                                            suitMode = c(0, 4, 4), # empirical + LN suitability
                                            suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                            suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                            initMode = 2,
                                            verbose = 1)

run_ms_CSL_Mest_test_W0$quantities$jnll # 2290.246
run_ms_CSL_Mest_test_W0$quantities$jnll_comp
run_ms_CSL_Mest_test_W0$quantities$M1 # 0.2603446
run_ms_CSL_Mest_test_W0$sdrep
run_ms_CSL_Mest_test_W0$estimated_params$log_phi
run_ms_CSL_Mest_test_W0$quantities$vulnerability #0.8450781 and 0.7672554
save(run_ms_CSL_Mest_test_W, file = "results/models_final/ms_LN_run_refit_weigth_test_OLD_GAMMA.Rdata")

### EXPLORE OUTPUTS ########################################################

load("results/models_final/ms_LN_run_refit_weigth_test.Rdata")
load("results/models_final/ss_run.Rdata")

run_ms_CSL_Mest_test_W$data_list$suit_styr
run_ms_CSL_Mest_test_W$data_list$suit_endyr

plot_diet_comp2(run_ms_CSL_Mest_test_W)
plot_diet_comp1(run_ms_CSL_Mest_test_W)
run_ms_CSL_Mest_test_W$estimated_params$diet_comp_weights # 1.0000000 0.8874628 2.7702246
run_ms_CSL_Mest_test_W$data_list$Diet_weights_mcallister #1.35885650 0.07358298 0.30996314
run_ms_CSL_Mest_test_W$sdrep
run_ms_CSL_Mest_test_W$estimated_params$log_phi
run_ms_CSL_Mest_test_W$quantities$diet_ESS  #0.00000 46.04381 65.87318

mod_list <- list(ss_run, run_ms_CSL_Mest_test_W)
mod_names <- c("ss model", "ms model")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1)

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)

#### get Neff ###
# Get all unique sample sizes from diet data
diet <- run_ms_CSL_Mest_test_W$data_list$diet_data
weights <- run_ms_CSL_Mest_test_W$estimated_params$diet_comp_weights

# For ATF (Pred = 2)
ATF_data <- diet %>%
  filter(Pred == 2) %>%
  select(Pred_age, Sample_size) %>%
  distinct()

theta_ATF <- exp(weights[2])
ATF_data$ESS <- ATF_data$Sample_size * (theta_ATF / (1 + theta_ATF))

# For SBF (Pred = 3)
SBF_data <- diet %>%
  filter(Pred == 3) %>%
  select(Pred_age, Sample_size) %>%
  distinct()

theta_SBF <- exp(weights[3])
SBF_data$ESS <- SBF_data$Sample_size * (theta_SBF / (1 + theta_SBF))

print(ATF_data)
print(SBF_data)

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

fit_out(run_ms_CSL_Mest_test)
fit_out(run_ms_CSL_Mest_test_W_DM_steady)

run_ms_CSL_Mest_test$quantities$M1#0.2638791
run_ms_CSL_Mest_test_W_DM_steady$quantities$M1 #0.2631309

# PLOT POPULATION DYNAMICS
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, incl_proj = TRUE)
plot_biomass(Rceattle = mod_list, model_names = mod_names, species =1) #Now biomass looks alike

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
plot_mortality(run_ms_CSL_Mest_test_W, type = 3, width = 8,
               height = 8, "figures/Manuscript1/M1_M2_plot_final")

plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species =1)
plot_depletion(Rceattle = mod_list, model_names = mod_names, species =1)

## DIAGNOSTICS
plot_diet_comp2(run_ms_CSL_Mest_test_W, file = "figures/Manuscript1/diagnostics/diet_fit")

plot_catch(Rceattle = mod_list, model_names = mod_names, width = 5,
           height = 3,
           lwd = 2,
           alpha = 0.3), file = "figures/Manuscript1/diagnostics/catch_plot")

plot_index(Rceattle = mod_list, model_names = mod_names, width = 5,
           height = 3), file = "figures/Manuscript1/diagnostics/survey_index_plot")

plot_selectivity(run_ms_CSL_Mest_test_W, width = 4,
                 height = 4, file = "figures/Manuscript1/diagnostics/selectivity_plot")

plot_comp(run_ms_CSL_Mest_test_W, cex = 3,
          lwd = 3), file = "figures/Manuscript1/diagnostics/comps_plot")

plot_comp(run_ms_CSL_Mest_test_W, species = 1) # Fitted composition data
plot_index(run_ms_CSL_Mest_test_W) # Fitted indices of abundance

plot_logindex(Rceattle = mod_list, model_names = mod_names, width = 5,
              height = 3, file = "figures/Manuscript1/diagnostics/survey_index_plot")  # Fitted log indices of abundance

plot_indexresidual(Rceattle = mod_list, model_names = mod_names, width = 5,
                   height = 3, top_adj = -0.8), file = "figures/Manuscript1/diagnostics/survey_index_residual_plot") # Residuals of log indices of abundance

plot_catch(run_ms_CSL_Mest_test_W) # Fitted catch series
plot_diet_comp2(run_ms_CSL_Mest_test_W)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Retrospective ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# we can also do retrospective peels and calculate Mohn's Rho on the CEATTLE
# NOTE: this is using mean historical F for projections as we changed it above
run_ms_CSL_Mest_test_W_retro <- retrospective(Rceattle = run_ms_CSL_Mest_test_W, peels = 5)

# Look at Mohns rho
run_ms_CSL_Mest_test_W_retro$mohns

retro_list <- run_ms_CSL_Mest_test_W_retro$Rceattle_list

# Overwrite endyr with endyr_peel for each peel so plot_timeseries clips correctly
retro_list_fixed <- lapply(retro_list, function(x) {
  if (!is.null(x$data_list$endyr_peel)) {
    x$data_list$endyr <- x$data_list$endyr_peel
  }
  x
})

plot_biomass(retro_list_fixed, species = 1, legend.pos = "topright",
             width = 5, height = 3,
             file = "figures/Manuscript1/diagnostics/retros_plot_FINAL")

# Plot retrospectives
plot_biomass(run_ms_CSL_Mest_test_W_retro$Rceattle_list, species = 1, legend.pos = "topright"), width = 5,
height = 3, file = "figures/Manuscript1/diagnostics/retros_plot")

# See how forecast changes
plot_biomass(run_ms_CSL_Mest_test_W_retro$Rceattle_list, incl_proj = TRUE, species = 1, legend.pos = "topright", width = 5,
             height = 3, file = "figures/Manuscript1/diagnostics/retros_forecast_plot")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Jitter ---- Not converging
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
jitters <- Rceattle::jitter(Rceattle = run_ms_CSL_Mest_test_W, njitter = 10, phase = TRUE)
hist(log(jitters$nll - min(jitters$nll)))

# Plot the jitters
plot_biomass(jitters$Rceattle_list)
length(jitters$Rceattle_list) # Some models did not converge, and are not returned


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Sensitivity ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#### 1) run the model with median gam_a and gam_b values
inits_W = run_ms_CSL_Mest_test_W$estimated_params
map_W = run_ms_CSL_Mest_test_W$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
#Median VALUES
inits_W$log_gam_a = c(0, log(3.35), log(2.94))  # Median
inits_W$log_gam_b = c(0, log(1.83), log(1.120))

run_ms_CSL_Mest_test_W_Sen1 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W$data_list,
                                                 inits = inits_W, # Initial parameters from single species ests
                                                 map = map_W,
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
                                                 suitMode = c(0, 4, 4), # empirical + LN suitability
                                                 suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                                 suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                                 initMode = 2,
                                                 verbose = 1)


run_ms_CSL_Mest_test_W_Sen1$quantities$jnll #median: 2368.572
run_ms_CSL_Mest_test_W_Sen1$quantities$jnll_comp
save(run_ms_CSL_Mest_test_W_Sen1, file = "results/models_final/ms_LN_run_refit_weigth_Sen1.Rdata")

##Broader selectivity
inits_W = run_ms_CSL_Mest_test_W$estimated_params
map_W = run_ms_CSL_Mest_test_W$map # gam_a, gam_b, and log_phi are turned off here

inits_W$log_gam_a = c(0, log(3.7), log(3.1)) #mean
inits_W$log_gam_b = c(0, log(1.1), log(0.5)) #broades

run_ms_CSL_Mest_test_W_Sen2 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W$data_list,
                                                 inits = inits_W, # Initial parameters from single species ests
                                                 map = map_W,
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
                                                 suitMode = c(0, 4, 4), # empirical + LN suitability
                                                 suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                                 suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                                 initMode = 2,
                                                 verbose = 1)

save(run_ms_CSL_Mest_test_W_Sen2, file = "results/models_final/ms_LN_run_refit_weigth_Sen2.Rdata")

mod_list <- list(run_ms_CSL_Mest_test_W, run_ms_CSL_Mest_test_W_Sen1, run_ms_CSL_Mest_test_W_Sen2)
mod_names <- c("ms model", "sen1", "sen2")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1)

plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species =1, width = 5,
                 height = 4,
                 lwd = 2,              # slightly thinner looks cleaner at 300 dpi
                 mod_cex = 1,        # slightly larger text
                 alpha = 0.3,
                 file = "figures/Manuscript1/sensitivity/Rec_plot_combined")

plot_diet_comp1(run_ms_CSL_Mest_test_W_Sen1, file = "figures/Manuscript1/sensitivity/diet_fit_gam_b_broader")

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, width = 7,
                  height = 7,
                  lwd = 2)

#### 2) Run with different values of "other food"
run_ms_CSL_Mest_test_W$data_list$other_food<- 2 * run_ms_CSL_Mest_test_W$data_list$other_food

run_ms_CSL_Mest_test_W_Sens2 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W$data_list,
                                                  inits = run_ms_CSL_Mest_test_W$estimated_params, # Initial parameters from single species ests
                                                  map = run_ms_CSL_Mest_test_W$map,
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
                                                  suitMode = c(0, 4, 4), # empirical + LN suitability
                                                  suit_styr = c(1980, 2013, 2005),  # hake, ATF, sablefish
                                                  suit_endyr = c(2019, 2018, 2008), # hake, ATF, sablefish
                                                  initMode = 2,
                                                  verbose = 1)

run_ms_CSL_Mest_test_W_Sens2$quantities$jnll #2362.957
run_ms_CSL_Mest_test_W_Sens2$quantities$jnll_comp

mod_list <- list(run_ms_CSL_Mest_test_W, run_ms_CSL_Mest_test_W_Sens2)
mod_names <- c("ms model", "ms model (x2 other food")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)

plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE), width = 7,
height = 7,
lwd = 2, file = "figures/Manuscript1/sensitivity/Bio_consumed_Other_foodx2")

plot_m_at_age(Rceattle = mod_list, model_names = mod_names, age =6 , species =1, file= "figures/Manuscript1/sensitivity/M_Age6_broad_gamb_plot")
plot_m2_at_age_prop(Rceattle = mod_list, model_names = mod_names, species =1, age = 1)# file= "figures/Manuscript1/sensitivity/M2_Age6_broad_gamb_plot") # Predation M by each species

