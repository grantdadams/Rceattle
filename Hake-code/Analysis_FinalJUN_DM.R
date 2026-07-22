######################################################
## Run the age-comps with DM instead of multinomal
######################################################
SBF_ATF_hakedata <- Rceattle::read_data(file = "300426_SBF_ATF_Hake_Final.xlsx")
SBF_ATF_hakedata_DM<- SBF_ATF_hakedata
head(SBF_ATF_hakedata_DM$fleet_control)

## DEFINE DM for compositional data
SBF_ATF_hakedata_DM$fleet_control$Comp_loglike <- c(1,1) # default: DM for age comps

# Manually add diet_ll_type -- length must equal nspp
SBF_ATF_hakedata_DM$diet_ll_type <- rep(1L, SBF_ATF_hakedata_DM$nspp)  # default: DM for diet comps

##Run single species model
ss_run_DM <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata_DM,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


ss_run_DM$quantities$jnll #2128.156
ss_run_DM$quantities$jnll_comp
ss_run_DM$quantities$M1 #Fixed at 0.21

save(ss_run_DM, file = "results/Models_July11/Comps_DM/ss_run_DM.Rdata")

ss_run_M_DM <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata_DM,
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

ss_run_M_DM$quantities$jnll #2183.483
ss_run_M_DM$quantities$jnll_comp
ss_run_M_DM$quantities$M1 #0.2284249
save(ss_run_M_DM, file = "results/Models_July11/Comps_DM/ss_run_M_DM.Rdata")

ss_run_M_free_DM <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata_DM,
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

ss_run_M_free_DM$quantities$jnll #2183.483
ss_run_M_free_DM$quantities$jnll_comp
ss_run_M_free_DM$quantities$M1 #0.2736712

load("results/Models_July11/Comps_DM/run_ms_CSL_Mest_test_W_DM_steady.Rdata")

mod_list <- list(ss_run_DM, ss_run_M_DM, ss_run_M_free_DM, run_ms_CSL_Mest_test_W_DM_steady)
mod_names <- c("ss_fix", "ss_Mprior", "ss_Mfree", "ms_Mprior")

# Plot biomass trajectory
plot_ssb(Rceattle = mod_list, model_names = mod_names, species = 1) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species = 1)


#############################################################################################
#### MULTISPECIES MODELS ###################################################################
## Run multispecies model estimating M1
ms_run_DM <- Rceattle::fit_mod(data_list = SBF_ATF_hakedata_DM,
                            inits = ss_run_DM$estimated_params, # Initial parameters from single species ests
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


ms_run_DM$quantities$jnll #2188
ms_run_DM$quantities$M1 #0.30
save(ms_run_DM, file = "results/Models_July11/Comps_DM/ms_run_DM.Rdata")

#############################################################################################################
## LOGNORMAL SUITABILITY MODEL + DM LIKELIHOOD DISTRIBUTION
#############################################################################################################

# Create initial parameter list:
test_data <- SBF_ATF_hakedata_DM

inits = ms_run_DM$estimated_params
map = ms_run_DM$map # gam_a, gam_b, and log_phi are turned off here

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
run_ms_CSL_Mest_prior_DM<- Rceattle::fit_mod(data_list = test_data,
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

run_ms_CSL_Mest_prior_DM$quantities$jnll #2259.102
run_ms_CSL_Mest_prior_DM$quantities$jnll_comp
run_ms_CSL_Mest_prior_DM$quantities$M1 #0.26
run_ms_CSL_Mest_prior_DM$sdrep
run_ms_CSL_Mest_prior_DM$estimated_params$log_phi
run_ms_CSL_Mest_prior_DM$quantities$vulnerability #0.8450781 and 0.7672554
plot_diet_comp2(run_ms_CSL_Mest_prior_DM)

## Here we see that ln_M1 is the parameter with higher gradient
gr <- run_ms_CSL_Mest_prior$obj$gr() # Get gradients (unnamed vector)
pars <- run_ms_CSL_Mest_prior$obj$par # Get parameters and names (same order as gradient)
grcheck2 <- data.frame(names = names(pars), par = pars, gr = abs(gr[1,]))

###########################################################################################
## Lets try to estabilize the model and reduce gradient for M1 by re-fitting the model again.
############################################################################################
## Diet_comps_Multinomial
run_ms_CSL_Mest_test_DM <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_prior_DM$data_list,
                                          inits = run_ms_CSL_Mest_prior_DM$estimated_params, # Initial parameters from single species ests
                                          map = run_ms_CSL_Mest_prior_DM$map,
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

run_ms_CSL_Mest_test_DM$quantities$jnll #2305.217
run_ms_CSL_Mest_test_DM$quantities$M1 #0.2661498
run_ms_CSL_Mest_test_DM$quantities$jnll_comp
run_ms_CSL_Mest_test_DM$sdrep

save(run_ms_CSL_Mest_test, file = "results/models_final/ms_LN_run_refit.Rdata")

gr <- run_ms_CSL_Mest_test_DM$obj$gr() # Get gradients (unnamed vector)
pars <- run_ms_CSL_Mest_test_DM$obj$par # Get parameters and names (same order as gradient)
grcheck <- data.frame(names = names(pars), par = pars, gr = abs(gr[1,]))


###########################################################################################
## test model weigthing (Diet comps Dirichlet multinomial)
###########################################################################################
#load("results/models_final/ms_LN_run_refit.Rdata")
map<- run_ms_CSL_Mest_test_DM$map
map$mapList$diet_comp_weights <- c(NA, 2, 3) # spp 1 fixed, spp 2&3 estimated
map$mapFactor$diet_comp_weights <- factor(map$mapList$diet_comp_weights)

run_ms_CSL_Mest_test_W_DM <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_DM$data_list,
                                            inits = run_ms_CSL_Mest_test_DM$estimated_params, # Initial parameters from single species ests
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


run_ms_CSL_Mest_test_W_DM$quantities$jnll #2293.756
run_ms_CSL_Mest_test_W_DM$quantities$jnll_comp
run_ms_CSL_Mest_test_W_DM$quantities$M1 # 0.2642167
run_ms_CSL_Mest_test_W_DM$sdrep
run_ms_CSL_Mest_test_W_DM$estimated_params$log_phi
run_ms_CSL_Mest_test_W_DM$quantities$vulnerability #0.8450781 and 0.7672554
run_ms_CSL_Mest_test_W_DM$quantities$unweighted_jnll_comp
plot_diet_comp2(run_ms_CSL_Mest_test_W_DM)
plot_diet_comp2(run_ms_CSL_Mest_test_W)

save(run_ms_CSL_Mest_test_W_DM, file = "results/Models_July11/Comps_DM/ms_LN_run_refit_weigth_DM.Rdata")

run_ms_CSL_Mest_test_W_DM_steady <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W_DM$data_list,
                                               inits = run_ms_CSL_Mest_test_W_DM$estimated_params, # Initial parameters from single species ests
                                               map = run_ms_CSL_Mest_test_W_DM$map,
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

run_ms_CSL_Mest_test_W_DM_steady$quantities$jnll
run_ms_CSL_Mest_test_W_DM_steady$quantities$jnll_comp
run_ms_CSL_Mest_test_W_DM_steady$sdrep
run_ms_CSL_Mest_test_W_DM_steady$quantities$vulnerability
run_ms_CSL_Mest_test_W_DM_steady$estimated_params$log_phi
save(run_ms_CSL_Mest_test_W_DM_steady, file = "results/Models_July11/Comps_DM/run_ms_CSL_Mest_test_W_DM_steady.Rdata")

mod_list <- list(run_ms_CSL_Mest_test_W_DM_steady, run_ms_CSL_Mest_prior)
mod_names <- c("MS_model_Est_DM", "MS_model_Est")

# Plot biomass trajectory
plot_ssb(Rceattle = mod_list, model_names = mod_names, species = 1) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species = 1)

plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species = 1)

plot_m2_at_age_prop(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names, species = 1)

plot_mortality(run_ms_CSL_Mest_test_W_DM_steady, type = 3, width = 8,
               height = 8, "figures/Manuscript1/M1_M2_plot_final2")

## DIAGNOSTICS
plot_diet_comp2(run_ms_CSL_Mest_test_W, file = "figures/Manuscript1/diagnostics/diet_fit")

plot_catch(Rceattle = mod_list, model_names = mod_names, width = 5,
           height = 4,
           lwd = 2,
           alpha = 0.3, file = "figures/Manuscript1/diagnostics/catch_plot3")

plot_index(Rceattle = mod_list, model_names = mod_names, width = 5,
           height = 4, file = "figures/Manuscript1/diagnostics/survey_index_plot3")

plot_selectivity(run_ms_CSL_Mest_test_W_DM_steady, width = 4,
                 height = 4, file = "figures/Manuscript1/diagnostics/selectivity_plot3")

plot_comp(run_ms_CSL_Mest_test_W_DM_steady, cex = 3,
          lwd = 3, file = "figures/Manuscript1/diagnostics/comps_plot3")

plot_comp(run_ms_CSL_Mest_test_W_DM_steady, species = 1) # Fitted composition data
plot_index(run_ms_CSL_Mest_test_W_DM_steady) # Fitted indices of abundance

plot_logindex(Rceattle = mod_list, model_names = mod_names, width = 5,
              height = 3, file = "figures/Manuscript1/diagnostics/survey_index_plot3")  # Fitted log indices of abundance

plot_indexresidual(Rceattle = mod_list, model_names = mod_names, width = 5.5,
                   height = 4, top_adj = -0.8, file = "figures/Manuscript1/diagnostics/survey_index_residual_plot3") # Residuals of log indices of abundance


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Sensitivity ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#### 1) run the model with median gam_a and gam_b values
inits_W = run_ms_CSL_Mest_test_W_DM_steady$estimated_params
map_W = run_ms_CSL_Mest_test_W_DM_steady$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
#Median VALUES
inits_W$log_gam_a = c(0, log(3.35), log(2.94))  # Median
inits_W$log_gam_b = c(0, log(1.83), log(1.120))

run_ms_CSL_Mest_test_W_Sen1 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W_DM_steady$data_list,
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


run_ms_CSL_Mest_test_W_Sen1$quantities$jnll #median: 2260.897
run_ms_CSL_Mest_test_W_Sen1$quantities$jnll_comp
save(run_ms_CSL_Mest_test_W_Sen1, file = "results/Models_July11/Comps_DM/ms_LN_run_refit_weigth_Sen1.Rdata")

##Broader selectivity
inits_W = run_ms_CSL_Mest_test_W_DM_steady$estimated_params
map_W = run_ms_CSL_Mest_test_W_DM_steady$map # gam_a, gam_b, and log_phi are turned off here

inits_W$log_gam_a = c(0, log(3.7), log(3.1)) #mean
inits_W$log_gam_b = c(0, log(1.1), log(0.5)) #broades

run_ms_CSL_Mest_test_W_Sen2 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W_DM_steady$data_list,
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

save(run_ms_CSL_Mest_test_W_Sen2, file = "results/Models_July11/Comps_DM/ms_LN_run_refit_weigth_Sen2.Rdata")

mod_list <- list(run_ms_CSL_Mest_test_W_DM_steady, run_ms_CSL_Mest_test_W_Sen1, run_ms_CSL_Mest_test_W_Sen2)
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
run_ms_CSL_Mest_test_W_DM_steady$data_list$other_food<- 2 * run_ms_CSL_Mest_test_W_DM_steady$data_list$other_food

run_ms_CSL_Mest_test_W_Sens2 <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_test_W_DM_steady$data_list,
                                                  inits = run_ms_CSL_Mest_test_W_DM_steady$estimated_params, # Initial parameters from single species ests
                                                  map = run_ms_CSL_Mest_test_W_DM_steady$map,
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
save(run_ms_CSL_Mest_test_W_Sens2, file = "results/Models_July11/Comps_DM/ms_LN_run_refit_weigth_Sen3_Bother.Rdata")

mod_list <- list(run_ms_CSL_Mest_test_W_DM_steady, run_ms_CSL_Mest_test_W_Sens2)
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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Retrospective ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# we can also do retrospective peels and calculate Mohn's Rho on the CEATTLE
# NOTE: this is using mean historical F for projections as we changed it above
run_ms_CSL_Mest_test_W_retro <- retrospective(Rceattle = run_ms_CSL_Mest_test_W_DM_steady, peels = 5)

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
             file = "figures/Manuscript1/diagnostics/retros_plot_FINAL2")
