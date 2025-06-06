library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data("BS2017SS") # Single-species data. ?BS2017SS for more information on the data
data("BS2017MS") # Multi-species data. Note: the only difference is the residual mortality (M1_base) is lower

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")


################################################
# Estimation
################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


# Single-species, but estimate M
ss_run_M <- Rceattle::fit_mod(data_list = mydata,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              M1Fun = build_M1(M1_model = 1,
                                               M1_use_prior = FALSE,
                                               M2_use_prior = FALSE),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

plot_biomass(ss_run_M, add_ci = TRUE)


# - Multi-species
# For the a multispecies model we from the single species parameters.
ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            M1Fun = build_M1(M1_model = 1,
                                             updateM1 = TRUE,
                                             M1_use_prior = FALSE,
                                             M2_use_prior = FALSE),
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

plot_biomass(ms_run, add_ci = TRUE)

################################################
# Plotting
################################################
# We can plot all runs
mod_list <- list(ss_run, ss_run_M, ms_run)
mod_names <- c("Single-species", "Single-species estimate M", "Multi-species")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

#####################################################################
# FIX SUITABILITY AND SUM across prey ages (NEW PART)
# Create initial parameter list:
test_data<-BS2017MS

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here
ms_run$map$mapList$log_gam_a

# Create a list prey size preference
# Set weight ratio parameters
inits$log_gam_a = c(-6.61, -3.91, -3.46)  # Mean log weight ratio
inits$log_gam_b = c(1.99, 1.69, 1.89)  # Standard deviation of log weight ratio f
test_data$diet_data$Prey_age <- -999 # Set prey age to negative values to sum across ages

a<- test_data$diet_data
head(a)

test_data_sum<- test_data$diet_data %>%
  group_by(Pred, Prey, Pred_sex, Prey_sex, Pred_age, Prey_age, Year, Sample_size) %>%
  summarise(Stomach_proportion_by_weight = sum(Stomach_proportion_by_weight))
head(test_data_sum)

test_data$diet_data<- test_data_sum #plug in the new data

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)

# Do this to estimate vulnerability (log_phi) :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapFactor$log_phi <- factor(map$mapList$log_phi)

#map$mapList$log_phi <- inits$log_phi # If we set this, the log_phi is not gonna be estimated. Should I do this or not?
#map$mapFactor$log_phi <- factor(map$mapList$log_phi) #If I do the above, then I should also do this.

#I did both and they gave different results. If we do not set #map$mapList$log_phi <- inits$log_phi , the estimated log_phi = 1.649402, and it shouldn't be above 1, right?

run_ms_LN <- Rceattle::fit_mod(data_list = test_data,
                               inits = inits, # Initial parameters from single species ests
                               map = map,
                               M1Fun = build_M1(M1_model = 1,
                                                updateM1 = TRUE,
                                                M1_use_prior = FALSE,
                                                M2_use_prior = FALSE),
                               file = NULL, # Don't save
                               estimateMode = 0, # Estimate
                               niter = 3, # 3 iterations around population and predation dynamics
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # MSVPA based
                               suitMode = 4, # empirical suitability
                               verbose = 1)



run_ms_LN$quantities$jnll_comp
run_ms_LN$data_list$Diet_weights_mcallister # Three values
run_ms_LN$quantities$vulnerability #Now it is 0.31 for all spss

## MODEL RE-Weigthing
# Access the calculated weights
run_ms_LN$data_list$Diet_comp_weights
mcallister_weights <- run_ms_LN$data_list$Diet_weights_mcallister

# Use the calculated weights for re-running
test_data$Diet_comp_weights <- run_ms_LN$data_list$Diet_weights_mcallister
inits$diet_comp_weights<- run_ms_LN$data_list$Diet_weights_mcallister #If I don't reset inits, the model give me the same mcallister weight, so I guess I have to do this so we initiate form previous mcallister weight instead of 1?


# Re-run the model with updated weights
run_ms_LN_rw <- Rceattle::fit_mod(data_list = test_data,
                                  inits = inits, # Initial parameters from multispecies ests
                                  map = map,
                                  M1Fun = build_M1(M1_model = 1,
                                                   updateM1 = TRUE,
                                                   M1_use_prior = FALSE,
                                                   M2_use_prior = FALSE),
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  niter = 3, # 3 iterations around population and predation dynamics
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 1, # MSVPA based
                                  suitMode = 4, # empirical suitability
                                  verbose = 1)
# JNLL
run_ms_LN_rw$quantities$jnll_comp ##Nothing changes
run_ms_LN_rw$quantities$vulnerability
run_ms_LN_rw$data_list$Diet_weights_mcallister #nothing changes.
run_ms_LN_rw$data_list$Diet_comp_weights #same as mcallister

# Plot
plot_biomass(list(run_ms_LN, run_ms_LN_rw), model_names = c("run_ms_LN", "run_ms_LN_rw")) ##Nothing changes
plot_b_eaten(list(run_ms_LN, run_ms_LN_rw), model_names = c("run_ms_LN", "run_ms_LN_rw")) ##Nothing changes

# We can plot all runs
mod_list <- list(ss_run, ss_run_M, ms_run, run_ms_LN_rw)
mod_names <- c("Single-species", "Single-species estimate M", "Multi-species", "run_ms_LN_rw")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

# Plot mortality and predation
plot_b_eaten(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey by each predator


# Compare negative log-likelihoods
print("Original model fit:")
print(ms_run$quantities$jnll_comp)
print(ms_run$quantities$jnll) #10116.54

print("LN model fit:")
print(run_ms_LN$quantities$jnll_comp)
print(run_ms_LN$quantities$jnll) #12308.74

print("Reweighted model fit:")
print(run_ms_LN_rw$quantities$jnll_comp)
print(run_ms_LN_rw$quantities$jnll) #10160.9

## MODEL RE-Weigthing 2 ##################################################
# Access the calculated weights
run_ms_LN_rw$data_list$Diet_comp_weights
mcallister_weights <- run_ms_LN_rw$data_list$Diet_weights_mcallister

# Use the calculated weights for re-running
test_data$Diet_comp_weights <- run_ms_LN_rw$data_list$Diet_weights_mcallister #Increase for testing effect
inits$diet_comp_weights<- run_ms_LN_rw$data_list$Diet_weights_mcallister


# Re-run the model with updated weights
run_ms_LN_rw2 <- Rceattle::fit_mod(data_list = test_data,
                                   inits = inits, # Initial parameters from multispecies ests
                                   map = map,
                                   M1Fun = build_M1(M1_model = 1,
                                                    updateM1 = TRUE,
                                                    M1_use_prior = FALSE,
                                                    M2_use_prior = FALSE),
                                   file = NULL, # Don't save
                                   estimateMode = 0, # Estimate
                                   niter = 3, # 3 iterations around population and predation dynamics
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 1, # MSVPA based
                                   suitMode = 4, # empirical suitability
                                   verbose = 1)
# JNLL
run_ms_LN_rw2$quantities$jnll_comp ##Nothing changes
run_ms_LN_rw2$quantities$vulnerability
run_ms_LN_rw2$data_list$Diet_weights_mcallister #nothing changes.
run_ms_LN_rw2$data_list$Diet_comp_weights #same as mcallister

# Plot
mod_list <- list(ms_run, run_ms_LN, run_ms_LN_rw, run_ms_LN_rw2)
mod_names <- c("Multi-species", "Multi-species_LN", "run_ms_LN_rw", "run_ms_LN_rw2")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)

# Plot mortality and predation
plot_b_eaten(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey by each predator

run_ms_LN$estimated_params$log_phi
run_ms_LN_rw$estimated_params$log_phi
