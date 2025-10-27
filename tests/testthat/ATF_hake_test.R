library(Rceattle)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)

#access toke: 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

ATF_hakedata <- Rceattle::read_data(file = "070725_ATF_Hake_model_Neff.xlsx")

ss_run <- Rceattle::fit_mod(data_list = ATF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run$quantities$jnll #2073.326
#save(ss_run, file = "results/models/ATF/ATF_ss_run.Rdata")

ms_run <- Rceattle::fit_mod(data_list = ATF_hakedata,
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

ms_run$quantities$jnll #2085.478
save(ms_run, file = "results/models/ATF/ATF_ms_run.Rdata")

plot_diet_comp(ms_run) #diet estimate hake canni is good, in this model we do not have ATF predation, hence = 0

# Correct way to access the data sent to TMB
ms_run$obj$env$data$n_stomach_obs
head(ms_run$obj$env$data$stomach_id)

#####################################################################
# FIX SUITABILITY AND SUM across prey ages (NEW PART)
# Create initial parameter list:
test_data <- ATF_hakedata

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
inits$log_gam_a = c(0, 3.7)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, 1.8)

#ORIGINAL (TRUE) VALUES
#inits$log_gam_a = c(0, 3.3)  # Median log weight ratio for ATF, 0 for other species (pred/prey)
#inits$log_gam_b = c(0, 1.8)

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on ATF
inits$log_phi[2,2] <- -999 # Set ATF do not feed on ATF
#inits$log_phi[2,1] <- 4

# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on atf
map$mapList$log_phi[2,2] <- NA #so we dont estimate atf on atf

map$mapFactor$log_phi <- factor(map$mapList$log_phi)

#In this model we are debugging estimateMode = 3
run_ms_LN_0 <- Rceattle::fit_mod(data_list = test_data,
                               inits = inits, # Initial parameters from single species ests
                               map = map,
                               M1Fun = build_M1(M1_model = 0,
                                                updateM1 = FALSE,
                                                M1_use_prior = FALSE,
                                                M2_use_prior = FALSE),
                               file = NULL, # Don't save
                               estimateMode = 3, # debug
                               niter = 3, # 3 iterations around population and predation dynamics
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # MSVPA based
                               loopnum = 5,
                               phase = TRUE,
                               suitMode = c(0, 4), # empirical + LN suitability
                               initMode = 2,
                               verbose = 1)

run_ms_LN_0$quantities$jnll_comp
run_ms_LN_0$quantities$diet_prop_hat
run_ms_LN_0$quantities$jnll
run_ms_LN_0$quantities$jnll_comp
run_ms_LN_0$estimated_params$log_phi
run_ms_LN_0$quantities$vulnerability
run_ms_LN_0$estimated_params$diet_comp_weights
run_ms_LN_0$data_list$Diet_weights_mcallister

# Correct way to access the data sent to TMB
run_ms_LN_0$obj$env$data$n_stomach_obs
head(run_ms_LN_0$obj$env$data$stomach_id)

plot_diet_comp(run_ms_LN_0) #diet estimate hake canni is goog but ATF on hake close to 0

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

run_ms_LN_1$quantities$jnll_comp
run_ms_LN_1$estimated_params$diet_comp_weights
run_ms_LN_1$data_list$Diet_weights_mcallister

plot_diet_comp2(run_ms_LN_1)
plot_b_eaten_prop(run_ms_LN_1)

#In this model we are estimating estimateMode = 0 and NO est M1 but updating it
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

#In this model we are estimating estimateMode = 0 and M1
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

run_ms_LN_3$quantities$jnll #2370.025
run_ms_LN_3$quantities$jnll_comp
run_ms_LN_4$quantities$jnll_comp
save(run_ms_LN_3, file = "results/models/ATF/ATF_run_ms_LN_estM.Rdata")

load("results/models/ATF/ATF_run_ms_LN_estM.Rdata")

run_ms_LN_2$quantities$diet_hat
run_ms_LN_2$quantities$diet_prop_hat
run_ms_LN_2$quantities$jnll
run_ms_LN_1$estimated_params$log_phi
run_ms_LN_3$quantities$vulnerability
run_ms_LN_2$estimated_params$diet_comp_weights
run_ms_LN_1$data_list$Diet_weights_mcallister

plot_diet_comp(run_ms_LN_3)
plot_diet_comp2(run_ms_LN_4)
plot_b_eaten_prop(run_ms_LN_3)


models <- list(ss_run, ms_run,  run_ms_LN_3)
model_names <- c("ss_run", "Cannib", "run_ms_LN_estM")

# Predator: 1=Hake_F, 2=ATF_F, 3=Hake_M, 4=ATF_M
# Prey:     1=Hake_F, 2=ATF_F, 3=Hake_M, 4=ATF_M
sum(run_ms_LN_1$quantities$suitability[2,1,1:5 ,1:5 ,1])
run_ms_LN_1$quantities$suit_other[2,1,1:5 ,1]
run_ms_LN_1$quantities$avail_food[2,1,1:5 ,1]
run_ms_LN_1$quantities$vulnerability
run_ms_LN_1$quantities$vulnerability_other
run_ms_LN_1$data_list$other_food
run_ms_LN_1$quantities$diet_hat

dim(run_ms_LN_3$quantities$M2_at_age)

sum(run_ms_LN_3$quantities$M2_at_age[1,1,1:5,1])

plot_mortality(run_ms_LN_3, type = 1)

# Index and catch fits can show multiple models
plot_index(mod_list, model_names= mod_names)
plot_logindex(mod_list, model_names= mod_names)
plot_catch(mod_list, model_names= mod_names)

# Composition plots can be only show one model at a time
plot_comp(ss_run) # - Produces many plots!
plot_comp(run_ms_LN_3)

# Plot mortality-at-age
library(ggplot2)
plot_mortality(ms_run)
plot_mortality(run_ms_LN_3)

# Plot ration
plot_ration(run_ms_LN_3)
plot_ration(model2)

plot_selectivity(run_ms_LN_3)
plot_selectivity(ss_run)

model1<- ss_run
# Look at jnll composition
# - pull “dev-name-change” branch to get unweighted likelihood components!
model1$quantities$jnll_comp
sum(model1$quantities$jnll_comp[-20,]) # - Sum across everything except diet likelihood
sum(model1$quantities$unweighted_jnll_comp[-20,]) # - Excludes weights

# Get residual mortality (time variant M1)
model1$quantities$M1_at_age[1,1,,]

###plots
##plot single species models
mod_list <- list(ss_run, ms_run, run_ms_LN_3)
mod_names <- c("ss_run", "Cannib", "run_ms_LN_estM")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_m_at_age(Rceattle = mod_list, model_names = mod_names)
plot_mortality(run_ms_LN_3)

dim(run_ms_LN_3$quantities$M1_at_age)
dim(run_ms_LN_3$quantities$M2_at_age)
dim(run_ms_LN_3$quantities$M_at_age)

run_ms_LN_3$quantities$M
#=================================================================================================
### increase weigth comps =============== (I THINK THIS MAKE THINGS WORST)
test_data <- ATF_hakedata

# Fixed weights, constrained preferences
test_data$Diet_comp_weights <- c(1, 4)

inits = ms_run$estimated_params
map = ms_run$map # gam_a, gam_b, and log_phi are turned off here
#ms_run$map$mapList$log_gam_a

inits$diet_comp_weights <- c(1, 4)

# Create a list prey size preference
# Set weight ratio parameters
inits$log_gam_a = c(0, 3.006)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, 1.887)  # Standard deviation of log weight ratio for ATF, 0 for other species

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on ATF
inits$log_phi[2,2] <- -999 # Set ATF do not feed on ATF
#inits$log_phi[2,1] <- 4  # set a better starting value

# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on atf
map$mapList$log_phi[2,2] <- NA #so we dont estimate atf on atf
#map$mapList$log_phi[2,1] <- NA

map$mapFactor$log_phi <- factor(map$mapList$log_phi)

# Run without reweighting
run_ms_LN_M_wg <- Rceattle::fit_mod(data_list = test_data,
                                    inits = inits, # Initial parameters from single species ests
                                    map = map,
                                    M1Fun = build_M1(M1_model = 0,
                                                     updateM1 = TRUE,
                                                     M1_use_prior = FALSE,
                                                     M2_use_prior = FALSE),
                                    file = NULL, # Don't save
                                    estimateMode = 0, # Estimate
                                    niter = 3, # 3 iterations around population and predation dynamics
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 1, # MSVPA based
                                    loopnum = 5,
                                    phase = TRUE,
                                    suitMode = c(0, 4), # empirical + LN suitability
                                    initMode = 2,
                                    verbose = 1)

run_ms_LN_M_wg$quantities$jnll #5333.507
run_ms_LN_M_wg$estimated_params$log_phi #0.5, 6.10
run_ms_LN_M_wg$quantities$vulnerability #0,997
run_ms_LN_M_wg$quantities$vulnerability_other #0.00221676
run_ms_LN_M_wg$estimated_params$diet_comp_weights #1,2
run_ms_LN_M_wg$data_list$Diet_weights_mcallister #0.961507022 0.002632437

###plots
##plot single species models
mod_list <- list(run_ms_LN_M_wg, ms_run, run_ms_LN_M)
mod_names <- c("run_ms_LN_M_wg", "ms_run", "run_ms_LN_M")

# Plot biomass trajectory (CHAOS!!!)
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_ration(run_ms_LN_M)
plot_ration(run_ms_LN_M_wg)

#=======================================================================================
## MODEL RE-Weigthing
# This function performs iterative reweighting of diet composition weights
perform_model_reweighting <- function(initial_model, test_data, inits, map,
                                      max_iterations = 10,
                                      verbose = TRUE) {

  # Initialize storage for results
  model_runs <- list()
  iteration_metrics <- data.frame(
    iteration = integer(),
    jnll = numeric(),
    weight_change = numeric()
  )

  # Store initial model
  current_model <- initial_model
  model_runs[[1]] <- current_model

  # Get initial weights
  previous_weights <- current_model$data_list$Diet_weights_mcallister

  if(verbose) {
    cat("Starting iterative reweighting process...\n")
    cat("Initial JNLL:", current_model$quantities$jnll_comp, "\n")
  }

  # Iterative reweighting loop - runs all iterations
  for(i in 1:max_iterations) {

    if(verbose) cat("\n--- Iteration", i, "---\n")

    # Update data and inits with new weights
    test_data$Diet_comp_weights <- current_model$data_list$Diet_weights_mcallister
    inits$diet_comp_weights <- current_model$data_list$Diet_weights_mcallister

    # Fit model with updated weights
    current_model <- Rceattle::fit_mod(data_list = test_data,
                                       inits = inits, # Initial parameters from single species ests
                                       map = map,
                                       M1Fun = build_M1(M1_model = 0,
                                                        updateM1 = TRUE,
                                                        M1_use_prior = FALSE,
                                                        M2_use_prior = FALSE),
                                       file = NULL, # Don't save
                                       estimateMode = 0, # Estimate
                                       niter = 3, # 3 iterations around population and predation dynamics
                                       random_rec = FALSE, # No random recruitment
                                       msmMode = 1, # MSVPA based
                                       loopnum = 5,
                                       phase = TRUE,
                                       suitMode = c(0, 4), # empirical + LN suitability
                                       initMode = 2,
                                       verbose = 1)

    # Store model
    model_runs[[i + 1]] <- current_model

    # Get current weights
    current_weights <- current_model$data_list$Diet_weights_mcallister

    # Calculate metrics for tracking
    weight_change <- max(abs(current_weights - previous_weights))
    current_jnll <- current_model$quantities$jnll_comp

    # Store metrics
    iteration_metrics <- rbind(iteration_metrics,
                               data.frame(
                                 iteration = i,
                                 jnll = current_jnll,
                                 weight_change = weight_change
                               ))

    if(verbose) {
      cat("JNLL:", current_jnll, "\n")
      cat("Max weight change:", weight_change, "\n")
      cat("Vulnerability:", current_model$quantities$vulnerability, "\n")
    }

    # Update previous weights for next iteration
    previous_weights <- current_weights
  }

  if(verbose) {
    cat("\nCompleted all", max_iterations, "iterations\n")
  }

  # Return results
  return(list(
    models = model_runs,
    iteration_metrics = iteration_metrics,
    final_model = current_model,
    n_iterations = max_iterations
  ))
}

# Usage example:
reweight_results <- perform_model_reweighting(
  initial_model = run_ms_LN_M, #run_ms_LN_M_wg,
  test_data = test_data,
  inits = inits,
  map = map,
  max_iterations = 5,
  verbose = TRUE
)

rw_model <- reweight_results$final_model
all_models <- reweight_results$models

reweight_results$iteration_metrics

rw_model$quantities$jnll #2204.869
rw_model$estimated_params$log_phi #0.5, 28.39
rw_model$quantities$vulnerability
rw_model$estimated_params$diet_comp_weights #1,1
rw_model$data_list$Diet_weights_mcallister #1.39011120, 0.05770769

###plots
##plot single species models
mod_list <- list(rw_model, run_ms_LN_M)
mod_names <- c("rw_model", "run_ms_LN_M")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names) #Now biomass looks alike
plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names) #this looks pretty different
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_ration(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_b_eaten(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_diet_comp_corrected(rw_model)

# Show first 5 values for ATF (species 2, sex 1, first 5 ages, year 1:5) ## SAME VALUES
print(run_ms_LN_M_wg$quantities$consumption_at_age[2, 1, 1:5, 1:5])
print(run_ms_LN_M$quantities$consumption_at_age[2, 1, 1:5, 1:5])

#Show suitability (ATF age 5 -> hake age 1) ## SAME VALUES
print(run_ms_LN_2$quantities$suitability[2, 1, 5, 1:5, 1:5])
print(run_ms_LN_1$quantities$suitability[2, 1, 5, 1:5, 1:5])

#Hake abundance (age 1, year 1) ## DIFFERENT AS EXPECTED BASE ON DIFF IN BIOMASS AND RATION
print(run_ms_LN_M_wg$quantities$avgN_at_age[1, 1, 1, 1])
print(run_ms_LN_M$quantities$avgN_at_age[1, 1, 1, 1])

#ATF abundance (age 5, year 1) ## (VALUES DO NOT CHANGE)
print(run_ms_LN_M_wg$quantities$avgN_at_age[2, 1, 5, 1])
print(run_ms_LN_M$quantities$avgN_at_age[2, 1, 5, 1])

## Vulnerability, Vulnerability-other, log_phi, diet_comps_weigh, diet_weights_mcallister DIFFER BETWEEN MODELS


# this gives us all biomas eaten (canibalism + ATF)
run_ms_LN_M$quantities$B_eaten_as_prey

#if we want to filer out only biomass eaten by ATF we need to:
run_ms_LN_M$quantities$B_eaten
plot_b_eaten_prop(run_ms_LN_M)
plot_b_eaten_prop(run_ms_LN_M_wg)
plot_b_eaten(run_ms_LN_M)

#filter ATF:(hake = 1, ATF = 2 and 4)
hake_eaten_by_atf <- c()
for(yr in 1:dim(run_ms_LN_M$quantities$B_eaten)[5]){
  hake_eaten_by_atf[yr] <- sum(run_ms_LN_M$quantities$B_eaten[c(2,4), 1, , , yr])
}

hake_eaten_by_atf

# 1. Get the sequence of years from your model object
start_year <- run_ms_LN_M$data_list$styr
all_years <- seq(from = start_year, length.out = length(hake_eaten_by_atf))

# 2. Combine your data into a data frame
biomass_data <- data.frame(
  Year = all_years,
  Biomass_Eaten = hake_eaten_by_atf /1000000 #change from tonnes to million tonnes.
)

# 3. Filter the data to include years up to 2020
biomass_to_plot <- biomass_data[biomass_data$Year <= 2019, ]

# 4. Create the plot
plot(
  x = biomass_to_plot$Year,
  y = biomass_to_plot$Biomass_Eaten,
  type = "l",  # Use "l" for a line plot
  lwd = 2,     # Make the line thicker
  col = "dodgerblue",
  xlab = "Year",
  ylab = "Biomass of Hake Eaten (tons)",
  main = "Hake Consumption by Arrowtooth Flounder",
  ylim = c(0, max(biomass_to_plot$Biomass_Eaten, na.rm = TRUE) * 1.1), # Ensure y-axis starts at 0
  las = 1 # Make y-axis labels horizontal
)

# Add grid lines for better readability
grid()

# Plot ration for female ATF (species 2, sex 1)
plot(model$quantities$ration[2, 1, , 1], type = 'b',
     xlab = "ATF Age", ylab = "Annual Ration (kg/predator)", main = "Input Ration by Predator Age")
