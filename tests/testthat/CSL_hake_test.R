library(Rceattle)

#CANN_hakedata <- Rceattle::read_data(file = "051225_Cannib_hake.xlsx")
CSL_hakedata <- Rceattle::read_data(file = "051225_CSL_Hake_model_FINAL.xlsx")

ss_run <- Rceattle::fit_mod(data_list = CSL_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run$quantities$jnll #1749.188
#save(ss_run, file = "results/models/ATF/ATF_ss_run.Rdata")

## MODEL with cannibalism and CSL
ms_run <- Rceattle::fit_mod(data_list = CSL_hakedata,
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

ms_run_Mest <- Rceattle::fit_mod(data_list = CSL_hakedata,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            M1Fun = build_M1(M1_model = 1,  #do not estimate mortality!
                                             updateM1 = TRUE,
                                             M1_use_prior = FALSE,
                                             M2_use_prior = FALSE),
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            phase = FALSE,
                            suitMode = 0, # empirical suitability
                            initMode = 2, # Fished start with init devs
                            verbose = 1)


ms_run$quantities$jnll #1762.849
plot_diet_comp3(ms_run)
ms_run_Mest$quantities$jnll #158173.3

# Check the M1 value in the crashed model
print(ms_run_Mest$quantities$M1[1, 1, 1:5, 1])

# Check the M1 value in the working model
print(ms_run$quantities$M1[1, 1, 1:5, 1])


models <- list(ss_run, ms_run)
model_names <- c("ss_run", "ms_run_Mnoets")

plot_biomass(models, model_names= model_names)
plot_recruitment(models, model_names= model_names)
plot_ssb(models, model_names= model_names)

plot_diet_comp(ms_run)
plot_diet_comp(ms_run_Mest)

plot_b_eaten(models, model_names= model_names)
plot_ration(models, model_names= model_names)
plot_b_eaten_prop(models, model_names= model_names)
plot_m2_at_age_prop(models, model_names= model_names)





ms_run$quantities$ration[2,1,1:5 ,1]

#save(ms_run, file = "results/models/ATF/ATF_ms_run.Rdata")

plot_diet_comp(ms_run) #diet estimate hake canni is good, in this model we do not have ATF predation, hence = 0
plot_diet_comp2(ms_run)
plot_b_eaten_prop(ms_run)
plot_m2_at_age_prop(ms_run)

plot_ration(ms_run)
plot_ration(ms_run_tonnes)
plot


sum(ms_run$quantities$suitability[2,1,1:5 ,1:5 ,1])
sum(ms_run_tonnes$quantities$suitability[2,1,1:5 ,1:5 ,1])
ms_run$quantities$suit_other[2,1,1:5 ,1]
ms_run$quantities$avail_food[2,1,1:5 ,1]
ms_run$quantities$vulnerability
ms_run$quantities$vulnerability_other
ms_run$data_list$other_food
ms_run$quantities$diet_hat
ms_run$quantities$B_eaten
ms_run$quantities$consumption_at_age[2,1,1:5 ,1]

sum(ms_run$quantities$avgN_at_age[2,2,1:5 ,1])
ms_run$quantities$biomass_at_age[2,2,1:5 ,1]
ms_run$quantities$M2_at_age[1,1,1:5 ,1]
ms_run$quantities$B_eaten

sum(ms_run$quantities$diet_prop[2, 1, 1:5, 1:5, 1])

# This is the intermediate calculation before normalization
sum(ms_run$quantities$stom_div_bio[2, 1, 1:5, 1:5, 1])

# Sum of diet proportions across all prey for each predator age
ms_run$quantities$diet_prop_hat[2, 2, 1:5, 1]




# Should be: diet_prop / (hake_N × hake_weight)
# If NaN, either division by zero or diet_prop is bad

summary(ms_run$data_list$diet_data)


# Correct way to access the data sent to TMB
ms_run$obj$env$data$n_stomach_obs
head(ms_run$obj$env$data$stomach_id)

ms_run$quantities$jnll_comp
ms_run$quantities$diet_prop_hat
ms_run$quantities$jnll
ms_run$quantities$jnll_comp
ms_run$estimated_params$log_phi
ms_run$quantities$vulnerability
ms_run$estimated_params$diet_comp_weights
ms_run$data_list$Diet_weights_mcallister

model1<- ms_run
# Look at jnll composition
# - pull “dev-name-change” branch to get unweighted likelihood components!
model1$quantities$jnll_comp
sum(model1$quantities$jnll_comp[-20,]) # - Sum across everything except diet likelihood
sum(model1$quantities$unweighted_jnll_comp[-20,]) # - Excludes weights



# Check age range
range(ms_run$data_list$diet_data$Pred_age)  # Should be 1 to 24
range(ms_run$data_list$diet_data$Prey_age)      # Should be 1 to 20 (or 0-19?)

# Check diet values
range(ms_run$data_list$diet_data$Stomach_proportion_by_weight)
# Should show 0.0695, 0.0022, etc. (not all zeros!)

# 1. Check what the Model thinks the dimensions are
print(paste("Model Ages:", paste(ms_run$data_list$nages, collapse=", ")))
print(paste("Model Sexes:", paste(ms_run$data_list$nsex, collapse=", ")))

# 2. Check what is in your Diet Data CSV
diet_data <- ms_run$data_list$diet_data

print(paste("Max Predator Age in Data:", max(diet_data$Pred_age)))
print(paste("Max Prey Age in Data:", max(diet_data$Prey_age)))

models <- list(ss_run, ms_run)
model_names <- c("ss_run", "Cannib+SeaLion")

plot_biomass(ss_run)
plot_biomass(models, model_names= model_names)
plot_recruitment(models, model_names= model_names)
plot_ssb(models, model_names= model_names)

plot_diet_comp(ms_run)

# Look at the first few rows of the data ACTUALLY used by the model
head(ms_run$data_list$diet_data)

# Look at a few values for Predator 2 (CSL) eating Prey 1 (Hake)
print(head(ms_run$quantities$suitability[2,1,1:5, 1:5, 1]))
# Sum B_eaten across all ages/years to get total tonnes
total_hake_eaten_tonnes <- sum(ms_run$quantities$B_eaten[2,1,,,])
print(total_hake_eaten_tonnes)
# Check Male (Sex 2) Age 20 consumption
# [Pred, Prey, PredAge, PreyAge, Year]
print(sum(ms_run$quantities$B_eaten[2, 1, 20, , 1]))
# Check M2 for Hake (Prey 1)
# Look at M2 for Age 1 Hake across the first few years
print(head(ms_run$quantities$M2_at_age[1, 1, 1, ]))

run_ms_LN_3<- ms_run

# 1. Get the Ration per Individual (Output from Model)
#    Indices: [Species, Sex, Age, Year]
#    For CSL (Sp 2), Female (Sex 1), Age 5, Year 1 (1980)
model_ration <- ms_run$quantities$ration[2, 1, 5, 1]

# 2. Get the Population N (Output from Model)
model_n <- ms_run$quantities$avgN_at_age[2, 1, 5, 1]

# 3. Calculate Total Model Consumption
model_total <- model_ration * model_n

model_total[2, 1, 5, 1]

# 4. Print Comparison
print(paste("Rceattle Ration (Per Indiv):", model_ration))
print(paste("Rceattle Population N:", model_n))
print(paste("Rceattle Total Consumption:", model_total))

# Calculate total CSL consumption of hake
total_hake_eaten <- sum(ms_run$quantities$B_eaten[2, 1, , , ])  # CSL eating Hake, all ages/years

# Compare to expected:
# From your bioenergetics:
# - Total CSL population × average ration × proportion hake in diet
# - Should be in millions of kg for the whole time series

# Per year:
hake_eaten_per_year <- apply(ms_run$quantities$B_eaten[2, 1, , , ], 3, sum)
plot(hake_eaten_per_year, type='l', main="Total Hake Consumed by CSL")
dev.off()

# Manual calculation for year 1:
csl_N <- ms_run$quantities$avgN_at_age[2, 1, , 1]  # Males, year 1
csl_ration <- ms_run$quantities$ration[2, 1, , 1]  # kg/year
csl_suitability <- ms_run$quantities$suitability[2, 1, , , 1]  # CSL eating Hake

ms_run$quantities$consumption_at_age[2, 1, , 1]

# Total consumption of hake age-1:
manual_calc <- sum(csl_N * csl_ration * csl_suitability)
model_calc <- sum(ms_run$quantities$B_eaten[2, 1, , , 1])


manual_calc <- sum(csl_N * csl_ration)

# Should be similar:
data.frame(
  Manual = manual_calc,
  Model = model_calc,
  Difference = manual_calc - model_calc
)

model_consumption <- apply(ms_run$quantities$B_eaten[2, 1, , , ], 3, sum)

# Compare N vs avgN
ms_run$quantities$N_at_age[2, 2, 5, 1]      # N at start of year
ms_run$quantities$avgN_at_age[2, 2, 5, 1]   # Average during year
# avgN should be smaller than N

ms_run$quantities$B_eaten[2, 1, 1:5,1:5, 1]
ss_run$quantities$biomass_at_age[1, 1, 1:5, 1]


devtools::document()
