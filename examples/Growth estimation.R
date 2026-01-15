# Example fitting internal growth via a von Bertalanffy function in Rceattle
# Includes CAAL data
# - this follows Giancarlo's example in WHAM (https://giancarlomcorrea.netlify.app/labs/OFI_WK_2023/examples/case3.R)
# - from his time-varying growth workshop (https://giancarlomcorrea.netlify.app/post/ofi-workshop/).
# - I modified the models so that all fleets use length-based logistic selectivity.
# - see "tests/comparison/WHAM growth comparison.R" for full comparison

library(Rceattle)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# To run model we use data from Giancarlo's class formatted for Rceattle:
# - see line 134 of "tests/comparison/WHAM growth comparison.R" for modification
data(whamGrowthData)

# The model includes conditional-age-at-length data for the survey
# - The lengths from the caal data will be used for growth estimation
# - if the number of lengths in the caal_data for each species does not match "data$nlengths" it will produce an error.
head(whamGrowthData$caal_data)

# - To convert length-at-age to weight-at-age Rceattle uses time- and sex-invariant power equation: W = a * L ^ b
# - with parameters specified in the data for each species
whamGrowthData$alpha_wt_len
whamGrowthData$beta_wt_len

# - The month of spawning and a survey/fishery occurs is also specified in the data
whamGrowthData$spawn_month         # Spawning
whamGrowthData$fleet_control$Month # Fleets


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Estimation ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# To estimate growth we can set "growthFun" using switches in "build_growth"
# See "?build_growth" for all details

vbgf_model <- Rceattle::fit_mod(data_list = whamGrowthData,
                            inits = NULL,       # Initial parameters at default
                            estimateMode = 0,   # Estimate
                            growthFun = build_growth(growth_model = 1), # Von Bert
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0,        # Single species mode
                            phase = FALSE,
                            initMode = 2,       # Unfished disequilibrium
                            verbose = 1)

richards_model <- Rceattle::fit_mod(data_list = whamGrowthData,
                                inits = NULL,       # Initial parameters at default
                                estimateMode = 0,   # Estimate
                                growthFun = build_growth(growth_model = 2), # Richards
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0,        # Single species mode
                                phase = FALSE,
                                initMode = 2,       # Unfished disequilibrium
                                verbose = 1)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Plot and compare ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
plot_biomass(list(vbgf_model, richards_model), model_names = c("VBGF", "Richards"), incl_proj = TRUE)

# - AIC
vbgf_model$opt$AIC
richards_model$opt$AIC
