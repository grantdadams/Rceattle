# Grant Adams, Kirstin Holsman, Andre Punt - April 2021
# Function to run Gulf of Alaska single-species CEATTLE models in TMB
# Citation:
# Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn, M. W., Ianelli, J. N., Spies, I., … Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fisheries Research, 251(October 2021), 106303. https://doi.org/10.1016/j.fishres.2022.106303

################################################
# Load packages
################################################

library(Rceattle)
library(readxl)

################################################
# Load data
################################################

# Pollock
data("GOApollock")
# Differences between CEATTLE and SAFE: the SAFE model penalizes the first 7 and last recruitment deviates, while I penalize them all. Besides that its parameterized the exact same.

# Cod
data("GOAcod")
mydata_pcod$pmature[1,2:13] <- 2 # Spawn wt from SS model includes sex-ratio and maturity already, so setting Pmature (age-at-maturity) to 2 to have CEATTLE calculations be the same
# The cod model was  a bit trickier to recreate because of the use of internal age-length estimation/conditional age-at-length and more flexible selectivity parameterization than I currently have set up in CEATTLE. I use the output weight-at-age, mortality, and terminal year age-length-key from SS as inputs into CEATTLE and set age-based selectivity in CEATTLE to be the form that is most similar to the terminal year length-based selectivity pattern in SS. I fit the same index, catch, and length-at-age data as well as marginal age- and length-comp data that the SAFE model doesn't fit. I ignore annual varying selectivity (you could do it in CEATTLE... see meta_data sheet in the Excel file, I just stopped trying after some convergence issues). I ignored the aging error (could add in, but was having convergence problems). I also do not use the Taylor and Methot 2013 recruitment bias correction thats in SS and did not have a 2016 bump in mortality from the blob.

# ATF
data("GOAatf")
# Differences between CEATTLE and SAFE: ATF composition sample sizes and non-parametric selectivity penalties are the same for each sex while they are different in the SAFE model (differences between sexes for both bits are < 5 in the SAFE model. This was a convenience for coding


################################################
# Fit models
################################################

# Pollock
GOApollock$styr = 1977 # The SAFE model starts at 1970, so change styr to 1970 to run the full time series model (data is in there). I start them all at 1977 because thats the years with overlap.
pollock_model <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Use default phasing


# Arrowtooth flounder
GOAatf$styr = 1977 # The SAFE model starts at 1961, so change styr to 1961 to run the full time series model (data is in there). I start them all at 1977 because thats the years with overlap.
atf_model <- Rceattle::fit_mod(
  data_list = GOAatf,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Use default phasing


# Cod
# -- Cod start year is 1977 so no need to change
cod_model <- Rceattle::fit_mod(
  data_list = GOAcod,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Use default phasing

# Reweighting the cod model helps using Macalliser Ianelli weights
mydata_pcod$fleet_control$Comp_weights <- cod_model$data_list$fleet_control$Est_weights_macallister

# -- Refit
cod_model <- Rceattle::fit_mod(
  data_list = mydata_pcod,
  inits = mydata_pcod$estimated_params, # Start from the unweighted model's MLEs
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Phase


################################################
# Compare with SAFE Models
################################################
# Columns = year, pollock, cod, atf - outputs here are 1977-2018
data("GOAsafe2018")
safe2018biomass <- GOAsafe2018$biomass
safe2018ssb <- GOAsafe2018$ssb
safe2018rec <- GOAsafe2018$recruitment

# Assign data to CEATTLE object
pollock_model_safe <- pollock_model
atf_model_safe <- atf_model
cod_model_safe <- cod_model

# - Pollock
pollock_model_safe$quantities$biomass[1,1:42] <- t(safe2018biomass[1:42,c(2)]) * 1000 # Rescaling - SAFE outputs in 1,000 mt
pollock_model_safe$quantities$biomassSSB[1,1:42] <- t(safe2018ssb[1:42,c(2)]) * 1000 # Rescaling - SAFE outputs in 1,000 mt

# - ATF
atf_model_safe$quantities$biomass[1,1:42] <- t(safe2018biomass[1:42,c(4)]) * 1000 # Rescaling - SAFE outputs in 1,000 mt
atf_model_safe$quantities$biomassSSB[1,1:42] <- t(safe2018ssb[1:42,c(4)]) * 1000 # Rescaling - SAFE outputs in 1,000 mt

# - Cod
cod_model_safe$quantities$biomass[1,1:42] <- t(safe2018biomass[1:42,c(3)])
cod_model_safe$quantities$biomassSSB[1,1:42] <- t(safe2018ssb[1:42,c(3)])

# Convert Pollock to age-3 biomass rather than age-1 biomass a-la SAFE model
pollock_model$quantities$biomass[1,1:49] <- colSums(pollock_model$quantities$biomassByage[1,3:10,1:49])

# Plot biomass trends
plot_biomass(list(pollock_model, pollock_model_safe), model_names = c("CEATTLE", "SAFE"))
plot_biomass(list(atf_model, atf_model_safe), model_names = c("CEATTLE", "SAFE"))
plot_biomass(list(cod_model, cod_model_safe), model_names = c("CEATTLE", "SAFE"))


################################################
# Example diagnostics
################################################
plot_index(pollock_model)
plot_catch(pollock_model)
plot_selectivity(pollock_model)
plot_indexresidual(pollock_model)
plot_logindex(pollock_model)
plot_recruitment(pollock_model, add_ci = TRUE)
plot_comp(pollock_model)

