# This script fits CEATTLE to data from Cod, Hake, and Herring off Georges Bank

# The model is a recreation of Curti et al 2013
# - Curti, K. L., Collie, J. S., Legault, C. M., & Link, J. S. (2013). Evaluating the performance of a multispecies statistical catch-at-age model. Canadian Journal of …, 484(January), 470–484. https://doi.org/10.1139/cjfas-2012-0229

# NOTE: The model uses input consumption-to-biomass data (Pyrs) and Ceq = 4
# The model also fixes suitability parameters to that estimated externally from the data


# - Load the data and package
library(Rceattle)
data(GeorgesBank3spp)
write_data(GeorgesBank3spp, file = "Georges_Bank_3spp.xlsx")

# - Review excel data and load
georges_bank <- read_data(file = "Georges_Bank_3spp.xlsx")
georges_bank$est_M1 <- c(0,0,0) # Fix M1 to input value (M1_base)


################################################
# Estimation
################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = georges_bank,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


# - Multi-species
# For the a multispecies model we start from the single species parameters.
georges_bank$est_M1 <- c(0,0,0) # Fix residual M

# - We also need to fix suitability and input values
georges_bank$msmMode = 1
georges_bank$suitMode = 4
inits <- ss_run$estimated_params

# - Update suitability parameters
inits$log_phi[1,] <- c(0.618,1.549,3.255) # Predator rows, prey columns (cod)
inits$log_phi[2,] <- c(-Inf,5.396,3.223) # Predator rows, prey columns (hake)
inits$log_phi[3,] <- -Inf # Herring

# Time varying weight based lognormal suitability
inits$log_gam_a <- log(c(4,3,0)) # Mean of lognormal
# - Update to be species specific
# Cod = 4.159, 4.833, 3.996
# Silver hake = NA, 3.946, 2.261

inits$log_gam_b <- log(c(2,2,0)) # Variance of lognormal
# - Update to be species specific
# Cod = 2.259, 1.875, 1.433
# Silver hake = NA, 2.979, 1.093

# - Map out suit params
map <- build_map(georges_bank, inits)
map$mapList$log_phi <- matrix(NA, 3, 3)
map$mapList$log_gam_a <- rep(NA, 3)
map$mapList$log_gam_b <- rep(NA, 3)

map$mapFactor$log_phi <- as.factor(map$mapList$log_phi)
map$mapFactor$log_gam_a <- as.factor(map$mapList$log_gam_a)
map$mapFactor$log_gam_b <- as.factor(map$mapList$log_gam_b)

ms_run <- Rceattle::fit_mod(data_list = georges_bank,
                            map = map,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 4, # Time-varying weight lognormal
                            verbose = 1)


################################################
# Plotting
################################################
# We can plot all runs
mod_list <- list(ss_run, ms_run)
mod_names <- c("Single-species", "Multi-species")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

# Plot mortality and predation
plot_b_eaten(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey by each predator
plot_mort(Rceattle = ms_run, type = 3) # Mortality-at-age time series

# Run diagnostics
plot_selectivity(Rceattle = ms_run)
plot_comp(ms_run) # Fitted survey composition data
plot_index(ms_run) # Fitted indices of abundance
plot_catch(ms_run) # Fitted catch series


# Numbers at age (thousands)
par(mfrow = c(3,1), mar = c(2,4,1,1))
plot(y = colSums(ms_run$quantities$NByage[1,1,,1:30])/1000, x = 1978:2007, type = "l", ylab = "N")
lines(y = colSums(ss_run$quantities$NByage[1,1,,1:30])/1000, x = 1978:2007, lty = 2)
legend("topright", "Cod", bty = "n")

plot(y = colSums(ms_run$quantities$NByage[2,1,,1:30])/1000, x = 1978:2007, type = "l", ylab = "N")
lines(y = colSums(ss_run$quantities$NByage[2,1,,1:30])/1000, x = 1978:2007, lty = 2)
legend("topright", "Hake", bty = "n")

plot(y = colSums(ms_run$quantities$NByage[3,1,,1:30])/1000, x = 1978:2007, type = "l", ylab = "N")
lines(y = colSums(ss_run$quantities$NByage[3,1,,1:30])/1000, x = 1978:2007, lty = 2)
legend("topright", "Herring", bty = "n")
