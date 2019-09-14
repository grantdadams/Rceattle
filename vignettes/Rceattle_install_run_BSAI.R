# Grant Adams, Kirstin Holsman, Andre Punt - April 2019
# Function to run CEATTLE model in TMB
# Citation:
# Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., and Moffitt, E. A. 2015. A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134: 360–378.

# Install devtools if you don't already have it
install.packages("devtools")
# Install TMB and rtools https://cran.r-project.org/bin/windows/Rtools/
# Instructions can be found here for non-pc: https://github.com/kaskr/adcomp/wiki/Download
install.packages('TMB', type = 'source')
# Try "TMB::runExample(all = TRUE)" to see if TMB works

# Install Rceattle
devtools::install_github("grantdadams/Rceattle", auth_token = "4925b42ac46f1e0aefd671e9dc0c1cf1b3157017")



library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data

# Write data to excel
# Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")
#
# # Change the data how you want in excel
# # Read the data back in
# mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")


################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            TMBfilename = "ceattle_v01_06",
                            cpp_directory = "inst/executables",
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = 1, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            silent = TRUE,
                            recompile = TRUE)
# Type ?fit_mod for more details
ests <- mod_objects$estimated_params
inits <- ss_run$estimated_params
inits$ln_mn_rec <- ests$ln_mn_rec

inits$ln_mean_F[1:3] <- ests$ln_mean_F

inits$log_srv_q[7] <- ests$log_srv_q[4]
inits$sel_inf[,4:7,1] <- ests$srv_sel_inf
inits$sel_slp[,4:7,1] <- ests$srv_sel_slp
inits$sel_coff[1:3,] <- ests$fsh_sel_coff
inits$init_dev <- ests$init_dev
inits$F_dev[1:3,] <- ests$F_dev
inits$rec_dev <- ests$rec_dev


ss_run2 <- Rceattle::fit_mod(data_list = BS2017SS,
                            TMBfilename = "ceattle_v01_06",
                            cpp_directory = "inst/executables",
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            silent = TRUE,
                            recompile = TRUE)


round(ss_run2$quantities$jnll_comp, 3)
round(mod_objects$quantities$jnll_comp, 3)

sum(ss_run2$quantities$jnll_comp[-4,])
sum(mod_objects$quantities$jnll_comp[-8,])


ss_run2$quantities$srv_bio_hat
mod_objects$quantities$srv_bio_hat


ss_run2$quantities$S[1,1,2,]
mod_objects$quantities$S[1,2,]

ss_run2$quantities$NByage[1,1,,10:13]
mod_objects$quantities$NByage[1,,10:13]

ss_run2$quantities$Zed[1,1,,10:13]
mod_objects$quantities$Zed[1,,10:13]

ss_run2$quantities$sel[7,1,,10:13]
mod_objects$quantities$srv_sel[4,,10:13]

ss_run2$quantities$sel[1,1,,10:13]
mod_objects$quantities$fsh_sel[1,,10:13]




ss_run2$quantities$flt_spp
ss_run2$quantities$flt_type
ss_run2$quantities$flt_sel_type


mod_objects$quantities$jnll_comp

# The you can plot the model results using using
plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run)


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
# Or we can use the previous data set

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run2$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            niter = 10, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            silent = TRUE,
                            recompile = FALSE)

# We can plot both runs as well:
mod_list <- list(ss_run2, ms_run)
mod_names <- c("SS", "MS")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

plot_selectivity(Rceattle = mod_list, model_names = mod_names)
plot_mort(Rceattle = mod_list, model_names = mod_names, age = 2)

# Run diagnostics
plot_fsh_comp(ms_run) # Fitted fishery composition data
plot_srv_comp(ms_run) # Fitted survey composition data
plot_index(ms_run) # Fitted indices of abundance
plot_catch(ms_run) # Fitted catch series
