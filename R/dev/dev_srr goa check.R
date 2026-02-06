# Code to run Gulf of Alaska CEATTLE model in TMB
# Citation:
# Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn, M. W., Ianelli, J. N., Spies, I., ... & Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fisheries Research, 251, 106303.

library(Rceattle)
library(dplyr)

################################################
# Data
################################################
# Example
# To run the 2018 single species assessment for the Gulf of Alaska, a data file must first be loaded:
data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data


################################################
# Estimation
################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
# GOA2018SS$fleet_control$proj_F_prop <- rep(1, nrow(GOA2018SS$fleet_control))
#GOA2018SS$fleet_control$Time_varying_sel[8] <- 0

GOA2018SS$fleet_control$Sel_norm_bin1[2] <- 10
ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


GOA2018SS$fleet_control$Time_varying_sel[8] <- 5
ss_run2 <- Rceattle::fit_mod(data_list = GOA2018SS,
                            inits = ss_run$estimated_params, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

load("~/Documents/GitHub/goa save.RData")
plot_biomass(list(mod_objects, ss_run))

inits <- mod_objects$estimated_params

names(ss_run$estimated_params)[!names(ss_run$estimated_params) %in% names(inits)]

inits$R_ln_sd <- inits$ln_rec_sigma
inits$ln_rec_sigma <- NULL

inits$index_ln_q <- inits$ln_srv_q
inits$ln_srv_q <- NULL

inits$index_q_dev <- inits$ln_srv_q_dev
inits$ln_srv_q_dev <- NULL

inits$index_q_ln_sd <- inits$ln_sigma_srv_q
inits$ln_sigma_srv_q <- NULL

inits$sel_dev_ln_sd <- inits$ln_sigma_sel
inits$ln_sigma_sel <- NULL

inits$index_q_dev_ln_sd <- inits$ln_sigma_time_varying_srv_q
inits$ln_sigma_time_varying_srv_q <- NULL

inits$index_ln_sd <- inits$ln_sigma_srv_index
inits$ln_sigma_srv_index <- NULL

inits$catch_ln_sd <- inits$ln_sigma_fsh_catch
inits$ln_sigma_fsh_catch <- NULL

inits$index_q_beta <- ss_run$estimated_params$index_q_beta

inits$beta_rec_pars <- ss_run$estimated_params$beta_rec_pars

inits$index_q_rho <- ss_run$estimated_params$index_q_rho

names(ss_run$estimated_params)[!names(ss_run$estimated_params) %in% names(inits)]
names(inits)[!names(inits) %in% names(ss_run$estimated_params)]

inits[c("ln_pop_scalar", "ln_sex_ratio_sigma", "logH_1", "logH_1a", "logH_1b", "logH_2", "logH_3", "H_4", "log_gam_a", "log_gam_b", "log_phi")] = NULL

map <- build_map(ss_run$data_list, inits)
names(map$mapList)[!names(map$mapList) %in% names(inits)]
names(inits)[!names(inits) %in% names(map$mapList)]

ss_run$data_list$fleet_control$Comp_loglike <- -1
ss_init <- Rceattle::fit_mod(data_list = ss_run$data_list,
                            inits = inits, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 3, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)
sum(ss_init$quantities$sel_at_age-mod_objects$quantities$sel_at_age)
sum(mod_objects$quantities$NByage-ss_init$quantities$N_at_age)
sum(mod_objects$quantities$comp_hat-ss_init$quantities$comp_hat)
sum(mod_objects$quantities$comp_n-ss_init$quantities$comp_n)
check <- (mod_objects$quantities$comp_obs-ss_init$quantities$comp_obs)

sum(ss_init$quantities$comp_hat)
round(ss_init$quantities$jnll_comp, 4)-round(mod_objects$quantities$jnll_comp, 4)

check <- mod_objects$data_list$comp_data %>%
  filter(Fleet_code == 10)
sum(check[,-1])

ss_init$quantities$ln_catch_sd-mod_objects$quantities$fsh_log_sd_hat
ss_init$quantities$catch_hat-mod_objects$quantities$fsh_bio_hat
ss_init$quantities$catch_hat[1]
mod_objects$quantities$fsh_bio_hat[1]


mod_objects$quantities$F_spp-ss_init$quantities$F_spp
sum(mod_objects$quantities$F_flt-ss_init$quantities$F_flt)
sum(mod_objects$quantities$F_flt_age-ss_init$quantities$F_flt_age)
sum(mod_objects$quantities$M-ss_init$quantities$M_at_age)
sum(mod_objects$quantities$Z-ss_init$quantities$Z_at_age)
(mod_objects$quantities$NByage[1,1,,1]-ss_init$quantities$N_at_age[1,1,,1])
sum(ss_init$data_list$wt[,-1]-mod_objects$data_list$wt[,-1], na.rm = T)

ss_init$quantities$N_at_age[1,1,,1]


ss_init$quantities$mort_sum[1,]
mod_objects$quantities$mort_sum[1,]

mod_objects$data_list$initMode
ss_init$data_list$initMode

ss_init$map$mapList$ln_Finit
mod_objects$map$mapList$ln_Finit

ss_init$estimated_params$ln_Finit
mod_objects$estimated_params$ln_Finit

ss_init$quantities$M1[1,,]
mod_objects$quantities$M1[1,,]

ss_init$quantities$init_dev[1,]-mod_objects$quantities$init_dev[1,]
ss_init$quantities$Finit
mod_objects$quantities$Finit


mod_objects$quantities$NByage[1,1,,1]
ss_init$quantities$N_at_age[1,1,,1]
ss_init$quantities$N_at_age[1,1,1,1] * exp(-ss_init$quantities$M_at_age[1,1,1,1] + ss_init$estimated_params$init_dev[1,1] - exp(-10))

plot_biomass(list(ss_init, mod_objects))
plot_recruitment(list(ss_init, mod_objects))

mod_objects2 <- ss_init
mod_objects2$quantities$catch_hat <- mod_objects$quantities$fsh_bio_hat
plot_catch(list(mod_objects2, ss_init))
