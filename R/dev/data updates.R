data("GOA2018SS")
data_list <- GOA2018SS

# if(any(data_list$fleet_control$Selectivity == 2 & data_list$fleet_control$Time_varying_sel > 1)){
#   data_list$fleet_control <- data_list$fleet_control %>%
#       dplyr::mutate(Sel_curve_pen1 = ifelse(Selectivity == 2 & Time_varying_sel > 1, Time_varying_sel, NA),
#                     Sel_curve_pen2 = ifelse(Selectivity == 2 & Time_varying_sel > 1, Sel_sd_prior, 0),
#                     Time_varying_sel = ifelse(Selectivity == 2 & Time_varying_sel > 1, 0, Time_varying_sel),
#                     Sel_sd_prior = ifelse(Selectivity == 2 & Time_varying_sel > 1, 0, Sel_sd_prior)
#
#       ) %>%
#     dplyr::relocate(Sel_curve_pen1, .after = Nselages) %>%
#     dplyr::relocate(Sel_curve_pen2, .after = Sel_curve_pen1)
#   print("Updating format where 'Selectivity == 2'. Moving non-parametric penalties to 'Sel_curve_pen1' and 'Sel_curve_pen2'.")
# }

# data_list$fleet_control <- data_list$fleet_control %>%
# dplyr::mutate(Time_varying_sel = 0, Sel_sd_prior = 0) %>%
#   dplyr::relocate(Sel_curve_pen1, .after = Nselages) %>%
#   dplyr::relocate(Sel_curve_pen2, .after = Sel_curve_pen1)

# data_list$R_sexr <- NULL
# data_list$est_sex_ratio <- NULL
# data_list$sex_ratio_sigma <- NULL
# data_list$aLW <- NULL
# data_list$maturity <- data_list$pmature
# data_list$pmature <- NULL
# data_list$weight <- data_list$wt
# data_list$wt <- NULL
#
# data_list$diet_data <- data_list$stom_prop_data
# data_list$stom_prop_data <- NULL
#
# ss_run <- Rceattle::fit_mod(data_list = data_list,
#                             inits = NULL, # Initial parameters = 0
#                             file = NULL, # Don't save
#                             estimateMode = 0, # Estimate
#                             random_rec = FALSE, # No random recruitment
#                             msmMode = 0, # Single species mode
#                             phase = TRUE,
#                             verbose = 1)
# plot_biomass(ss_run)
data_list$ration_data <- data_list$Pyrs
data_list$Pyrs <- NULL
GOA2018SS <- data_list

usethis::use_data(GOA2018SS, overwrite = TRUE)
