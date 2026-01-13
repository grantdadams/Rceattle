data("GOApollock")
data_list <- GOApollock
# if(is.null(data_list$fleet_control$Sel_curve_pen1)){ data_list$fleet_control$Sel_curve_pen1 = NA}
# if(is.null(data_list$fleet_control$Comp_loglike)){ data_list$fleet_control$Comp_loglike = 0}
# if(is.null(data_list$fleet_control$Sel_curve_pen2)){ data_list$fleet_control$Sel_curve_pen2 = NA}
# data_list$fleet_control <- data_list$fleet_control %>%
#   dplyr::mutate(CAAL_loglike = 0) %>%
#   dplyr::mutate(CAAL_weights = 1) %>%
#   dplyr::rename(Bin_first_selected = Age_first_selected,
#                 Catchability = Estimate_q,
#                 Time_varying_sel_sd_prior = Sel_sd_prior) %>%
#   dplyr::select(Fleet_name, Fleet_code, Fleet_type, Species, Month, Selectivity_index, Selectivity, N_sel_bins, Sel_curve_pen1, Sel_curve_pen2, Time_varying_sel, Time_varying_sel_sd_prior, Bin_first_selected, Sel_norm_bin1, Sel_norm_bin2, Comp_loglike, Comp_weights, CAAL_loglike, Weight1_Numbers2, Weight_index, Age_transition_index, Q_index, Catchability, Q_prior, Q_sd_prior, Time_varying_q, Time_varying_q_sd_prior, Estimate_index_sd, Index_sd_prior, Estimate_catch_sd, Catch_sd_prior, proj_F_prop)

# Rename "Estimate_q" to "Catchability"
# Rename "Sel_sd_prior" to "Time_varying_sel_sd_prior"
# data_list$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 7))
# colnames(data_list$caal_data) <- c("Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", "CAAL_1")


colnames(data_list$maturity) <- c("Species", paste0("Age", 1:max(data_list$nages)))
colnames(data_list$sex_ratio) <- c("Species", paste0("Age", 1:max(data_list$nages)))

GOApollock <- data_list

usethis::use_data(GOApollock, overwrite = TRUE)
