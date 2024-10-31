library(dplyr)
data("GOA2018SS")

# Transpose fleet_control if long format
if(sum(colnames(GOA2018SS$fleet_control)[1:2] == c("Fleet_name", "Fleet_code")) != 2){ #, "Fleet_type", "Species", "Selectivity_index", "Selectivity")) != 6){
  GOA2018SS$fleet_control <- as.data.frame(t(GOA2018SS$fleet_control))
  colnames(GOA2018SS$fleet_control) <- GOA2018SS$fleet_control[1,]
  GOA2018SS$fleet_control <- GOA2018SS$fleet_control[-1,]
  GOA2018SS$fleet_control <- cbind(data.frame(Fleet_name = rownames(GOA2018SS$fleet_control)),
                                   GOA2018SS$fleet_control)
  rownames(GOA2018SS$fleet_control) = NULL
  GOA2018SS$fleet_control[,-which(colnames(GOA2018SS$fleet_control) %in% c("Fleet_name", "Time_varying_q"))] <- apply(
    GOA2018SS$fleet_control[,-which(colnames(GOA2018SS$fleet_control) %in% c("Fleet_name", "Time_varying_q"))], 2, as.numeric)
}

GOA2018SS$fleet_control <- GOA2018SS$fleet_control %>%
  dplyr::select(-c(Accumatation_age_lower, Accumatation_age_upper)) %>%
  mutate(Comp_loglike = 0) %>%
  relocate(Comp_loglike, .after = Age_first_selected) %>%
  mutate(Age_max_selected = NA) %>%
  relocate(Age_max_selected, .after = Age_first_selected)

usethis::use_data(GOA2018SS, overwrite = T)
