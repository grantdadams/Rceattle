#' Combine data sets. Will use the env_data data set from data_set1 and diet data will have to be updated.
#'
#' @param data_list1 Rceattle data_list 1
#' @param data_list2 Rceattle data_list 2
#'
#' @export
#'
combine_data <- function(data_list1 = NULL, data_list2 = NULL){

  data_list_new <- data_list1

  dat_names <- names(data_list_new)

  names_not_used <- c("nspp", "styr", "endyr", "projyr")

  vec_names <- c("spnames", "nsex", "spawn_month", "nages", "minage", "nlengths","pop_wt_index", "ssb_wt_index", "est_M1", "pop_alk_index", "sigma_rec_prior", "other_food", "estDynamics", "Ceq", "Cindex","Pvalue", "fday", "CA","CB", "Qc", "Tco",  "Tcm",  "Tcl",  "CK1", "CK4", "beta_wt_len", "alpha_wt_len") # Object names of vectors

  mat_names <- c("fleet_control", "index_data", "catch_data", "comp_data", "caal_data", "emp_sel", "NByageFixed", "age_trans_matrix", "age_error", "weight",   "maturity", "sex_ratio", "M1_base", "ration_data", "diet_data") # Object names of matrices

  # Get index from data_set1 of the 4 indices
  fleet_index1 <- max(data_list1$fleet_control$Fleet_code, na.rm = TRUE)
  weight_index1 <- max(data_list1$weight$Wt_index, na.rm = TRUE)
  alk_index1 <- max(data_list1$age_trans_matrix$Age_transition_index, na.rm = TRUE)
  nspp1 <- data_list1$nspp
  q_index1 <- max(data_list1$fleet_control$Q_index, na.rm = TRUE)
  sel_index1 <- max(data_list1$fleet_control$Selectivity_index, na.rm = TRUE)

  # Update vector indices
  data_list2[["pop_wt_index"]] <- data_list2[["pop_wt_index"]] + weight_index1
  data_list2[["ssb_wt_index"]] <- data_list2[["ssb_wt_index"]] + weight_index1
  data_list2[["pop_age_transition_index"]] <- data_list2[["pop_age_transition_index"]] + alk_index1

  # Update fleet control
  data_list2$fleet_control$Age_transition_index <- data_list2$fleet_control$Age_transition_index + alk_index1
  data_list2$fleet_control$Selectivity_index <- data_list2$fleet_control$Selectivity_index + sel_index1
  data_list2$fleet_control$Q_index <- data_list2$fleet_control$Q_index + q_index1
  data_list2$fleet_control$Species <- data_list2$fleet_control$Species + nspp1
  data_list2$fleet_control$Weight_index <- data_list2$fleet_control$Weight_index + weight_index1
  data_list2$fleet_control$Fleet_code <- data_list2$fleet_control$Fleet_code + fleet_index1

  # Update fleet and spp indices of matrices
  for(i in 2:15){
    if(!is.null(data_list2[[mat_names[i]]]$Species)){
      data_list2[[mat_names[i]]]$Species <- data_list2[[mat_names[i]]]$Species + nspp1
    }
    if(!is.null(data_list2[[mat_names[i]]]$Fleet_code)){
      data_list2[[mat_names[i]]]$Fleet_code <- data_list2[[mat_names[i]]]$Fleet_code + fleet_index1
    }
  }

  # Update stomach pred sp and prey sp
  data_list2$diet_data$Pred <- data_list2$diet_data$Pred + nspp1
  data_list2$diet_data$Prey <- data_list2$diet_data$Prey + nspp1

  # Update weight and alk indices
  data_list2$weight$Wt_index <- data_list2$weight$Wt_index + weight_index1
  data_list2$age_trans_matrix$Age_transition_index <- data_list2$age_trans_matrix$Age_transition_index + alk_index1



  # Combine vectors
  for(i in vec_names){
    data_list_new[[i]] <- c(data_list1[[i]], data_list2[[i]])
  }

  # Combine matrices
  for(i in mat_names[1:17]){
    data_list_new[[i]] <- plyr::rbind.fill(data_list1[[i]], data_list2[[i]])
  }

  data_list1$env_data <- merge(data_list1$env_data, data_list2$env_data, by = "Year", all = TRUE)


  # Add new species
  data_list_new$nspp <- data_list1$nspp + data_list2$nspp

  return(data_list_new)
}
