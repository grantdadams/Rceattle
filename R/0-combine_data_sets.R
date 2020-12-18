#' Combine data sets. Will use the env_data data set from data_set1 and diet data will have to be updated.
#'
#' @param data_list1
#' @param data_list2
#'
#' @export
#'
combine_data <- function(data_list1 = NULL, data_list2 = NULL){

  data_list_new <- data_list1

  dat_names <- names(data_list_new)

  names_not_used <- c("nspp", "styr", "endyr", "projyr")

  vec_names <- c("spnames", "nsex", "spawn_month", "R_sexr", "nages", "minage", "nlengths","pop_wt_index", "ssb_wt_index", "pop_alk_index", "sigma_rec_prior", "other_food", "estDynamics", "proj_F", "est_sex_ratio", "sex_ratio_sigma","Ceq", "Cindex","Pvalue", "fday", "CA","CB", "Qc", "Tco",  "Tcm",  "Tcl",  "CK1", "CK4") # Object names of vectors

  mat_names <- c("fleet_control", "srv_biom", "fsh_biom", "comp_data", "emp_sel", "NByageFixed", "age_trans_matrix", "age_error", "wt",   "pmature", "sex_ratio", "M1_base", "Mn_LatAge", "aLW", "Pyrs", "UobsAge", "UobsWtAge", "env_data") # Object names of matrices

  # Get index from data_set1 of the 4 indices
  fleet_index1 <- max(data_list1$fleet_control$Fleet_code, na.rm = TRUE)
  weight_index1 <- max(data_list1$wt$Wt_index, na.rm = TRUE)
  alk_index1 <- max(data_list1$age_trans_matrix$ALK_index, na.rm = TRUE)
  nspp1 <- data_list1$nspp
  q_index1 <- max(data_list1$fleet_control$Q_index, na.rm = TRUE)
  sel_index1 <- max(data_list1$fleet_control$Selectivity_index, na.rm = TRUE)
  nspp1 <- data_list1$nspp

  # Update vector indices
  data_list2[["pop_wt_index"]] <- data_list2[["pop_wt_index"]] + weight_index1
  data_list2[["ssb_wt_index"]] <- data_list2[["ssb_wt_index"]] + weight_index1
  data_list2[["pop_alk_index"]] <- data_list2[["pop_alk_index"]] + alk_index1

  # Update fleet control
  data_list2$fleet_control$ALK_index <- data_list2$fleet_control$ALK_index + alk_index1
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
  for(i in 16:17){
    data_list2[[mat_names[i]]]$Pred <- data_list2[[mat_names[i]]]$Pred + nspp1
    data_list2[[mat_names[i]]]$Prey <- data_list2[[mat_names[i]]]$Prey + nspp1
  }

  # Update wt and alk indices
  data_list2$wt$Wt_index <- data_list2$wt$Wt_index + weight_index1
  data_list2$age_trans_matrix$ALK_index <- data_list2$age_trans_matrix$ALK_index + alk_index1



  # Combine vectors
  for(i in vec_names){
    data_list_new[[i]] <- c(data_list1[[i]], data_list2[[i]])
  }

  # Combine matrices
  for(i in mat_names[1:17]){
    data_list_new[[i]] <- plyr::rbind.fill(data_list1[[i]], data_list2[[i]])
  }


  # Add new species
  data_list_new$nspp <- data_list2$nspp + data_list2$nspp

  return(data_list_new)
}
