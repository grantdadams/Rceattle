#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB CEATTLE
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @export
rearrange_dat <- function(data_list){

  # Step 1 - remove numeric objects from control
  data_list$fleet_control <- data_list$fleet_control[,-which(colnames(data_list$fleet_control) %in% c("Sel_sd_prior", "Log_q_prior", "Q_sd_prior", "Survey_sd_prior", "Sel_sd_prior", "proj_F", "Catch_sd_prior"))]

  # Step 2 -  Seperate survey biomass info from observation
  data_list$srv_biom_ctl <- data_list$srv_biom[,c("Fleet_code", "Species", "Year")]
  data_list$srv_biom_n <- as.matrix(data_list$srv_biom[,c("Month")])
  data_list$srv_biom_obs <- data_list$srv_biom[,c("Observation", "CV")]

  # Step 3 -  Seperate catch biomass info from observation
  data_list$fsh_biom_ctl <- data_list$fsh_biom[,c("Fleet_code", "Species", "Year")]
  data_list$fsh_biom_n <- as.matrix(data_list$fsh_biom[,c("Month")])
  data_list$fsh_biom_obs <- data_list$fsh_biom[,c("Catch", "CV")]

  # Step 4 -  Seperate survey comp info from observation
  data_list$comp_ctl <- data_list$comp_data[,c("Fleet_code", "Species", "Sex", "Age0_Length1", "Year")]
  data_list$comp_n <- data_list$comp_data[,c("Month", "Sample_size")]
  data_list$comp_obs <- data_list$comp_data[,grep("Comp_", colnames(data_list$comp_data))]

  # Step 5 -  Seperate uobs info from observation
  data_list$UobsWtAge_ctl <- data_list$UobsWtAge[,c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age", "Year")]
  data_list$UobsWtAge <- data_list$UobsWtAge[,c("Sample_size", "Stomach_proportion_by_weight")]

  data_list$UobsAge_ctl <- data_list$UobsAge[,c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age", "Year")]
  data_list$UobsAge <- data_list$UobsAge[,c("Sample_size", "Stomach_proportion_by_number")]


  # Step 6 -  Seperate survey empirical selectivity info from observation
  data_list$emp_sel_ctl <- as.matrix(data_list$emp_sel[,c("Fleet_code", "Species", "Sex", "Year")])
  data_list$emp_sel_obs <- as.matrix(data_list$emp_sel[,grep("Comp_", colnames(data_list$emp_sel))])

  # Make data_list names different
  data_list$fleet_control$Fleet_name <- suppressWarnings(as.numeric(as.character(data_list$fleet_control$Fleet_name)))

  # Species names
  data_list$spnames <- NULL

  # Rearrange age-transition matrix
  age_trans_matrix <- data_list$age_trans_matrix
  unique_alk <- unique(as.character(age_trans_matrix$ALK_index))
  alk <- array(0, dim = c(length(unique_alk), 2, max(data_list$nages, na.rm = T), max(data_list$nlengths, na.rm = T)))


  for (i in 1:nrow(age_trans_matrix)) {
    alk_ind <- as.numeric(as.character(age_trans_matrix$ALK_index[i]))
    sex <- as.numeric(as.character(age_trans_matrix$Sex[i]))
    sp <- as.numeric(as.character(age_trans_matrix$Species[i]))
    age <- as.numeric(as.character(age_trans_matrix$Age[i])) - data_list$minage[sp] + 1

    if (age > data_list$nages[sp]) {
      message(paste0("Error: number of ages in age_trans_matrix for species: ", sp, " and index: ", alk_ind))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    # Assign
    if(sex == 0){ sex = c(1, 2)}
    for(j in 1:length(sex)){
      alk[alk_ind, sex[j], age, 1:data_list$nlengths[sp]] <- as.numeric(as.character(age_trans_matrix[i, (1:data_list$nlengths[sp]) + 5]))

      # Normalize
      alk[alk_ind, sex[j], age, 1:data_list$nlengths[sp]] <- alk[alk_ind, sex[j], age, 1:data_list$nlengths[sp]] / sum(alk[alk_ind, sex[j], age, 1:data_list$nlengths[sp]], na.rm = TRUE)
    }
  }
  data_list$age_trans_matrix <- alk



  # Normalize comp data
  for(i in 1:nrow(data_list$comp_obs)){
    data_list$comp_obs[i,] = data_list$comp_obs[i,] / sum(data_list$comp_obs[i,], na.rm = TRUE)
  }

  # Set up wt array
  wt_matrix <- data_list$wt
  unique_wt <- unique(as.character(wt_matrix$Wt_index))
  wt <- array(0, dim = c(length(unique_wt), 2, max(data_list$nages, na.rm = T), length(data_list$styr:data_list$endyr)))

  for (i in 1:nrow(wt_matrix)) {

    wt_ind <- as.numeric(as.character(wt_matrix$Wt_index[i]))
    sp <- as.numeric(as.character(wt_matrix$Species[i]))
    sex <- as.numeric(as.character(wt_matrix$Sex[i]))
    yr <- as.numeric(as.character(wt_matrix$Year[i])) - data_list$styr + 1
    if(sex == 0){ sex = c(1, 2)}
    for(j in 1:length(sex)){
      if(sum(grepl("[[:space:]]", as.character(wt_matrix[i, (1:data_list$nages[sp]) + 5])))){
        stop(paste("Space found in wt data: row", i))
      }
      wt[wt_ind, sex[j], 1:data_list$nages[sp], yr] <- as.numeric(as.character(wt_matrix[i, (1:data_list$nages[sp]) + 5]))
    }
  }
  data_list$wt <- wt

  # Set up M1 array
  m1 <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T)))

  for (i in 1:nrow(data_list$M1_base)) {
    sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
    sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))
    if(sex == 0){ sex = c(1, 2)}
    for(j in 1:length(sex)){
      m1[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$M1_base[i,-c(1,2)])
    }
  }
  data_list$M1_base <- m1


  # Set up Mn_LatAge array
  Mn_LatAge <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T)))

  for (i in 1:nrow(data_list$Mn_LatAge)) {
    sp <- as.numeric(as.character(data_list$Mn_LatAge$Species[i]))
    sex <- as.numeric(as.character(data_list$Mn_LatAge$Sex[i]))
    if(sex == 0){ sex = c(1, 2)}
    for(j in 1:length(sex)){
      Mn_LatAge[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$Mn_LatAge[i,-c(1,2)])
    }
  }
  data_list$Mn_LatAge <- Mn_LatAge

  # Make data.frames into matrices
  for(i in 1:length(data_list)){
    if(class(data_list[[i]]) == "data.frame"){
      data_list[[i]] <- as.matrix(data_list[[i]])
    }
  }

  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "fsh_biom",    "srv_biom", "comp_data")
  for(i in 1:length(items_to_remove)){
    data_list[[items_to_remove[i]]] <- NULL
  }

  return(data_list)
}
