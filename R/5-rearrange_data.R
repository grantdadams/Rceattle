#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB CEATTLE
#'
#' @param dat_list a data_list created from \code{\link{build_dat}}.
#' @export
rearrange_dat <- function(dat_list){

  # Step 1 - remove numeric objects from control
  dat_list$fleet_control <- dat_list$fleet_control[,-which(colnames(dat_list$fleet_control) %in% c("Sel_sd_prior", "Log_q_prior", "Q_sd_prior", "Survey_sd_prior", "Sel_sd_prior", "proj_F", "Catch_sd_prior"))]

  # Step 2 -  Seperate survey biomass info from observation
  dat_list$srv_biom_ctl <- dat_list$srv_biom[,c("Fleet_code", "Species", "Year")]
  dat_list$srv_biom_n <- as.matrix(dat_list$srv_biom[,c("Month")])
  dat_list$srv_biom_obs <- dat_list$srv_biom[,c("Observation", "CV")]

  # Step 3 -  Seperate catch biomass info from observation
  dat_list$fsh_biom_ctl <- dat_list$fsh_biom[,c("Fleet_code", "Species", "Year")]
  dat_list$fsh_biom_n <- as.matrix(dat_list$fsh_biom[,c("Month")])
  dat_list$fsh_biom_obs <- dat_list$fsh_biom[,c("Catch", "CV")]

  # Step 4 -  Seperate survey comp info from observation
  dat_list$comp_ctl <- dat_list$comp_data[,c("Fleet_code", "Species", "Sex", "Age0_Length1", "Year")]
  dat_list$comp_n <- dat_list$comp_data[,c("Month", "Sample_size")]
  dat_list$comp_obs <- dat_list$comp_data[,grep("Comp_", colnames(dat_list$comp_data))]

  # Step 6 -  Seperate survey empirical selectivity info from observation
  dat_list$emp_sel_ctl <- as.matrix(dat_list$emp_sel[,c("Fleet_code", "Species", "Year", "Sex")])
  dat_list$emp_sel_obs <- as.matrix(dat_list$emp_sel[,grep("Comp_", colnames(dat_list$emp_sel))])

  # Make data_list names different
  dat_list$fleet_control$Fleet_name <- suppressWarnings(as.numeric((dat_list$fleet_control$Fleet_name)))

  # Species names
  dat_list$spnames <- NULL

  # Make data.frames into matrices
  for(i in 1:length(dat_list)){
    if(class(dat_list[[i]]) == "data.frame"){
      dat_list[[i]] <- as.matrix(dat_list[[i]])
    }
  }

  items_to_remove <- c("emp_sel",    "fsh_comp",    "srv_comp",    "fsh_biom",    "srv_biom")
  for(i in 1:length(items_to_remove)){
    dat_list[[items_to_remove[i]]] <- NULL
  }

  # Normalize age-transition matrix
  for(i in 1:dim(dat_list$age_trans_matrix)[3]){
    dat_list$age_trans_matrix[,,i] = dat_list$age_trans_matrix[,,i] / rowSums(dat_list$age_trans_matrix[,,i], na.rm = T)
  }


  # Normalize srv comp
  for(i in nrow(dat_list$srv_comp_obs)){
    dat_list$srv_comp_obs[i,] = dat_list$srv_comp_obs[i,] / sum(dat_list$srv_comp_obs[i,], na.rm = TRUE)
  }

  # Normalize fsh comp
  for(i in nrow(dat_list$fsh_comp_obs)){
    dat_list$fsh_comp_obs[i,] = dat_list$fsh_comp_obs[i,] / sum(dat_list$fsh_comp_obs[i,], na.rm = TRUE)
  }

  # Set up wt array
  wt_matrix <- dat_list$wt
  unique_wt <- unique(as.character(wt_matrix$Wt_index))
  wt <- array(0, dim = c(length(unique_wt), 2, max(data_list$nages, na.rm = T), length(data_list$styr:data_list$endyr)))

  for (i in 1:nrow(wt_matrix)) {

    wt_ind <- as.numeric(as.character(wt_matrix$Wt_index[i]))
    sp <- as.numeric(as.character(wt_matrix$Species[i]))
    sex <- as.numeric(as.character(wt_matrix$Sex[i]))
    yr <- as.numeric(as.character(wt_matrix$Year[i])) - data_list$styr + 1
    if(sex == 0){ sex = c(1, 2)}

    wt[wt_ind, sex, 1:data_list$nages[sp], yr] <- as.numeric(as.character(wt_matrix[i, (1:data_list$nages[sp]) + 5]))
  }
  dat_list$wt <- wt

  # Set up M1 array

  return(dat_list)
}
