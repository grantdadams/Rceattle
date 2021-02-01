#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB CEATTLE
#'
#' @param dat_list a data_list created from \code{\link{build_dat}}.
#' @export
rearrange_dat <- function(dat_list){

  # Step 1 - remove q priors
  dat_list$srv_control <- dat_list$srv_control[,-which(colnames(dat_list$srv_control) %in% c("Sel_sd_prior", "Log_q_prior", "Q_sd_prior", "Survey_sd_prior"))]

  # Step 2 -  Seperate survey biomass info from observation
  dat_list$srv_biom_ctl <- dat_list$srv_biom[,c("Survey_code", "Species", "Sex", "Year")]
  dat_list$srv_biom_n <- as.matrix(dat_list$srv_biom[,c("Month")])
  dat_list$srv_biom_obs <- dat_list$srv_biom[,c("Observation", "CV")]

  # Step 3 -  Seperate catch biomass info from observation
  dat_list$fsh_biom_ctl <- dat_list$fsh_biom[,c("Fishery_code", "Species", "Sex", "Year")]
  dat_list$fsh_biom_n <- as.matrix(dat_list$fsh_biom[,c("Month")])
  dat_list$fsh_biom_obs <- dat_list$fsh_biom[,c("Catch", "CV")]

  # Step 4 -  Seperate survey comp info from observation
  dat_list$srv_comp_ctl <- dat_list$srv_comp[,c("Survey_code", "Species", "Sex", "Age0_Length1", "Year")]
  dat_list$srv_comp_n <- dat_list$srv_comp[,c("Month", "Sample_size")]
  dat_list$srv_comp_obs <- dat_list$srv_comp[,grep("Comp_", colnames(dat_list$srv_comp))]

  # Step 5 -  Seperate catch comp info from observation
  dat_list$fsh_comp_ctl <- dat_list$fsh_comp[,c("Fishery_code", "Species","Sex", "Age0_Length1", "Year")]
  dat_list$fsh_comp_n <- dat_list$fsh_comp[,c("Month", "Sample_size")]
  dat_list$fsh_comp_obs <- dat_list$fsh_comp[,grep("Comp_", colnames(dat_list$fsh_comp))]

  # Remove first row of empirical selectivity if all NAs

  # Step 6 -  Seperate survey empirical selectivity info from observation
  dat_list$srv_emp_sel_ctl <- as.matrix(dat_list$srv_emp_sel[,c("Survey_code", "Species", "Year")])
  dat_list$srv_emp_sel_obs <- as.matrix(dat_list$srv_emp_sel[,grep("Comp_", colnames(dat_list$srv_emp_sel))])

  # Step 7 -  Seperate fishery empirical selectivity from observation
  dat_list$fsh_emp_sel_ctl <- as.matrix(dat_list$fsh_emp_sel[,c("Fishery_code", "Species", "Year")])
  dat_list$fsh_emp_sel_obs <- as.matrix(dat_list$fsh_emp_sel[,grep("Comp_", colnames(dat_list$fsh_emp_sel))])

  # Make data_list names different
  dat_list$fsh_control$Fishery_name <- suppressWarnings(as.numeric((dat_list$fsh_control$Fishery_name)))
  dat_list$srv_control$Survey_name <- suppressWarnings(as.numeric((dat_list$srv_control$Survey_name)))

  # projected fishery
  # dat_list$proj_F <- dat_list$fsh_control$proj_F
  dat_list$fsh_control <- dat_list$fsh_control[,-which(colnames(dat_list$fsh_control) %in% c("Sel_sd_prior", "proj_F", "Catch_sd_prior"))]

  # Species names
  dat_list$spnames <- NULL

  # Make data.frames into matrices
  for(i in 1:length(dat_list)){
    if(class(dat_list[[i]]) == "data.frame"){
      dat_list[[i]] <- as.matrix(dat_list[[i]])
    }
  }

  items_to_remove <- c("fsh_emp_sel", "srv_emp_sel",    "fsh_comp",    "srv_comp",    "fsh_biom",    "srv_biom")
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

  # Set up index matrix for survey selectivity

  return(dat_list)
}
