#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB CEATTLE
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @export
rearrange_dat <- function(data_list){
  '%!in%' <- function(x,y)!('%in%'(x,y))

  # Step 1 - remove numeric objects from control
  data_list$ln_srv_q_prior <- log(data_list$fleet_control$Q_prior)

  data_list$fleet_control <- data_list$fleet_control[,-which(colnames(data_list$fleet_control) %in% c("Sel_sd_prior", "Q_prior", "Q_sd_prior", "Time_varying_q_sd_prior", "Survey_sd_prior", "proj_F", "Catch_sd_prior"))]

  data_list$fleet_control$Time_varying_sel <- round(data_list$fleet_control$Time_varying_sel)

  # Step 2 -  Seperate survey biomass info from observation
  data_list$srv_biom_ctl <- data_list$srv_biom[,c("Fleet_code", "Species", "Year")]
  data_list$srv_biom_n <- as.matrix(data_list$srv_biom[,c("Month")])
  data_list$srv_biom_obs <- data_list$srv_biom[,c("Observation", "Log_sd")]

  # Step 3 -  Seperate catch biomass info from observation
  data_list$fsh_biom_ctl <- data_list$fsh_biom[,c("Fleet_code", "Species", "Year")]
  data_list$fsh_biom_n <- as.matrix(data_list$fsh_biom[,c("Month")])
  data_list$fsh_biom_obs <- data_list$fsh_biom[,c("Catch", "Log_sd")]

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
  yrs <- data_list$styr:data_list$endyr
  if(nrow(data_list$emp_sel) > 0 ){
    for(i in 1:nrow(data_list$emp_sel)){
      # Fill in years
      if(data_list$emp_sel$Year[i] == 0){

        # Change first year
        data_list$emp_sel$Year[i] <- yrs[1]

        # Change the rest
        emp_sel_tmp <- data_list$emp_sel[i,]
        emp_sel_tmp$Year <- yrs[2]
        emp_sel <- emp_sel_tmp

        for(yr in 3:length(yrs)){
          emp_sel_tmp$Year <- yrs[yr]
          emp_sel <- rbind(emp_sel, emp_sel_tmp)
        }

        data_list$emp_sel <- rbind(data_list$emp_sel, emp_sel)
      }
    }
  }

  data_list$emp_sel_ctl <- as.matrix(data_list$emp_sel[,c("Fleet_code", "Species", "Sex", "Year")])
  data_list$emp_sel_obs <- matrix(as.numeric(unlist(data_list$emp_sel[,grep("Comp_", colnames(data_list$emp_sel))])), nrow = nrow(data_list$emp_sel_ctl))

  # Make data_list names different
  data_list$fleet_control$Fleet_name <- suppressWarnings(as.numeric(as.character(data_list$fleet_control$Fleet_name)))

  # Species names
  data_list$spnames <- NULL

  # Input missing acucmulation ages - defaults to age range
  for(i in 1:nrow(data_list$fleet_control)){

    if(is.na(data_list$fleet_control$Nselages[i])){
      data_list$fleet_control$Nselages[i] = -999
    }

    # Lower
    if(is.na(data_list$fleet_control$Accumatation_age_lower[i])){
      data_list$fleet_control$Accumatation_age_lower[i] <- data_list$minage[data_list$fleet_control$Species[i]]
    }

    # Upper
    if(is.na(data_list$fleet_control$Accumatation_age_upper[i])){
      data_list$fleet_control$Accumatation_age_upper[i] <- data_list$nages[data_list$fleet_control$Species[i]]
    }

    # Selected age
    if(is.na(data_list$fleet_control$Age_first_selected[i])){
      data_list$fleet_control$Age_first_selected[i] <- data_list$minage[data_list$fleet_control$Species[i]]
    }
  }

  # Rearrange age-transition matrix
  age_trans_matrix <- data_list$age_trans_matrix
  unique_age_transition <- unique(as.character(age_trans_matrix$Age_transition_index))
  for(i in 1:length(data_list$pop_age_transition_index)){
    if(data_list$pop_age_transition_index[i] %!in% unique_age_transition){
      stop("Check population age_transition index, not in age_transition file")
    }
  }
  age_transition <- array(0, dim = c(length(unique_age_transition), 2, max(data_list$nages, na.rm = T), max(data_list$nlengths, na.rm = T)))


  for (i in 1:nrow(age_trans_matrix)) {
    age_transition_ind <- as.numeric(as.character(age_trans_matrix$Age_transition_index[i]))
    sex <- as.numeric(as.character(age_trans_matrix$Sex[i]))
    sp <- as.numeric(as.character(age_trans_matrix$Species[i]))
    age <- as.numeric(as.character(age_trans_matrix$Age[i])) - data_list$minage[sp] + 1

    if (age > data_list$nages[sp]) {
      message(paste0("Error: number of ages in age_trans_matrix for species: ", sp, " and index: ", age_transition_ind))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    # Assign
    if(sex == 0){ sex = c(1, 2)}
    for(j in 1:length(sex)){
      age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]] <- as.numeric(as.character(age_trans_matrix[i, (1:data_list$nlengths[sp]) + 5]))

      # Normalize
      age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]] <- age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]] / sum(age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]], na.rm = TRUE)
    }
  }
  data_list$age_trans_matrix <- age_transition


  # Rearrange age_error matrices
  arm <- array(0, dim = c(data_list$nspp, max(data_list$nages, na.rm = T), max(data_list$nages, na.rm = T)))

  data_list$age_error <- as.data.frame(data_list$age_error) # FIXME: somethin is up with data.frames
  for (i in 1:nrow(data_list$age_error)) {
    sp <- as.numeric(as.character(data_list$age_error$Species[i]))
    true_age <- as.numeric(as.character(data_list$age_error$True_age[i])) - data_list$minage[sp] + 1

    if (true_age > data_list$nages[sp]) {
      message(paste0("Error: number of true ages specified in age_error for species: ", sp))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    arm[sp, true_age, 1:data_list$nages[sp]] <- as.numeric(as.character(data_list$age_error[i, (1:data_list$nages[sp]) + 2]))

    # Normalize
    arm[sp, true_age, 1:data_list$nages[sp]] <- arm[sp, true_age, 1:data_list$nages[sp]] / sum(arm[sp, true_age, 1:data_list$nages[sp]], na.rm = TRUE)
  }
  data_list$age_error <- arm



  # Normalize comp data
  for(i in 1:nrow(data_list$comp_obs)){
    data_list$comp_obs[i,] = as.numeric(data_list$comp_obs[i,]) / sum(as.numeric(data_list$comp_obs[i,]), na.rm = TRUE)
  }

  # Set up wt array
  wt_matrix <- data_list$wt
  unique_wt <- unique(as.character(wt_matrix$Wt_index))
  for(i in 1:length(data_list$pop_wt_index)){
    if(data_list$pop_wt_index[i] %!in% unique_wt){
      stop("Check population weight index, not in weight file")
    }
  }

  for(i in 1:length(data_list$ssb_wt_index)){
    if(data_list$ssb_wt_index[i] %!in% unique_wt){
      stop("Check SSB weight index, not in weight file")
    }
  }

  wt <- array(0, dim = c(length(unique_wt), 2, max(data_list$nages, na.rm = T), length(data_list$styr:data_list$endyr)))

  for (i in 1:nrow(wt_matrix)) {

    wt_ind <- as.numeric(as.character(wt_matrix$Wt_index[i]))
    sp <- as.numeric(as.character(wt_matrix$Species[i]))
    sex <- as.numeric(as.character(wt_matrix$Sex[i]))
    yr <- as.numeric(as.character(wt_matrix$Year[i])) - data_list$styr + 1

    if(yr <= data_list$endyr - data_list$styr + 1){

      # If year == 0, set weight to all years
      if(yr == (-data_list$styr + 1)){
        yr = 1:length(data_list$styr:data_list$endyr)
      }

      if(sex == 0){ sex = c(1, 2)}
      for(j in 1:length(sex)){
        if(sum(grepl("[[:space:]]", as.character(wt_matrix[i, (1:data_list$nages[sp]) + 5])))){
          stop(paste("Space found in wt data: row", i))
        }
        wt[wt_ind, sex[j], 1:data_list$nages[sp], yr] <- as.numeric(as.character(as.matrix(unlist(wt_matrix[i, (1:data_list$nages[sp]) + 5]))))
      }
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
      m1[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = T)) + 2])
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
      Mn_LatAge[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$Mn_LatAge[i,(1:max(data_list$nages, na.rm = T)) + 2])
    }
  }
  data_list$Mn_LatAge <- Mn_LatAge


  # Set up NByageFixed array
  NByageFixed <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T), length(data_list$styr:data_list$projyr)))

  if(nrow(data_list$NByageFixed) > 0){
    for (i in 1:nrow(data_list$NByageFixed)) {

      sp <- as.numeric(as.character(data_list$NByageFixed$Species[i]))
      sex <- as.numeric(as.character(data_list$NByageFixed$Sex[i]))
      yr <- as.numeric(as.character(data_list$NByageFixed$Year[i])) - data_list$styr + 1


      if(sex == 0){ sex = c(1, 2)}

      for(j in 1:length(sex)){
        if(yr > 0){
          NByageFixed[sp, sex[j], 1:max(data_list$nages, na.rm = T), yr] <- as.numeric(data_list$NByageFixed[i,c(1:max(data_list$nages, na.rm = T))+4])
        }
      }
    }
  }

  data_list$NByageFixed <- NByageFixed


  # Set up environmental indices
  data_list$env_yrs <-  data_list$env_data$Year
  env_cols <- ncol(data_list$env_data) - 1
  data_list$env_index <-  as.matrix(data_list$env_data[,2:ncol(data_list$env_data)],ncol = env_cols)


  # Set up pyrs array
  Pyrs <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages), length(data_list$styr:data_list$endyr))) # FIXME: Change for forecast

  for (i in 1:nrow(data_list$Pyrs)) {
    sp <- as.numeric(as.character(data_list$Pyrs$Species[i]))
    sex <- as.numeric(as.character(data_list$Pyrs$Sex[i]))
    if(sex == 0){
      sex = c(1,2)
    }
    yr <- as.numeric(as.character(data_list$Pyrs$Year[i])) - data_list$styr + 1

    if(yr <= (data_list$endyr - data_list$styr + 1)){
      Pyrs[sp, sex, 1:data_list$nages[sp], yr] <- as.numeric(as.character(as.matrix(unlist(data_list$Pyrs[i, (1:data_list$nages[sp]) + 3]))))
    }
  }
  data_list$Pyrs <- Pyrs

  # Remove species column from alw, pmature, sex_ratio
  data_list$sex_ratio <- data_list$sex_ratio[,-1]
  data_list$pmature <- data_list$pmature[,-1]
  data_list$aLW <- data_list$aLW[,-1]


  # Make data.frames into matrices
  for(i in 1:length(data_list)){
    if(class(data_list[[i]])[1] == "data.frame"){
      data_list[[i]] <- as.matrix(data_list[[i]])
    }
  }

  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "fsh_biom",    "srv_biom", "comp_data", "env_data")
  for(i in 1:length(items_to_remove)){
    data_list[[items_to_remove[i]]] <- NULL
  }

  return(data_list)
}
