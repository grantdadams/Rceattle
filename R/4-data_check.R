#' Function to check data for errors. Does not update the data set!
#'
#' @param data_list Rceattle data list
#'
data_check <- function(data_list) {

  # --- 1. Base Species checks ----
  if(data_list$nspp != max(data_list$weight$Species)) stop("`nspp` does not match the number of species in the weight data. Check `nspp` or `weight`")
  if(length(data_list$spnames) != data_list$nspp) stop("species names not included for all species")
  if(length(data_list$spawn_month) != data_list$nspp) stop("'spawn_month' not included for all species")
  if(length(data_list$nages) != data_list$nspp) stop("'nages' not included for all species")
  if(length(data_list$nlengths) != data_list$nspp) stop("'nlengths' not included for all species")
  if(length(data_list$M1_base) == 1) stop("M1 is a single value, please make it age/species specific")
  if(length(data_list$other_food) != data_list$nspp) stop("'other_food' not included for all species")
  if(sum(data_list$other_food < 0) > 0) stop("Other food for one species is negative")

  # Species checks ----
  for(sp in 1:data_list$nspp){
    if(sum(data_list$nages[sp] < data_list$fleet_control$N_sel_bins[which(data_list$fleet_control$Species == sp)], na.rm = TRUE) > 1){
      stop(paste("N_sel_bins is greater than nages for species", sp))
    }
  }

  # Fleet checks ----
  for(flt in 1:nrow(data_list$fleet_control)){
    if(!is.na(data_list$fleet_control$Catchability[flt])){
      if((data_list$fleet_control$Catchability[flt] == "AR1" & data_list$fleet_control$Time_varying_q[flt] > (ncol(data_list$env_data) - 1))|
         (data_list$fleet_control$Catchability[flt] == "AR1" & data_list$fleet_control$Time_varying_q[flt] < 1)){
        stop("For catchability type 6 environmental index specified in 'Time_varying_q' is greater than number of indices in 'env_data'")
      }
    }

    # Time-varying sel not possible
    if(data_list$fleet_control$Time_varying_sel[flt] == 2 & data_list$fleet_control$Selectivity[flt] %in% c("NonParametric", "Hake")){
      stop("For non-parametric selectivities, 'Time_varying_sel' cant not be a random walk")
    }

    # Max sel age > nages
    data_list$fleet_control$Sel_norm_bin1[flt] <- ifelse(data_list$fleet_control$Sel_norm_bin1[flt] > data_list$nages[data_list$fleet_control$Species[flt]], data_list$nages[data_list$fleet_control$Species[flt]], data_list$fleet_control$Sel_norm_bin1[flt])
  }

  # - Mirroring warnings
  mirror_sel <- data_list$fleet_control %>%
    dplyr::group_by(Selectivity_index) %>%
    dplyr::filter(n() > 1 ) %>%
    dplyr::ungroup()
  if(nrow(mirror_sel) > 0){
    message(paste0("Selectivity for ", paste(mirror_sel$Fleet_name, collapse = ", "), " is mirrored with another fleet"))
  }

  mirror_q <- data_list$fleet_control %>%
    dplyr::filter(!is.na(Catchability)) %>%
    dplyr::group_by(Q_index) %>%
    dplyr::filter(n() > 1 ) %>%
    dplyr::ungroup()
  if(nrow(mirror_q) > 0){
    message(paste0("Catchability for ", paste(mirror_q$Fleet_name, collapse = ", "), " is mirrored with another fleet"))
  }


  # * Validate selectivity and catchability inputs ----
  # Allowed inputs include both the old integer codes and the new string names
  valid_sel <- c(0:7, "Fixed", "Logistic", "NonParametric", "DoubleLogistic", "DescendingLogistic", "Hake", "2DAR1", "3DAR1")
  valid_q <- c(0:6, "Fixed", "Estimated", "Estimated-with-prior", "Analytical", "Environmental", "AR1")

  invalid_sel <- data_list$fleet_control %>%
    dplyr::filter(!Selectivity %in% valid_sel)

  invalid_q <- data_list$fleet_control %>%
    dplyr::filter(!Catchability %in% valid_q)

  # Throw clear errors to guide the user
  if(nrow(invalid_sel) > 0) {
    stop(paste("Invalid Selectivity specified for fleets:",
               paste(invalid_sel$Fleet_name, collapse = ", "),
               ".\nPlease use an integer code 0:7 or one of:",
               paste(c("Empirical", "Logistic", "Non-parametric", "Double-Logistic", "Descending-Logistic", "Hake-Non-parametric"), collapse = ", ")))
  }

  if(nrow(invalid_q) > 0) {
    stop(paste("Invalid Catchability specified for fleets:",
               paste(invalid_q$Fleet_name, collapse = ", "),
               ".\nPlease use an integer code 0:6 or one of:",
               paste(c("Fixed", "Estimated", "Estimated-with-prior", "Analytical", "Environmental", "AR1-Deviates"), collapse = ", ")))
  }



  # Weight-at-age ----
  # * Year range ----
  wt_yr <- data_list$weight %>%
    dplyr::group_by(Wt_index, Sex) %>%
    dplyr::distinct(Year) %>%
    dplyr::mutate(Tmp_ind = paste0("index = ", Wt_index," & sex = ", Sex))

  for(ind in unique(wt_yr$Tmp_ind)){

    tmp_wt <- wt_yr %>%
      dplyr::filter(Tmp_ind == ind) %>%
      dplyr::distinct(Year) %>%
      dplyr::pull(Year)

    # If not time-varying
    if(length(tmp_wt) > 1){

      # If time-varying, check weight spans range
      if(any(!(data_list$styr:data_list$endyr) %in% tmp_wt)){
        stop(paste0("Weight data for ", ind, " does not span all hindcast years"))
      }
    }
  }


  # * Index checks ----
  wt_index <- wt_index <- data_list$weight %>%
    dplyr::distinct(Wt_index, Species, Sex)

  # - Data checks ----
  if(any(!data_list$pop_wt_index %in% wt_index$Wt_index)){
    stop("Check population weight index, not in weight file")
  }

  if(any(!data_list$ssb_wt_index %in% wt_index$Wt_index)){
    stop("Check SSB weight index, not in weight file")
  }

  for(sp in 1:data_list$nspp){
    wt_no <- data_list$weight %>%
      dplyr::filter(Species != sp) %>% # Not species
      dplyr::distinct(Wt_index) %>%
      pull(Wt_index)

    wt_yes <- data_list$weight %>%
      dplyr::filter(Species == sp) %>% # Is species
      dplyr::distinct(Wt_index) %>%
      pull(Wt_index)

    if(any( wt_no %in% wt_yes )){
      stop("Check weight indices (Wt_index), the same weight index was used for multiple species")
    }
  }

  if(any(data_list$weight %>%
         dplyr::select(-c(Wt_name, Wt_index, Species, Sex, Year)) %>%
         ncol() < data_list$nages)){
    stop("Weight data does not span range of ages")
  }


  # Biological data ----
  if(ncol(data_list$maturity) < max(data_list$nages)){
    stop("Maturity-at-age (maturity) does not span all ages")
  }

  if(ncol(data_list$sex_ratio) < max(data_list$nages)){
    stop("Sex ratio does not span all ages")
  }


  # ration_data ----
  if(nrow(data_list$ration_data) > 0){
    if(any(data_list$ration_data %>%
           dplyr::select(-c(Species, Sex, Year)) %>%
           ncol() < data_list$nages)){
      stop("'ration_data' data does not span range of ages")
    }
  }



  # Age transition matrix ----
  if(any(data_list$age_trans_matrix %>%
         dplyr::select(-c(Age_transition_name, Age_transition_index, Species, Sex, Age)) %>%
         ncol() < data_list$nlengths)){
    stop("`age_trans_matrix` data does not span range of lengths")
  }

  for(sp in 1:data_list$nspp){
    ages_tmp <- data_list$age_trans_matrix %>%
      as.data.frame() %>%
      dplyr::filter(Species == sp) %>%
      dplyr::pull(Age)
    if(!all(data_list$minage[sp]:data_list$nages[sp] %in% ages_tmp)){
      message(paste("`age_trans_matrix` data does not span range of age for species", sp, "will fill with 0s"))
    }
  }


  # Age error matrix ----
  if(any(data_list$age_error %>%
         as.data.frame() %>%
         dplyr::select(-c(Species, True_age)) %>%
         ncol() < data_list$nages)){
    stop("`age_error` observed ages do not span range of ages")
  }


  for(sp in 1:data_list$nspp){
    ages_tmp <- data_list$age_error %>%
      as.data.frame() %>%
      dplyr::filter(Species == sp) %>%
      dplyr::pull(True_age)
    if(!all(data_list$minage[sp]:data_list$nages[sp] %in% ages_tmp)){
      message(paste("`age_error` data does not span range of true ages for species", sp, "will fill with 0s"))
    }
  }

  # # Age matrix
  #
  # if(ncol(data_list$NByageFixed) != max(data_list$nages, na.rm = T)+4){
  #   message(paste0("NByageFixed does not include all ages"))
  # }
  #
  # if(ncol(data_list$weight) != max(data_list$nages, na.rm = T)+4){
  #   stop(paste0("Weight-at-age (weight) does not include all ages"))
  # }

  # Switches ---
  if(any(data_list$suitMode %in% c(1, 3))){
    stop("Length based suitability not yet implemented")
  }

  if(sum(data_list$fleet_control$proj_F_prop, na.rm = TRUE) == 0 & data_list$HCR > 0){
    stop("HCR is > 0 and 'proj_F_prop' is 0")
  }


  # Diet data ----
  # - Pred age
  if(nrow(data_list$diet_data) > 0){
    Max_age = data_list$diet_data %>%
      dplyr::group_by(Pred) %>%
      dplyr::summarise(Max_age = max(Pred_age)) %>%
      dplyr::arrange(Pred)

    if(any(Max_age$Max_age > data_list$nages)){
      stop("Pred ages in 'diet_data' > 'nages'")
    }

    if(sum(duplicated(data_list$diet_data)) > 0){
      stop("Diet data includes duplicated rows")
    }

    # - Prey age
    Max_age = data_list$diet_data %>%
      dplyr::group_by(Prey) %>%
      dplyr::summarise(Max_age = max(Prey_age)) %>%
      dplyr::arrange(Prey)

    if(any(Max_age$Max_age > data_list$nages)){
      stop("Prey ages in 'diet_data' > 'nages'")
    }

    # - Stomach proportion > 1
    diet_sum = data_list$diet_data %>%
      dplyr::group_by(Pred, Pred_age, Pred_sex, Year) %>%
      dplyr::summarise(diet_sum = sum(Stomach_proportion_by_weight))

    if(any(diet_sum$diet_sum > 1)){
      stop("Stomach proportion in `diet_data` for some predators-at-age/sex/year is > 1")
    }
  } else {
    if(data_list$msmMode > 0){
      stop("No diet data included")
    }
  }

  # Diet composition weights
  if(any(is.na(data_list$Diet_comp_weights) & data_list$suitMode > 0)){
    stop("Diet composition likelihood weight for a species with estimated suitability is NA")
  }

  if(any(!(data_list$Diet_loglike %in% c(0, 1)) & data_list$suitMode > 0)){
    stop("Diet composition likelihood for a species with estimated suitability is not 0 or 1")
  }

  # Sexes ----
  m1_sex <- data_list$M1_base %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(max_sex = max(Sex)) %>%
    dplyr::arrange(Species)

  if(any(m1_sex$max_sex > data_list$nsex)){
    stop("'M1_base' has more sexes than specified in 'nsex'")
  }

  wt_sex <- data_list$weight %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(max_sex = max(Sex)) %>%
    dplyr::arrange(Species)

  if(any(wt_sex$max_sex > data_list$nsex)){
    stop("'weight' has more sexes than specified in 'nsex'")
  }

  if(nrow(data_list$ration_data) > 0){
    ration_data_sex <- data_list$ration_data %>%
      dplyr::group_by(Species) %>%
      dplyr::summarise(max_sex = max(Sex)) %>%
      dplyr::arrange(Species)

    if(any(ration_data_sex$max_sex > data_list$nsex)){
      stop("'ration_data' has more sexes than specified in 'nsex'")
    }
  } else{
    message("No ration_data (ration) data")
  }

  # Environmental data----
  if(any(data_list$srr_indices > ncol(data_list$env_index))){stop("'srr_indices' greater than the number of indices included")}
  if(any(data_list$M1_indices > ncol(data_list$env_index))){stop("'M1_indices' greater than the number of indices included")}
}
