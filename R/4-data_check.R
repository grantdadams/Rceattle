#' Function to clean data for Rceattle runs
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
clean_data <- function(data_list){

  # Transpose fleet control
  data_list <- Rceattle::transpose_fleet_control(data_list)

  # Remove years of data previous to start year and after end year
  # - Data in likelihood
  data_list$index_data <- data_list$index_data %>%
    dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
  data_list$catch_data <- data_list$catch_data %>%
    dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
  data_list$comp_data <- data_list$comp_data %>%
    dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
  data_list$diet_data <- data_list$diet_data %>%
    dplyr::filter(Year >= data_list$styr & Year <= data_list$projyr | Year == 0) %>%
    dplyr::arrange(Pred, Pred_sex, Pred_age, Prey, Prey_sex, Prey_age, Year)

  # - Fixed data
  data_list$weight <- data_list$weight %>%
    dplyr::filter(Year >= data_list$styr & Year <= data_list$projyr | Year == 0)
  data_list$diet_data <- as.data.frame(data_list$diet_data) %>%
    dplyr::filter(Year >= data_list$styr & Year <= data_list$projyr | Year == 0)
  data_list$emp_sel <- data_list$emp_sel %>%
    dplyr::filter(Year >= data_list$styr & Year <= data_list$projyr | Year == 0)
  data_list$NByageFixed <- data_list$NByageFixed %>%
    dplyr::filter(Year >= data_list$styr & Year <= data_list$projyr | Year == 0)
  data_list$ration_data <- data_list$ration_data %>%
    dplyr::filter(Year >= data_list$styr & Year <= data_list$projyr | Year == 0)


  # Add temp multi-species SB0
  if(is.null(data_list$MSSB0)){
    data_list$MSSB0 <- rep(999, data_list$nspp)
    data_list$MSB0 <- rep(999, data_list$nspp)
  }


  # Extend catch data to proj year for projections (where data are missing)
  if(data_list$projyr > data_list$endyr){
    for(flt in (unique(data_list$catch_data$Fleet_code))){
      catch_data_sub <- data_list$catch_data %>%
        dplyr::filter(Fleet_code == flt)
      yrs_proj <- (data_list$endyr + 1):data_list$projyr
      yrs_proj <- yrs_proj[which(!yrs_proj %in% catch_data_sub$Year)]
      nyrs_proj <- length(yrs_proj)
      proj_catch_data <- data.frame(Fleet_name = rep(catch_data_sub$Fleet_name[1], nyrs_proj),
                                    Fleet_code = rep(flt, nyrs_proj),
                                    Species = rep(catch_data_sub$Species[1], nyrs_proj),
                                    Year = yrs_proj,
                                    Month = rep(catch_data_sub$Month[length(catch_data_sub$Month)], nyrs_proj),
                                    Selectivity_block = rep(catch_data_sub$Selectivity_block[length(catch_data_sub$Selectivity_block)], nyrs_proj),
                                    Catch = rep(NA, nyrs_proj),
                                    Log_sd = rep(catch_data_sub$Log_sd[length(catch_data_sub$Log_sd)], nyrs_proj))
      data_list$catch_data <- rbind(data_list$catch_data, proj_catch_data)
    }
  }

  data_list$catch_data <- data_list$catch_data[
    with(data_list$catch_data, order(Fleet_code, Year)),]

  return(data_list)
}


#' Function to check data for errors and update formatting where necessary
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
data_check <- function(data_list) {

  # Species checks
  if(data_list$nspp != max(data_list$weight$Species)){
    stop("`nspp` does not match the number of species in the weight data. Check `nspp` or `weight`")
  }

  if(length(data_list$spnames) != data_list$nspp){
    stop("species names not included for all species")
  }

  if (length(data_list$M1_base) == 1) {
    stop("M1 is a single value, please make it age/species specific")
  }

  if (sum(data_list$other_food < 0) > 0) {
    stop("Other food for one species is negative")
  }

  # Species checks ----
  for(sp in 1:data_list$nspp){
    if(sum(data_list$nages[sp] < data_list$fleet_control$Nselages[which(data_list$fleet_control$Species == sp)], na.rm = TRUE) > 1){
      stop(paste("Nselages is greater than nages for species", sp))
    }
  }

  # Fleet checks ----
  for(flt in 1:nrow(data_list$fleet_control)){
    if(!is.na(data_list$fleet_control$Estimate_q[flt])){
      if((data_list$fleet_control$Estimate_q[flt] == 6 & data_list$fleet_control$Time_varying_q[flt] > (ncol(data_list$env_data) - 1))|
         (data_list$fleet_control$Estimate_q[flt] == 6 & data_list$fleet_control$Time_varying_q[flt] < 1)){
        stop("For catchability type 6 environmental index specified in 'Time_varying_q' is greater than number of indices in 'env_data'")
      }
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
    warning(paste0("Selectivity for ", paste(mirror_sel$Fleet_name, collapse = ", "), " is mirrored with another fleet"))
  }

  mirror_q <- data_list$fleet_control %>%
    dplyr::filter(!is.na(Estimate_q)) %>%
    dplyr::group_by(Q_index) %>%
    dplyr::filter(n() > 1 ) %>%
    dplyr::ungroup()
  if(nrow(mirror_q) > 0){
    warning(paste0("Catchability for ", paste(mirror_q$Fleet_name, collapse = ", "), " is mirrored with another fleet"))
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
      warning(paste("`age_trans_matrix` data does not span range of age for species", sp, "will fill with 0s"))
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
      warning(paste("`age_error` data does not span range of true ages for species", sp, "will fill with 0s"))
    }
  }

  # # Age matrix
  #
  # if(ncol(data_list$NByageFixed) != max(data_list$nages, na.rm = T)+4){
  #   print(paste0("NByageFixed does not include all ages"))
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
    warning("No ration_data (ration) data")
  }

  # Environmental data----
  if(any(data_list$srr_indices > ncol(data_list$env_index))){stop("'srr_indices' greater than the number of indices included")}
  if(any(data_list$M1_indices > ncol(data_list$env_index))){stop("'M1_indices' greater than the number of indices included")}
}


#' Function to check for missing switches for map and parameter functions
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
switch_check <- function(data_list){

  # Estimation
  if(is.null(data_list$estDynamics)){
    data_list$estDynamics <- rep(0, data_list$nspp)
    print("'estDynamics' are not included in data, assuming 0")
  }

  # Stock-recruitment
  if(is.null(data_list$srr_fun)){
    data_list$srr_fun <- 0
    data_list$srr_pred_fun <- 0
    data_list$srr_est_mode <- 0
    print("'srr_fun' are not included in data, assuming 0")
  }

  # Diet data weights
  if(is.null(data_list$Diet_comp_weights)){
    data_list$Diet_comp_weights <- rep(1, data_list$nspp)
    print("'Diet_comp_weights' are not included in data, assuming 1")
  }

  # Normalization age
  if(is.null(data_list$fleet_control$Sel_norm_bin1)){
    data_list$fleet_control$Sel_norm_bin1 <- NA
    print("'Sel_norm_bin1' not specified in 'fleet_control', assuming 'NA'")
  }


  if(is.null(data_list$fleet_control$Sel_norm_bin2)){
    data_list$fleet_control$Sel_norm_bin2 <- NA
    print("'Sel_norm_bin2' not specified in 'fleet_control', assuming 'NA'")
  }

  # Sel curve penalties
  if(is.null(data_list$fleet_control$Sel_curve_pen1)){
    data_list$fleet_control$Sel_curve_pen1 <- 0
    print("'Sel_curve_pen1' not specified in 'fleet_control', assuming '0'")
  }

  if(is.null(data_list$fleet_control$Sel_curve_pen2)){
    data_list$fleet_control$Sel_curve_pen2 <- 0
    print("'Sel_curve_pen2' not specified in 'fleet_control', assuming '0'")
  }

  if(any(data_list$fleet_control$Selectivity == 2 & !is.na(data_list$fleet_control$Time_varying_sel) & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1)))){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::mutate(Sel_curve_pen1 = ifelse(Selectivity == 2 & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1)), Time_varying_sel, NA),
                    Sel_curve_pen2 = ifelse(Selectivity == 2 & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1)), Sel_sd_prior, NA),
                    Time_varying_sel = ifelse(Selectivity == 2 & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1)), 0, Time_varying_sel),
                    Sel_sd_prior = ifelse(Selectivity == 2 & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1)), 0, Sel_sd_prior)

      )
    print("Updating format where 'Selectivity == 2'. Moving non-parametric penalties to 'Sel_curve_pen1' and 'Sel_curve_pen2'.")
  }


  if(any(data_list$fleet_control$Selectivity == 2 & is.na(data_list$fleet_control$Sel_curve_pen1))){
    stop("'Sel_curve_pen1' is NA in 'fleet_control' for fleet with non-parametric selectivity ('Selectivity = 2')")
  }

  if(any(data_list$fleet_control$Selectivity == 2 & is.na(data_list$fleet_control$Sel_curve_pen2))){
    stop("'Sel_curve_pen2' is NA in 'fleet_control' for fleet with non-parametric selectivity ('Selectivity = 2')")
  }

  # Comp weights
  if(is.null(data_list$fleet_control$Comp_loglike)){
    data_list$fleet_control$Comp_loglike <- -1
    print("'Comp_loglike' not specified in 'fleet_control', assuming multinomial")
  }

  # Mortality
  if(is.null(data_list$M1_model)){
    data_list$M1_model <- rep(0, data_list$nspp)
    print("'M1_model' are not included in data, assuming 0")
  }

  if(is.null(data_list$msmMode)){
    data_list$msmMode <- 0
    print("'msmMode' are not included in data, assuming 0")
  }

  if(is.null(data_list$M1_re)){
    data_list$M1_re <- rep(0, data_list$nspp)
    print("'M1_re' is not in data, assuming 0 for all species")
  }

  # Init mode
  if(is.null(data_list$initMode)){
    data_list$initMode = 2
    print("'initMode' not input, setting to 2 (default)")
  }

  return(data_list)
}
