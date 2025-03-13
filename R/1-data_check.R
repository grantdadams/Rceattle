data_check <- function(data_list) {

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
    data_list$fleet_control$Age_max_selected[flt] <- ifelse(data_list$fleet_control$Age_max_selected[flt] > data_list$nages[data_list$fleet_control$Species[flt]], data_list$nages[data_list$fleet_control$Species[flt]], data_list$fleet_control$Age_max_selected[flt])
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
  if(any(data_list$pop_wt_index %!in% wt_index$Wt_index)){
    stop("Check population weight index, not in weight file")
  }

  if(any(data_list$ssb_wt_index %!in% wt_index$Wt_index)){
    stop("Check SSB weight index, not in weight file")
  }

  if(any(wt_index %>% dplyr::count(Wt_index) %>% dplyr::pull(n) > 1)){
    stop("Check weight indices (Wt_index), the same weight index was used for multiple species")
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

  # # Age matrix
  #
  # if(ncol(data_list$NByageFixed) != max(data_list$nages, na.rm = T)+4){
  #   print(paste0("NByageFixed does not include all ages"))
  # }
  #
  # if(ncol(data_list$weight) != max(data_list$nages, na.rm = T)+4){
  #   stop(paste0("Weight-at-age (weight) does not include all ages"))
  # }

}
