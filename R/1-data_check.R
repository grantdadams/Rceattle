data_check <- function(data_list) {

  if (length(data_list$M1_base) == 1) {
    stop("M1 is a single value, please make it age/species specific")
  }

  if (sum(data_list$other_food < 0) > 0) {
    stop("Other food for one species is negative")
  }

  # Species checks
  for(sp in 1:data_list$nspp){
    if(sum(data_list$nages[sp] < data_list$fleet_control$Nselages[which(data_list$fleet_control$Species == sp)], na.rm = TRUE) > 1){
      stop(paste("Nselages is greater than nages for species", sp))
    }
  }

  # Fleet checks
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

  # Weight at age check
  if(data_list$styr < min(data_list$wt$Year, na.rm = TRUE)){
    stop("Weight data does not go back to the start year")
  }

  if(data_list$endyr > max(data_list$wt$Year, na.rm = TRUE)){
    stop("Weight data does not to the end year")
  }
  

  # Wt data
  data_list$wt <- data_list$wt %>%
    mutate(
      Wt_index = as.numeric(as.character(Wt_index)),
      Species = as.numeric(as.character(Species)),
      Sex = as.numeric(as.character(Sex)),
      Year = as.numeric(as.character(Year)) - data_list$styr + 1)

  unique_wt <- unique(as.numeric(data_list$wt$Wt_index))

  # - Data checks ----
  if(any(data_list$pop_wt_index %!in% unique_wt)){
    stop("Check population weight index, not in weight file")
  }

  if(any(data_list$ssb_wt_index %!in% unique_wt)){
    stop("Check SSB weight index, not in weight file")
  }

  if(any(data_list$wt %>%
         dplyr::select(-c(Wt_name, Wt_index, Species, Sex, Year)) %>%
         ncol() < data_list$nages)){
    stop("Weight data does not span range of ages")
  }

  # # Age matrix
  #
  # if(ncol(data_list$NByageFixed) != max(data_list$nages, na.rm = T)+4){
  #   print(paste0("NByageFixed does not include all ages"))
  # }
  #
  # if(ncol(data_list$wt) != max(data_list$nages, na.rm = T)+4){
  #   stop(paste0("Weight-at-age (wt) does not include all ages"))
  # }

}
