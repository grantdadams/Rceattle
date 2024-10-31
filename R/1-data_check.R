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

  # # Age matrix
  #
  # if(ncol(data_list$NByageFixed) != max(data_list$nages, na.rm = T)+4){
  #   print(paste0("NByageFixed does not include all ages"))
  # }
  #
  # if(ncol(data_list$wt) != max(data_list$nages, na.rm = T)+4){
  #   stop(paste0("Weight-at-age (wt) does not include all ages"))
  # }

  # Check catch units are all the same if HCR is CMSY
  if(data_list$HCR == 1){
    if(length(unique(
      data_list$fleet_control$Weight1_Numbers2[which(data_list$fleet_control$Fleet_type == 1)])) > 1
    ){
      stop("HCR == 1 (CMSY), but catch units are not the same")
    }
  }
}
