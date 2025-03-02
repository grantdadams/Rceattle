data_check <- function(data_list) {

  if (length(data_list$BTempC) != length(data_list$Tyrs)) {
    stop("Length of temperature series is not the same as the number of temperature years")
  }

  if (length(data_list$M1_base) == 1) {
    stop("M1 is a single value, please make it age/species specific")
  }

  if (sum(data_list$other_food < 0) > 0) {
    stop("Other food for one species is negative")
  }

  for(sp in 1:data_list$nspp){
    if(sum(data_list$nages[sp] < data_list$fleet_control$Nselages[which(data_list$fleet_control$Species == sp)], na.rm = TRUE) > 1){
      stop(paste("Nselages is greater than nages for species", sp))
    }
  }

  if(data_list$styr < min(data_list$wt$Year, na.rm = TRUE)){
    stop("Weight data does not go back to the start year")
  }

  if(data_list$endyr > max(data_list$wt$Year, na.rm = TRUE)){
    stop("Weight data does not to the end year")
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
