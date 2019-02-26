#' Function to construct the TMB map argument for CEATTLE
#'
#' @description Reads a parameter list to construct map
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param params a parameter list created from \code{\link{build_param}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA except dummy
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE, which sets the map for ln_rec_sigma to NA
#'
#' @return a list of map arguments for each parameter
#' @export
build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE) {

  map_list <- params

  # STEP 1 -- Convert map_list to seperate parameters
  for(i in 1:length(map_list)){
    map_list[[i]] <- replace(map_list[[i]], values = c(1:length(map_list[[i]])))
  }


  # STEP 2 -- NA out parameters not to be estimated
  # Initial population deviations - map out last age and ages not seen
  for(i in 1:nrow(map_list$init_dev)){
    if((data_list$nages[i]-1) < ncol(map_list$init_dev)){
      map_list$init_dev[i,  (data_list$nages[i]) : ncol(map_list$init_dev) ] <- NA
    }
  }

  # Survey selectivity coefficients
  for( i in 1: nrow(data_list$srv_control)){
    if(data_list$srv_control$Selectivity[i] == 0){ # Empirical

      # Map out non-parametric
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))

      # Map out logistic and double logistic
      map_list$srv_sel_slp[1:2, i] <- NA
      map_list$srv_sel_inf[1:2, i] <- NA
    }
    if(data_list$srv_control$Selectivity[i] == 1){ # Logitistic

      # Map out non-parametric
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))

      # Map out double logistic
      map_list$srv_sel_slp[2, i] <- NA
      map_list$srv_sel_inf[2,i] <- NA
    }
    if(data_list$srv_control$Selectivity[i] == 2){ # Non-parametric at age
      map_list$srv_sel_slp[1:2, i] <- NA
      map_list$srv_sel_inf[1:2, i] <- NA

      # If nselages is  < max(nselages)
      if(data_list$nselages[i] < max(data_list$nselages, na.rm = TRUE)){
        map_list$srv_sel_coff[i, (data_list$nselages[i] + 1):max(data_list$nselages)]  <- replace(map_list$srv_sel_coff[i, (data_list$nselages[i] + 1):max(data_list$nselages)], values = rep(NA, length(map_list$srv_sel_coff[i, (data_list$nselages[i] + 1):max(data_list$nselages)])))
      }
    }
    if(data_list$srv_control$Selectivity[i] == 3){ # Double logistic
      # Map out non-parametric
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))
    }
  }

  # Catchability of surveys
  for( i in 1: nrow(data_list$srv_control)){
    # If not estimating turn of
    if(data_list$srv_control$Estimate_q[i] == 0){
      map_list$log_srv_q[i] <- NA
    }
  }


  # Recruitment deviation sigmas - turn off if not estimating
  if(random_rec == FALSE){
    map_list$ln_rec_sigma <- replace(map_list$ln_rec_sigma, values = rep(NA, length(map_list$ln_rec_sigma)))
  }

  ######################################################
  ####### Predation bits
  ######################################################
  # 1. Turn off all predation parameters for single species
  if(data_list$msmMode == 0){

    # Suitability parameters
    map_list$log_gam_a <- replace(map_list$log_gam_a, values = rep(NA, length(map_list$log_gam_a)))
    map_list$log_gam_b <- replace(map_list$log_gam_b, values = rep(NA, length(map_list$log_gam_b)))
    map_list$phi <- replace(map_list$phi, values = rep(NA, length(map_list$phi)))

    # Multispecies
    map_list$logH_1 <- replace(map_list$logH_1, values = rep(NA, length(map_list$logH_1)))
    map_list$logH_1a <- replace(map_list$logH_1a, values = rep(NA, length(map_list$logH_1a)))
    map_list$logH_1b <- replace(map_list$logH_1b, values = rep(NA, length(map_list$logH_1b)))

    map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
    map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
    map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))

  }

  # 2. MSVPA based predation
  if(data_list$msmMode == 1){
    # Multispecies
    map_list$logH_1 <- replace(map_list$logH_1, values = rep(NA, length(map_list$logH_1)))
    map_list$logH_1a <- replace(map_list$logH_1a, values = rep(NA, length(map_list$logH_1a)))
    map_list$logH_1b <- replace(map_list$logH_1b, values = rep(NA, length(map_list$logH_1b)))

    map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
    map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
    map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))

  }

  # 3. Kinzey and Punt predation equations
  if(data_list$msmMode > 1){

    # Holling Type 1
    if(data_list$msmMode == 2){
      map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
      map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }

    # Holling Type 2
    if(data_list$msmMode == 3){
      map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }

    # Holling Type 3
    if(data_list$msmMode == 4){
      map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }

    # Predator interference
    if(data_list$msmMode == 5){
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }

    # Predator preemption
    if(data_list$msmMode == 6){
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }

    # Hassell-Varley
    if(data_list$msmMode == 7){
      map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
    }

    # Ecosim
    if(data_list$msmMode == 8){
      map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }
  }


  ######################################################
  ####### Suitability bits
  ######################################################
  if(data_list$msmMode > 0){

    # 2.1. Empirical suitability
    if(data_list$suitMode == 0){
      map_list$log_gam_a <- replace(map_list$log_gam_a, values = rep(NA, length(map_list$log_gam_a)))
      map_list$log_gam_b <- replace(map_list$log_gam_b, values = rep(NA, length(map_list$log_gam_b)))
      map_list$phi <- replace(map_list$phi, values = rep(NA, length(map_list$phi)))
    }

    # 2.2. GAMMA suitability
    if(data_list$suitMode %in% c(1:3)){
      map_list$phi <- replace(map_list$phi, values = rep(NA, length(map_list$phi)))
    }

    # 2.3. and 2.4 Lognormal
    if(data_list$suitMode %in% c(4:5)){
      # Use all the parameters
    }
  }


  # STEP 3 - set up debug - I.E. turn off all parameters besides dummy
  map_list$dummy <- NA
  if(debug == TRUE){
    for(i in 1:length(map_list)){
      map_list[[i]] <- replace(map_list[[i]], values = rep(NA, length(map_list[[i]])))
    }
    map_list$dummy = 1
  }


  # STEP 4 -- Conver to factor
  for(i in 1:length(map_list)){
    map_list[[i]] <- factor(map_list[[i]])
  }


  return(map_list)
}
