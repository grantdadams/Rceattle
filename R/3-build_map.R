#' Function to construct the TMB map argument for CEATTLE
#'
#' @description Reads a parameter list to construct map
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param params a parameter list created from \code{\link{build_params}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA except dummy
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE, which sets the map for ln_rec_sigma to NA
#'
#' @return a list of map arguments for each parameter
#' @export
build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE) {

  # Get year objects
  nyrs <- data_list$endyr - data_list$styr + 1
  nyrs_proj <- data_list$projyr - data_list$styr + 1
  yrs_proj <- (nyrs + 1):nyrs_proj

  # Get params to map
  map_list <- params

  # STEP 1 -- Convert map_list to seperate parameters
  for(i in 1:length(map_list)){
    map_list[[i]] <- replace(map_list[[i]], values = c(1:length(map_list[[i]])))
  }


  # -- 1 Map out future fishing mortality
  map_list$proj_F <- as.numeric(replace(map_list$proj_F, values = rep(NA, length(map_list$proj_F))))

  # -- 2 Map out future recruitment deviations
  map_list$rec_dev[,yrs_proj] <- as.numeric(replace(map_list$rec_dev[,yrs_proj], values = rep(NA, length(map_list$rec_dev[,yrs_proj]))))


  # STEP 2 -- NA out parameters not to be estimated

  # -- 3 Initial population deviations - map out last age and ages not seen
  for(i in 1:nrow(map_list$init_dev)){
    if((data_list$nages[i]-1) < ncol(map_list$init_dev)){
      map_list$init_dev[i,  (data_list$nages[i]) : ncol(map_list$init_dev) ] <- NA
    }
  }

  # -- 4 Fishery selectivity coefficients
  for( i in 1: nrow(data_list$fsh_control)){

    # -- 4.0. Empirical - sel_type = 0
    if(data_list$fsh_control$Selectivity[i] == 0){

      # Map out non-parametric
      map_list$fsh_sel_coff[i,] <- replace(map_list$fsh_sel_coff[i,], values = rep(NA, length(map_list$fsh_sel_coff[i,])))

      # Map out logistic and double logistic
      map_list$fsh_sel_slp[1:2, i] <- NA
      map_list$fsh_sel_inf[1:2, i] <- NA

      # Map out logistic and double logistic deviates
      map_list$fsh_sel_slp_dev[1:2, i,] <- NA
      map_list$fsh_sel_inf_dev[1:2, i,] <- NA

      # Map out selectivity var
      map_list$ln_sigma_fsh_sel[i] <- NA
    }

    # -- 4.1. Logitistic - sel_type = 1
    if(data_list$fsh_control$Selectivity[i] == 1){

      # Map out non-parametric
      map_list$fsh_sel_coff[i,] <- replace(map_list$fsh_sel_coff[i,], values = rep(NA, length(map_list$fsh_sel_coff[i,])))

      # Map out double logistic
      map_list$fsh_sel_slp[2, i] <- NA
      map_list$fsh_sel_inf[2, i] <- NA

      # Map out logistic and double logistic deviates
      map_list$fsh_sel_slp_dev[2, i,] <- NA
      map_list$fsh_sel_inf_dev[2, i,] <- NA

      # Map out time varying parameters if not used
      if(data_list$srv_control$Time_varying_sel[i] == 0){
        map_list$fsh_sel_slp_dev[1, i,] <- NA
        map_list$fsh_sel_inf_dev[1, i,] <- NA
      }

      # Map out selectivity var if not using time-varying or using random walk
      if(data_list$srv_control$Time_varying_sel[i] %in% c(0, 1)){
        map_list$ln_sigma_fsh_sel[i] <- NA
      }
    }

    # -- 4.2. Non-parametric - sel_type = 2
    if(data_list$fsh_control$Selectivity[i] == 2){ # Non-parametric at age

      # If nselages is  < max(nselages)
      if(data_list$fsh_control$Nselages[i] < max(data_list$fsh_control$Nselages, na.rm = TRUE)){
        mapped_ages <- (data_list$fsh_control$Nselages[i] + 1):max(data_list$fsh_control$Nselages, na.rm = T)
        map_list$fsh_sel_coff[i, mapped_ages]  <- replace(map_list$fsh_sel_coff[i, mapped_ages], values = rep(NA, length(map_list$fsh_sel_coff[i, mapped_ages])))
      }

      # Map out logistic and double logistic
      map_list$fsh_sel_slp[1:2, i] <- NA
      map_list$fsh_sel_inf[1:2, i] <- NA

      # Map out logistic and double logistic deviates
      map_list$fsh_sel_slp_dev[1:2, i,] <- NA
      map_list$fsh_sel_inf_dev[1:2, i,] <- NA

      # Map out selectivity var
      map_list$ln_sigma_fsh_sel[i] <- NA
    }

    # -- 4.3. Double logistic - sel_type = 3
    if(data_list$fsh_control$Selectivity[i] == 3){ # Double logistic
      # Map out non-parametric
      map_list$fsh_sel_coff[i,] <- replace(map_list$fsh_sel_coff[i,], values = rep(NA, length(map_list$fsh_sel_coff[i,])))

      # Map out time varying parameters if not used
      if(data_list$srv_control$Time_varying_sel[i] == 0){
        map_list$fsh_sel_slp_dev[1:2, i,] <- NA
        map_list$fsh_sel_inf_dev[1:2, i,] <- NA
      }

      # Map out selectivity var if not using time-varying or using random walk
      if(data_list$srv_control$Time_varying_sel[i] %in% c(0, 1, 3)){
        map_list$ln_sigma_fsh_sel[i] <- NA
      }
    }
  }


  # -- 5 Survey selectivity coefficients
  for( i in 1: nrow(data_list$srv_control)){

    # -- 5.0. Empirical - sel_type = 0
    if(data_list$srv_control$Selectivity[i] == 0){

      # Map out non-parametric
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))

      # Map out logistic and double logistic
      map_list$srv_sel_slp[1:2, i] <- NA
      map_list$srv_sel_inf[1:2, i] <- NA

      # Map out logistic and double logistic deviates
      map_list$srv_sel_slp_dev[1:2, i,] <- NA
      map_list$srv_sel_inf_dev[1:2, i,] <- NA

      # Map out selectivity var
      map_list$ln_sigma_srv_sel[i] <- NA
    }

    # -- 5.1. Logitistic - sel_type = 1
    if(data_list$srv_control$Selectivity[i] == 1){

      # Map out non-parametric
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))

      # Map out double logistic
      map_list$srv_sel_slp[2, i] <- NA
      map_list$srv_sel_inf[2, i] <- NA

      # Map out logistic and double logistic deviates
      map_list$srv_sel_slp_dev[2, i,] <- NA
      map_list$srv_sel_inf_dev[2, i,] <- NA

      # Map out time varying parameters if not used
      if(data_list$srv_control$Time_varying_sel[i] == 0){
        map_list$srv_sel_slp_dev[1, i,] <- NA
        map_list$srv_sel_inf_dev[1, i,] <- NA
      }

      # Map out selectivity var if not using time-varying or using random walk
      if(data_list$srv_control$Time_varying_sel[i] %in% c(0, 1)){
        map_list$ln_sigma_srv_sel[i] <- NA
      }
    }

    # -- 5.2. Non-parametric - sel_type = 2
    if(data_list$srv_control$Selectivity[i] == 2){ # Non-parametric at age

      # If nselages is  < max(nselages)
      if(data_list$srv_control$Nselages[i] < max(data_list$srv_control$Nselages, na.rm = TRUE)){
        mapped_ages <- (data_list$srv_control$Nselages[i] + 1):max(data_list$srv_control$Nselages, na.rm = T)
        map_list$srv_sel_coff[i, mapped_ages]  <- replace(map_list$srv_sel_coff[i, mapped_ages], values = rep(NA, length(map_list$srv_sel_coff[i, mapped_ages])))
      }

      # Map out logistic and double logistic
      map_list$srv_sel_slp[1:2, i] <- NA
      map_list$srv_sel_inf[1:2, i] <- NA

      # Map out logistic and double logistic deviates
      map_list$srv_sel_slp_dev[1:2, i,] <- NA
      map_list$srv_sel_inf_dev[1:2, i,] <- NA

      # Map out selectivity var
      map_list$ln_sigma_srv_sel[i] <- NA
    }

    # -- 5.3. Double logistic - sel_type = 3
    if(data_list$srv_control$Selectivity[i] == 3){ # Double logistic
      # Map out non-parametric
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))

      # Map out time varying parameters if not used
      if(data_list$srv_control$Time_varying_sel[i] == 0){
        map_list$srv_sel_slp_dev[1:2, i,] <- NA
        map_list$srv_sel_inf_dev[1:2, i,] <- NA
      }

      # Map out selectivity var if not using time-varying or using random walk or time block
      if(data_list$srv_control$Time_varying_sel[i] %in% c(0, 1, 3)){
        map_list$ln_sigma_srv_sel[i] <- NA
      }
    }
  }




  # -- 6. Survey control
  for( i in 1: nrow(data_list$srv_control)){
    # Catchability of surveys
    # If not estimating turn of
    if(data_list$srv_control$Estimate_q[i] %in% c(0,2)){
      map_list$log_srv_q[i] <- NA
      map_list$ln_srv_q_dev[i,] <- NA
      map_list$ln_sigma_srv_q[i] <- NA
    }

    # Time-varying catchability of surveys
    # If not estimating turn of
    if(data_list$srv_control$Time_varying_q[i] == 0){
      map_list$ln_srv_q_dev[i,] <- NA
      map_list$ln_sigma_srv_q[i] <- NA
    }

    # Random walk - map out q_sd
    if(data_list$srv_control$Time_varying_q[i] == 1){
      map_list$ln_sigma_srv_q[i] <- NA
    }

    # Time block - map out q_sd
    if(data_list$srv_control$Time_varying_q[i] == 3){
      map_list$ln_sigma_srv_q[i] <- NA
    }

    # Standard deviation of surveys index
    # If not estimating turn of
    if(data_list$srv_control$Estimate_survey_sd[i] %in% c(0,2)){
      map_list$ln_sigma_srv_index[i] <- NA
    }
  }


  # -- 7. Share survey q and selectivity
  sel_index <- data_list$srv_control$Selectivity_index
  sel_index_tested <- c()

  q_index <- data_list$srv_control$Q_index
  q_index_tested <- c()
  rows_tests <- c()
  for(i in 1: nrow(data_list$srv_control)){
    sel_test <- sel_index[i] %in% sel_index_tested
    q_test <- q_index[i] %in% q_index_tested

    # If selectivity is the same as a previous index
    if(sel_test){
      sel_duplicate <- which(sel_index_tested == sel_index[i])[1]
      sel_duplicate_vec <- c(which(sel_index_tested == sel_index[i]), i)

      # Error check selectivity type
      if(length(unique(data_list$srv_control$Selectivity[sel_duplicate_vec])) > 1){
        warning("Survey selectivity of surveys with same Selectivity_index is not the same")
        warning(paste0("Double check Selectivity in srv_control of surveys:", paste(data_list$srv_control$Survey_name[sel_duplicate_vec])))
      }


      # Error check time-varying selectivity type
      if(length(unique(data_list$srv_control$Time_varying_sel[sel_duplicate_vec])) > 1){
        warning("Time varying survey selectivity of surveys with same Selectivity_index is not the same")
        warning(paste0("Double check Time_varying_sel in srv_control of surveys:", paste(data_list$srv_control$Survey_name[sel_duplicate_vec])))
      }

      # FIXME add checks for surveys sel sigma

      # Make selectivity maps the same if selectivity is the same
      map_list$srv_sel_slp[1:2, i] <- map_list$srv_sel_slp[1:2, sel_duplicate]
      map_list$srv_sel_inf[1:2, i] <- map_list$srv_sel_inf[1:2, sel_duplicate]
      map_list$srv_sel_coff[i,] <- map_list$srv_sel_coff[sel_duplicate,]
      map_list$srv_sel_slp_dev[1:2, i,] <- map_list$srv_sel_slp_dev[1:2, sel_duplicate,]
      map_list$srv_sel_inf_dev[1:2, i,] <- map_list$srv_sel_inf_dev[1:2, sel_duplicate,]
      map_list$ln_sigma_srv_sel[i] <- map_list$ln_sigma_srv_sel[sel_duplicate]
    }


    # If catchability is the same as a previous index
    if(q_test){
      q_duplicate <- which(q_index_tested == q_index[i])[1]
      q_duplicate_vec <- c(which(q_index_tested == q_index[i]), i)

      # Error check selectivity type
      if(length(unique(data_list$srv_control$Estimate_q[q_duplicate_vec])) > 1){
        warning("Survey catchability of surveys with same Q_index is not the same")
        warning(paste0("Double check Estimate_q in srv_control of surveys:", paste(data_list$srv_control$Survey_name[q_duplicate_vec])))
      }


      # Error check time-varying selectivity type
      if(length(unique(data_list$srv_control$Time_varying_q[q_duplicate_vec])) > 1){
        warning("Time varying survey catchability of surveys with same Q_index is not the same")
        warning(paste0("Double check Time_varying_q in srv_control of surveys:", paste(data_list$srv_control$Survey_name[q_duplicate_vec])))
      }

      # FIXME add checks for surveys q sigma

      # Make catchability maps the same if selectivity is the same
      map_list$ln_srv_q_dev[i,] <- map_list$ln_srv_q_dev[sel_duplicate,]
      map_list$ln_sigma_srv_q[i] <- map_list$ln_sigma_srv_q[sel_duplicate]
    }


    # Add index
    sel_index_tested <- c(sel_index_tested, sel_index[i])
    q_index_tested <- c(q_index_tested, q_index[i])
  }


  # -- 8. Fishery control
  for( i in 1: nrow(data_list$fsh_control)){
    # Standard deviation of fishery time series
    # If not estimating turn of
    if(data_list$fsh_control$Estimate_catch_sd[i] %in% c(0,2)){
      map_list$ln_sigma_fsh_catch[i] <- NA
    }
  }


  # -- 9. Recruitment deviation sigmas - turn off if not estimating
  if(random_rec == FALSE){
    map_list$ln_rec_sigma <- replace(map_list$ln_rec_sigma, values = rep(NA, length(map_list$ln_rec_sigma)))
  }



  # -- 10. Map out Fdev for years with 0 catch to very low number
  fsh_biom <- data_list$fsh_biom
  fsh_ind <- fsh_biom$Fishery_code[which(fsh_biom$Catch_kg == 0)]
  yr_ind <- fsh_biom$Year[which(fsh_biom$Catch_kg == 0)] - data_list$styr + 1

  map_list$F_dev[fsh_ind, yr_ind] <- NA
  map_list$fsh_sel_slp_dev[1:2, i, yr_ind] <- NA
  map_list$fsh_sel_inf_dev[1:2, i, yr_ind] <- NA

  ######################################################
  ####### Predation bits
  ######################################################
  # 1. Turn off all predation parameters for single species
  if(data_list$msmMode == 0){

    # Suitability parameters
    map_list$log_gam_a <- replace(map_list$log_gam_a, values = rep(NA, length(map_list$log_gam_a)))
    map_list$log_gam_b <- replace(map_list$log_gam_b, values = rep(NA, length(map_list$log_gam_b)))
    map_list$log_phi <- replace(map_list$log_phi, values = rep(NA, length(map_list$log_phi)))

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
      map_list$log_phi <- replace(map_list$log_phi, values = rep(NA, length(map_list$log_phi)))
    }

    # 2.2. GAMMA suitability
    if(data_list$suitMode %in% c(1:3)){
      map_list$log_phi <- replace(map_list$log_phi, values = rep(NA, length(map_list$log_phi)))
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
  map_list_grande <- list()
  map_list_grande[[1]] <- map_list
  map_list_grande[[2]] <- map_list

  for(i in 1:length(map_list_grande[[1]])){
    map_list_grande[[1]][[i]] <- factor(map_list_grande[[1]][[i]])
  }


  return(map_list_grande)
}
