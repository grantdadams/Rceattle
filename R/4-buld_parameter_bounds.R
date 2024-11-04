#' Build parameter bounds
#'
#' Function to build parameter bounds based on Holsman et al 2015 and Kinzey and Punt 2010
#'
#' @param param_list Parameter list object built from \code{\link{build_params}}
#'
#' @return List of upper and lower bounds
#' @export
#'
build_bounds <- function(param_list = NULL, data_list) {

  upper_bnd <- param_list
  lower_bnd <- param_list

  # General bounds
  for (i in 1:length(param_list)) {
    upper_bnd[[i]] <- replace(upper_bnd[[i]], values = rep(Inf, length(upper_bnd[[i]])))
    lower_bnd[[i]] <- replace(lower_bnd[[i]], values = rep(-Inf, length(lower_bnd[[i]])))
  }

  # # Predator selectivity Bounds for gamma suitability
  # if (data_list$suitMode %in% c(1:2)) {
  #     lower_bnd$log_gam_a <- replace(lower_bnd$log_gam_a, values = rep(1e-10, length(lower_bnd$log_gam_a)))
  #     upper_bnd$log_gam_a <- replace(upper_bnd$log_gam_a, values = rep(19.9, length(upper_bnd$log_gam_a)))
  # } else {
  #     lower_bnd$log_gam_b <- replace(lower_bnd$log_gam_b, values = rep(-10, length(lower_bnd$log_gam_b)))
  #     upper_bnd$log_gam_b <- replace(upper_bnd$log_gam_b, values = rep(20, length(upper_bnd$log_gam_b)))
  # }
  #
  # # Kinzey functional form
  # lower_bnd$logH_3 <- replace(lower_bnd$logH_3, values = rep(-30, length(lower_bnd$logH_3)))
  # upper_bnd$logH_3 <- replace(upper_bnd$logH_3, values = rep(-1e-06, length(upper_bnd$logH_3)))
  #
  # lower_bnd$H_4 <- replace(lower_bnd$H_4, values = rep(-0.1, length(lower_bnd$H_4)))
  # upper_bnd$H_4 <- replace(upper_bnd$H_4, values = rep(20, length(upper_bnd$H_4)))

  # Recruitment
  lower_bnd$rec_dev <- replace(lower_bnd$rec_dev, values = rep(-15, length(lower_bnd$rec_dev)))
  upper_bnd$rec_dev <- replace(upper_bnd$rec_dev, values = rep(15, length(upper_bnd$rec_dev)))

  # Selectivity
  lower_bnd$ln_sel_slp <- replace(lower_bnd$ln_sel_slp, values = rep(0.01, length(lower_bnd$ln_sel_slp)))
  for(flt in 1:nrow(data_list$fleet_control)){
    upper_bnd$ln_sel_slp[,flt,] <- replace(upper_bnd$ln_sel_slp[,flt,], values = rep(data_list$nages[data_list$fleet_control$Species[flt]]+0.5, length(upper_bnd$ln_sel_slp[,flt,])))
  }


  # Selectivity deviates
  # - Slope
  for(i in 1:nrow(data_list$fleet_control)){
    # If using blocks don't put bounds on deviates, as these are estimated
    if(data_list$fleet_control$Time_varying_sel[i] != 3){
      lower_bnd$ln_sel_slp_dev[,i,,] <- replace(lower_bnd$ln_sel_slp_dev[,i,,], values = rep(-5, length(lower_bnd$ln_sel_slp_dev[,i,,])))
      upper_bnd$ln_sel_slp_dev[,i,,] <- replace(upper_bnd$ln_sel_slp_dev[,i,,], values = rep(5, length(upper_bnd$ln_sel_slp_dev[,i,,])))
    }
  }

  # - Asymptotic
  for(i in 1:nrow(data_list$fleet_control)){
    # If using blocks don't put bounds on deviates, as these are estimated
    if(data_list$fleet_control$Time_varying_sel[i] != 3){
      lower_bnd$sel_inf_dev[,i,,] <- replace(lower_bnd$sel_inf_dev[,i,,], values = rep(-5, length(lower_bnd$sel_inf_dev[,i,,])))
      upper_bnd$sel_inf_dev[,i,,] <- replace(upper_bnd$sel_inf_dev[,i,,], values = rep(5, length(upper_bnd$sel_inf_dev[,i,,])))
    }
  }


  # N0
  lower_bnd$init_dev <- replace(lower_bnd$init_dev, values = rep(-1000, length(lower_bnd$init_dev)))
  upper_bnd$init_dev <- replace(upper_bnd$init_dev, values = rep(23, length(upper_bnd$init_dev)))

  # Survey variance
  lower_bnd$ln_sigma_srv_index <- replace(lower_bnd$ln_sigma_srv_index, values = rep(-10, length(lower_bnd$ln_sigma_srv_index)))
  upper_bnd$ln_sigma_srv_index <- replace(upper_bnd$ln_sigma_srv_index, values = rep(10, length(upper_bnd$ln_sigma_srv_index)))

  # Fishery variance
  lower_bnd$ln_sigma_fsh_catch <- replace(lower_bnd$ln_sigma_fsh_catch, values = rep(-10, length(lower_bnd$ln_sigma_fsh_catch)))
  upper_bnd$ln_sigma_fsh_catch <- replace(upper_bnd$ln_sigma_fsh_catch, values = rep(3, length(upper_bnd$ln_sigma_fsh_catch)))

  # F
  lower_bnd$F_dev <- replace(lower_bnd$F_dev, values = rep(-1000, length(lower_bnd$F_dev)))
  upper_bnd$F_dev <- replace(upper_bnd$F_dev, values = rep(10, length(upper_bnd$F_dev)))

  lower_bnd$ln_M1 <- replace(lower_bnd$ln_M1, values = rep(log(0.001), length(lower_bnd$ln_M1)))
  upper_bnd$ln_M1 <- replace(upper_bnd$ln_M1, values = rep(log(2), length(upper_bnd$ln_M1)))

  # Combine bounds in list
  bounds <- list(upper = upper_bnd, lower = lower_bnd)


  # Make sure inits are within bounds
  if (sum(unlist(bounds$upper) < as.numeric(unlist(param_list)), na.rm = TRUE) > 0 | sum(as.numeric(unlist(param_list)) <
                                                                                         unlist(bounds$lower), na.rm = TRUE) > 0) {
    lower_check <- param_list
    upper_check <- param_list
    param_check <- data.frame(matrix(NA, nrow = length(param_list), ncol = 3))
    colnames(param_check) <- c("Parameter", "Lower", "Upper")
    param_check$Parameter <- names(param_list)

    for (i in 1:length(param_list)) {
      lower_check[[i]] <- param_list[[i]] < lower_bnd[[i]]
      upper_check[[i]] <- param_list[[i]] > upper_bnd[[i]]
      param_check$Lower[i] <- sum(lower_check[[i]], na.rm = TRUE)
      param_check$Upper[i] <- sum(upper_check[[i]], na.rm = TRUE)
    }

    print("Non-zero value indicates error in initial value")
    print(param_check)
    stop("Initial parameter values are not within bounds")
  }

  return(bounds)
}
