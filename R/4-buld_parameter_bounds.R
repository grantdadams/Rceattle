#' Build parameter bounds
#'
#' Function to build parameter bounds based on Holsman et al 2015 and Kinzey and Punt 2010
#'
#' @param param_list Parameter list object built from \code{\link{build_params}}
#'
#' @return List of upper and lower bounds
#' @export
#'
build_bounds <- function( param_list = NULL){

  upper_bnd <- param_list
  lower_bnd <- param_list

  # General bounds
  for(i in 1:length(param_list)){
    upper_bnd[[i]] <- replace(upper_bnd[[i]], values = rep(Inf, length(upper_bnd[[i]])))
    lower_bnd[[i]] <- replace(lower_bnd[[i]], values = rep(-Inf, length(lower_bnd[[i]])))
  }

  # Predator selectivity
  lower_bnd$log_gam_a <- replace(lower_bnd$log_gam_a, values = rep(1.0e-10, length(lower_bnd$log_gam_a)))
  upper_bnd$log_gam_a <- replace(upper_bnd$log_gam_a, values = rep(19.9, length(upper_bnd$log_gam_a)))

  lower_bnd$log_gam_b <- replace(lower_bnd$log_gam_b, values = rep(-5.2, length(lower_bnd$log_gam_b)))
  upper_bnd$log_gam_b <- replace(upper_bnd$log_gam_b, values = rep(10, length(upper_bnd$log_gam_b)))

  # Functional form
  lower_bnd$logH_3 <- replace(lower_bnd$logH_3, values = rep(-30, length(lower_bnd$logH_3)))
  upper_bnd$logH_3 <- replace(upper_bnd$logH_3, values = rep(-0.000001, length(upper_bnd$logH_3)))

  lower_bnd$H_4 <- replace(lower_bnd$H_4, values = rep(-0.1, length(lower_bnd$H_4)))
  upper_bnd$H_4 <- replace(upper_bnd$H_4, values = rep(20, length(upper_bnd$H_4)))

  # Recruitment
  lower_bnd$rec_dev <- replace(lower_bnd$rec_dev, values = rep(-10, length(lower_bnd$rec_dev)))
  upper_bnd$rec_dev <- replace(upper_bnd$rec_dev, values = rep(10, length(upper_bnd$rec_dev)))

  # N0
  lower_bnd$init_dev <- replace(lower_bnd$init_dev, values = rep(-10, length(lower_bnd$init_dev)))
  upper_bnd$init_dev <- replace(upper_bnd$init_dev, values = rep(10, length(upper_bnd$init_dev)))

  # F
  lower_bnd$F_dev <- replace(lower_bnd$F_dev, values = rep(-10, length(lower_bnd$F_dev)))
  upper_bnd$F_dev <- replace(upper_bnd$F_dev, values = rep(10, length(upper_bnd$F_dev)))

  bounds <- list(upper= upper_bnd, lower = lower_bnd)


  # Make sure inits are within bounds
  if( sum(unlist(bounds$upper) < as.numeric(unlist(param_list))) > 0 | sum(as.numeric(unlist(param_list)) < unlist(bounds$lower)) > 0 ){
    lower_check <- param_list
    upper_check <- param_list
    param_check <- data.frame(matrix(NA, nrow = length(param_list), ncol = 3))
    colnames(param_check) <- c("Parameter", "Lower", "Upper")
    param_check$Parameter <- names(param_list)

    for(i in 1:length(param_list)){
      lower_check[[i]] <- param_list[[i]] < lower_bnd[[i]]
      upper_check[[i]] <- param_list[[i]] > upper_bnd[[i]]
      param_check$Lower[i] <- sum(lower_check[[i]])
      param_check$Upper[i] <- sum(upper_check[[i]])
    }

    print("Non-zero value indicates error in initial value")
    print(param_check)
    stop("Initial parameter values are not within bounds")
  }

  return(bounds)
}
