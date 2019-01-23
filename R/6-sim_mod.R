#' Simulate Rceattle data
#'
#' @description  Simulates data used in Rceattle from the expected values etimated from Rceattle. The variances and uncertainty are the same as used in the operating model. The function currently simulates (assumed distribution) the following: survey biomass (log-normal), survey catch-at-length/age (multinomial), EIT biomass (log-normal), EIT catch-at-length/age (multinomial), total catch (kg) (log-normal), and catch-at-length/age.
#'
#' @param ceattle_obj CEATTLE model object exported from \code{\link{Rceattle}}
#' @export
sim_mod <- function( ceattle_obj ){
  #TODO Options for simulation diet data: multinomial, sqrt-normal, dirichlet, multinomial
  dat_sim <- ceattle_obj$data_list

  for(sp in 1:ceattle_obj$data_list$nspp){ #FIXME: may need to adjust for multiple surveys

    # Slot 0 -- BT survey biomass -- NFMS annual BT survey
    srv_biom_lse = dat_sim$srv_biom_se / dat_sim$srv_biom
    srv_biom_lse = sqrt( log( (  srv_biom_lse^2 ) + 1));

    dat_sim$srv_biom[sp,] = replace(dat_sim$srv_biom[sp,], values = exp(rnorm(length(dat_sim$srv_biom[sp,]), mean = log(ceattle_obj$quantities$srv_bio_hat[sp,1:dat_sim$nyrs_srv_age[sp]]), sd = srv_biom_lse[sp,1:dat_sim$nyrs_srv_age[sp]])))


    # Slot 1 -- BT survey age composition -- NFMS annual BT survey
    for(yr in 1:dat_sim$nyrs_srv_age[sp]){
      dat_sim$srv_age_obs[yr,,sp] = replace( dat_sim$srv_age_obs[yr,,sp], values = rmultinom(n = 1, size = dat_sim$srv_age_n[sp, yr], prob = ceattle_obj$quantities$srv_age_hat[sp, , yr]))
    }


    if(sp == 1){
      # Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
      dat_sim$obs_eit = replace(dat_sim$obs_eit, values = exp(rnorm(length(dat_sim$obs_eit), mean = log(ceattle_obj$quantities$eit_hat[1:dat_sim$n_eit]), sd = .2)))


      # Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
      for(yr in 1:dat_sim$n_eit){
        dat_sim$obs_eit_age[yr,] = replace( dat_sim$obs_eit_age[yr,], values = rmultinom(n = 1, size = dat_sim$eit_age_n[yr], prob = ceattle_obj$quantities$eit_age_comp_hat[, yr]))
      }
    }


    # Slot 4 -- Total catch -- Fishery observer data
    dat_sim$tcb_obs[sp,] = replace(dat_sim$tcb_obs[sp,], values = exp(rnorm(length(dat_sim$tcb_obs[sp,]), mean = log(ceattle_obj$quantities$tc_biom_hat[sp,1:dat_sim$nyrs_tc_biom[sp]]), sd = 0.05)))


    # Slot 5 -- Fishery age composition -- Fishery observer data
    for(yr in 1:dat_sim$nyrs_fsh_comp[sp]){
      dat_sim$obs_catch[yr,,sp] = replace( dat_sim$obs_catch[yr,,sp], values = rmultinom(n = 1, size = 200, prob = ceattle_obj$quantities$fsh_age_hat[sp, , yr]))
    }

    #FIXME include diet data
  }


  return(dat_sim)
}



#' Evaluate simulation performance
#'
#' @description Function to evaluate the simulation performance with regard to bias using the median relative error (MRE) and precision using the coefficient of variation.
#'
#' @param operating_mod CEATTLE model object exported from \code{\link{Rceattle}} to be used as the operating model
#' @param simulation_mods List of CEATTLE model objects exported from \code{\link{Rceattle}} fit to simulated data
compare_sim <- function( operating_mod, simulation_mods, object = "quantities"){

  # Get differences
  sim_mre <- list()
  sim_mse <- list()
  sim_mean <- list()
  sim_median <- list()
  sim_sd <- list()
  sim_cv <- list()
  sim_params <- list()

    for( j in 1:length(names(operating_mod[[object]]))){
      param <- names(operating_mod[[object]])[j]

      sim_mre[[param]] <- list()
      sim_mse[[param]] <- list()
      sim_mean[[param]] <- list()
      sim_sd[[param]] <- list()
      sim_cv[[param]] <- list()
      sim_params[[param]] <- list()

      om_params <- operating_mod[[object]][[param]]

      for(i in 1:length(simulation_mods)){

        sm_params <- simulation_mods[[i]][[object]][[param]]

        sim_params[[param]][[i]] <- sm_params
        sim_mre[[param]][[i]] <- (sm_params - om_params)/om_params
        sim_mse[[param]][[i]] <- (sm_params - om_params)^2
      }

      param_dim <- length(dim(om_params))

      # If 1 value
      if(param_dim == 0){
        sim_mean[[param]] <- mean(unlist(sim_params[[param]]))
        sim_sd[[param]] <- sd(unlist(sim_params[[param]]))
        sim_cv[[param]] <- sim_sd[[param]] / sim_mean[[param]]

        sim_mre[[param]] <- median(unlist(sim_mre[[param]]))
        sim_mse[[param]] <- mean(unlist(sim_mse[[param]]))
      }

      # If multiple values
      if(param_dim > 0){
        # Get mean, sd, and CV
        sim_mean[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, mean)
        sim_median[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, median)
        sim_sd[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, sd)
        sim_cv[[param]] <- sim_sd[[param]] / sim_mean[[param]]

        sim_mre[[param]] <- apply(simplify2array(sim_mre[[param]]), 1:param_dim, median)
        sim_mse[[param]] <- apply(simplify2array(sim_mse[[param]]), 1:param_dim, mean)
      }
    }


  result_list <- list(Mean = sim_mean, Median = sim_median, SD = sim_sd, CV = sim_cv, MRE = sim_mre, MSE = sim_mse)
  return(result_list)
}

