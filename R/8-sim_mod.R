#' Simulate Rceattle data
#'
#' @description  Simulates data used in Rceattle from the expected values etimated from Rceattle. The variances and uncertainty are the same as used in the operating model. The function currently simulates (assumed distribution) the following: survey biomass (log-normal), survey catch-at-length/age (multinomial), EIT biomass (log-normal), EIT catch-at-length/age (multinomial), total catch (kg) (log-normal), and catch-at-length/age.
#'
#' @param ceattle_obj CEATTLE model object exported from \code{\link{Rceattle}}
#' @expected TRUE/FALSE, whether to simulate the data or export the expected value
#' @export
sim_mod <- function( ceattle_obj, expected = FALSE ){
  #TODO Options for simulation diet data: multinomial, sqrt-normal, dirichlet, multinomial
  dat_sim <- ceattle_obj$data_list

  for(sp in 1:ceattle_obj$data_list$nspp){ #FIXME: may need to adjust for multiple surveys

    # Slot 0 -- BT survey biomass -- NFMS annual BT survey
    srv_biom_lse = dat_sim$srv_biom_se / dat_sim$srv_biom
    srv_biom_lse = sqrt( log( (  srv_biom_lse^2 ) + 1));

    if(expected){# Simulate
      values <- exp(rnorm(length(dat_sim$srv_biom[sp,]), mean = log(ceattle_obj$quantities$srv_bio_hat[sp,1:dat_sim$nyrs_srv_age[sp]]), sd = srv_biom_lse[sp,1:dat_sim$nyrs_srv_age[sp]]))
    } else{ # Expected value
      values <- ceattle_obj$quantities$srv_bio_hat[sp,1:dat_sim$nyrs_srv_age[sp]]
    }

    dat_sim$srv_biom[sp,] = replace(dat_sim$srv_biom[sp,], values = values)


    # Slot 1 -- BT survey age composition -- NFMS annual BT survey
    for(yr in 1:dat_sim$nyrs_srv_age[sp]){

      if(expected){# Simulate
        values <- rmultinom(n = 1, size = dat_sim$srv_age_n[sp, yr], prob = ceattle_obj$quantities$srv_age_hat[sp, , yr])
      } else{ # Expected value
        values <- ceattle_obj$quantities$srv_age_hat[sp, , yr]
      }

      dat_sim$srv_age_obs[yr,,sp] = replace( dat_sim$srv_age_obs[yr,,sp], values = values)
    }


    if(sp == 1){
      # Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
      if(expected){# Simulate
        values <- exp(rnorm(length(dat_sim$obs_eit), mean = log(ceattle_obj$quantities$eit_hat[1:dat_sim$n_eit]), sd = .2))
      } else{ # Expected value
        values <- ceattle_obj$quantities$eit_hat[1:dat_sim$n_eit]
      }

      dat_sim$obs_eit = replace(dat_sim$obs_eit, values = values)


      # Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
      for(yr in 1:dat_sim$n_eit){
        if(expected){# Simulate
          values <- rmultinom(n = 1, size = dat_sim$eit_age_n[yr], prob = ceattle_obj$quantities$eit_age_comp_hat[, yr])
        } else{ # Expected value
          values <- ceattle_obj$quantities$eit_age_comp_hat[, yr]
        }
        dat_sim$obs_eit_age[yr,] = replace( dat_sim$obs_eit_age[yr,], values = values)
      }
    }


    # Slot 4 -- Total catch -- Fishery observer data
    if(expected){# Simulate
      values <- exp(rnorm(length(dat_sim$tcb_obs[sp,]), mean = log(ceattle_obj$quantities$tc_biom_hat[sp,1:dat_sim$nyrs_tc_biom[sp]]), sd = 0.05))
    } else{ # Expected value
      values <- ceattle_obj$quantities$tc_biom_hat[sp,1:dat_sim$nyrs_tc_biom[sp]]
    }
    dat_sim$tcb_obs[sp,] = replace(dat_sim$tcb_obs[sp,], values = values)


    # Slot 5 -- Fishery age composition -- Fishery observer data
    for(yr in 1:dat_sim$nyrs_fsh_comp[sp]){
      if(expected){# Simulate
        values <- rmultinom(n = 1, size = 200, prob = ceattle_obj$quantities$fsh_age_hat[sp, , yr])
      } else{ # Expected value
        values <- ceattle_obj$quantities$fsh_age_hat[sp, , yr]
      }
      dat_sim$obs_catch[yr,,sp] = replace( dat_sim$obs_catch[yr,,sp], values = values)
    }

    # Slot 5 -- Diet composition from lognormal suitability
    # 4D
    if(length(dim(dat_sim$UobsWtAge)) == 4){
      for(r_age in 1:dat_sim$nages[sp]){
        if(ceattle_obj$data_list$suitMode %in% c(2:3) & expected & sum(ceattle_obj$quantities$mn_UobsWtAge_hat[sp,,r_age,] > 0) > 0){
          values <- rmultinom(n = 1, size = 20, prob = ceattle_obj$quantities$mn_UobsWtAge_hat[sp,,r_age,]) #FIXME change sample size
        } else {
          values <- ceattle_obj$quantities$mn_UobsWtAge_hat[sp,,r_age,]
        }
        dat_sim$UobsWtAge[sp,,r_age,] <- replace(dat_sim$UobsWtAge[sp,,r_age,], values = values)
      }
    }

    # 5D
    if(length(dim(dat_sim$UobsWtAge)) == 5){
      for(yr in 1:dat_sim$nyrs_fsh_comp[sp]){
        for(r_age in 1:dat_sim$nages[sp]){
          if(ceattle_obj$data_list$suitMode %in% c(2:3) & expected& sum(ceattle_obj$quantities$UobsWtAge_hat[sp,,r_age, , yr] > 0) > 0){
            values <- rmultinom(n = 1, size = 20, prob = ceattle_obj$quantities$UobsWtAge_hat[sp,,r_age,,yr]) #FIXME change sample size
          } else {
            values <- ceattle_obj$quantities$UobsWtAge_hat[sp,,r_age,,yr]
          }
          dat_sim$UobsWtAge[sp,,r_age,,yr] <- replace(dat_sim$UobsWtAge[sp,,r_age,,yr], values = values)
        }
      }
    }

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

