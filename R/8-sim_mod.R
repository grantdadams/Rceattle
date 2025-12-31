#' Simulate Rceattle data
#'
#' @description  Simulates data used in Rceattle from the expected values etimated from Rceattle. The variances and uncertainty are the same as used in the operating model. The function currently simulates (assumed distribution) the following: survey biomass (log-normal), survey catch-at-length/age (multinomial), EIT biomass (log-normal), EIT catch-at-length/age (multinomial), total catch (kg) (log-normal), and catch-at-length/age.
#'
#' @param Rceattle CEATTLE model object exported from \code{\link{Rceattle}}
#' @param simulate TRUE/FALSE, whether to simulate the data or export the expected value
#'
#' @export
sim_mod <- function(Rceattle, simulate = FALSE) {
  # TODO Options for simulation diet data: multinomial, sqrt-normal, dirichlet, multinomial
  dat_sim <- Rceattle$data_list


  # Slot 0 -- BT survey biomass -- NFMS annual BT survey
  ln_index_sd = Rceattle$quantities$ln_index_sd

  if (simulate) {
    # Simulate
    values <- exp(rnorm(length(dat_sim$index_data$Observation), mean = log(Rceattle$quantities$index_hat) - (ln_index_sd^2)/2, sd = ln_index_sd))
  } else {
    # Estimated value
    values <- Rceattle$quantities$index_hat
  }
  dat_sim$index_data$Observation = values


  # Slot 1 -- Age composition
  for (obs in 1:nrow(dat_sim$comp_data)) {
    if (simulate & (sum(Rceattle$quantities$comp_hat[obs,], na.rm = TRUE) > 0)) {
      # Simulate
      #FIXME add dirichlet multinomial
      values <- rmultinom(n = 1, size = dat_sim$comp_data$Sample_size[obs], prob = Rceattle$quantities$comp_hat[obs,])
    } else {
      # Expected value
      values <- Rceattle$quantities$comp_hat[obs, ]
    }
    dat_sim$comp_data[obs, 9:ncol(dat_sim$comp_data)] = values
  }




  # Slot 2 -- Total catch -- Fishery observer data
  fsh_biom_lse = Rceattle$quantities$ln_catch_sd

  if (simulate) {
    # Simulate
    values <- exp(rnorm(length(dat_sim$catch_data$Catch), mean = log(Rceattle$quantities$catch_hat) - (fsh_biom_lse^2)/2, sd = fsh_biom_lse))
  } else {
    # simulate value
    values <- Rceattle$quantities$catch_hat
  }

  dat_sim$catch_data$Catch = values


  #TODO
  # # Slot 5 -- Diet composition from lognormal suitability 4D
  # if (length(dim(dat_sim$diet_data)) == 4) {
  #     for (sp in 1:dat_sim$nspp) {
  #         for (r_age in 1:dat_sim$nages[sp]) {
  #             if (Rceattle$data_list$suitMode > 0 & simulate & sum(Rceattle$quantities$mn_UobsWtAge_hat[sp, , r_age,
  #                                                                                                       ] > 0) > 0) {
  #                 values <- rmultinom(n = 1, size = dat_sim$stom_tau[sp], prob = Rceattle$quantities$mn_UobsWtAge_hat[sp,
  #                                                                                                                     , r_age, ])  #FIXME change sample size
  #             } else {
  #                 values <- Rceattle$quantities$mn_UobsWtAge_hat[sp, , r_age, ]
  #             }
  #             dat_sim$diet_data[sp, , r_age, ] <- replace(dat_sim$diet_data[sp, , r_age, ], values = values)
  #         }
  #     }
  # }
  #
  # # 5D
  # if (length(dim(dat_sim$diet_data)) == 5) {
  #     for (sp in 1:dat_sim$nspp) {
  #         for (yr in 1:dat_sim$nyrs_fsh_comp[sp]) {
  #             for (r_age in 1:dat_sim$nages[sp]) {
  #                 if (Rceattle$data_list$suitMode > 0 & simulate & sum(Rceattle$quantities$UobsWtAge_hat[sp, , r_age,
  #                                                                                                        , yr] > 0) > 0) {
  #                     values <- rmultinom(n = 1, size = dat_sim$stom_tau[sp], prob = Rceattle$quantities$UobsWtAge_hat[sp,
  #                                                                                                                      , r_age, , yr])  #FIXME change sample size
  #                 } else {
  #                     values <- Rceattle$quantities$UobsWtAge_hat[sp, , r_age, , yr]
  #                 }
  #                 dat_sim$diet_data[sp, , r_age, , yr] <- replace(dat_sim$diet_data[sp, , r_age, , yr], values = values)
  #             }
  #         }
  #     }
  # }

  return(dat_sim)
}

#' Sample historical recruitment deviates and place in the projection
#'
#' @param Rceattle CEATTLE model object exported from \code{\link{Rceattle}}
#' @param sample_rec Include resampled recruitment deviates from the"hindcast" in the projection of the OM. Resampled deviates are used rather than sampling from N(0, sigmaR) because initial deviates bias R0 low. If false, uses mean of recruitment deviates.
#' @param update_model Update model dynamics. Default = TRUE
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier \code{mean rec * (1 + (rec_trend/projection years) * 1:projection years)}. Can be of length 1 or of length nspp. If length 1, all species get the same trend.
#'
#' @returns Rceattle model
#' @export
#'
sample_rec <- function(Rceattle, sample_rec = TRUE, update_model = TRUE, rec_trend = 0){

  # Years for simulations
  hind_yrs <- (Rceattle$data_list$styr) : Rceattle$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  proj_yrs <- (Rceattle$data_list$endyr + 1) : Rceattle$data_list$projyr
  proj_nyrs <- length(proj_yrs)

  # - Adjust rec trend
  if(length(rec_trend)==1){
    rec_trend = rep(rec_trend, Rceattle$data_list$nspp)
  }

  # Replace future rec devs ----
  #FIXME - update non-sample rec for stock recruit relationship
  for(sp in 1:Rceattle$data_list$nspp){

    # -- where SR curve is estimated directly
    if(Rceattle$data_list$srr_fun == Rceattle$data_list$srr_pred_fun){
      if(sample_rec){ # Sample devs from hindcast
        rec_dev <- sample(x = Rceattle$estimated_params$rec_dev[sp, 1:hind_nyrs], size = proj_nyrs, replace = TRUE) + log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
      } else{ # Set to mean rec otherwise
        rec_dev <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs]) * (1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs))  - log(Rceattle$quantities$R0[sp]) # - Scale mean rec for rec trend
      }
    }

    # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
    if(Rceattle$data_list$srr_fun != Rceattle$data_list$srr_pred_fun){
      if(sample_rec){ # Sample devs from hindcast
        rec_dev <- sample(x = (log(Rceattle$quantities$R) - log(Rceattle$quantities$R_hat))[sp, 1:hind_nyrs],
                          size = proj_nyrs, replace = TRUE) + log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
      } else{ # Set to mean rec otherwise
        rec_dev <- log(mean((log(Rceattle$quantities$R) - log(Rceattle$quantities$R_hat))[sp, 1:hind_nyrs]) * (1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
      }
    }

    # - Update OM with devs
    Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1] <- replace(
      Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1],
      values =  rec_dev)
  }

  if(update_model){
    # * Update fit ----
    estMode <- Rceattle$data_list$estimateMode
    Rceattle <-
      suppressWarnings(
        suppressMessages(
          fit_mod(
            data_list = Rceattle$data_list,
            inits = Rceattle$estimated_params,
            map =  NULL,
            bounds = NULL,
            file = NULL,
            estimateMode = 3,
            HCR = build_hcr(HCR = Rceattle$data_list$HCR,
                            DynamicHCR = Rceattle$data_list$DynamicHCR,
                            Ftarget = Rceattle$data_list$Ftarget,
                            Flimit = Rceattle$data_list$Flimit,
                            Ptarget = Rceattle$data_list$Ptarget,
                            Plimit = Rceattle$data_list$Plimit,
                            Alpha = Rceattle$data_list$Alpha,
                            Pstar = Rceattle$data_list$Pstar,
                            Sigma = Rceattle$data_list$Sigma,
                            Fmult = Rceattle$data_list$Fmult,
                            HCRorder = Rceattle$data_list$HCRorder
            ),
            recFun = build_srr(srr_fun = Rceattle$data_list$srr_fun,
                               srr_pred_fun  = Rceattle$data_list$srr_pred_fun,
                               proj_mean_rec  = Rceattle$data_list$proj_mean_rec,
                               srr_meanyr = Rceattle$data_list$srr_meanyr,
                               srr_hat_styr = Rceattle$data_list$srr_hat_styr,
                               srr_hat_endyr = Rceattle$data_list$srr_hat_endyr,
                               srr_est_mode  = Rceattle$data_list$srr_est_mode ,
                               srr_prior  = Rceattle$data_list$srr_prior,
                               srr_prior_sd   = Rceattle$data_list$srr_prior_sd,
                               Bmsy_lim = Rceattle$data_list$Bmsy_lim,
                               srr_indices = Rceattle$data_list$srr_indices),
            M1Fun =     build_M1(M1_model = Rceattle$data_list$M1_model,
                                 M1_re = Rceattle$data_list$M1_re,
                                 updateM1 = FALSE,  # Dont update M1 from data, fix at previous parameters
                                 M1_use_prior = Rceattle$data_list$M1_use_prior,
                                 M2_use_prior = Rceattle$data_list$M2_use_prior,
                                 M_prior = Rceattle$data_list$M_prior,
                                 M_prior_sd = Rceattle$data_list$M_prior_sd,
                                 M1_indices = Rceattle$data_list$M1_indices),
            growthFun = build_growth(growth_model = Rceattle$data_list$growth_model,
                                     growth_re = Rceattle$data_list$growth_re,
                                     growth_indices = Rceattle$data_list$growth_indices),
            random_rec = Rceattle$data_list$random_rec,
            niter = Rceattle$data_list$niter,
            msmMode = Rceattle$data_list$msmMode,
            avgnMode = Rceattle$data_list$avgnMode,
            suitMode = Rceattle$data_list$suitMode,
            suit_styr = Rceattle$data_list$suit_styr,
            suit_endyr = Rceattle$data_list$suit_endyr,
            initMode = Rceattle$data_list$initMode,
            phase = FALSE,
            loopnum = Rceattle$data_list$loopnum,
            getsd = TRUE,
            verbose = 0)
        )
      )
    Rceattle$data_list$estimateMode <- estMode
  }

  return(Rceattle)
}

#' Evaluate simulation performance
#'
#' @description Function to evaluate the simulation performance with regard to bias using the median relative error (MRE) and precision using the coefficient of variation.
#'
#' @param operating_mod CEATTLE model object exported from \code{\link{Rceattle}} to be used as the operating model
#' @param simulation_mods List of CEATTLE model objects exported from \code{\link{Rceattle}} fit to simulated data
compare_sim <- function(operating_mod, simulation_mods, object = "quantities") {
  # TODO update

  # Get differences
  sim_mre <- list()
  sim_mse <- list()
  sim_mean <- list()
  sim_median <- list()
  sim_sd <- list()
  sim_cv <- list()
  sim_params <- list()

  for (j in 1:length(names(operating_mod[[object]]))) {
    param <- names(operating_mod[[object]])[j]

    sim_mre[[param]] <- list()
    sim_mse[[param]] <- list()
    sim_mean[[param]] <- list()
    sim_sd[[param]] <- list()
    sim_cv[[param]] <- list()
    sim_params[[param]] <- list()

    om_params <- operating_mod[[object]][[param]]

    for (i in 1:length(simulation_mods)) {

      sm_params <- simulation_mods[[i]][[object]][[param]]

      sim_params[[param]][[i]] <- sm_params
      sim_mre[[param]][[i]] <- (sm_params - om_params)/om_params
      sim_mse[[param]][[i]] <- (sm_params - om_params)^2
    }

    param_dim <- length(dim(om_params))

    # If 1 value
    if (param_dim == 0) {
      sim_mean[[param]] <- mean(unlist(sim_params[[param]]))
      sim_sd[[param]] <- sd(unlist(sim_params[[param]]))
      sim_cv[[param]] <- sim_sd[[param]]/sim_mean[[param]]

      sim_mre[[param]] <- median(unlist(sim_mre[[param]]))
      sim_mse[[param]] <- mean(unlist(sim_mse[[param]]))
    }

    # If multiple values
    if (param_dim > 0) {
      # Get mean, sd, and CV
      sim_mean[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, mean)
      sim_median[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, median)
      sim_sd[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, sd)
      sim_cv[[param]] <- sim_sd[[param]]/sim_mean[[param]]

      sim_mre[[param]] <- apply(simplify2array(sim_mre[[param]]), 1:param_dim, median)
      sim_mse[[param]] <- apply(simplify2array(sim_mse[[param]]), 1:param_dim, mean)
    }
  }


  result_list <- list(Mean = sim_mean, Median = sim_median, SD = sim_sd, CV = sim_cv, MRE = sim_mre, MSE = sim_mse)
  return(result_list)
}

