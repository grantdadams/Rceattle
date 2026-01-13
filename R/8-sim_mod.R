#' Simulate Rceattle data
#'
#' @description Simulates data used in Rceattle from the expected values estimated
#' from an existing Rceattle model. The variances and uncertainty are consistent
#' with those used in the operating model. The function simulates: survey biomass
#' (log-normal), catch-at-age/length composition (multinomial), conditional-age-at-length (CAAL; multinomial),
#' and total catch (log-normal).
#'
#' @param Rceattle A CEATTLE model object exported from \code{\link{Rceattle}}.
#' @param simulate Logical. If \code{TRUE}, simulates data from distributions.
#'   If \code{FALSE}, returns the expected values (hats).
#'
#' @return A \code{data_list} object containing the simulated or expected data
#'   values, formatted for use in \code{Rceattle}.
#' @export
#'
sim_mod <- function(Rceattle, simulate = FALSE) {
  # TODO Options for simulation diet data: multinomial, sqrt-normal, dirichlet, multinomial
  dat_sim <- Rceattle$data_list
  quantities <- Rceattle$quantities


  # Indices of abundance/biomass ----
  ln_index_sd <- quantities$ln_index_sd
  index_hat <- quantities$index_hat

  if (simulate) {
    # Log-normal simulation with bias correction
    dat_sim$index_data$Observation <- exp(rnorm(
      n = length(index_hat),
      mean = log(index_hat) - (ln_index_sd^2) / 2,
      sd = ln_index_sd
    ))
  } else {
    # Expected value
    dat_sim$index_data$Observation <- index_hat
  }


  # Age/Length composition ----
  if(nrow(dat_sim$comp_data) > 0){
    comp_hat <- quantities$comp_hat
    for (obs in 1:nrow(dat_sim$comp_data)) {
      prob_vec <- comp_hat[obs, ]
      sum_prob <- sum(prob_vec, na.rm = TRUE)

      if (simulate && sum_prob > 0) {
        # --- Multinomial ---
        # FIXME: add Dirichlet-multinomial option
        sim_comp <- rmultinom(n = 1, size = dat_sim$comp_data$Sample_size[obs], prob = prob_vec)
        dat_sim$comp_data[obs, 9:ncol(dat_sim$comp_data)] <- as.vector(sim_comp)
      } else {
        dat_sim$comp_data[obs, 9:ncol(dat_sim$comp_data)] <- prob_vec
      }
    }
  }


  # CAAL ----
  caal_hat <- quantities$caal_hat
  if(nrow(dat_sim$caal_data) > 0){
    for (obs in 1:nrow(dat_sim$caal_data)) {
      prob_vec <- caal_hat[obs, ]
      sum_prob <- sum(prob_vec, na.rm = TRUE)

      if (simulate && sum_prob > 0) {
        # FIXME: add Dirichlet-multinomial option
        sim_comp <- rmultinom(n = 1, size = dat_sim$caal_data$Sample_size[obs], prob = prob_vec)
        dat_sim$caal_data[obs, 7:ncol(dat_sim$caal_data)] <- as.vector(sim_comp)
      } else {
        dat_sim$caal_data[obs, 7:ncol(dat_sim$caal_data)] <- prob_vec
      }
    }
  }



  # Catch ----
  ln_catch_sd <- quantities$ln_catch_sd
  catch_hat <- quantities$catch_hat

  if (simulate) {
    # Log-normal simulation with bias correction
    dat_sim$catch_data$Catch <- exp(rnorm(
      n = length(dat_sim$catch_data$Catch),
      mean = log(catch_hat) - (ln_catch_sd^2) / 2,
      sd = ln_catch_sd
    ))
  } else {
    # Expected values
    dat_sim$catch_data$Catch <- catch_hat
  }


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


#' Generate Length-at-Age Transition Matrix
#'
#' This function calculates a probability transition matrix that defines the
#' probability of a fish of a given age belonging to specific length bins.
#' It supports Von Bertalanffy and Richards growth models and includes
#' a Stock Synthesis (SS) style plus-group correction.
#'
#' @param fracyr Numeric. Fraction of the year (0 = Jan 1st).
#' @param nsex_sp Integer. Number of sexes for the species.
#' @param nages_sp Integer. Number of age classes.
#' @param nlengths_sp Integer. Number of length bins.
#' @param nyrs Integer. Number of years in the simulation.
#' @param lengths_sp Vector. Boundaries of the length bins.
#' @param minage_sp Numeric. The reference age (L1) for growth estimation.
#' @param maxage_sp Numeric. The age at which growth enters the asymptotic phase.
#' @param growth_params_sp Array. Dimensions (sex, yr, 4).
#'   Params: K, L1, Linf, Richards m.
#' @param growth_ln_sd_sp Array. Dimensions (sex, 2).
#'   Log-SD of length: 1st param is SD at minage, 2nd param is SD at maxage.
#' @param growth_model_sp Integer. 1 = Von Bertalanffy, 2 = Richards.
#'
#' @return A 4D array of probabilities with dimensions (sex, age, length, year).
get_growth_matrix_r <- function(fracyr, nsex_sp, nages_sp, nlengths_sp, nyrs,
                                lengths_sp, minage_sp, maxage_sp,
                                growth_params_sp, growth_ln_sd_sp, growth_model_sp) {

  # Define names for the dimensions
  dim_names <- list(
    sex    = paste0("Sex_", 1:nsex_sp),
    age    = paste0("Age_", 1:nages_sp),
    length = paste0("Len_", lengths_sp),
    year   = paste0("Year_", 1:nyrs)
  )

  # Initialize Output: (sex, age, ln, yr)
  growth_matrix <- array(0, dim = c(nsex_sp, nages_sp, nlengths_sp, nyrs),
                         dimnames = dim_names)
  length_at_age <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
                         dimnames =   list(
                           sex    = paste0("Sex_", 1:nsex_sp),
                           age    = paste0("Age_", 0:(nages_sp - 1)),
                           year   = paste0("Year_", 1:nyrs)
                         ))
  length_sd     <- array(0, dim = c(nsex_sp, nages_sp, nyrs))

  l_min <- lengths_sp[1]
  l_max <- lengths_sp[nlengths_sp]

  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      # --- 1. Calculate Mean Length at Age ---
      # Params: 1:K, 2:L1, 3:Linf, 4:Richards_m
      k    <- growth_params_sp[s, y, 1]
      l1   <- growth_params_sp[s, y, 2]
      linf <- growth_params_sp[s, y, 3]

      b_len <- (l1 - l_min) / minage_sp

      for(a in 1:nages_sp) {
        current_age <- a + fracyr

        if (growth_model_sp == 1) { # VB
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            if(y == 1){
              length_at_age[s, a, y] <- linf + (l1 - linf) * (exp(-k * (current_age - minage_sp)))
            } else{
              if(a == nages_sp){ # linear growth + growth equation
                last_linear = l_min + b_len * minage_sp # last age (cont) with linear growth

                length_at_age[s, a, y] = last_linear + (last_linear - linf) * (exp(-k * (current_age - minage_sp)) - 1.0)
              }else{
                length_at_age[s, a, y] <- length_at_age[s, a-1, y-1] + (length_at_age[s, a-1, y-1] - growth_params_sp[s, y-1, 3]) * (exp(-growth_params_sp[s, y-1, 1]) - 1)
              }
            }
          }
        } else if (growth_model_sp == 2) { # Richards
          m <- growth_params_sp[s, y, 4]
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            if(y == 1){
              length_at_age[s, a, y] <- (linf^m + (l1^m - linf^m) * (exp(-k * (current_age - minage_sp))))^(1/m)
            } else{
              if(a == nages_sp){ # linear growth + growth equation
                last_linear = l_min + b_len * minage_sp # last age (cont) with linear growth

                length_at_age[s, a, y] = (last_linear^m + (last_linear^m - linf^m) * (exp(-k * (current_age - minage_sp)) - 1.0))^(1/m)
              }else{
                lagk <- growth_params_sp[s, y-1, 2]
                lagm <- growth_params_sp[s, y-1, 4]
                laglinf <-  growth_params_sp[s, y-1, 3]
                length_at_age[s, a, y] <- (length_at_age[s, a-1, y-1]^lagm + (length_at_age[s, a-1, y-1]^lagm - laglinf^lagm) * (exp(-lagk) - 1))^1/lagm
              }
            }
          }
        }

        # --- 2. Plus Group Correction (SS Style) ---
        if(a == nages_sp) {
          diff <- growth_params_sp[s, y, 3] - length_at_age[s, a, y]
          ages <- 0:(nages_sp-1)
          weight_a <- exp(-0.2 * ages)
          vals <- length_at_age[s, a, y] + (ages / nages_sp) * diff
          length_at_age[s, a, y] <- sum(vals * weight_a) / sum(weight_a)
        }

        # --- 3. SD Calculation ---
        sd1 <- exp(growth_ln_sd_sp[s, 1])
        sda <- exp(growth_ln_sd_sp[s, 2])

        if(current_age < minage_sp) {
          length_sd[s, a, y] <- sd1
        } else if(a == nages_sp) {
          length_sd[s, a, y] <- sda
        } else {
          slope <- (sda - sd1) / (linf - l1) # Match C++ interpolation
          length_sd[s, a, y] <- sd1 + slope * (length_at_age[s, a, y] - l1)
        }

        # --- 4. Matrix Distribution ---
        for(l in 1:nlengths_sp) {
          if(l == 1) {
            fac1 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- pnorm(fac1)
          } else if(l == nlengths_sp) {
            fac1 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- 1 - pnorm(fac1)
          } else {
            fac1 <- (lengths_sp[l+1] - length_at_age[s, a, y]) / length_sd[s, a, y]
            fac2 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- pnorm(fac1) - pnorm(fac2)
          }
        }
      }
    }
  }

  return(list(length_at_age = length_at_age, growth_matrix = growth_matrix))
}


#' Calculate Predicted Weight-at-Age
#'
#' Converts a growth matrix (length-at-age probabilities) into mean weight-at-age
#' using a length-weight relationship (W = a * L^b).
#'
#' @param nsex_sp Integer. Number of sexes.
#' @param nages_sp Integer. Number of age classes.
#' @param nlengths_sp Integer. Number of length bins.
#' @param nyrs Integer. Number of years.
#' @param lengths_sp Vector. Boundaries of the length bins.
#' @param growth_matrix Array. 4D array (sex, age, length, year) from get_growth_matrix_r.
#' @param lw_params Array. Dimensions (sex, yr, 2).
#'   Params: 1st is alpha (a), 2nd is beta (b).
#'
#' @details The function calculates midpoints for length bins to avoid bias.
#' For the first bin, it assumes the width is equal to the second bin's width.
#' The final weight-at-age is the expected value across all length bins for that age.
#'
#' @return A 3D array of mean weights with dimensions (sex, age, year).
get_weight_at_age_r <- function(nsex_sp, nages_sp, nlengths_sp, nyrs,
                                lengths_sp, length_at_age, growth_matrix, lw_params) {
  # Define names for the dimensions
  dim_names <- list(
    sex  = paste0("Sex_", 1:nsex_sp),
    age  = paste0("Age_", 1:nages_sp),
    year = paste0("Year_", 1:nyrs)
  )

  # Output: (sex, age, yr)
  waa <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
               dimnames = dim_names)


  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      alpha <- lw_params[s, y, 1]
      beta  <- lw_params[s, y, 2]

      # Weight at length for all bins
      wal <- alpha * (lengths_sp + (lengths_sp[2] - lengths_sp[1])/2) ^ beta

      for(a in 1:nages_sp) {
        # Matrix multiply: Prob(length | age) * Weight(length)
        waa[s, a, y] <- sum(growth_matrix[s, a, , y] * wal)
      }
    }
  }
  return(waa)
}


