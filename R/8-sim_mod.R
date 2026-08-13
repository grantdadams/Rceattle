#' Resolve `Index_distribution` to its integer family code
#'
#' Accepts either spelling a `fleet_control` column can hold -- the name
#' (`"MVN"`) or the code (`1`) -- because `index_ll_type` is built inside
#' `rearrange_data()` and is not carried on a fitted model's `data_list`.
#' Anything unrecognized falls back to lognormal, matching the column default.
#'
#' @param x The `Index_distribution` column.
#' @return Integer vector of codes, one per fleet. See `index_distribution_map`.
#' @noRd
.index_family_codes <- function(x) {
  if (is.null(x)) return(integer(0))
  chr <- trimws(as.character(x))
  num <- suppressWarnings(as.integer(chr))
  nmd <- as.integer(index_distribution_map[chr])
  out <- ifelse(!is.na(num), num, nmd)
  out[is.na(out)] <- 0L
  as.integer(out)
}


#' Draw survey-index observations under each fleet's own likelihood
#'
#' A simulator has to draw from the observation model the likelihood assumes, or
#' `self_test()` measures recovery under a data-generating process the estimation
#' model never saw -- and reports it as if it had. The four families need three
#' different draws (`ceattle.cpp`, section 8.2):
#'
#' * `Lognormal` (0): log scale, `exp(N(log(hat) - ba * sd^2 / 2, sd))`.
#' * `Normal` (3): natural scale with an ABSOLUTE sd, `N(hat, sd)`. No lognormal
#'   bias term -- the likelihood has none.
#' * `MVN` (1) / `MVNORM` (2): natural scale and correlated,
#'   `hat + chol(Sigma)' z`. The two differ only by a constant in the likelihood,
#'   so they simulate identically.
#'
#' Sigma is positional: its rows follow the fleet's *fitted* observations
#' (`Year` in `(0, endyr]`, `Observation > 0`) in `index_data` order -- the same
#' predicate `.align_index_cov()` and `rearrange_data()` use. Rows outside that
#' set are left as they are, because Sigma says nothing about them.
#'
#' A natural-scale draw can come out non-positive, which no index can be. Those
#' rows fail the `Observation > 0` gate on the refit and drop out of the fitted
#' set (`.align_index_cov()` re-keys Sigma around them, so nothing errors), which
#' would quietly shrink the data a self-test is scored on. Warn instead.
#'
#' @param data_list Data list being simulated into.
#' @param index_hat Predicted index, one per `index_data` row.
#' @param log_index_sd Per-row observation sd, as `ceattle.cpp` reports it: a
#'   log-scale sd for `Lognormal`, an absolute one for `Normal`.
#' @param ba_obs Observation bias-adjustment multiplier.
#' @return Numeric vector of simulated observations, in `index_data` row order.
#' @noRd
.sim_index_obs <- function(data_list, index_hat, log_index_sd, ba_obs) {
  idx <- data_list$index_data
  fc  <- data_list$fleet_control
  obs <- as.numeric(index_hat)
  if (is.null(fc) || is.null(idx) || !nrow(idx)) return(obs)

  fam <- .index_family_codes(fc$Index_distribution)
  if (!length(fam)) fam <- rep(0L, nrow(fc))

  # Every fleet lognormal -- the column default, and the only family that existed
  # before MVN/MVNORM/Normal. Draw all rows in one call, exactly as this did
  # before the per-family dispatch, so a seeded self_test(), jitter() or
  # run_mse(simulate_data = TRUE) on an existing model reproduces bit for bit.
  # The per-fleet loop below consumes the RNG stream in a different order, which
  # would silently change every simulated data set with more than one index fleet
  # while leaving single-survey models looking unaffected.
  if (all(fam == 0L)) {
    return(exp(stats::rnorm(
      n    = length(obs),
      mean = log(index_hat) - ba_obs * (log_index_sd^2) / 2,
      sd   = log_index_sd)))
  }

  nonpos <- character(0)

  for (i in seq_len(nrow(fc))) {
    code <- fc$Fleet_code[i]
    rows <- which(idx$Fleet_code == code)
    if (!length(rows)) next
    f <- fam[i]

    if (f %in% c(1L, 2L)) {
      fit_rows <- rows[idx$Year[rows] > 0 &
                         idx$Year[rows] <= data_list$endyr &
                         idx$Observation[rows] > 0]
      if (!length(fit_rows)) next
      nm    <- as.character(fc$Fleet_name[i])
      Sigma <- data_list$index_cov[[nm]]
      if (is.null(Sigma)) Sigma <- data_list$index_cov[[as.character(code)]]
      if (is.null(Sigma)) {
        stop(sprintf(paste0("Fleet '%s' uses an %s index likelihood but no covariance ",
                            "matrix was supplied in index_cov, so its observations ",
                            "cannot be simulated."),
                     nm, c("MVN", "MVNORM")[f]), call. = FALSE)
      }
      Sigma <- as.matrix(Sigma)
      Sigma <- (Sigma + t(Sigma)) / 2   # as rearrange_data() does
      if (nrow(Sigma) != length(fit_rows)) {
        stop(sprintf(paste0("Index covariance matrix for fleet '%s' is %d x %d but the ",
                            "fleet has %d fitted survey observations to simulate."),
                     nm, nrow(Sigma), ncol(Sigma), length(fit_rows)), call. = FALSE)
      }
      L <- tryCatch(t(chol(Sigma)), error = function(e) {
        stop(sprintf(paste0("Index covariance matrix for fleet '%s' is not positive ",
                            "definite (%s), so it cannot be simulated from."),
                     nm, conditionMessage(e)), call. = FALSE)
      })
      obs[fit_rows] <- index_hat[fit_rows] +
        as.numeric(L %*% stats::rnorm(length(fit_rows)))
      if (any(obs[fit_rows] <= 0)) nonpos <- c(nonpos, nm)

    } else if (f == 3L) {
      obs[rows] <- stats::rnorm(length(rows), mean = index_hat[rows],
                                sd = log_index_sd[rows])
      if (any(obs[rows] <= 0)) nonpos <- c(nonpos, as.character(fc$Fleet_name[i]))

    } else {
      obs[rows] <- exp(stats::rnorm(
        n    = length(rows),
        mean = log(index_hat[rows]) - ba_obs * (log_index_sd[rows]^2) / 2,
        sd   = log_index_sd[rows]))
    }
  }

  if (length(nonpos)) {
    warning("Simulated a non-positive survey index for fleet(s) ",
            paste(unique(nonpos), collapse = ", "),
            ". Those observations are dropped when the model is refitted, so a ",
            "self-test is scored on less data than the original fit. The ",
            "observation error is large relative to the index.", call. = FALSE)
  }
  obs
}


#' Simulate Rceattle data
#'
#' @description Simulates data used in Rceattle from the expected values estimated
#' from an existing Rceattle model. The variances and uncertainty are consistent
#' with those used in the operating model. The function simulates: survey biomass
#' (under the fleet's own \code{Index_distribution} -- lognormal, natural-scale
#' normal, or the correlated MVN/MVNORM draw from the supplied covariance),
#' catch-at-age/length composition (multinomial or dirichlet-multinomial), conditional-age-at-length (CAAL; multinomial or dirichlet-multinomial),
#' and total catch (log-normal).
#'
#' @param Rceattle A CEATTLE model object exported from \code{Rceattle}.
#' @param simulate Logical. If \code{TRUE}, simulates data from distributions.
#'   If \code{FALSE}, returns the expected values (hats).
#'
#' @return A \code{data_list} object containing the simulated or expected data
#'   values, formatted for use in \code{Rceattle}.
#' @export
#'
sim_mod <- function(Rceattle, simulate = FALSE) {
  dat_sim <- Rceattle$data_list
  quantities <- Rceattle$quantities

  # Observation bias-adjustment multiplier (default 1). The index/catch
  # lognormal likelihood fits to mean log(hat) - bias_adjust_obs * sigma^2/2
  # (ceattle.cpp, JNLL_INDEX / JNLL_CATCH), so the simulator has to apply the
  # SAME offset or the estimation model is fitted to data drawn from a different
  # mean than its own likelihood assumes -- a systematic bias in scale (and so in
  # catchability), not noise, which no number of simulations averages away.
  # Mirrors residuals.Rceattle(), which resolves the flag the same way.
  ba_obs <- dat_sim$bias_adjust_obs
  if (is.null(ba_obs) && !is.null(Rceattle$obj)) {
    ba_obs <- Rceattle$obj$env$data$bias_adjust_obs
  }
  if (is.null(ba_obs)) ba_obs <- 1
  ba_obs <- as.numeric(ba_obs)[1]


  # Indices of abundance/biomass ----
  log_index_sd <- quantities$log_index_sd
  index_hat <- quantities$index_hat

  if (simulate) {
    # Per fleet, under its own Index_distribution: lognormal, natural-scale
    # normal, or a correlated MVN/MVNORM draw from the supplied covariance.
    # Drawing every fleet as lognormal (as this did) simulates from an
    # observation model the likelihood does not assume, which a self-test then
    # reports as if it had.
    dat_sim$index_data$Observation <-
      .sim_index_obs(dat_sim, index_hat, log_index_sd, ba_obs)
  } else {
    # Expected value
    dat_sim$index_data$Observation <- index_hat
  }


  # Age/Length composition ----
  if(nrow(dat_sim$comp_data) > 0){
    comp_hat <- quantities$comp_hat
    # Proportion columns found by name; the number of identifying columns
    # ahead of them is not fixed.
    comp_cols <- .composition_cols(dat_sim$comp_data, "Comp_")
    for (obs in 1:nrow(dat_sim$comp_data)) {
      prob_vec <- comp_hat[obs, ]
      sum_prob <- sum(prob_vec, na.rm = TRUE)
      n_eff = dat_sim$comp_data$Sample_size[obs]
      flt <- dat_sim$comp_data$Fleet_code[obs]

      if (simulate && sum_prob > 0) {
        .ll <- as.character(dat_sim$fleet_control$Comp_distribution[flt])

        # --- Multinomial ---
        # "MultinomialAFSC" differs from "Multinomial" only in the likelihood's
        # normalisation, not in the sampling distribution.
        if(.ll %in% c("Multinomial", "MultinomialAFSC")){
          sim_comp <- rmultinom(n = 1, size = dat_sim$comp_data$Sample_size[obs], prob = prob_vec)

        # --- Dirichlet-multinomial ---
        } else if(.ll == "DirichletMultinomial"){
          # Theta is the overdispersion/precision parameter.
          theta <- exp(Rceattle$estimated_params$comp_weights[flt])

          # Dirichlet parameters alpha = theta * expected_proportions
          alpha <- prob_vec * theta

          # 1. Draw from Dirichlet
          # Using a Gamma distribution trick to get Dirichlet draws
          dir_draw <- rgamma(length(alpha), shape = alpha, rate = 1)
          dir_draw <- dir_draw / sum(dir_draw)

          # 2. Draw from Multinomial using the Dirichlet-adjusted probabilities
          sim_comp <- rmultinom(n = 1, size = n_eff, prob = dir_draw)

        } else {
          stop(sprintf(
            "sim_mod(): unsupported Comp_distribution '%s' for fleet %s.",
            .ll, dat_sim$fleet_control$Fleet_name[flt]), call. = FALSE)
        }

        dat_sim$comp_data[obs, comp_cols] <- as.vector(sim_comp)
      } else {
        dat_sim$comp_data[obs, comp_cols] <- prob_vec
      }
    }
  }


  # CAAL ----
  caal_hat <- quantities$caal_hat
  if(nrow(dat_sim$caal_data) > 0){
    caal_cols <- .composition_cols(dat_sim$caal_data, "CAAL_")
    for (obs in 1:nrow(dat_sim$caal_data)) {
      prob_vec <- caal_hat[obs, ]
      sum_prob <- sum(prob_vec, na.rm = TRUE)
      n_eff = dat_sim$caal_data$Sample_size[obs]
      flt <- dat_sim$caal_data$Fleet_code[obs]

      if (simulate && sum_prob > 0) {
        .ll <- as.character(dat_sim$fleet_control$CAAL_distribution[flt])

        # --- Multinomial ---
        if(.ll %in% c("Multinomial", "MultinomialAFSC")){
          sim_caal <- rmultinom(n = 1, size = dat_sim$caal_data$Sample_size[obs], prob = prob_vec)

        # --- Dirichlet-multinomial ---
        } else if(.ll == "DirichletMultinomial"){
          # Theta is the overdispersion/precision parameter.
          theta <- exp(Rceattle$estimated_params$caal_weights[flt])

          # Dirichlet parameters alpha = theta * expected_proportions
          alpha <- prob_vec * theta

          # 1. Draw from Dirichlet
          # Using a Gamma distribution trick to get Dirichlet draws
          dir_draw <- rgamma(length(alpha), shape = alpha, rate = 1)
          dir_draw <- dir_draw / sum(dir_draw)

          # 2. Draw from Multinomial using the Dirichlet-adjusted probabilities
          sim_caal <- rmultinom(n = 1, size = n_eff, prob = dir_draw)

        } else {
          stop(sprintf(
            "sim_mod(): unsupported CAAL_distribution '%s' for fleet %s.",
            .ll, dat_sim$fleet_control$Fleet_name[flt]), call. = FALSE)
        }

        dat_sim$caal_data[obs, caal_cols] <- as.vector(sim_caal)

      } else {
        dat_sim$caal_data[obs, caal_cols] <- prob_vec
      }
    }
  }



  # Catch ----
  log_catch_sd <- quantities$log_catch_sd
  catch_hat <- quantities$catch_hat

  if (simulate) {
    # Log-normal simulation with bias correction
    dat_sim$catch_data$Catch <- exp(stats::rnorm(
      n = length(dat_sim$catch_data$Catch),
      mean = log(catch_hat) - ba_obs * (log_catch_sd^2) / 2,
      sd = log_catch_sd
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

#' Sample historical recruitment deviations into the projection
#'
#' @param Rceattle CEATTLE model object exported from \code{Rceattle}
#' @param sample_rec Include resampled recruitment deviations from the hindcast in the OM projection. Resampled deviations are used rather than drawing from N(0, sigmaR) because the initial deviations bias R0 low. If FALSE, uses the mean recruitment deviation.
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
        # `log(R) - log(R_hat)` is already a log-scale deviation centred near
        # zero, so its mean is routinely negative and the outer log() this
        # used to take returned NaN. Take the mean deviation directly and add
        # the log trend, mirroring the sampling branch above.
        rec_dev <- mean((log(Rceattle$quantities$R) - log(Rceattle$quantities$R_hat))[sp, 1:hind_nyrs]) +
          log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
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
          # Rebuild so the resampled projection recruitment propagates into the
          # reported dynamics. .refit_like() carries the source model's whole
          # specification across the refit -- HCR, stock-recruit, M, growth, and
          # the catchability / selectivity / composition linkages.
          .refit_like(
            data_list    = Rceattle$data_list,
            inits        = Rceattle$estimated_params,
            estimateMode = 3,
            getsd        = TRUE)
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
#' @param operating_mod CEATTLE model object exported from \code{Rceattle} to be used as the operating model
#' @param simulation_mods List of CEATTLE model objects exported from \code{Rceattle} fit to simulated data
#' @param object character string specifying which part of the model to compare (default = "quantities")
#' @return A data frame summarizing simulation performance metrics
#' @export
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
#' @param growth_log_sd_sp Array. Dimensions (sex, 2).
#'   Log-SD of length: 1st param is SD at minage, 2nd param is SD at maxage.
#' @param growth_model_sp Integer. 1 = Von Bertalanffy, 2 = Richards.
#'
#' @return A 4D array of probabilities with dimensions (sex, age, length, year).
get_growth_matrix_r <- function(fracyr, nsex_sp, nages_sp, nlengths_sp, nyrs,
                                lengths_sp, minage_sp, maxage_sp,
                                growth_params_sp, growth_log_sd_sp, growth_model_sp) {

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
          diff <- growth_params_sp[s, y, 3] - length_at_age[s, a, y] # Linf - current size
          ages <- 0:(nages_sp)
          weight_a <- exp(-0.2 * ages)
          vals <- length_at_age[s, a, y] + (ages / nages_sp) * diff
          length_at_age[s, a, y] <- sum(vals * weight_a) / sum(weight_a)
        }

        # --- 3. SD Calculation ---
        sd1 <- exp(growth_log_sd_sp[s, 1])
        sda <- exp(growth_log_sd_sp[s, 2])

        if(current_age <= minage_sp) {
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
            fac1 <- (lengths_sp[2] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- stats::pnorm(fac1)
          } else if(l == nlengths_sp) {
            fac1 <- (lengths_sp[nlengths_sp] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- 1 - stats::pnorm(fac1)
          } else {
            fac1 <- (lengths_sp[l+1] - length_at_age[s, a, y]) / length_sd[s, a, y]
            fac2 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- stats::pnorm(fac1) - stats::pnorm(fac2)
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
#' @param length_at_age Array. Mean length at age from get_growth_matrix_r.
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


