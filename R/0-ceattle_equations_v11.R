#' Calculate Fishing Mortality and FSPRs
#'
#' This function calculates fishing mortality and fishing selectivity for multiple fleets and species over a specified number of years.
#'
#' @param n_flt Integer. The number of fleets.
#' @param nspp Integer. The number of species.
#' @param max_nsex Integer. The number of sexes considered for each species.
#' @param max_nages Integer. The number of age classes considered for each species.
#' @param nyrs Integer. The total number of years for the analysis.
#' @param nyrs_hind Integer. The number of years for hindcasting.
#' @param flt_spp Integer vector. A vector indicating the species associated with each fleet.
#' @param flt_type Integer vector. A vector indicating the type of each fleet (1 for active fleets).
#' @param sel Array. An array of selectivity values with dimensions (n_flt, max_nsex, max_nages, nyrs).
#' @param ln_mean_F Numeric vector. The natural logarithm of mean fishing mortality for each fleet.
#' @param F_dev Array. An array of fishing mortality deviations with dimensions (n_flt, nyrs).
#' @param Ftarget Numeric vector. The target fishing mortality for each species.
#' @param Fmult Numeric vector. Multipliers for fishing mortality for each species.
#' @param Flimit Numeric vector. Limit fishing mortality for each species.
#' @param QnormHCR Numeric vector. HCR normalization values for each species.
#' @param forecast Integer vector. A vector indicating whether a forecast is being run for each species (0 or 1).
#' @param proj_F_prop Numeric vector. Projected fishing proportions for each fleet.
#' @param nages
#' @param nsex
#' @param HCR
#'
#' @return A list containing:
#' \item{F_spp}{Array of total fishing mortality by species and year.}
#' \item{F_flt_age}{Array of fishing mortality by fleet, sex, age, and year.}
#' \item{F_spp_age}{Array of total fishing mortality by species, sex, age, and year.}
#' \item{Ftarget_age_spp}{Array of target fishing mortality by species, sex, age, and year.}
#' \item{Flimit_age_spp}{Array of limit fishing mortality by species, sex, age, and year.}
#'
#' @examples
#' # Example usage of the function
#' results <- calculate_fishing_mortality(n_flt = 5, nspp = 3, nyrs = 10,
#'                                          max_nsex = 2, max_nages = 5, nyrs = 20,
#'                                          nyrs_hind = 10, flt_spp = c(1, 2, 1, 3, 2),
#'                                          flt_type = c(1, 1, 1, 1, 1),
#'                                          sel = array(runif(5 * 2 * 5 * 10), dim = c(5, 2, 5, 10)),
#'                                          ln_mean_F = runif(5),
#'                                          F_dev = array(runif(5 * 10), dim = c(5, 10)),
#'                                          Ftarget = runif(3),
#'                                          Fmult = runif(3),
#'                                          Flimit = runif(3),
#'                                          QnormHCR = runif(3),
#'                                          forecast = c(1, 1, 0),
#'                                          proj_F_prop = runif(5))
#'
calculate_fishing_mortality <- function(n_flt, nspp, nages, max_nages, nsex, max_nsex, nyrs, nyrs_hind,
                                        flt_spp, flt_type, sel, ln_mean_F, F_dev, Ftarget, Fmult,
                                        Flimit, QnormHCR, forecast, proj_F_prop, HCR) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize arrays
  F_flt <- array(0, dim = c(n_flt, nyrs))
  F_spp <- array(0, dim = c(nspp, nyrs))
  F_flt_age <- array(0, dim = c(n_flt, max_nsex, max_nages, nyrs))
  F_spp_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  Ftarget_age_spp <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  Flimit_age_spp <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  proj_F <- array(0, dim = c(nspp, nyrs)) # Initialize proj_F array

  for (flt in 1:n_flt) {
    sp <- flt_spp[flt] # Temporary index of fishery survey

    if (flt_type[flt] == 1) {
      for (age in 1:nages[sp]) {
        for (sex in 1:nsex[sp]) {
          for (yr in 1:nyrs) {

            # Hindcast
            if (yr <= nyrs_hind) {
              F_flt_age[flt, sex, age, yr] <- sel[flt, sex, age, yr] * exp(ln_mean_F[flt] + F_dev[flt, yr])
            }

            # Forecast
            if (yr > nyrs_hind) {
              # Apply HCRs
              proj_F[sp, yr] <- switch(HCR + 1,
                                       0, # No fishing
                                       Ftarget[sp], # CMSY
                                       Ftarget[sp], # Constant F
                                       Ftarget[sp], # Constant F that achieves X% of SSB0
                                       Ftarget[sp] * Fmult[sp], # Constant Fspr
                                       Ftarget[sp], # NPFMC Tier 3 HCR
                                       Flimit[sp] + QnormHCR[sp], # PFMC Category 1 HCR
                                       Ftarget[sp] # SESSF Tier 1 HCR
              )

              # Set F to zero if not running forecast
              if (forecast[sp] == 0) {
                proj_F[sp, yr] <- 0
              }

              F_flt_age[flt, sex, age, yr] <- sel[flt, sex, age, yr] * proj_F_prop[flt] * proj_F[sp, yr]
            }

            # Sum F across fleets
            F_spp_age[sp, sex, age, yr] <- F_spp_age[sp, sex, age, yr] + F_flt_age[flt, sex, age, yr]
          }
        }
      }

      # F across fleets or species
      for (yr in 1:nyrs) {
        # Hindcast
        if (yr <= nyrs_hind) {
          F_flt[flt, yr] <- exp(ln_mean_F[flt] + F_dev[flt, yr])
          F_spp[sp, yr] <- F_spp[sp, yr] + exp(ln_mean_F[flt] + F_dev[flt, yr])
        }

        # Forecast
        if (yr > nyrs_hind) {
          F_flt[flt, yr] <- proj_F_prop[flt] * proj_F[sp, yr]
          F_spp[sp, yr] <- F_spp[sp, yr] + proj_F_prop[flt] * proj_F[sp, yr]
        }
      }


      # Calculate F target by age and sex for reference points
      Flimit_age_spp[sp,,,] <- Flimit_age_spp[sp,,,] +
        sel[flt,,,] * proj_F_prop[flt] * Flimit[sp]
      Ftarget_age_spp[sp,,,] <- Ftarget_age_spp[sp,,,] +
        sel[flt,,,] * proj_F_prop[flt] * Ftarget[sp]
    }
  }

  # Report
  RTMB::REPORT(F_spp)
  RTMB::REPORT(F_spp_age)
  RTMB::REPORT(F_flt)
  RTMB::REPORT(F_flt_age)
  RTMB::REPORT(Ftarget_age_spp)
  RTMB::REPORT(Flimit_age_spp)

  # Return results as a list
  return(list(F_spp = F_spp, F_spp_age = F_spp_age, F_flt = F_flt, F_flt_age = F_flt_age,
              Ftarget_age_spp = Ftarget_age_spp, Flimit_age_spp = Flimit_age_spp,
              proj_F = proj_F))
}


#' Calculate SPR-Based Reference Points
#'
#' This function calculates the spawning potential ratio (SPR) based reference points for multiple species.
#' It computes total mortality-at-age and SPR values based on the provided parameters and data.
#'
#' @param nspp Number of species.
#' @param nages A vector containing the number of ages for each species.
#' @param nyrs Number of years.
#' @param nyrs_hind Number of years for hindcasting.
#' @param initMode Initialization mode for the population.
#' @param R_sexr A vector of recruitment values for each species.
#' @param M1_at_age A function that returns natural mortality at age.
#' @param M2_at_age A function that returns additional mortality at age.
#' @param Flimit_age_spp A matrix of limit fishing mortality at age for each species.
#' @param Ftarget_age_spp A matrix of target fishing mortality at age for each species.
#' @param wt A function that returns weight at age.
#' @param ssb_wt_index An index for spawning stock biomass weight.
#' @param pmature A function that returns the proportion mature at age.
#' @param spawn_month A vector of spawning months for each species.
#' @param max_nages
#' @param Finit
#'
#' @return A list containing the calculated SPR values: SPR0, SPRFinit, SPRlimit, and SPRtarget.
#'
#' @examples
#' # Example usage
#' result <- calculate_spr_reference_points(nspp, nages, nyrs, nyrs_hind,
#'                                           initMode, R_sexr, M1_at_age, M2_at_age,
#'                                           Flimit_age_spp, Ftarget_age_spp,
#'                                           wt, ssb_wt_index, pmature, spawn_month,
#'                                           Finit)
#'
calculate_spr_reference_points <- function(nspp, nages, max_nages, nyrs, nyrs_hind,
                                           initMode, R_sexr, M1_at_age, M2_at_age,
                                           Flimit_age_spp, Ftarget_age_spp,
                                           wt, ssb_wt_index, pmature, spawn_month,
                                           Finit) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize NbyageSPR array with appropriate dimensions
  NbyageSPR <- array(0, dim = c(4, nspp, max_nages))  # Assuming max_nages is the maximum age across species
  SPR0 <- rep(0, nspp)
  SPRFinit <- rep(0, nspp)
  SPRlimit <- rep(0, nspp)
  SPRtarget <- rep(0, nspp)

  for (sp in 1:nspp) {

    if(initMode != 3) {
      Finit[sp] <- 0  # If population starts out at equilibrium set Finit to 0 (R_init and R0 will be the same)
    }

    NbyageSPR[1, sp, 1] <- R_sexr[sp]  # F = 0
    NbyageSPR[2, sp, 1] <- R_sexr[sp]  # F = Flimit
    NbyageSPR[3, sp, 1] <- R_sexr[sp]  # F = Ftarget
    NbyageSPR[4, sp, 1] <- R_sexr[sp]  # F = Finit

    for (age in 2:(nages[sp]-1)) {
      NbyageSPR[1, sp, age] <- NbyageSPR[1, sp, age-1] * exp(-M1_at_age[sp, 1, age-1] - M2_at_age[sp, 1, age-1, nyrs_hind])
      NbyageSPR[2, sp, age] <- NbyageSPR[2, sp, age-1] * exp(-(M1_at_age[sp, 1, age-1] + Flimit_age_spp[sp, 1, age-1, nyrs_hind]))
      NbyageSPR[3, sp, age] <- NbyageSPR[3, sp, age-1] * exp(-(M1_at_age[sp, 1, age-1] + Ftarget_age_spp[sp, 1, age-1, nyrs_hind]))
      NbyageSPR[4, sp, age] <- NbyageSPR[4, sp, age-1] * exp(-(M1_at_age[sp, 1, age-1] + Finit[sp]))
    }

    # Plus group
    NbyageSPR[1, sp, nages[sp]] <- NbyageSPR[1, sp, nages[sp]-1] * exp(-M1_at_age[sp, 1, nages[sp]]) /
      (1 - exp(-M1_at_age[sp, 1, nages[sp]]))

    NbyageSPR[2, sp, nages[sp]] <- NbyageSPR[2, sp, nages[sp]-1] *
      exp(-(M1_at_age[sp, 1, nages[sp]] + Flimit_age_spp[sp, 1, nages[sp]-1, nyrs_hind])) /
      (1 - exp(-(M1_at_age[sp, 1, nages[sp]] + Flimit_age_spp[sp, 1, nages[sp], nyrs_hind])))

    NbyageSPR[3, sp, nages[sp]] <- NbyageSPR[3, sp, nages[sp]-1] *
      exp(-(M1_at_age[sp, 1, nages[sp]] + Ftarget_age_spp[sp, 1, nages[sp]-1, nyrs_hind])) /
      (1 - exp(-(M1_at_age[sp, 1, nages[sp]] + Ftarget_age_spp[sp, 1, nages[sp], nyrs_hind])))

    NbyageSPR[4, sp, nages[sp]] <- NbyageSPR[4, sp, nages[sp]-1] *
      exp(-(M1_at_age[sp, 1, nages[sp]] + Finit[sp])) /
      (1 - exp(-(M1_at_age[sp, 1, nages[sp]] + Finit[sp])))

    # Calculate SPR
    for (age in 1:nages[sp]) {
      SPR0[sp] <- SPR0[sp] + NbyageSPR[1, sp, age] * wt[ssb_wt_index[sp], 1, age, nyrs_hind] *
        pmature[sp, age] * exp(-M1_at_age[sp, 1, age] * spawn_month[sp]/12)

      SPRlimit[sp] <- SPRlimit[sp] + NbyageSPR[2, sp, age] * wt[ssb_wt_index[sp], 1, age, nyrs_hind] *
        pmature[sp, age] * exp(-(M1_at_age[sp, 1, age] + Flimit_age_spp[sp, 1, age, nyrs_hind]) * spawn_month[sp]/12)

      SPRtarget[sp] <- SPRtarget[sp] + NbyageSPR[3, sp, age] * wt[ssb_wt_index[sp], 1, age, nyrs_hind] *
        pmature[sp, age] * exp(-(M1_at_age[sp, 1, age] + Ftarget_age_spp[sp, 1, age, nyrs_hind]) * spawn_month[sp]/12)

      SPRFinit[sp] <- SPRFinit[sp] + NbyageSPR[4, sp, age] * wt[ssb_wt_index[sp], 1, age, 1] *
        pmature[sp, age] * exp(-(M1_at_age[sp, 1, age] + Finit[sp]) * spawn_month[sp]/12)
    }
  }

  # Report
  RTMB::REPORT(SPR0)
  RTMB::REPORT(SPRFinit)
  RTMB::REPORT(SPRlimit)
  RTMB::REPORT(SPRtarget)

  # Return results as a list
  return(list(SPR0 = SPR0, # Estimated spawning biomass per recruit at F = 0
              SPRFinit = SPRFinit, SPRlimit = SPRlimit, SPRtarget = SPRtarget))
}


#' Calculate Stock-Recruit Parameters
#'
#' @description
#' Calculates steepness and R0 parameters for different stock-recruit relationships
#'
#' @param nspp Number of species
#' @param srr_fun Stock-recruit relationship function type (0-5)
#' @param rec_pars Matrix of recruitment parameters
#' @param beta_rec_pars Matrix of beta parameters for environmental effects
#' @param env_index_srr Matrix of environmental indices
#' @param SPR0 Vector of unfished spawning potential ratio
#' @param SPRFinit Vector of initial spawning potential ratio
#'
#' @return List containing:
#' \itemize{
#'   \item steepness - Vector of steepness parameters
#'   \item R0 - Vector of unfished recruitment
#'   \item R_init - Vector of initial recruitment
#'   \item zero_N_pen - Vector of penalties
#' }
#'
#' @details
#' Implements 6 different stock-recruit relationships:
#' \itemize{
#'   \item 0: Random about mean
#'   \item 1: Random about mean with environmental linkage
#'   \item 2: Beverton-Holt
#'   \item 3: Beverton-Holt with environmental impacts
#'   \item 4: Ricker
#'   \item 5: Ricker with environmental impacts
#' }
calculate_sr_parameters <- function(nspp, srr_fun, rec_pars, beta_rec_pars,
                                    env_index_srr, SPR0, SPRFinit) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize outputs
  steepness <- numeric(nspp)
  R0 <- numeric(nspp)
  R_init <- numeric(nspp)
  zero_N_pen <- numeric(nspp)
  penalty <- 0.0

  for(sp in 1:nspp) {

    switch(srr_fun + 1, # R switch is 1-based

           { # Random about mean
             steepness[sp] <- 0.99
             R_init[sp] <- R0[sp] <- exp(rec_pars[sp, 1])
           },

           { # Random about mean with environmental linkage
             steepness[sp] <- 0.99
             beta_rec_tmp <- beta_rec_pars[sp,]
             env_rec_tmp <- env_index_srr[1,]
             srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
             R_init[sp] <- R0[sp] <- exp(rec_pars[sp, 1] + srr_mult)
           },

           { # Beverton-Holt
             steepness[sp] <- exp(rec_pars[sp, 2]) * SPR0[sp]/(4.0 + exp(rec_pars[sp, 2]) * SPR0[sp])
             R0[sp] <- (exp(rec_pars[sp, 2])-1.0/SPR0[sp]) / exp(rec_pars[sp, 3])
             R_init[sp] <- (exp(rec_pars[sp, 2])-1/SPRFinit[sp]) / exp(rec_pars[sp, 3])
           },

           { # Beverton-Holt with environmental impacts
             beta_rec_tmp <- beta_rec_pars[sp,]
             env_rec_tmp <- env_index_srr[1,]
             srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
             srr_alpha <- exp(rec_pars[sp, 2] + srr_mult)
             steepness[sp] <- srr_alpha * SPR0[sp]/(4.0 + srr_alpha * SPR0[sp])
             R0[sp] <- (srr_alpha-1.0/SPR0[sp]) / exp(rec_pars[sp, 3])
             R_init[sp] <- (srr_alpha-1.0/SPRFinit[sp]) / exp(rec_pars[sp, 3])
           },

           { # Ricker
             steepness[sp] <- 0.2 * exp(0.8*log(exp(rec_pars[sp, 2]) * SPR0[sp]))

             # R at F0
             ricker_intercept <- exp(rec_pars[sp, 2]) * SPR0[sp] - 1.0
             # pos_tmp <- posfun(ricker_intercept)
             # ricker_intercept <- pos_tmp$ans + 1.0
             # zero_N_pen[sp] <- zero_N_pen[sp] + pos_tmp$penalty
             R0[sp] <- log(ricker_intercept)/(exp(rec_pars[sp, 3]) * SPR0[sp]/1000000.0)

             # R at equilibrium F
             ricker_intercept <- exp(rec_pars[sp, 2]) * SPRFinit[sp] - 1.0
             # pos_tmp <- posfun(ricker_intercept)
             # ricker_intercept <- pos_tmp$ans + 1.0
             # zero_N_pen[sp] <- zero_N_pen[sp] + pos_tmp$penalty
             R_init[sp] <- log(ricker_intercept)/(exp(rec_pars[sp, 3]) * SPRFinit[sp]/1000000.0)

           },

           { # Ricker with environmental impacts
             beta_rec_tmp <- beta_rec_pars[sp,]
             env_rec_tmp <- env_index_srr[1,]
             srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
             srr_alpha <- exp(rec_pars[sp, 2] + srr_mult)
             steepness[sp] <- 0.2 * exp(0.8*log(srr_alpha * SPR0[sp]))

             ricker_intercept <- srr_alpha * SPR0[sp] - 1.0
             # pos_tmp <- posfun(ricker_intercept)
             # ricker_intercept <- pos_tmp$ans + 1.0
             # zero_N_pen[sp] <- zero_N_pen[sp] + pos_tmp$penalty
             R0[sp] <- log(ricker_intercept)/(exp(rec_pars[sp, 3]) * SPR0[sp]/1000000.0)

             ricker_intercept <- srr_alpha * SPRFinit[sp] - 1.0
             # pos_tmp <- posfun(ricker_intercept)
             # ricker_intercept <- pos_tmp$ans + 1.0
             # zero_N_pen[sp] <- zero_N_pen[sp] + pos_tmp$penalty
             R_init[sp] <- log(ricker_intercept)/(exp(rec_pars[sp, 3]) * SPRFinit[sp]/1000000.0)
           },

           stop("Invalid srr_fun value")
    )
  }

  # Report
  RTMB::REPORT(steepness)
  RTMB::REPORT(R0)
  RTMB::REPORT(R_init)
  RTMB::REPORT(zero_N_pen)

  # Return results as a list
  return(list(
    steepness = steepness,
    R0 = R0,
    R_init = R_init, # Equilibrium recruitment at F = Finit (non-equilibrium).
    zero_N_pen = zero_N_pen
  ))
}


#' Calculate Expected Recruitment
#'
#' This function calculates the expected recruitment for a given number of species
#' based on specified recruitment functions and environmental effects.
#'
#' @param nspp Number of species.
#' @param nyrs Number of years.
#' @param srr_pred_fun Recruitment function type (0-5).
#' @param R0 Initial recruitment values for each species.
#' @param SPRFinit Initial spawning potential ratio for each species.
#' @param beta_rec_pars Matrix of beta recruitment parameters.
#' @param env_index_srr Matrix of environmental indices for recruitment.
#' @param rec_pars Matrix of recruitment parameters.
#' @param ssb Matrix of spawning stock biomass values.
#' @param minage Vector of minimum ages for each species.
#'
#' @return A matrix of expected recruitment values.
calculate_expected_recruitment <- function(nspp, nyrs, srr_pred_fun, R0, SPRFinit, beta_rec_pars, env_index_srr, rec_pars, ssb, minage) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize the R_hat matrix
  R_hat <- matrix(0, nrow = nspp, ncol = nyrs)

  for (sp in 1:nspp) {

    # Year 1 (aren't fit in likelihood)
    R_hat[sp, 1] <- switch(
      srr_pred_fun + 1,
      # Random about mean (e.g. Alaska)
      R0[sp],

      # Random about mean with environmental effects
      {
        beta_rec_tmp <- beta_rec_pars[sp,]
        env_rec_tmp <- env_index_srr[1,]
        srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
        R0[sp] * exp(srr_mult)
      },

      # Beverton-Holt
      (exp(rec_pars[sp, 1]) - 1/SPRFinit[sp]) / exp(rec_pars[sp, 2]), # (Alpha-1/SPR0)/beta

      # Beverton-Holt with environmental impacts on alpha
      {
        beta_rec_tmp <- beta_rec_pars[sp,]
        env_rec_tmp <- env_index_srr[1,]
        srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
        srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
        (srr_alpha - 1/SPRFinit[sp]) / exp(rec_pars[sp, 2]) # (Alpha-1/SPR0)/beta
      },

      # Ricker
      log(exp(rec_pars[sp, 1]) * SPRFinit[sp]) / (exp(rec_pars[sp, 2]) * SPRFinit[sp] / 1000000), # FIXME - make time-varying

      # Ricker with environmental impacts on alpha
      {
        beta_rec_tmp <- beta_rec_pars[sp,]
        env_rec_tmp <- env_index_srr[1,]
        srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
        srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
        log(srr_alpha * SPRFinit[sp]) / (exp(rec_pars[sp, 2]) * SPRFinit[sp] / 1000000) # FIXME - make time-varying
      },

      stop("Invalid 'srr_pred_fun'")
    )

    # Year 1+
    for (yr in 2:nyrs) {
      R_hat[sp, yr] <- switch(
        srr_pred_fun + 1,
        # Random about mean (e.g. Alaska)
        R0[sp],

        # Random about mean with environmental effects
        {
          beta_rec_tmp <- beta_rec_pars[sp,]
          env_rec_tmp <- env_index_srr[yr,]
          srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
          R0[sp] * exp(srr_mult)
        },

        # Beverton-Holt
        exp(rec_pars[sp, 1]) * ssb[sp, yr - minage[sp]] / (1 + exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]]),

        # Beverton-Holt with environmental impacts on alpha
        {
          beta_rec_tmp <- beta_rec_pars[sp,]
          env_rec_tmp <- env_index_srr[yr,]
          srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
          srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
          srr_alpha * ssb[sp, yr - minage[sp]] / (1 + exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]])
        },

        # Ricker
        exp(rec_pars[sp, 1]) * ssb[sp, yr - minage[sp]] * exp(-exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]] / 1000000), # Divide by 1e6 for numerical stability

        # Ricker with environmental impacts on alpha
        {
          beta_rec_tmp <- beta_rec_pars[sp,]
          env_rec_tmp <- env_index_srr[yr,]
          srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
          srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
          srr_alpha * ssb[sp, yr - minage[sp]] * exp(-exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]] / 1000000)
        },

        stop("Invalid 'srr_pred_fun'")
      )
    }
  }

  # Report
  RTMB::REPORT(R_hat)

  return(R_hat)
}


#' Calculate Population Dynamics
#'
#' This function calculates the population dynamics for multiple species over a specified number of years,
#' including initial abundance, recruitment, and spawning biomass.
#'
#' @param nspp Number of species.
#' @param nyrs Number of years.
#' @param nsex A vector containing the number of sexes for each species.
#' @param nages A vector containing the number of ages for each species.
#' @param nyrs_hind Number of hindcast years.
#' @param proj_mean_rec A flag indicating whether to use mean recruitment.
#' @param srr_pred_fun A function identifier for the stock-recruitment relationship.
#' @param R A matrix of recruitment values.
#' @param avg_R A vector of average recruitment values.
#' @param R0 A vector of initial recruitment values.
#' @param rec_dev A matrix of recruitment deviations.
#' @param beta_rec_pars A matrix of recruitment parameters.
#' @param env_index_srr A matrix of environmental indices.
#' @param M_at_age A 4D array of natural mortality rates at age.
#' @param Ftarget_age_spp A 4D array of fishing mortality rates at age.
#' @param pmature A matrix of maturity values.
#' @param spawn_month A vector of spawning months.
#' @param wt A 4D array of weights.
#' @param ssb_wt_index A vector of spawning biomass weight indices.
#' @param R_sexr A vector of sex ratios.
#' @param minage A vector of minimum ages.
#' @param MSSB0 A vector of multi-species spawning biomass.
#' @param MSB0 A vector of multi-species biomass.
#' @param msmMode A flag indicating whether to run in multi-species mode.
#' @param max_nsex
#' @param max_nages
#' @param N_at_age
#'
#' @return A list containing the calculated population dynamics:
#'   - NByage0: Abundance at age for the initial year.
#'   - NByageF: Abundance at age for the forecast year.
#'   - DynamicNByage0: Dynamic abundance at age for the initial year.
#'   - DynamicNByageF: Dynamic abundance at age for the forecast year.
#'   - B0: Biomass at the initial year.
#'   - SB0: Spawning biomass at the initial year.
#'   - SBF: Spawning biomass at the forecast year.
#'   - DynamicB0: Dynamic biomass at the initial year.
#'   - DynamicSB0: Dynamic spawning biomass at the initial year.
#'   - DynamicSBF: Dynamic spawning biomass at the forecast year.
#'
calculate_depletion_reference_points <- function(nspp, nyrs, nsex, nages, max_nsex, max_nages, nyrs_hind, proj_mean_rec,
                                                 srr_pred_fun, R, avg_R, R0, rec_dev,
                                                 beta_rec_pars, env_index_srr, N_at_age, M_at_age,
                                                 Ftarget_age_spp, pmature, spawn_month,
                                                 wt, ssb_wt_index, R_sexr, minage,
                                                 MSSB0, MSB0, msmMode) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize arrays to zero with appropriate dimensions
  NByage0 <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  NByageF <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  DynamicNByageF <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  DynamicNByage0 <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  B0 <- array(0, dim = c(nspp, nyrs))
  SB0 <- array(0, dim = c(nspp, nyrs))
  SBF <- array(0, dim = c(nspp, nyrs))
  DynamicB0 <- array(0, dim = c(nspp, nyrs))
  DynamicSB0 <- array(0, dim = c(nspp, nyrs))
  DynamicSBF <- array(0, dim = c(nspp, nyrs))

  # Loop through species and years
  for(sp in 1:nspp) {
    for(yr in 1:nyrs) {

      # Year 1 initialization
      if(yr == 1) {
        for(sex in 1:nsex[sp]) {
          for(age in 1:nages[sp]) {
            NByage0[sp, sex, age, yr] <- NByageF[sp, sex, age, yr] <-
              DynamicNByageF[sp, sex, age, yr] <- DynamicNByage0[sp, sex, age, yr] <-
              N_at_age[sp, sex, age, 1]
          }
        }
      }

      # Recruitment Year > 1
      if(yr > 1) {
        # Option 1a: Use mean rec
        if(proj_mean_rec == 1 && srr_pred_fun < 2) {
          # Equilibrium RPs
          NByage0[sp, 1, 1, yr] <- avg_R[sp]
          NByageF[sp, 1, 1, yr] <- avg_R[sp]

          # Dynamic RPs
          if(yr <= nyrs_hind) {
            DynamicNByageF[sp, 1, 1, yr] <- DynamicNByage0[sp, 1, 1, yr] <- R[sp, yr]
          }
          if(yr > nyrs_hind){
            DynamicNByageF[sp, 1, 1, yr] <- DynamicNByage0[sp, 1, 1, yr] <-
              exp(log(avg_R[sp]) + rec_dev[sp, yr])
          }
        }

        # Option 1b: Use mean rec and env
        if(proj_mean_rec == 1 && srr_pred_fun > 1) {
          beta_rec_tmp <- beta_rec_pars[sp, ]
          env_rec_tmp <- env_index_srr[yr, ]
          srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
          NByage0[sp, 1, 1, yr] <- avg_R[sp] * exp(srr_mult)

          if(yr <= nyrs_hind) {
            DynamicNByageF[sp, 1, 1, yr] <- DynamicNByage0[sp, 1, 1, yr] <- R[sp, yr]
          }
          if(yr > nyrs_hind){
            DynamicNByageF[sp, 1, 1, yr] <- DynamicNByage0[sp, 1, 1, yr] <-
              exp(log(avg_R[sp]) + rec_dev[sp, yr] + srr_mult)
          }
        }

        # Option 2: Use SRR
        if(proj_mean_rec == 0) {
          # Switch for different recruitment functions
          switch(srr_pred_fun,
                 # Case 1: Random about mean
                 {
                   NByage0[sp, 1, 1, yr] <- R0[sp]
                   NByageF[sp, 1, 1, yr] <- R0[sp]
                   DynamicNByage0[sp, 1, 1, yr] <- R0[sp] * exp(rec_dev[sp, yr])
                   DynamicNByageF[sp, 1, 1, yr] <- R0[sp] * exp(rec_dev[sp, yr])
                 },
                 # Case 2: Random about mean with environmental covariates
                 {
                   NByage0[sp, 1, 1, yr] <- R0[sp]
                   NByageF[sp, 1, 1, yr] <- R0[sp]
                   beta_rec_tmp <- beta_rec_pars[sp, ]
                   env_rec_tmp <- env_index_srr[yr, ]
                   srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
                   DynamicNByage0[sp, 1, 1, yr] <- R0[sp] * exp(rec_dev[sp, yr] + srr_mult)
                   DynamicNByageF[sp, 1, 1, yr] <- R0[sp] * exp(rec_dev[sp, yr] + srr_mult)
                 },
                 # Additional cases for other recruitment functions...
                 stop("Invalid srr_fun")
          )
        }

        # Account for sex ratio
        # - Females
        NByage0[sp, 1, 1, yr] <- NByage0[sp, 1, 1, yr] * R_sexr[sp]
        NByageF[sp, 1, 1, yr] <- NByageF[sp, 1, 1, yr] * R_sexr[sp]

        DynamicNByage0[sp, 1, 1, yr] <- DynamicNByage0[sp, 1, 1, yr] * R_sexr[sp]
        DynamicNByageF[sp, 1, 1, yr] <- DynamicNByageF[sp, 1, 1, yr] * R_sexr[sp]

        # N-at-age calculations for year > 1
        for(sex in 1:nsex[sp]) {

          # - Males
          if(sex == 2){
            NByage0[sp, 2, 1, yr] <- NByage0[sp, 1, 1, yr] / R_sexr[sp] * (1 - R_sexr[sp])
            NByageF[sp, 2, 1, yr] <- NByageF[sp, 1, 1, yr] / R_sexr[sp] * (1 - R_sexr[sp])

            DynamicNByage0[sp, 2, 1, yr] <- DynamicNByage0[sp, 1, 1, yr] / R_sexr[sp] * (1 - R_sexr[sp])
            DynamicNByageF[sp, 2, 1, yr] <- DynamicNByageF[sp, 1, 1, yr] / R_sexr[sp] * (1 - R_sexr[sp])
          }

          for(age in 2:(nages[sp]-1)) {
            NByage0[sp, sex, age, yr] <- NByage0[sp, sex, age-1, yr-1] *
              exp(-M_at_age[sp, sex, age-1, yr-1])

            NByageF[sp, sex, age, yr] <- NByageF[sp, sex, age-1, yr-1] *
              exp(-M_at_age[sp, sex, age-1, yr-1] - Ftarget_age_spp[sp, sex, age-1, yr-1])

            DynamicNByage0[sp, sex, age, yr] <- DynamicNByage0[sp, sex, age-1, yr-1] *
              exp(-M_at_age[sp, sex, age-1, yr-1])

            DynamicNByageF[sp, sex, age, yr] <- DynamicNByageF[sp, sex, age-1, yr-1] *
              exp(-M_at_age[sp, sex, age-1, yr-1] - Ftarget_age_spp[sp, sex, age-1, yr-1])
          }

          # Plus group calculations
          NByage0[sp, sex, nages[sp], yr] <-
            NByage0[sp, sex, nages[sp]-1, yr-1] * exp(-M_at_age[sp, sex, nages[sp]-1, yr-1]) +
            NByage0[sp, sex, nages[sp], yr-1] * exp(-M_at_age[sp, sex, nages[sp], yr-1])

          NByageF[sp, sex, nages[sp], yr] <-
            NByageF[sp, sex, nages[sp]-1, yr-1] *
            exp(-M_at_age[sp, sex, nages[sp]-1, yr-1] - Ftarget_age_spp[sp, sex, nages[sp]-1, yr-1]) +
            NByageF[sp, sex, nages[sp], yr-1] *
            exp(-M_at_age[sp, sex, nages[sp], yr-1] - Ftarget_age_spp[sp, sex, nages[sp], yr-1])

          DynamicNByage0[sp, sex, nages[sp], yr] <-
            DynamicNByage0[sp, sex, nages[sp]-1, yr-1] * exp(-M_at_age[sp, sex, nages[sp]-1, yr-1]) +
            DynamicNByage0[sp, sex, nages[sp], yr-1] * exp(-M_at_age[sp, sex, nages[sp], yr-1])

          DynamicNByageF[sp, sex, nages[sp], yr] <-
            DynamicNByageF[sp, sex, nages[sp]-1, yr-1] *
            exp(-M_at_age[sp, sex, nages[sp]-1, yr-1] - Ftarget_age_spp[sp, sex, nages[sp]-1, yr-1]) +
            DynamicNByageF[sp, sex, nages[sp], yr-1] *
            exp(-M_at_age[sp, sex, nages[sp], yr-1] - Ftarget_age_spp[sp, sex, nages[sp], yr-1])
        }
      }

      # Calculate Dynamic SB0 and SB at F target
      if(yr < nyrs_hind){
        yr_ind <- yr
      }
      if(yr >= nyrs_hind){
        yr_ind <- nyrs_hind
      }

      for(age in 1:nages[sp]) {
        SB0[sp, yr] <- SB0[sp, yr] +
          NByage0[sp, 1, age, yr] * wt[ssb_wt_index[sp], 1, age, nyrs_hind] *
          pmature[sp, age] * exp(-M_at_age[sp, 1, age, yr] * spawn_month[sp]/12)

        SBF[sp, yr] <- SBF[sp, yr] +
          NByageF[sp, 1, age, yr] * wt[ssb_wt_index[sp], 1, age, nyrs_hind] *
          pmature[sp, age] * exp(-(M_at_age[sp, 1, age, yr] +
                                     Ftarget_age_spp[sp, 1, age, yr]) * spawn_month[sp]/12)

        DynamicSB0[sp, yr] <- DynamicSB0[sp, yr] +
          DynamicNByage0[sp, 1, age, yr] * wt[ssb_wt_index[sp], 1, age, yr_ind] *
          pmature[sp, age] * exp(-M_at_age[sp, 1, age, yr] * spawn_month[sp]/12)

        DynamicSBF[sp, yr] <- DynamicSBF[sp, yr] +
          DynamicNByageF[sp, 1, age, yr] * wt[ssb_wt_index[sp], 1, age, yr_ind] *
          pmature[sp, age] * exp(-(M_at_age[sp, 1, age, yr] +
                                     Ftarget_age_spp[sp, 1, age, yr]) * spawn_month[sp]/12)

        for(sex in 1:nsex[sp]) {
          B0[sp, yr] <- B0[sp, yr] +
            NByage0[sp, sex, age, yr] * wt[ssb_wt_index[sp], sex, age, nyrs_hind]

          DynamicB0[sp, yr] <- DynamicB0[sp, yr] +
            DynamicNByage0[sp, sex, age, yr] * wt[ssb_wt_index[sp], sex, age, yr_ind]
        }
      }

      # Input SB0 for multi-species mode
      if(msmMode > 0) {
        SB0[sp, yr] <- MSSB0[sp]
        B0[sp, yr] <- MSB0[sp]
      }
    }
  }

  # Return results as a list
  return(list(NByage0 = NByage0, NByageF = NByageF, DynamicNByage0 = DynamicNByage0,
              DynamicNByageF = DynamicNByageF, B0 = B0, SB0 = SB0, SBF = SBF,
              DynamicB0 = DynamicB0, DynamicSB0 = DynamicSB0, DynamicSBF = DynamicSBF))
}



#' Calculate Temperature-Dependent Consumption and Ration
#'
#' @description
#' Calculates temperature-dependent consumption and ration for multiple species across years,
#' implementing various temperature dependence functions and consumption models.
#'
#' @param nspp Integer, number of species
#' @param nyrs Integer, number of years
#' @param nyrs_hind Integer, number of hindcast years
#' @param Ceq Vector, consumption equation type for each species
#' @param Qc Vector, Q coefficient for each species
#' @param Tcm Vector, maximum temperature for each species
#' @param Tco Vector, optimum temperature for each species
#' @param Tcl Vector, limiting temperature for each species
#' @param CK1 Vector, first temperature coefficient for each species
#' @param CK4 Vector, fourth temperature coefficient for each species
#' @param CA Vector, intercept of allometric mass function
#' @param CB Vector, slope of allometric mass function
#' @param fday Vector, number of foraging days per year
#' @param env_index Matrix, environmental index by year and index type
#' @param Cindex Vector, index for environmental variables
#' @param wt Array, weight at age by population, sex, and year
#' @param pop_wt_index Vector, population weight indices
#' @param Pvalue Vector, proportion of maximum consumption
#' @param Pyrs Array, proportion of year specific values
#' @param nsex Vector, number of sexes per species
#' @param nages Vector, number of ages per species
#' @param max_nsex
#' @param max_nages
#'
#' @return List containing:
#' \itemize{
#'   \item fT: Matrix of temperature functions by species and year
#'   \item consumption_at_age: Array of consumption at age
#'   \item ration: Array of annual ration in kg/year
#' }
#'
#' @details
#' Implements multiple temperature dependence functions:
#' 1. Exponential (Stewart et al. 1983)
#' 2. Warm-water species (Kitchell et al. 1977)
#' 3. Cool and cold-water species (Thornton and Lessem 1979)
#' 4. No temperature dependence
#'
calculate_ration <- function(nspp, nyrs, max_nsex, max_nages, nyrs_hind, Ceq, Qc, Tcm, Tco, Tcl,
                             CK1, CK4, CA, CB, fday, env_index, Cindex,
                             wt, pop_wt_index, Pvalue, Pyrs, nsex, nages) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize outputs
  fT <- matrix(0, nrow = nspp, ncol = nyrs)
  consumption_at_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  ration <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))

  # Calculate temperature function of consumption
  for(sp in 1:nspp) {
    for(yr in 1:nyrs) {

      switch(Ceq[sp],
             # Case 1: Exponential function (Stewart et al. 1983)
             {
               fT[sp, yr] <- exp(Qc[sp] * env_index[yr, Cindex[sp]])
             },

             # Case 2: Warm-water species (Kitchell et al. 1977)
             {
               Yc <- log(Qc[sp]) * (Tcm[sp] - Tco[sp] + 2.0)
               Zc <- log(Qc[sp]) * (Tcm[sp] - Tco[sp])
               Vc <- (Tcm[sp] - env_index[yr, Cindex[sp]]) / (Tcm[sp] - Tco[sp])
               Xc <- Zc^2 * (1 + sqrt(1 + 40/Yc))^2 / 400
               fT[sp, yr] <- Vc^Xc * exp(Xc * (1 - Vc))
             },

             # Case 3: Cool and cold-water species (Thornton and Lessem 1979)
             {
               G2 <- (1/(Tcl[sp] - Tcm[sp])) * log((0.98 * (1 - CK4[sp])) / (CK4[sp] * 0.02))
               L2 <- exp(G2 * (Tcl[sp] - env_index[yr, Cindex[sp]]))
               Kb <- (CK4[sp] * L2) / (1 + CK4[sp] * (L2 - 1))
               G1 <- (1/(Tco[sp] - Qc[sp])) * log((0.98 * (1 - CK1[sp])) / (CK1[sp] * 0.02))
               L1 <- exp(G1 * (env_index[yr, Cindex[sp]] - Qc[sp]))
               Ka <- (CK1[sp] * L1) / (1 + CK1[sp] * (L1 - 1))
               fT[sp, yr] <- Ka * Kb
             },

             # Case 4: No temperature dependence
             {
               fT[sp, yr] <- 1.0
             }
      )
    }
  }

  # Calculate historic ration
  for(sp in 1:nspp) {
    for(sex in 1:nsex[sp]) {
      for(age in 1:nages[sp]) {
        for(yr in 1:nyrs) {
          if(yr <= nyrs_hind) {
            # Hindcast calculation
            consumption_at_age[sp, sex, age, yr] <- CA[sp] *
              (wt[pop_wt_index[sp], sex, age, yr] * 1000)^(1 + CB[sp]) *
              fT[sp, yr] * fday[sp]

            consumption_at_age[sp, sex, age, yr] <- consumption_at_age[sp, sex, age, yr] *
              Pvalue[sp] * Pyrs[sp, sex, age, yr]

          }
          if(yr > nyrs_hind){
            # Projection calculation
            consumption_at_age[sp, sex, age, yr] <- CA[sp] *
              (wt[pop_wt_index[sp], sex, age, nyrs_hind] * 1000)^(1 + CB[sp]) *
              fT[sp, yr] * fday[sp]

            consumption_at_age[sp, sex, age, yr] <- consumption_at_age[sp, sex, age, yr] *
              Pvalue[sp] * Pyrs[sp, sex, age, nyrs_hind]
          }

          # Calculate annual ration in kg/year
          ration[sp, sex, age, yr] <- consumption_at_age[sp, sex, age, yr] / 1000
        }
      }
    }
  }

  # Report
  RTMB::REPORT(consumption_at_age)
  RTMB::REPORT(ration)

  return(ration)
}


#' Reorganize Stomach Content Proportion Observations
#'
#' @description
#' This function reorganizes stomach content proportion observations into a diet proportion matrix,
#' accounting for predator-prey sex combinations and averaging across years when specified.
#'
#' @param stom_prop_obs Matrix of stomach proportion observations
#' @param stom_prop_ctl Matrix of control data for stomach proportions
#' @param minage Vector of minimum ages by species
#' @param nspp Number of species
#' @param nyrs Number of years
#' @param nyrs_hind Number of hindcast years
#' @param styr Start year
#' @param nsex Vector, number of sexes per species
#' @param nages Vector, number of ages per species
#' @param max_nsex
#' @param max_nages
#'
#' @return A list containing:
#' \itemize{
#'   \item diet_prop: Array of diet proportions
#'   \item r_sexes: Matrix of predator sex indices
#'   \item k_sexes: Matrix of prey sex indices
#' }
#'
#' @details
#' The function processes stomach content observations and organizes them into a comprehensive
#' diet proportion matrix. It handles both single-sex and two-sex models, and can process
#' both year-specific and year-averaged data.
reorganize_stomach_content <- function(stom_prop_obs, stom_prop_ctl, minage, nspp,
                                       nyrs, nsex, nages, max_nsex, max_nages, nyrs_hind, styr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize outputs
  n_obs <- nrow(stom_prop_obs)
  r_sexes <- matrix(as.integer(1), nrow=n_obs, ncol=2)
  k_sexes <- matrix(as.integer(1), nrow=n_obs, ncol=2)

  # Correct dimensions for diet_prop
  diet_prop <- array(0, dim=c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))

  # Precompute indices for predator and prey
  rsp <- as.integer(stom_prop_ctl[, 1])      # Index of predator
  ksp <- as.integer(stom_prop_ctl[, 2])      # Index of prey
  r_sex <- as.integer(stom_prop_ctl[, 3])    # Index of predator sex
  k_sex <- as.integer(stom_prop_ctl[, 4])    # Index of prey sex
  r_age <- as.integer(stom_prop_ctl[, 5] - minage[rsp] + 1)  # Index of predator age
  k_age <- as.integer(stom_prop_ctl[, 6] - minage[ksp] + 1)  # Index of prey age
  flt_yr <- as.integer(stom_prop_ctl[, 7])   # Index of year

  # Set predator and prey sex indices for joint sex data
  # - Sex = 0 by nsex = 2
  if(max_nsex > 2){
    r_sexes[, ] <- cbind(ifelse(r_sex > 0, r_sex, as.integer(1)), ifelse(r_sex > 0, r_sex, as.integer(2)))
    k_sexes[, ] <- cbind(ifelse(k_sex > 0, k_sex, as.integer(1)), ifelse(k_sex > 0, k_sex, as.integer(2)))
  }

  # Process for each observation
  for(stom_ind in 1:n_obs) {
    # Extract current indices
    current_rsp <- as.integer(rsp[stom_ind])
    current_ksp <- as.integer(ksp[stom_ind])
    current_r_age <- as.integer(r_age[stom_ind])
    current_k_age <- as.integer(k_age[stom_ind])
    current_diet_yr <- as.integer(flt_yr[stom_ind])

    if(current_diet_yr > 0) { # Annual diet data
      yr <- current_diet_yr - styr + 1
      if(yr < nyrs_hind) {
        diet_prop[current_rsp, as.integer(r_sexes[stom_ind, ]), current_r_age,
                  current_ksp, as.integer(k_sexes[stom_ind, ]), current_k_age, yr] <- as.numeric(stom_prop_obs[stom_ind, 2])
      }
    }

    if(current_diet_yr == 0) { # Average diet data
      diet_prop[current_rsp, as.integer(r_sexes[stom_ind, ]), current_r_age,
                current_ksp, as.integer(k_sexes[stom_ind, ]), current_k_age, ] <- stom_prop_obs[stom_ind, 2]
    }
  }

  # Report
  # RTMB::REPORT(r_sexes)
  # RTMB::REPORT(k_sexes)

  # Return
  return(diet_prop = diet_prop)
}


#' Calculate Other Food Stomach Content
#'
#' @description
#' Calculates the proportion of other food in stomach content for each predator species,
#' sex, age and year combination.
#'
#' @param nyrs Number of years
#' @param nspp Number of species
#' @param nsex Vector containing number of sexes for each species
#' @param nages Vector containing number of ages for each species
#' @param diet_prop Array of diet proportions [pred_sp*sex, prey_sp*sex, pred_age, prey_age, year]
#' @param other_food Vector of other food values for each species
#' @param max_nsex
#' @param max_nages
#'
#' @return Array of other food diet proportions [pred_sp, pred_sex, pred_age, year]
#'
#' @details
#' The function initializes other food diet proportion as 1 and subtracts the proportions
#' of all other prey species. If other food value exists for a predator species (>0),
#' the proportion is penalized by dividing by the other food value.
#'
calculate_other_food_diet_prop <- function(nyrs, nspp, nsex, nages, max_nsex, max_nages, diet_prop, other_food) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize array for other food diet proportions
  other_food_diet_prop <- array(1, dim = c(nspp, max_nsex, max_nages, nyrs))

  # Loop through years
  for(yr in 1:nyrs) {
    for(rsp in 1:nspp) { # Loop through predator species
      for(r_sex in 1:nsex[rsp]) { # Loop through predator sex
        for(r_age in 1:nages[rsp]) { # Loop through predator age

          # Subtract diet proportions for all prey items
          other_food_diet_prop[rsp, r_sex, r_age, yr] = 1 -
            sum(diet_prop[rsp, r_sex, r_age, , , , yr])

          # Apply other food penalties
          if(other_food[rsp] > 0) {
            other_food_diet_prop[rsp, r_sex, r_age, yr] <-
              other_food_diet_prop[rsp, r_sex, r_age, yr] / other_food[rsp]
          }
          if(other_food[rsp] == 0) {
            other_food_diet_prop[rsp, r_sex, r_age, yr] <- 0
          }
        }
      }
    }
  }

  return(other_food_diet_prop)
}


#' Calculate MSVPA based Suitability for Predator-Prey Interactions
#'
#' @description
#' This function calculates suitability coefficients for predator-prey interactions
#' based on stomach content data and population parameters.
#'
#' @param diet_prop Array of diet proportions [pred_sp x prey_sp x pred_age x prey_age x year]
#' @param avgN_at_age Array of average numbers at age [species x sex x age x year]
#' @param wt Array of weights [index x sex x age x year]
#' @param pop_wt_index Vector of population weight indices
#' @param other_food_diet_prop Array of other food proportions [pred_sp x pred_sex x pred_age x year]
#' @param nspp Integer number of species
#' @param nsex Vector of number of sexes per species
#' @param nages Vector of number of ages per species
#' @param nyrs Integer number of years
#' @param nyrs_hind Integer number of hindcast years
#' @param suit_styr Integer start year for suitability calculation
#' @param suit_endyr Integer end year for suitability calculation
#' @param nyrs_suit Integer number of years for suitability calculation
#' @param msmMode Integer indicating MSVPA mode (1 or 2)
#'
#' @return List containing:
#' \itemize{
#'   \item suit_main: Main suitability coefficients
#'   \item suit_other: Other food suitability coefficients
#' }
#'
calculate_MSVPA_suitability <- function(diet_prop, avgN_at_age, wt, pop_wt_index,
                                        other_food_diet_prop, nspp, nsex, max_nsex, nages, max_nages, nyrs,
                                        nyrs_hind, suit_styr, suit_endyr, nyrs_suit,
                                        msmMode) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize arrays
  stom_div_bio <- array(0, dim = c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))
  suma_suit <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))
  suit_main <- array(0, dim = c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))
  suit_other <- array(1, dim = c(nspp, max_nsex, max_nages, nyrs))

  # Calculate stomach proportion over biomass
  for(yr in 1:nyrs) {
    for(ksp in 1:nspp) {
      for(rsp in 1:nspp) {
        for(k_sex in 1:nsex[ksp]) {
          for(r_sex in 1:nsex[rsp]) {
            for(r_age in 1:nages[rsp]) {
              for(k_age in 1:nages[ksp]) {

                if(yr < nyrs_hind){
                  yr_ind <- yr
                }
                if(yr >= nyrs_hind){
                  yr_ind <- nyrs_hind
                }

                #if(avgN_at_age[ksp, k_sex, k_age, yr] > 0) {
                stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] <-
                  diet_prop[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] /
                  avgN_at_age[ksp, k_sex, k_age, yr]

                if(msmMode == 2) {
                  stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] <-
                    stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] /
                    avgN_at_age[ksp, k_sex, k_age, yr]
                }
                # }

                # Handle non-finite values
                # if(!is.finite(stom_div_bio[as.integer(rsp + (nspp * (r_sex-1))),
                #                            as.integer(ksp + (nspp * (k_sex-1))),
                #                            r_age, k_age, yr])) {
                #   stom_div_bio[as.integer(rsp + (nspp * (r_sex-1))),
                #                as.integer(ksp + (nspp * (k_sex-1))),
                #                r_age, k_age, yr] <- 0
                # }

                if(wt[pop_wt_index[ksp], k_sex, k_age, yr_ind] != 0) {
                  stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] <-
                    stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] /
                    wt[pop_wt_index[ksp], k_sex, k_age, yr_ind]

                  suma_suit[rsp, r_sex, r_age, yr] <-
                    suma_suit[rsp, r_sex, r_age, yr] +
                    stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr]
                }
              }
            }
          }
        }
      }
    }
  }

  # Calculate suitability
  for(rsp in 1:nspp) {
    for(r_sex in 1:nsex[rsp]) {
      for(r_age in 1:nages[rsp]) {
        for(ksp in 1:nspp) {
          for(k_sex in 1:nsex[ksp]) {
            for(k_age in 1:nages[ksp]) {

              # Calculate average suitability
              for(yr in suit_styr:suit_endyr) {
                # if(suma_suit[rsp, r_sex, r_age, yr] +
                #    other_food_diet_prop[rsp, r_sex, r_age, yr] > 0) {
                suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, 1] <-
                  suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, 1] +
                  stom_div_bio[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] /
                  (suma_suit[rsp, r_sex, r_age, yr] +
                     other_food_diet_prop[rsp, r_sex, r_age, yr])
                # }
              }

              # Average across years
              suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, 1] <-
                suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, 1] / nyrs_suit

              # # Handle non-finite values
              # if(!is.finite(suit_main[as.integer(rsp + (nspp * (r_sex-1))),
              #                         as.integer(ksp + (nspp * (k_sex-1))),
              #                         r_age, k_age, 1])) {
              #   suit_main[as.integer(rsp + (nspp * (r_sex-1))),
              #             as.integer(ksp + (nspp * (k_sex-1))),
              #             r_age, k_age, 1] <- 0
              # }

              # Fill in remaining years
              for(yr in 2:nyrs) {
                suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] <-
                  suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, 1]
              }

              # Calculate other suitability
              for(yr in 1:nyrs) {
                suit_other[rsp, r_sex, r_age, yr] <-
                  suit_other[rsp, r_sex, r_age, yr] -
                  suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, yr]
              }
            }
          }
        }
      }
    }
  }

  return(list(suit_main = suit_main,
              suit_other = suit_other))
}


#' Calculate GAMMA Size-Based Suitability Matrix
#'
#' @description
#' Calculates predator-prey suitability based on either length or weight ratios using
#' gamma distribution.
#'
#' @param nspp Integer, number of species
#' @param nages Vector of integers, number of ages for each species
#' @param nsex Vector of integers, number of sexes for each species
#' @param nyrs Integer, number of years
#' @param nyrs_hind Integer, number of hindcast years
#' @param laa Array [species, sex, age, year] of length-at-age values
#' @param wt Array [species, sex, age, year] of weight values
#' @param pop_wt_index Vector of indices for population weights
#' @param vulnerability Matrix [pred, prey] of vulnerability coefficients
#' @param vulnerability_other Vector of other food vulnerability
#' @param gam_a Vector of gamma distribution shape parameters
#' @param gam_b Vector of gamma distribution rate parameters
#' @param suitMode Integer, 1 for length-based or 2 for weight-based suitability
#' @param max_nages
#' @param max_nsex
#'
#' @return List containing:
#' \itemize{
#'   \item suit_main: Array [pred_sp, pred_sex, pred_age, prey_sp, prey_sex, prey_age, year] of suitability values
#'   \item suit_other: Array [pred_sp, pred_sex, pred_age, year] of other food suitability
#' }
#'
#' @details
#' The function calculates size-based suitability coefficients using either length or weight
#' ratios between predators and prey. Suitability is calculated using a gamma distribution
#' and scaled to be between 0 and 1.
#'
calculate_gamma_suitability <- function(nspp, nages, nsex, nyrs, nyrs_hind,
                                        max_nages, max_nsex,
                                        laa, wt, pop_wt_index, vulnerability,
                                        vulnerability_other, gam_a, gam_b, suitMode) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize arrays
  suit_main <- array(0, dim = c(nspp, max_nsex, max_nages,
                                nspp, max_nsex, max_nages, nyrs))
  suit_other <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))

  # Loop through all dimensions
  for(rsp in 1:nspp) {
    for(r_age in 1:nages[rsp]) {
      for(r_sex in 1:nsex[rsp]) {
        # Set other food suitability
        suit_other[rsp, r_sex, r_age, ] <- vulnerability_other[rsp]

        for(ksp in 1:nspp) {
          for(k_sex in 1:nsex[ksp]) {
            for(k_age in 1:nages[ksp]) {
              # Calculate size ratio based on mode
              if(suitMode == 1) {
                log_size_ratio <- laa[rsp, r_sex, r_age, ] / laa[ksp, k_sex, k_age, ]
              }
              if(suitMode == 2) {
                log_size_ratio <- wt[pop_wt_index[rsp], r_sex, r_age, ] / wt[pop_wt_index[ksp], k_sex, k_age, ]
              }

              log_size_ratio <- log(log_size_ratio)

              # Calculate suitability for each year
              for(yr in 1:nyrs) {
                # Determine year index
                if(yr < nyrs_hind){
                  yr_ind <- yr
                }
                if(yr >= nyrs_hind){
                  yr_ind <- nyrs_hind
                }

                if(log_size_ratio[yr] > 0) {

                  # https:#github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp
                  suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] <-
                    vulnerability[rsp, ksp] *
                    dgamma(log_size_ratio[yr], gam_a[rsp], gam_b[rsp]) /
                    dgamma((gam_a[rsp]-1) * gam_b[rsp], gam_a[rsp], gam_b[rsp])
                }
              }
            }
          }
        }
      }
    }
  }

  return(list(suit_main = suit_main,
              suit_other = suit_other))
}


#' Calculate Lognormal Suitability Matrix
#'
#' @description
#' Calculates the suitability matrix for predator-prey interactions using either length-based
#' or weight-based lognormal suitability functions.
#'
#' @param nspp Integer, number of species
#' @param nsex Numeric vector, number of sexes for each species
#' @param nages Numeric vector, number of ages for each species
#' @param nyrs Integer, number of years
#' @param nyrs_hind Integer, number of hindcast years
#' @param laa Array, length-at-age matrix [species, sex, age, year]
#' @param wt Array, weight matrix [species, sex, age, year]
#' @param vulnerability Matrix, predator-prey vulnerability matrix
#' @param vulnerability_other Vector, vulnerability to other food
#' @param gam_a Vector, first parameter of normal distribution for each predator
#' @param gam_b Vector, second parameter of normal distribution for each predator
#' @param pop_wt_index Vector, indices for population weights
#' @param suitMode Integer, 3 for length-based or 4 for weight-based suitability
#' @param max_nages
#' @param max_nsex
#'
#' @return List containing:
#' \itemize{
#'   \item suit_main: Array of predator-prey suitability coefficients
#'   \item suit_other: Array of other food suitability coefficients
#' }
#'
#' @details
#' The function calculates suitability coefficients based on either length or weight ratios
#' between predators and prey. It uses a lognormal function scaled by vulnerability.
#'
calculate_lognormal_suitability <- function(nspp, nages, nsex, nyrs, nyrs_hind,
                                            max_nages, max_nsex,
                                            laa, wt, pop_wt_index, vulnerability,
                                            vulnerability_other, gam_a, gam_b, suitMode) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize arrays
  suit_main <- array(0, dim = c(nspp, max_nsex, max_nages,
                                nspp, max_nsex, max_nages, nyrs))
  suit_other <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))

  # Loop through all dimensions
  for(rsp in 1:nspp) {
    for(r_age in 1:nages[rsp]) {
      for(r_sex in 1:nsex[rsp]) {
        # Set other food suitability
        suit_other[rsp, r_sex, r_age, ] <- vulnerability_other[rsp]

        for(ksp in 1:nspp) {
          for(k_sex in 1:nsex[ksp]) {
            for(k_age in 1:nages[ksp]) {
              # Calculate size ratio based on mode
              if(suitMode == 3) {
                log_size_ratio <- laa[rsp, r_sex, r_age, ] / laa[ksp, k_sex, k_age, ]
              }
              if(suitMode == 4) {
                log_size_ratio <- wt[pop_wt_index[rsp], r_sex, r_age, ] / wt[pop_wt_index[ksp], k_sex, k_age, ]
              }

              log_size_ratio <- log(log_size_ratio)

              # Calculate suitability for each year
              for(yr in 1:nyrs) {
                # Determine year index
                if(yr < nyrs_hind){
                  yr_ind <- yr
                }
                if(yr >= nyrs_hind){
                  yr_ind <- nyrs_hind
                }

                if(log_size_ratio[yr] > 0) {
                  # Calculate suitability

                  suit_main[rsp, r_sex, r_age, ksp, k_sex, k_age, yr] <-
                    vulnerability[rsp, ksp] *
                    dnorm(log_size_ratio[yr_ind], gam_a[rsp], gam_b[rsp]) /
                    dnorm(gam_a[rsp], gam_a[rsp], gam_b[rsp]) # Divide by mode to scale to 1
                }
              }
            }
          }
        }
      }
    }
  }

  return(list(suit_main = suit_main,
              suit_other = suit_other))
}



#' Calculate Predation Mortality and Diet Components
#'
#' @description
#' Calculates available food, predation mortality (M2), diet proportions, and biomass eaten
#' in a multi-species model framework. Implements both Type 2 and Type 3 MSVPA approaches.
#'
#' @param nspp Integer, number of species
#' @param nsex Numeric vector, number of sexes per species
#' @param nages Numeric vector, number of ages per species
#' @param nyrs Integer, total number of years
#' @param nyrs_hind Integer, number of hindcast years
#' @param avgN_at_age Array[spp,sex,age,year], average numbers at age
#' @param suit_main Array, suitability matrix
#' @param suit_other Array, other food suitability
#' @param other_food Numeric vector, other food availability by species
#' @param wt Array, weight at age matrix
#' @param pop_wt_index Numeric vector, population weight indices
#' @param ration Array, consumption ration
#' @param msmMode Integer, MSVPA type (1=Type 2, 2=Type 3)
#' @param max_nages
#' @param max_nsex
#'
#' @return List containing:
#' \itemize{
#'   \item avail_food: Available food matrix
#'   \item M2_at_age: Predation mortality at age
#'   \item B_eaten_as_prey: Biomass eaten as prey
#'   \item B_eaten: Biomass eaten matrix
#'   \item diet_prop_hat: Predicted diet proportions
#'   \item avg_diet_prop_hat: Average predicted diet proportions
#' }
#'
calculate_predation <- function(nspp, nsex, nages, nyrs, nyrs_hind,
                                max_nages, max_nsex,
                                avgN_at_age, suit_main, suit_other, other_food,
                                wt, pop_wt_index, ration, msmMode) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize arrays
  avail_food <- array(0, dim=c(nspp, max_nsex, max_nages, nyrs))
  M2_at_age <- array(0, dim=c(nspp, max_nsex, max_nages, nyrs))
  M2_prop <- array(0, dim=c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))
  B_eaten_as_prey <- array(0, dim=c(nspp, max_nsex, max_nages, nyrs))
  B_eaten <- array(0, dim=c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))
  diet_prop_hat <- array(0, dim=c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))
  avg_diet_prop_hat <- array(0, dim=c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages))

  # Calculate available food
  for(rsp in 1:nspp) {
    for(r_sex in 1:nsex[rsp]) {
      for(r_age in 1:nages[rsp]) {
        for(yr in 1:nyrs) {

          yr_ind <- yr
          if (yr > nyrs_hind) yr_ind = nyrs_hind

          # Sum across prey species
          avail_food[rsp, r_sex, r_age, yr] <- sum(suit_main[rsp, r_sex, r_age, , , , yr] *
                                                     avgN_at_age[, , , yr]^msmMode *
                                                     wt[pop_wt_index, , , yr_ind])

          # Add other food
          avail_food[rsp, r_sex, r_age, yr] <- avail_food[rsp, r_sex, r_age, yr] +
            other_food[rsp] * suit_other[rsp, r_sex, r_age, yr]
        }
      }
    }
  }

  # Calculate predation mortality and related quantities
  for(ksp in 1:nspp) {
    for(k_sex in 1:nsex[ksp]) {
      for(k_age in 1:nages[ksp]) {
        for(yr in 1:nyrs) {

          yr_ind <- yr
          if (yr > nyrs_hind) yr_ind = nyrs_hind

          if(msmMode == 1) { # Type 2 MSVPA

            # Calculate predation mortality proportion
            M2_prop[, , , ksp, k_sex, k_age, yr] <-
              (avgN_at_age[, , , yr] *
                 ration[, , , yr] *
                 suit_main[, , , ksp, k_sex, k_age, yr]) /
              avail_food[, , , yr]

            # Calculate predation mortality
            M2_at_age[ksp, k_sex, k_age, yr] <- sum(
              M2_prop[, , , ksp, k_sex, k_age, yr], na.rm = TRUE
            )

            # Calculate biomass eaten by predator-prey combination
            B_eaten[, , , ksp, k_sex, k_age, yr] <-
              avgN_at_age[ksp, k_sex, k_age, yr] *
              wt[pop_wt_index[ksp], k_sex, k_age, yr_ind] *
              M2_prop[, , , ksp, k_sex, k_age, yr]

            # Calculate biomass eaten as prey
            B_eaten_as_prey[ksp, k_sex, k_age, yr] <- sum(
              B_eaten[, , , ksp, k_sex, k_age, yr]
            )

            # Calculate diet proportion
            diet_prop_hat[, , , ksp, k_sex, k_age, yr] <-
              avgN_at_age[ksp, k_sex, k_age, yr] *
              wt[pop_wt_index[ksp], k_sex, k_age, yr_ind] *
              suit_main[, , , ksp, k_sex, k_age, yr] /
              avail_food[, , , yr]

          }
          if(msmMode == 2) { # Type 3 MSVPA
            # TODO
          }


          # Calculate average diet proportion
          if(yr <= nyrs_hind) {
            avg_diet_prop_hat[, , , ksp, k_sex, k_age] <-
              avg_diet_prop_hat[, , , ksp, k_sex, k_age] +
              diet_prop_hat[, , , ksp, k_sex, k_age, yr] / nyrs_hind
          }
        }
      }
    }
  }


  # Report
  RTMB::REPORT(B_eaten_as_prey)
  RTMB::REPORT(M2_at_age)

  return(list(
    avail_food = avail_food,
    M2_at_age = M2_at_age,
    B_eaten_as_prey = B_eaten_as_prey,
    B_eaten = B_eaten,
    diet_prop_hat = diet_prop_hat,
    avg_diet_prop_hat = avg_diet_prop_hat
  ))
}

#' Calculate Index of Abundance/Biomass
#'
#' @description
#' Calculates predicted indices of abundance or biomass for each survey/species combination
#' based on population numbers-at-age and various adjustment factors.
#'
#' @param index_ctl Matrix containing control information for indices
#' @param index_n Matrix containing month information for indices
#' @param N_at_age Array of numbers at age [species, sex, age, year]
#' @param Z_at_age Array of total mortality at age [species, sex, age, year]
#' @param sel Array of selectivity [fleet, sex, age, year]
#' @param wt Array of weights [fleet_wt_index, sex, age, year]
#' @param flt_units Vector indicating units for each fleet (1=weight, 2=numbers)
#' @param flt_wt_index Vector indicating which weight index to use for each fleet
#' @param nages Vector of number of ages for each species
#' @param nsex Vector of number of sexes for each species
#' @param nyrs_hind Number of hindcast years
#' @param styr Start year of model
#'
#' @return Vector of predicted index values
#'
#' @details
#' The function calculates predicted survey indices by combining numbers-at-age with
#' mortality, selectivity and optional weight information. Indices can be in either
#' biomass or numbers. The calculation accounts for the timing of the survey within
#' the year through the month parameter.
#'
calculate_abundance_index <- function(index_ctl, index_n, N_at_age, Z_at_age, sel,
                                      wt, flt_units, flt_wt_index, nages, nsex,
                                      nyrs_hind, styr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  n_indices <- nrow(index_ctl)
  index_hat <- numeric(n_indices)

  for(index_ind in 1:n_indices) {
    # Get index parameters
    index <- index_ctl[index_ind, 1]  # Survey index (1-based)
    sp <- index_ctl[index_ind, 2]      # Species index (1-based)
    flt_yr <- index_ctl[index_ind, 3]  # Year
    mo <- index_n[index_ind, 1]        # Month

    # Initialize prediction
    index_hat[index_ind] <- 0

    # Adjust year index
    if(flt_yr > 0) {
      flt_yr <- flt_yr - styr + 1
    }
    if(flt_yr < 0) {
      flt_yr <- -flt_yr - styr + 1
    }

    # Set year index for selectivity/weight
    yr_ind <- flt_yr
    if (flt_yr > nyrs_hind) yr_ind = nyrs_hind

    # Loop over ages and sexes
    for(age in 1:nages[sp]) {  # Adjusted to use nages[sp]
      for(sex in 1:nsex[sp]) {  # Adjusted to use nsex[sp]

        # Base calculation
        base_calc <- N_at_age[sp, sex, age, flt_yr] *
          exp(-(mo / 12) * Z_at_age[sp, sex, age, flt_yr]) *
          sel[index, sex, age, yr_ind]

        # Add to index based on units
        if(flt_units[index] == 1) {
          # Weight
          index_hat[index_ind] <- index_hat[index_ind] +
            base_calc * wt[flt_wt_index[index], sex, age, yr_ind]
        }

        if(flt_units[index] == 2) {
          # Numbers
          index_hat[index_ind] <- index_hat[index_ind] + base_calc
        }
      }
    }
  }

  RTMB::REPORT(index_hat)
  return(index_hat)
}


#' Calculate Analytical Survey Catchability (q) Following Ludwig and Martell 1994
#'
#' @description
#' This function calculates the analytical survey catchability coefficient (q) using the method
#' described in Ludwig and Martell (1994). It handles both time-varying and non-time-varying
#' standard deviations.
#'
#' @param index_ctl Matrix containing control information for indices
#' @param index_n Matrix containing temporal information for indices
#' @param index_obs Matrix containing observed index values
#' @param index_hat Vector of predicted index values
#' @param est_sigma_index Vector indicating if sigma is estimated (>0) or time-varying (0)
#' @param est_index_q Vector indicating if analytical q should be used (3)
#' @param n_flt Integer number of fleets
#' @param nyrs_hind Integer number of hindcast years
#' @param styr Integer start year
#'
#' @return A list containing:
#' \itemize{
#'   \item index_q_analytical: Vector of analytical q values
#'   \item index_q: Matrix of q values by fleet and year
#' }
#'
#' @details
#' The function implements the following steps:
#' 1. Calculates observation counts and sums log ratios of observed to predicted values
#' 2. Handles both estimated and time-varying standard deviations
#' 3. Takes averages and exponentiates to get final q values
#' 4. Applies analytical q values when specified
#'
#' @references
#' Ludwig, D. and Martell, S. (1994). A method for calculating catchability coefficients
#'
calculate_analytical_q <- function(index_ctl, index_n, index_obs, index_hat,
                                   est_sigma_index, est_index_q, n_flt, nyrs_hind, styr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize vectors
  index_n_obs <- numeric(n_flt)
  index_q_analytical <- numeric(n_flt)

  # Calculate sums for q
  for(index_ind in 1:nrow(index_ctl)) {
    index <- index_ctl[index_ind, 1]
    flt_yr <- index_ctl[index_ind, 3]

    if(flt_yr > 0) {
      flt_yr <- flt_yr - styr
      mo <- index_n[index_ind, 1]  # Assuming this is the correct column

      if(flt_yr < nyrs_hind) {
        # Handle estimated standard deviation
        if(est_sigma_index[index] > 0) {
          index_n_obs[index] <- index_n_obs[index] + 1
          index_q_analytical[index] <- index_q_analytical[index] +
            log(index_obs[index_ind, 1] / index_hat[index_ind])
        }

        # Handle time-varying sigma
        if(est_sigma_index[index] == 0) {
          index_n_obs[index] <- index_n_obs[index] +
            1 / (index_obs[index_ind, 2])^2
          index_q_analytical[index] <- index_q_analytical[index] +
            log(index_obs[index_ind, 1] / index_hat[index_ind]) /
            (index_obs[index_ind, 2])^2
        }
      }
    }
  }

  # Take averages and set analytical q values
  for(index in 1:n_flt) {
    if(index_n_obs[index] > 0) {  # Avoid division by zero
      index_q_analytical[index] <- exp(index_q_analytical[index] / index_n_obs[index])
    }
    if(index_n_obs[index] <= 0){
      index_q_analytical[index] <- NA  # Handle case where no observations
    }
  }

  # Report
  RTMB::REPORT(index_q_analytical)

  return(index_q_analytical)
}


#' Calculate Analytical Standard Deviation Following Ludwig and Walters 1994
#'
#' @description
#' Calculates the analytical standard deviation of survey indices using the method described
#' in Ludwig and Walters (1994).
#'
#' @param index_ctl Matrix containing index control data with columns for index number and year
#' @param index_obs Matrix containing observed index values
#' @param index_hat Vector of predicted index values
#' @param n_flt Integer indicating number of fleets/surveys
#' @param nyrs_hind Integer indicating number of hindcast years
#' @param styr Integer indicating start year
#'
#' @return Vector containing the analytical standard deviations for each fleet/survey
#'
#' @details
#' The function calculates the standard deviation by:
#' 1. Counting number of observations per index
#' 2. Calculating squared log differences between observed and predicted values
#' 3. Taking the square root of the mean squared differences
#'
#' @references
#' Ludwig, D. and Walters, C.J. (1994)
#'
calculate_analytical_sd <- function(index_ctl, index_obs, index_hat, n_flt, nyrs_hind, styr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize vectors
  index_n_obs <- numeric(n_flt)
  ln_index_analytical_sd <- numeric(n_flt)

  # Loop through index control rows
  for(index_ind in 1:nrow(index_ctl)) {
    index <- index_ctl[index_ind, 1] # Survey index
    flt_yr <- index_ctl[index_ind, 2] # Year

    if(flt_yr > 0) {
      flt_yr <- flt_yr - styr

      if(flt_yr < nyrs_hind && index <= n_flt) {
        # Count observations
        index_n_obs[index] <- index_n_obs[index] + 1

        # Calculate squared log difference
        ln_index_analytical_sd[index] <- ln_index_analytical_sd[index] +
          (log(index_obs[index_ind, 1]) - log(index_hat[index_ind]))^2
      }
    }
  }

  # Calculate final standard deviations
  for(index in 1:n_flt) {
    if (index_n_obs[index] > 0) {
      ln_index_analytical_sd[index] <- sqrt(ln_index_analytical_sd[index] / index_n_obs[index])
    }
    if(index_n_obs[index] <= 0) {
      ln_index_analytical_sd[index] <- NA  # Handle division by zero
    }
  }

  return(ln_index_analytical_sd)
}


#' Estimate Catch at Age and Total Yield
#'
#' @description
#' Calculates catch at age and total yield (in kg or numbers) for each fishery index
#'
#' @param catch_ctl Matrix containing control information for catch calculations
#' @param catch_n Matrix containing catch numbers data
#' @param F_flt_age Array of fishing mortality at fleet/sex/age/year
#' @param Z_at_age Array of total mortality at species/sex/age/year
#' @param N_at_age Array of numbers at age
#' @param wt Array of weights
#' @param sel Array of selectivity values
#' @param flt_wt_index Vector of fleet weight indices
#' @param proj_F_prop Vector of projected F proportions
#' @param flt_units Vector indicating fleet units (1=weight, 2=numbers)
#' @param nsex Vector of number of sexes by species
#' @param nages Vector of number of ages by species
#' @param styr Start year
#' @param nyrs_hind Number of hindcast years
#'
#' @return List containing:
#' \itemize{
#'   \item catch_hat - Estimated catch by fishery index
#'   \item max_catch_hat - Maximum possible catch by fishery index
#' }
#'
#' @details
#' The function estimates catch and maximum possible catch for each fishery index.
#' Calculations are done either by weight or by numbers depending on fleet units.
#' Separate calculations are made for hindcast and projection periods.
#'
estimate_catch <- function(catch_ctl, catch_n, F_flt_age, Z_at_age, N_at_age,
                           wt, sel, flt_wt_index, proj_F_prop, flt_units,
                           nsex, nages, styr, nyrs_hind) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  n_fsh <- nrow(catch_ctl)
  catch_hat <- numeric(n_fsh)
  max_catch_hat <- numeric(n_fsh)

  for(fsh_ind in 1:n_fsh) {
    # Get indices
    flt <- catch_ctl[fsh_ind, 1]      # Fishery index
    sp <- catch_ctl[fsh_ind, 2]       # Species index
    flt_yr <- catch_ctl[fsh_ind, 3]   # Year index
    mo <- catch_n[fsh_ind, 1]         # Month index

    # Initialize
    catch_hat[fsh_ind] <- 0
    max_catch_hat[fsh_ind] <- 0

    # Adjust fleet year
    if(flt_yr > 0) {
      flt_yr <- flt_yr - styr + 1
    }
    if(flt_yr < 0) {
      flt_yr <- -flt_yr - styr + 1
    }

    # Set year index for calculations
    yr_ind <- flt_yr
    if(flt_yr > nyrs_hind) {yr_ind <- nyrs_hind}

    # Calculate catch by sex and age
    for(sex in 1:nsex[sp]) {
      for(age in 1:nages[sp]) {

        if(flt_units[flt] == 1) {  # By weight
          catch_hat[fsh_ind] <- catch_hat[fsh_ind] +
            F_flt_age[flt, sex, age, flt_yr] / Z_at_age[sp, sex, age, flt_yr] *
            (1 - exp(-Z_at_age[sp, sex, age, flt_yr])) *
            N_at_age[sp, sex, age, flt_yr] *
            wt[flt_wt_index[flt], sex, age, yr_ind]

          max_catch_hat[fsh_ind] <- max_catch_hat[fsh_ind] +
            N_at_age[sp, sex, age, flt_yr] *
            wt[flt_wt_index[flt], sex, age, yr_ind] *
            sel[flt, sex, age, yr_ind] *
            proj_F_prop[flt]

        }

        if(flt_units[flt] == 2) {  # By numbers
          catch_hat[fsh_ind] <- catch_hat[fsh_ind] +
            F_flt_age[flt, sex, age, flt_yr] / Z_at_age[sp, sex, age, flt_yr] *
            (1 - exp(-Z_at_age[sp, sex, age, flt_yr])) *
            N_at_age[sp, sex, age, flt_yr]

          max_catch_hat[fsh_ind] <- max_catch_hat[fsh_ind] +
            N_at_age[sp, sex, age, flt_yr] *
            sel[flt, sex, age, yr_ind] *
            proj_F_prop[flt]
        }
      }
    }
  }

  RTMB::REPORT(max_catch_hat)
  return(catch_hat = catch_hat)
}


#' Calculate Exploitable Biomass
#'
#' @description
#' Calculates the exploitable biomass for each species and year based on numbers at age,
#' weight at age, selectivity, and projection F proportions.
#'
#' @param n_flt Integer. Number of fleets
#' @param flt_spp Vector. Species index for each fleet
#' @param flt_type Vector. Type of fleet (1 = fishery)
#' @param nyrs Integer. Total number of years
#' @param nyrs_hind Integer. Number of hindcast years
#' @param nages Vector. Number of ages for each species
#' @param nsex Vector. Number of sexes for each species
#' @param N_at_age Array. Numbers at age [species, sex, age, year]
#' @param wt Array. Weight at age [fleet weight index, sex, age, year]
#' @param sel Array. Selectivity [fleet, sex, age, year]
#' @param flt_wt_index Vector. Weight index for each fleet
#' @param proj_F_prop Vector. Projection F proportions by fleet
#'
#' @return Matrix of exploitable biomass [species, year]
#'
calculate_exploitable_biomass <- function(n_flt, flt_spp, flt_type, nyrs, nyrs_hind,
                                          nages, nsex, N_at_age, wt, sel,
                                          flt_wt_index, proj_F_prop) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  # Initialize exploitable biomass matrix
  exploitable_biomass <- matrix(0, nrow = length(flt_spp), ncol = nyrs)

  # Loop through fleets
  for(flt in 1:n_flt) {

    sp <- flt_spp[flt]

    # Only calculate for fishery fleets
    if(flt_type[flt] == 1) {

      for(yr in 1:nyrs) {

        # Set year index for hindcast vs projection
        yr_ind <- yr
        if(yr > nyrs_hind) {yr_ind <- nyrs_hind}

        # Sum across ages and sexes
        for(age in 1:nages[sp]) {
          for(sex in 1:nsex[sp]) {
            exploitable_biomass[sp, yr] <- exploitable_biomass[sp, yr] +
              N_at_age[sp, sex, age, yr] *
              wt[flt_wt_index[flt], sex, age, yr_ind] *
              sel[flt, sex, age, yr_ind] *
              proj_F_prop[flt]
          }
        }
      }
    }
  }

  RTMB::REPORT(exploitable_biomass)
  return(exploitable_biomass)
}



#' Calculate Age and Length Composition
#'
#' @description
#' Calculates predicted age and length compositions for fishery and survey data,
#' accounting for aging error and sex-specific patterns.
#'
#' @param comp_ctl Matrix of composition control data with columns for fleet, species, sex, type and year
#' @param comp_n Matrix of composition sample sizes and month information
#' @param F_flt_age Array of fishing mortality at age by fleet
#' @param Z_at_age Array of total mortality at age
#' @param N_at_age Array of numbers at age
#' @param sel Array of selectivity at age
#' @param index_q Array of catchability coefficients
#' @param age_error Matrix of aging error
#' @param age_trans_matrix Array of age-length transition matrices
#' @param flt_type Vector of fleet types (1=fishery, 2=survey)
#' @param nages Vector of number of ages by species
#' @param nlengths Vector of number of lengths by species
#' @param nsex Vector of number of sexes by species
#' @param styr Start year
#' @param nyrs_hind Number of hindcast years
#' @param flt_age_transition_index Vector of indices for age transition matrices
#' @param comp_obs
#'
#' @return List containing:
#' \itemize{
#'   \item comp_hat: Matrix of predicted compositions
#'   \item age_hat: Matrix of predicted age compositions before aging error
#'   \item age_obs_hat: Matrix of predicted age compositions after aging error
#'   \item true_age_comp_hat: Matrix of true age compositions
#' }
#'
#' @details
#' This function implements the age and length composition calculations for both fishery
#' and survey data. It handles:
#' - Sex-specific and combined sex compositions
#' - Age and length compositions
#' - Aging error
#' - Age-length conversion
#' - Joint sex composition data
#'
#' @examples
#' \dontrun{
#' results <- estimate_comp(
#'   comp_ctl = comp_control_matrix,
#'   comp_n = comp_sample_sizes,
#'   F_flt_age = fishing_mortality,
#'   Z_at_age = total_mortality,
#'   N_at_age = numbers_at_age
#' )
#' }
#'
estimate_comp <- function(comp_ctl, comp_n, comp_obs, F_flt_age, Z_at_age, N_at_age,
                          sel, index_q, age_error, age_trans_matrix,
                          flt_type, nages, nlengths, nsex, styr,
                          nyrs_hind, flt_age_transition_index) {

  # Unsure what is tripping the function
  # "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  "diag<-" <- ADoverload("diag<-")

  # Initialize outputs
  age_obs_hat <- matrix(0, nrow=nrow(comp_obs), ncol=ncol(comp_obs))            # Estimated catch at observed age (accounts for ageing error)
  comp_hat <- matrix(0, nrow=nrow(comp_obs), ncol=ncol(comp_obs))               # Estimated comp
  age_hat <- matrix(0, nrow=nrow(comp_obs), ncol=ncol(comp_obs))                # Estimated catch at true age
  true_age_comp_hat  <- matrix(0, nrow=nrow(comp_obs), ncol=ncol(comp_obs))     # True estimated age composition

  # Loop through composition data
  for(comp_ind in 1:nrow(comp_ctl)) {

    # Extract indices
    flt <- comp_ctl[comp_ind, 1]
    sp <- comp_ctl[comp_ind, 2]
    flt_sex <- comp_ctl[comp_ind, 3]
    comp_type <- comp_ctl[comp_ind, 4]
    yr <- comp_ctl[comp_ind, 5]
    mo <- comp_n[comp_ind, 1]

    # Year calculations
    if(yr > 0) yr <- yr - styr + 1
    if(yr < 0) yr <- -yr - styr + 1

    # Determine year index for projections
    yr_ind <- yr
    if(yr > nyrs_hind) {yr_ind <- nyrs_hind}

    # Calculate catch at age
    for(age in 1:nages[sp]) {

      if(flt_type[flt] == 1) { # Fishery

        # Handle different sex cases
        if(flt_sex == 0) { # Combined sexes
          if(nsex[sp] == 1){
            age_hat[comp_ind, age] <- F_flt_age[flt, 1, age, yr] / Z_at_age[sp, 1, age, yr] *
              (1 - exp(-Z_at_age[sp, 1, age, yr])) * N_at_age[sp, 1, age, yr]
          }
          if(nsex[sp] == 2){
            age_hat[comp_ind, age] <- F_flt_age[flt, 1, age, yr] / Z_at_age[sp, 1, age, yr] *
              (1 - exp(-Z_at_age[sp, 1, age, yr])) * N_at_age[sp, 1, age, yr] +
              F_flt_age[flt, 2, age, yr] / Z_at_age[sp, 2, age, yr] *
              (1 - exp(-Z_at_age[sp, 2, age, yr])) * N_at_age[sp, 2, age, yr]
          }
          # age_hat[comp_ind, age] <- sum(
          #   F_flt_age[flt, 1:nsex[sp], age, yr] / Z_at_age[sp, 1:nsex[sp], age, yr] *
          #     (1 - exp(-Z_at_age[sp, 1:nsex[sp], age, yr])) * N_at_age[sp, 1:nsex[sp], age, yr])
        }

        if (flt_sex >= 1 && flt_sex <= 2) { # Sex-specific composition data
          sex <- flt_sex
          age_hat[comp_ind, age] <- F_flt_age[flt, sex, age, yr] / Z_at_age[sp, sex, age, yr] *
            (1 - exp(-Z_at_age[sp, sex, age, yr])) * N_at_age[sp, sex, age, yr]
        }

        if (flt_sex == 3) { # Joint composition data
          for(sex in 1:nsex[sp]) {
            age_hat[comp_ind, age + nages[sp] * (sex - 1)] <- F_flt_age[flt, sex, age, yr] / Z_at_age[sp, sex, age, yr] *
              (1 - exp(-Z_at_age[sp, sex, age, yr])) * N_at_age[sp, sex, age, yr]
          }
        }
      }

      if (flt_type[flt] == 2) { # Survey
        if(flt_sex == 0) { # Combined sexes
          if(nsex[sp] == 1){
            age_hat[comp_ind, age] <- N_at_age[sp, 1, age, yr] * sel[flt, 1, age, yr_ind] *
              index_q[flt, yr_ind] * exp(-(mo / 12.0) * Z_at_age[sp, 1, age, yr])
          }

          if(nsex[sp] == 2){
            age_hat[comp_ind, age] <- N_at_age[sp, 1, age, yr] * sel[flt, 1, age, yr_ind] *
              index_q[flt, yr_ind] * exp(-(mo / 12.0) * Z_at_age[sp, 1, age, yr]) +
              N_at_age[sp, 2, age, yr] * sel[flt, 2, age, yr_ind] *
              index_q[flt, yr_ind] * exp(-(mo / 12.0) * Z_at_age[sp, 2, age, yr])
          }
          # age_hat[comp_ind, age] <- sum(N_at_age[sp, 1:nsex[sp], age, yr] * sel[flt, 1:nsex[sp], age, yr_ind] *
          #                                 index_q[flt, yr_ind] * exp(-(mo / 12.0) * Z_at_age[sp, 1:nsex[sp], age, yr]))

        }

        if (flt_sex >= 1 && flt_sex <= 2) { # Sex-specific composition data
          sex <- flt_sex
          age_hat[comp_ind, age] <- N_at_age[sp, sex, age, yr] * sel[flt, sex, age, yr_ind] *
            index_q[flt, yr_ind] * exp(-(mo / 12.0) * Z_at_age[sp, sex, age, yr])
        }

        if (flt_sex == 3) { # Joint composition data
          for(sex in 1:nsex[sp]) {
            age_hat[comp_ind, age + nages[sp] * (sex - 1)] <- N_at_age[sp, sex, age, yr] * sel[flt, sex, age, yr_ind] *
              index_q[flt, yr_ind] * exp(-(mo / 12.0) * Z_at_age[sp, sex, age, yr])
          }
        }
      }
    }

    # Adjustment for joint sex composition data
    joint_adjust <- 1
    if(flt_sex == 3){joint_adjust <- 2}


    # Get true age comp
    true_age_comp_hat[comp_ind, 1:(nages[sp] * joint_adjust)] <- age_hat[comp_ind, 1:(nages[sp] * joint_adjust)] / sum(age_hat[comp_ind, 1:(nages[sp] * joint_adjust)])

    # Adjust for aging error
    # - Combined or single-sex
    age_obs_hat[comp_ind, 1:nages[sp]] <- age_hat[comp_ind, 1:nages[sp]] %*% age_error[sp, 1:nages[sp], 1:nages[sp]]

    # # Adjust for aging error for
    # - joint data
    if(flt_sex == 3) {
      age_obs_hat[comp_ind, (nages[sp] + 1):(nages[sp] * 2)] <- age_hat[comp_ind, (nages[sp] + 1):(nages[sp] * 2)] %*% age_error[sp, 1:nages[sp], 1:nages[sp]]
    }

    # Survey catch-at-age - standardize to sum to 1
    if (comp_type == 0) {
      comp_hat[comp_ind, ] <- age_obs_hat[comp_ind, ] / sum(age_obs_hat[comp_ind, ])
    }

    # Catch-at-length
    if (comp_type == 1) {
      sex <- 1
      if(flt_sex > 0 & flt_sex < 3) {sex = flt_sex} # Adjust sex for males/females

      # Convert from catch-at-age to catch-at-length
      # - Combined or single-sex
      comp_hat[comp_ind, 1:nlengths[sp]] <- age_obs_hat[comp_ind, 1:nages[sp]] %*% age_trans_matrix[flt_age_transition_index[flt], sex, 1:nages[sp], 1:nlengths[sp]]

      # Convert from catch-at-age to catch-at-length for
      # - joint comp data
      if (flt_sex == 3) {
        sex = 2
        comp_hat[comp_ind, (nlengths[sp] + 1):(nlengths[sp] * 2)] <- age_obs_hat[comp_ind, (nages[sp] + 1):(nages[sp] * 2)] %*% age_trans_matrix[flt_age_transition_index[flt], sex, 1:nages[sp], 1:nlengths[sp]]
      }

      # Standardize to sum to 1
      comp_hat[comp_ind, ] <- comp_hat[comp_ind, ] / sum(comp_hat[comp_ind, ])
    }
  }

  # Report
  RTMB::REPORT(age_hat)
  RTMB::REPORT(age_obs_hat)
  RTMB::REPORT(true_age_comp_hat)

  return(comp_hat = comp_hat) # matrix(as.numeric(comp_hat), nrow = nrow(comp_hat), ncol = ncol(comp_hat)))
}


#' Calculate Index Likelihood
#'
#' This function calculates the likelihood for index data
#'
#' @param index_obs A matrix containing observed index data.
#' @param index_ctl A matrix containing control indices for the survey.
#' @param est_sigma_index A vector indicating the estimation method for the standard deviation.
#' @param index_ln_sd A vector containing the natural logarithm of the standard deviation estimates.
#' @param ln_index_analytical_sd A vector containing analytical standard deviation estimates.
#' @param jnll_comp A matrix to store the joint negative log-likelihood components.
#' @param index_hat A vector containing the estimated indices.
#' @param flt_type A vector indicating the type of fleet.
#' @param endyr An integer indicating the end year for the analysis.
#'
#' @return A modified jnll_comp matrix with updated values based on the calculations.
#'
calculate_index_nll <- function(index_obs, index_ctl, est_sigma_index, index_ln_sd,
                                ln_index_analytical_sd, jnll_comp,
                                index_hat, flt_type, endyr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  index_std_dev <- 0
  ln_index_sd <- rep(0, nrow(index_obs)) # Estimated/fixed log index sd (kg)


  for(index_ind in 1:nrow(index_obs)) {

    index <- index_ctl[index_ind, 1]      # Temporary survey index
    flt_yr <- index_ctl[index_ind, 3]     # Temporary index for years of data

    if(flt_yr > 0) {
      # Set up variance
      index_std_dev <- switch(est_sigma_index[index] + 1,
                              index_obs[index_ind, 2],               # Case 0: Provided standard deviation
                              exp(index_ln_sd[index]),               # Case 1: Estimated standard deviation
                              ln_index_analytical_sd[index],         # Case 2: Analytical
                              stop("Invalid 'Estimate_sigma_index'"))     # Default case

      ln_index_sd[index_ind] <- index_std_dev

      # Only include years from hindcast
      if(flt_type[index] > 0) {
        if(flt_yr <= endyr) {
          if(index_obs[index_ind, 1] > 0) {
            jnll_comp[1, index] <- jnll_comp[1, index] -
              dnorm(log(index_obs[index_ind, 1]),
                    log(index_hat[index_ind]) - (index_std_dev^2) / 2,
                    index_std_dev,
                    log = TRUE)

            # Martin's version (commented out):
            # jnll_comp[1, index] <- jnll_comp[1, index] +
            #   0.5 * ((log(index_obs[index_ind, 1]) - log(index_hat[index_ind]) +
            #   (index_std_dev^2) / 2) / index_std_dev)^2
          }
        }
      }
    }
  }

  return(jnll_comp)
}


#' Calculate Negative Log Likelihood for Catch Data
#'
#' This function computes the negative log likelihood contributions for fishery observations
#' based on the provided catch data, standard deviations, and other parameters.
#'
#' @param catch_obs A matrix where each row corresponds to an observation and columns contain
#'                  catch data (including observed catches and standard deviations).
#' @param catch_ctl A matrix that contains control data for fisheries, including indices for
#'                  fishery types, species, and years.
#' @param catch_hat A vector of estimated catches.
#' @param catch_ln_sd A vector of log standard deviations for catches.
#' @param est_sigma_fsh A vector indicating the estimation method for standard deviations.
#' @param jnll_comp A matrix to accumulate negative log likelihood contributions.
#' @param F_dev A matrix of fishing mortality deviations.
#' @param flt_type A vector indicating the type of fishery.
#' @param styr The starting year for the analysis.
#' @param endyr The ending year for the analysis.
#'
#' @return The updated jnll_comp matrix with contributions from the fishery observations.
#'
calculate_catch_nll <- function(catch_obs, catch_ctl, catch_hat, catch_ln_sd,
                                est_sigma_fsh, jnll_comp,
                                F_dev, flt_type, styr, endyr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  ln_catch_sd <- rep(0, nrow(catch_obs)) # Estimated/fixed fishery log_sd (kg)


  # Slot 1 -- Total catch -- Fishery observer data
  for(fsh_ind in 1:nrow(catch_obs)) {

    flt <- catch_ctl[fsh_ind, 1]      # Temporary fishery index
    flt_yr <- catch_ctl[fsh_ind, 3]   # Temporary index for years of data
    yr <- flt_yr - styr + 1           # Temporary index of years. Start at 0

    if(flt_yr > 0) {

      # Set up variance
      fsh_std_dev <- switch(est_sigma_fsh[flt] + 1,
                            catch_obs[fsh_ind, 2],  # Provided standard deviation
                            exp(catch_ln_sd[flt]),  # Estimated standard deviation
                            stop("Invalid 'Estimate_sigma_catch'"))

      ln_catch_sd[fsh_ind] <- fsh_std_dev  # Save estimated log_sd

      # Add only years from hindcast
      if(flt_type[flt] == 1) {
        if(flt_yr <= endyr) {
          if(catch_obs[fsh_ind, 1] > 0) {

            # Negative log likelihood contribution
            jnll_comp[2, flt] <- jnll_comp[2, flt] -
              dnorm(log(catch_obs[fsh_ind, 1]),
                    log(catch_hat[fsh_ind]) - fsh_std_dev^2 / 2,
                    fsh_std_dev, log=TRUE)

            # Martin's version (commented out)
            # jnll_comp[2, flt] <- jnll_comp[2, flt] +
            #   0.5 * ((log(catch_obs[fsh_ind, 1]) - log(catch_hat[fsh_ind])) / fsh_std_dev)^2

            # Slot 12 -- Epsilon -- Annual fishing mortality deviation
            jnll_comp[13, flt] <- jnll_comp[13, flt] + F_dev[flt, yr]^2 # FIXME - move to single F per year
          }
        }
      }
    }
  }

  return(jnll_comp)
}


#' Calculate Negative Log-Likelihood for Age/Length Composition
#'
#' This function calculates the joint negative log-likelihood for age/length composition data based on observed and expected proportions.
#'
#' @param comp_obs A matrix of observed proportions, where rows correspond to different observations.
#' @param comp_hat A matrix of expected proportions, where rows correspond to different observations.
#' @param comp_ctl A matrix of control parameters, including fleet index, species index, sex index, composition type, and year.
#' @param comp_n A matrix of sample sizes corresponding to each observation.
#' @param nages A vector containing the number of ages for each species.
#' @param nlengths A vector containing the number of lengths for each species.
#' @param DM_pars A vector of parameters for the Dirichlet-Multinomial distribution.
#' @param flt_type A vector indicating the type of fleet.
#' @param comp_ll_type A vector indicating the likelihood type for composition.
#' @param comp_weights A vector of weights for each fleet.
#' @param jnll_comp A matrix to store the joint negative log-likelihood values.
#' @param endyr An integer indicating the end year for the data.
#'
#' @return A matrix of joint negative log-likelihood values updated in-place.
#'
#' @examples
#' # Example usage:
#' result <- calculate_comp_nll(comp_obs, comp_hat, comp_ctl, comp_n, nages, nlengths, DM_pars, flt_type, comp_ll_type, comp_weights, jnll_comp, endyr)
#'
calculate_comp_nll <- function(comp_obs, comp_hat, comp_ctl, comp_n, nages, nlengths, DM_pars, flt_type, comp_ll_type, comp_weights, jnll_comp, endyr) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")

  jnll_comp[3, ] <- 0
  # comp_nll = matrix(0, nrow(comp_obs), ncol(comp_obs))

  for (comp_ind in 1:nrow(comp_obs)) {

    flt <- comp_ctl[comp_ind, 1]            # Temporary fleet index
    sp <- comp_ctl[comp_ind, 2]             # Temporary index of species
    flt_sex <- comp_ctl[comp_ind, 3]        # Temporary index for comp sex (0 = combined, 1 = female, 2 = male, 3 = joint)
    comp_type <- comp_ctl[comp_ind, 4]      # Temporary index for comp type (0 = age, 1 = length)
    yr <- comp_ctl[comp_ind, 5]             # Temporary index for years of data

    # Adjustment for joint sex composition data
    joint_adjust <- 1
    if(flt_sex == 3){joint_adjust <- 2}

    # Number of ages/lengths
    if(comp_type == 0){
      n_comp <- nages[sp] * joint_adjust
    }

    if(comp_type == 1){
      n_comp <- nlengths[sp] * joint_adjust
    }

    # Select sections
    # - as.numeric solves error in AD for some-reason
    #FIXME
    comp_obs_tmp <- (comp_obs[comp_ind, 1:n_comp]) + 0.00001  # Observed proportion with offset
    comp_hat_tmp <- (comp_hat[comp_ind, 1:n_comp]) + 0.00001  # Expected proportion with offset

    # Convert observed prop to observed numbers
    comp_obs_tmp <- comp_obs_tmp * comp_n[comp_ind, 2]
    alphas <- sum(comp_obs_tmp) * comp_hat_tmp * DM_pars[flt] # DM alpha

    if (!is.na(sum(comp_obs_tmp)) && sum(comp_obs_tmp) > 0) { # Error checking
      # Only use years wanted
      if (yr <= endyr && yr > 0 && flt_type[flt] > 0) {

        switch(comp_ll_type[flt] + 2,

               # case -1 (now case 1)
               {
                 for(ln in 1:n_comp) {
                   # Martin's
                   jnll_comp[3, flt] <- jnll_comp[3, flt] -
                     (comp_weights[flt] * comp_n[comp_ind, 2] *
                        (comp_obs[comp_ind, ln] + 0.00001) *
                        log((comp_hat[comp_ind, ln] + 0.00001) / (comp_obs[comp_ind, ln] + 0.00001)))
                   # comp_nll[comp_ind, ln] <- comp_weights[flt] * comp_n[comp_ind, 2] *
                   #   (comp_obs[comp_ind, ln] + 0.00001) *
                   #   log((comp_hat[comp_ind, ln] + 0.00001) / (comp_obs[comp_ind, ln] + 0.00001))

                 }
               },

               # case 0 (now case 2) -- Full multinomial
               {
                 jnll_comp[3, flt] <- jnll_comp[3, flt] -
                   comp_weights[flt] * dmultinom(comp_obs_tmp, prob = comp_hat_tmp, log = TRUE)
               },

               # case 1 (now case 3) -- Dirichlet-multinomial
               {
                 jnll_comp[3, flt] <- jnll_comp[3, flt] -
                   ddirmultinom(comp_obs_tmp, alphas, log = TRUE)
               },

               # default
               stop("Invalid 'comp_ll_type'")
        )
      }
    }
  }

  return(jnll_comp)

  # list(jnll_comp = jnll_comp,
  #             comp_nll = comp_nll))  # Return the updated jnll_comp matrix
}


#' Calculate Survey Catchability Penalties for NLL
#'
#' This function computes the survey catchability deviates based on various
#' parameters and conditions. It updates the joint negative log-likelihood
#' components for each fleet.
#'
#' @param n_flt Integer. The number of fleets.
#' @param nyrs_hind Integer. The number of years for hindcasting.
#' @param est_index_q Numeric vector. Estimated index of catchability for each fleet.
#' @param index_varying_q Numeric vector. Indicates if catchability is time-varying for each fleet.
#' @param flt_type Numeric vector. Type of each fleet.
#' @param index_ln_q Numeric vector. Logarithm of catchability index.
#' @param index_ln_q_prior Numeric vector. Prior values for the catchability index.
#' @param index_q_sd Numeric vector. Standard deviation for the catchability index.
#' @param index_q_dev Matrix. Deviates for catchability index.
#' @param index_q_dev_sd Numeric vector. Standard deviation for the deviates.
#' @param env_index Matrix. Environmental index values.
#' @param rho_trans Function. A function to transform rho values.
#' @param jnll_comp Matrix. Joint negative log-likelihood components to be updated.
#' @param AR1 Function. A function to compute the AR1 process.
#' @param SCALE Function. A function to scale the AR1 process.
#'
#' @return Updated joint negative log-likelihood components matrix.
#'
calculate_catchability_nll <- function(n_flt, nyrs_hind, est_index_q, index_varying_q,
                                       flt_type, index_ln_q, index_ln_q_prior,
                                       index_q_sd, index_q_dev, index_q_dev_sd,
                                       env_index, rho_trans, jnll_comp, AR1, SCALE) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  for(flt in 1:n_flt) {

    est_q <- est_index_q[flt]
    vary_q <- index_varying_q[flt]
    q_dev_sd <- index_q_dev_sd[flt]

    if(flt_type[flt] == 2){

      # Prior on catchability
      if(est_q == 2) {
        jnll_comp[7, flt] <- jnll_comp[7, flt] - dnorm(index_ln_q[flt],
                                                       mean = index_ln_q_prior[flt],
                                                       sd = index_q_sd[flt],
                                                       log = TRUE)
      }

      # QAR1 deviates fit to environmental index
      if(est_q == 6) {
        rho <- rho_trans(index_q_rho[flt])
        index_q_dev_tmp <- index_q_dev[flt, ]

        # Using TMB::SCALE and TMB::AR1
        jnll_comp[7, flt] <- SCALE(AR1(rho), q_dev_sd)(index_q_dev_tmp)

        # Observation error - Fit to environmental index
        q_index <- index_varying_q[flt] - 1
        jnll_comp[8, flt] <- jnll_comp[8, flt] - sum(dnorm(env_index[1:nyrs_hind, q_index],
                                                           mean = index_q_dev[flt, 1:nyrs_hind],
                                                           sd = q_dev_sd,
                                                           log = TRUE))
      }

      # Penalized/random deviate likelihood
      if(((vary_q == 1) | (vary_q == 2)) & (est_q == 1 | est_q == 2)) {
        jnll_comp[8, flt] <- jnll_comp[8, flt] - sum(dnorm(index_q_dev[flt, 1:nyrs_hind],
                                                           mean = 0,
                                                           sd = q_dev_sd,
                                                           log = TRUE))
      }

      # Random walk
      if((vary_q == 4) & (est_q == 1 | est_q == 2)) {
        jnll_comp[8, flt] <- jnll_comp[8, flt] - sum(dnorm(index_q_dev[flt, 2:nyrs_hind] -
                                                             index_q_dev[flt, 1:(nyrs_hind - 1)],
                                                           mean = 0,
                                                           sd = q_dev_sd,
                                                           log = TRUE))
      }
    }
  }

  return(jnll_comp)
}


#' Calculate Selectivity Penalties and Deviates
#'
#' This function computes the selectivity penalties and likelihood components for various selectivity types
#' based on survey data. It handles both non-parametric and parametric selectivity models.
#'
#' @param n_flt Integer. The number of surveys.
#' @param flt_spp A vector of species associated with each survey.
#' @param flt_type A vector indicating the type of each survey.
#' @param flt_sel_type A vector indicating the selectivity type for each survey.
#' @param flt_varying_sel A vector indicating whether the selectivity is varying for each survey.
#' @param non_par_sel A 3D array of non-parametric selectivities.
#' @param sel_curve_pen A 2D array of selectivity curve penalties.
#' @param avg_sel A 2D array of average selectivities.
#' @param sel_inf_dev A 4D array of selectivity influence deviations.
#' @param ln_sel_slp_dev A 4D array of log selectivity slope deviations.
#' @param sel_dev_sd A vector of standard deviations for selectivity deviations.
#' @param flt_nselages A vector indicating the number of selectivity ages for each survey.
#' @param nsex A vector indicating the number of sexes for each species.
#' @param nages A vector indicating the number of ages for each species.
#' @param nyrs_hind Integer. The number of years for hindcasting.
#' @param jnll_comp A matrix to store the joint negative log-likelihood components.
#'
#' @return A matrix of joint negative log-likelihood components updated with selectivity penalties and deviates.
#'
calculate_selectivity_nll <- function(n_flt, flt_spp, flt_type, flt_sel_type, flt_varying_sel,
                                      non_par_sel, sel_curve_pen, avg_sel, sel_inf_dev,
                                      ln_sel_slp_dev, sel_dev_sd, flt_nselages, nsex, nages,
                                      nyrs_hind, jnll_comp) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  for(flt in 1:n_flt) { # Loop around surveys
    jnll_comp[4, flt] <- 0
    jnll_comp[5, flt] <- 0
    sp <- flt_spp[flt]

    # Ianelli non-parametric selectivity penalties
    if(flt_type[flt] > 0 && flt_sel_type[flt] == 2) {

      for(sex in 1:nsex[sp]) {
        for(age in 1:(nages[sp] - 1)) {
          # if( non_par_sel[flt, sex, age] > non_par_sel[flt, sex, age+1]) { # not differentiable so using max2
          # Non-parametric selectivity penalties
          jnll_comp[5, flt] = jnll_comp[5, flt] + sel_curve_pen[flt, 1] *
            max2( log(non_par_sel[flt, sex, age] / non_par_sel[flt, sex, age+1]), 1e-16)^2
          # }
        }

        jnll_comp[5, flt] = jnll_comp[5, flt] + sum(sel_curve_pen[flt, 2] * diff( diff( log( (non_par_sel[flt, sex, 1:nages[sp]])))) ^ 2)
      }



      # Survey selectivity normalization (non-parametric)
      jnll_comp[5, flt] <- jnll_comp[5, flt] + sum(avg_sel[flt, 1:nsex[sp]]^2)
    }

    # Penalized/random effect likelihood time-varying logistic/double-logistic selectivity deviates
    if((flt_varying_sel[flt] %in% c(1, 2)) && (flt_sel_type[flt] != 2) && (flt_sel_type[flt] != 5) && (flt_type[flt] > 0)) {

      for(sex in 1:nsex[sp]) {
        for(yr in 1:nyrs_hind) {
          jnll_comp[5, flt] <- jnll_comp[5, flt] -
            dnorm(sel_inf_dev[1, flt, sex, yr], 0, sel_dev_sd[flt], log = TRUE) -
            dnorm(ln_sel_slp_dev[1, flt, sex, yr], 0, 4 * sel_dev_sd[flt], log = TRUE)

          # Double logistic deviates
          if(flt_sel_type[flt] == 3) {
            jnll_comp[5, flt] <- jnll_comp[5, flt] -
              dnorm(sel_inf_dev[2, flt, sex, yr], 0, sel_dev_sd[flt], log = TRUE) -
              dnorm(ln_sel_slp_dev[2, flt, sex, yr], 0, 4 * sel_dev_sd[flt], log = TRUE)
          }
        }
      }
    }

    # Penalized/random effect likelihood time-varying non-parametric (Taylor et al 2014) selectivity deviates
    if((flt_varying_sel[flt] %in% c(1, 2)) && (flt_sel_type[flt] == 5) && (flt_type[flt] > 0)) {

      for(age in 1:flt_nselages[flt]) {
        for(sex in 1:nsex[sp]) {
          for(yr in 1:nyrs_hind) {
            jnll_comp[5, flt] <- jnll_comp[5, flt] -
              dnorm(sel_coff_dev[flt, sex, age, yr], 0, sel_dev_sd[flt], log = TRUE)
          }
        }
      }
    }

    # Random walk: Type 4 = random walk on ascending and descending for double logistic;
    # Type 5 = ascending only for double logistics
    if((flt_varying_sel[flt] %in% c(4, 5)) && (flt_sel_type[flt] != 2) && (flt_sel_type[flt] != 5) && (flt_type[flt] > 0)) {

      for(sex in 1:nsex[sp]) {
        for(yr in 2:nyrs_hind) { # Start at second year
          jnll_comp[5, flt] <- jnll_comp[5, flt] -
            dnorm(ln_sel_slp_dev[1, flt, sex, yr] - ln_sel_slp_dev[1, flt, sex, yr-1], 0, sel_dev_sd[flt], log = TRUE) -
            dnorm(sel_inf_dev[1, flt, sex, yr] - sel_inf_dev[1, flt, sex, yr-1], 0, 4 * sel_dev_sd[flt], log = TRUE)

          # Double logistic deviates
          if((flt_sel_type[flt] == 3) && (flt_varying_sel[flt] == 4)) {
            jnll_comp[5, flt] <- jnll_comp[5, flt] -
              dnorm(sel_inf_dev[2, flt, sex, yr] - sel_inf_dev[2, flt, sex, yr-1], 0, sel_dev_sd[flt], log = TRUE) -
              dnorm(ln_sel_slp_dev[2, flt, sex, yr] - ln_sel_slp_dev[2, flt, sex, yr-1], 0, 4 * sel_dev_sd[flt], log = TRUE)
          }
        }
      }
    }
  } # End selectivity loop

  return(jnll_comp)
}



#' Calculate Joint Negative Log-Likelihood Components
#'
#' This function calculates the joint negative log-likelihood components for recruitment parameters
#' based on stock-recruitment relationships and other model parameters.
#'
#' @param nspp Number of species.
#' @param srr_est_mode Stock-recruitment estimation mode.
#' @param srr_pred_fun Stock-recruitment prediction function.
#' @param steepness A numeric vector of steepness values for each species.
#' @param srr_prior A numeric vector of prior values for stock-recruitment.
#' @param srr_prior_sd A numeric vector of standard deviations for the priors.
#' @param rec_pars A matrix of recruitment parameters.
#' @param Bmsy_lim A numeric vector of Bmsy limits.
#' @param initMode Initialization mode.
#' @param nages A numeric vector of the number of ages for each species.
#' @param init_dev A matrix of initial abundance-at-age values.
#' @param rec_dev A matrix of annual recruitment deviations.
#' @param R_sd A numeric vector of standard deviations for recruitment.
#' @param R A matrix of recruitment values.
#' @param R_hat A matrix of estimated recruitment values.
#' @param srr_fun Stock-recruitment function type.
#' @param srr_hat_styr Starting year for SRR hat.
#' @param srr_hat_endyr Ending year for SRR hat.
#' @param nyrs_hind
#' @param jnll_comp A matrix to store the joint negative log-likelihood components.
#'
#' @return Updated jnll_comp matrix.
#'
calculate_recruitment_nll <- function(nspp, srr_est_mode, srr_pred_fun, steepness, srr_prior,
                                      srr_prior_sd, rec_pars, Bmsy_lim, initMode, nages,
                                      init_dev, rec_dev, R_sd, R, R_hat, srr_fun,
                                      srr_hat_styr, srr_hat_endyr, nyrs_hind, jnll_comp) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  for (sp in 1:nspp) {

    # Slot 9 -- stock-recruit prior for Beverton -- Lognormal
    if ((srr_est_mode == 2) && (srr_pred_fun == 2 || srr_pred_fun == 3)) {
      jnll_comp[9, sp] <- jnll_comp[9, sp] - dnorm(log(steepness[sp]),
                                                   log(srr_prior[sp]) + (srr_prior_sd[sp])^2 / 2.0,
                                                   srr_prior_sd[sp], log = TRUE)
    }

    # -- Beta
    if ((srr_est_mode == 3) && (srr_pred_fun == 2 || srr_pred_fun == 3)) {
      beta_alpha <- ((1 - srr_prior[sp]) / (srr_prior_sd[sp])^2 - 1 / srr_prior[sp]) * (srr_prior[sp])^2
      beta_beta <- beta_alpha * (1 / srr_prior[sp] - 1)
      jnll_comp[9, sp] <- jnll_comp[9, sp] - dbeta(steepness[sp], beta_alpha, beta_beta, log = TRUE)
    }

    # Slot 9 -- stock-recruit prior for Ricker
    if ((srr_est_mode == 2) && (srr_pred_fun == 4 || srr_pred_fun == 5)) {
      jnll_comp[9, sp] <- jnll_comp[9, sp] - dnorm(rec_pars[sp, 1],
                                                   log(srr_prior[sp]),
                                                   srr_prior_sd[sp],
                                                   log = TRUE)
    }

    # Slot 9 -- penalty for Bmsy > Bmsy_lim for Ricker
    if ((Bmsy_lim[sp] > 0) && ((srr_pred_fun == 4) || (srr_pred_fun == 5))) {
      bmsy <- 1.0 / exp(rec_pars[sp, 2])
      # pos_tmp <- posfun(Bmsy_lim[sp] / 1000000.0 - bmsy, 0.001)
      # bmsy <- pos_tmp$ans + 1.0
      # jnll_comp[9, sp] <- jnll_comp[9, sp] + 100 *  pos_tmp$penalty
    }

    # Slot 10 -- init_dev -- Initial abundance-at-age
    if (initMode > 0) {
      for (age in 2:nages[sp]) {
        jnll_comp[12, sp] <- jnll_comp[12, sp] - dnorm(init_dev[sp, age - 1],
                                                       (R_sd[sp])^2 / 2.0,
                                                       R_sd[sp],
                                                       log = TRUE)
      }
    }

    # Slot 11 -- Tau -- Annual recruitment deviation
    for (yr in 1:nyrs_hind) {
      jnll_comp[11, sp] <- jnll_comp[11, sp] - dnorm(rec_dev[sp, yr],
                                                     (R_sd[sp])^2 / 2.0,
                                                     R_sd[sp],
                                                     log = TRUE)
    }

    # Additional penalty for SRR curve (sensu AMAK/Ianelli)
    if ((srr_fun == 0) && (srr_pred_fun > 0)) {
      for (yr in srr_hat_styr:srr_hat_endyr) {
        jnll_comp[9, sp] <- jnll_comp[9, sp] - dnorm(log(R[sp, yr]),
                                                     log(R_hat[sp, yr]),
                                                     R_sd[sp],
                                                     log = TRUE)
      }
    }
  }

  return(jnll_comp)
}



#' Calculate the Joint Negative Log-Likelihood Component for M_at_age
#'
#' This function computes the joint negative log-likelihood component for M_at_age
#' based on the specified models and priors for multiple species, sexes, ages, and years.
#'
#' @param nspp Number of species.
#' @param nsex A vector containing the number of sexes for each species.
#' @param nages A vector containing the number of ages for each species.
#' @param nyrs Number of years.
#' @param M1_model A vector indicating the model type for each species.
#' @param M1_use_prior A vector indicating whether to use the prior for M1 for each species.
#' @param M2_use_prior A vector indicating whether to use the prior for M2 for each species.
#' @param M1_at_age A 3D array containing M1_at_age values for species, sexes, and ages.
#' @param M_at_age A 4D array containing total M_at_age values for species, sexes, ages, and years.
#' @param M_prior A vector containing the prior means for each species.
#' @param M_prior_sd A vector containing the prior standard deviations for each species.
#' @param jnll_comp A matrix to store the joint negative log-likelihood components.
#'
#' @return A matrix with updated joint negative log-likelihood components.
#'
calculate_mortality_nll <- function(nspp, nsex, nages, nyrs, M1_model, M1_use_prior, M2_use_prior,
                                    M1_at_age, M_at_age, M_prior, M_prior_sd, jnll_comp) {

  # Unsure what is tripping the function
  "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g

  for(sp in 1:nspp) {
    # Prior on M1_at_age only and using species specific M1_at_age
    if(M1_model[sp] == 1 && M1_use_prior[sp] == 1 && M2_use_prior[sp] == 0) {
      jnll_comp[15, sp] <- jnll_comp[15, sp] - dnorm(
        log(M1_at_age[sp, 1, 1]),
        log(M_prior[sp]) + (M_prior_sd[sp]^2) / 2,
        M_prior_sd[sp]
      )
    }

    for(sex in 1:nsex[sp]) {
      # Prior on M1_at_age only and using species and sex specific M1_at_age
      if(M1_model[sp] == 2 && M1_use_prior[sp] == 1 && M2_use_prior[sp] == 0) {
        jnll_comp[15, sp] <- jnll_comp[15, sp] - dnorm(
          log(M1_at_age[sp, sex, 1]),
          log(M_prior[sp]) + (M_prior_sd[sp]^2) / 2,
          M_prior_sd[sp]
        )
      }

      for(age in 1:nages[sp]) {
        # Prior on M1_at_age only and using species and sex specific M1_at_age
        if(M1_model[sp] == 3 && M1_use_prior[sp] == 1 && M2_use_prior[sp] == 0) {
          jnll_comp[15, sp] <- jnll_comp[15, sp] - dnorm(
            log(M1_at_age[sp, sex, age]),
            log(M_prior[sp]) + (M_prior_sd[sp]^2) / 2,
            M_prior_sd[sp]
          )
        }

        for(yr in 1:nyrs) {
          # Prior on total M_at_age (M1_at_age and M2_at_age)
          if(M1_use_prior[sp] == 1 && M2_use_prior[sp] == 1) {
            jnll_comp[15, sp] <- jnll_comp[15, sp] - dnorm(
              log(M_at_age[sp, sex, age, yr]),
              log(M_prior[sp]) + (M_prior_sd[sp]^2) / 2,
              M_prior_sd[sp]
            )
          }
        }
      }
    }
  }

  return(jnll_comp)
}


#' #' Maintain a value above a limit with penalty adjustment
#' #'
#' #' @param x A numeric value.
#' #' @param eps A numeric threshold.
#' #' @param penalty A numeric value representing the current penalty.
#' #' @return A list containing the adjusted value and updated penalty.
#' posfun <- function(x, eps = 0.001) {
#'   denom <- 2 - x / eps
#'   ans <- ifelse(x >= eps, x, eps / denom)
#'   penalty <- ifelse(x < eps, 0.01 * (x - eps)^2, 0)
#'   return(list(ans = ans, penalty = penalty))
#' }

#' Calculate the Dirichlet multinomial log-likelihood
#'
#' @param obs A numeric vector of observed counts.
#' @param alpha A numeric vector of Dirichlet parameters.
#' @param do_log A logical indicating whether to return the log-likelihood.
#' @return The log-likelihood or the exponentiated likelihood.
ddirmultinom <- function(obs, alpha, do_log = FALSE) {
  N <- sum(obs)
  phi <- sum(alpha)
  ll <- lgamma(N + 1) + lgamma(phi) - lgamma(N + phi)

  for (a in seq_along(obs)) {
    ll <- ll - lgamma(obs[a] + 1) + lgamma(obs[a] + alpha[a]) - lgamma(alpha[a])
  }

  return(if (do_log) ll else exp(ll))
}

#' Transform a value to be between -1 and 1
#'
#' @param x A numeric value.
#' @return A transformed numeric value between -1 and 1.
rho_trans <- function(x) {
  return(2 / (1 + exp(-2 * x)) - 1)
}


# load("~/Documents/GitHub/Rceattle/ss ebs comp minus1.RData")
# params <- mod_objects$estimated_params
# data_list <- rearrange_dat(mod_objects$data_list)
# data_list$forecast <- c(0,0,0)
# data_list$Ceq = rep(1,3)
# data_list$avgnMode = 0
# # rtmb_ceattle(params, data_list)
# mod_objects$quantities$jnll_comp

max2 <- function(x,y){
  return(0.5*(abs(x-y)+x+y))
}

#' Title
#'
#' @param params
#' @param data_list
#'
#' @return

#'
#' @examples
#'
rtmb_ceattle <- function(start_par, data_list_reorganized){
  require(RTMB)
  require(Rceattle)
  "[<-" <- RTMB::ADoverload("[<-") # AD overload https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g?pli=1
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")

  RTMB::getAll(start_par, data_list_reorganized)


  # ------------------------------------------------------------------------- #
  # 1. DATA TRANSFORMATION ----
  # ------------------------------------------------------------------------- #
  nyrs = projyr - styr + 1
  nyrs_hind = endyr - styr + 1

  suit_endyr = suit_endyr - styr + 1
  suit_styr = suit_styr - styr + 1
  nyrs_suit = suit_endyr - suit_styr + 1
  nyrs_srrmean = srr_meanyr - styr + 1

  srr_hat_styr = srr_hat_styr - styr + 1
  srr_hat_endyr = srr_hat_endyr - styr + 1
  if(nyrs_srrmean > nyrs_hind){nyrs_srrmean = nyrs_hind}
  if(srr_hat_styr == 1){srr_hat_styr = 2} # R_hat starts at > year 1.

  mo = 0;                                                              # Month float
  if (msmMode == 0) { niter = 1 }                                      # Number of iterations for SS mode
  n_flt <- nrow(fleet_control)
  # comp_obs <- matrix(as.numeric(comp_obs), nrow = nrow(comp_obs), ncol = ncol(comp_obs))


  # ------------------------------------------------------------------------- #
  # 2. FIXED INPUTS ----
  # ------------------------------------------------------------------------- #

  # 2.1. FIXED VALUES
  sd_ration = 0.05                  # SD of ration likelihood
  stom_tau = 20                     # Stomach sample size: FIXME - have as input


  # ------------------------------------------------------------------------- #
  # 3. PARAMETER TRANSFORMATION  ----
  # ------------------------------------------------------------------------- #
  DM_pars = exp(comp_weights) # Dirichlet-multinomial scalars
  # pop_scalar <- ln_pop_scalar; pop_scalar <- exp(ln_pop_scalar)             # Fixed n-at-age scaling coefficient; n = [nspp, nages]
  M1_at_age <- exp(ln_M1)                                                       # Residual or total natural mortality at age

  # - F reference points
  Flimit <- exp(ln_Flimit)                                                     # Target F parameter on natural scale
  Ftarget <- exp(ln_Ftarget)                                                   # Limit F parameter on natural scale
  Finit <- exp(ln_Finit)                                                       # Initial F for non-equilibrium age-structure

  # - Suitability
  # gam_a <- exp(log_gam_a)                                                    # Predator size-selectivity: shape parameter for gamma suitability, mean for normal of logs
  # gam_b <- exp(log_gam_b)                                                    # Predator size-selectivity: scale parameter for gamma suitability, sd for normal of logs

  # - Variance terms
  R_sd <- exp(R_ln_sd)  # Convert log sd to natural scale
  sel_dev_sd <- exp(sel_dev_ln_sd)
  index_q_sd <- exp(index_q_ln_sd)
  index_q_dev_sd <- exp(index_q_dev_ln_sd)
  # sex_ratio_sd = exp(sex_ratio_ln_sd)

  # -- Kinzey Functional response parameters
  # H_1 <- exp(logH_1)                                                        # Predation functional form; n = [nspp, nspp2];
  # H_1a <- exp(logH_1a)                                                      # Age adjustment to H_1; n = [1, nspp];
  # H_1b <- exp(logH_1b)                                                      # Age adjustment to H_1; n = [1, nspp];
  # H_2 <- exp(logH_2)                                                        # Predation functional form; n = [nspp, nspp]
  # H_3 <- exp(logH_3)                                                        # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;

  # ------------------------------------------------------------------------- #
  # 4. DERIVED QUANTITIES SECTION  ----
  # ------------------------------------------------------------------------- #

  # -- 4.2. Estimated population quantities
  M2_at_age <- array(0, dim=c(nspp, max_nsex, max_nages, nyrs))                    # Total predation mortality at age
  N_at_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                       # Numbers at age
  index_q <- matrix(0, n_flt, nyrs_hind)                                       # Estimated survey catchability
  M_at_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                       # Total natural mortality at age
  Z_at_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                       # Total mortality at age
  R <- matrix(0, nspp, nyrs)                                                    # Estimated recruitment (n)

  # R_hat <- matrix(0, nspp, nyrs)                                               # Expected recruitment given SR curve
  # mort_sum <- matrix(0, nspp, max_nages)
  # R0 <- rep(0, nspp)                                                             # Equilibrium recruitment at F = 0.                                                    # Equilibrium recruitment at F = Finit (non-equilibrium).
  # srr_alpha <- 0.0
  # steepness <- rep(0, nspp)                                                     # Expected % of R0 at 20% SSB0.
  # exploitable_biomass <- matrix(0, nspp, nyrs)                                 # Estimated exploitable biomass (kg)
  # biomass_depletion <- matrix(0, nspp, nyrs)                                   # Estimated biomass depletion
  # ssb_depletion <- matrix(0, nspp, nyrs)                                       # Estimated biomass depletion of spawning stock biomass
  #
  #
  # # sex_ratio_hat <- array(0, dim = c(nspp, max_nages, nyrs))                   # Estimated age-specific sex ratio
  # # avg_sex_ratio_hat <- matrix(0, nspp, nyrs)                                 # Estimated sex ratio across all ages
  # # S_at_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                     # Survival at age
  # R_sd <- rep(0, nspp)                                                           # Standard deviation of recruitment variation
  # zero_N_pen <- rep(0, nspp)                                                   # Additional penalty to add to likelihood if n-at-age goes < 0
  # # sex_ratio_sd <- rep(0, nspp)                                               # Variance of sex ratio
  #
  # # -- 4.3. Selectivity parameters
  # sel_dev_sd <- rep(0, n_flt)                                                  # Standard deviation of selectivity deviates
  #
  # # -- 4.4. Fishery components
  # # catch_hat <- rep(0, nrow(catch_obs))                                         # Estimated fishery yield/numbers (kg)
  # # max_catch_hat <- rep(0, nrow(catch_obs))                                     # Estimated exploitable biomass/numbers by fleet (kg)
  # ln_catch_sd <- rep(0, nrow(catch_obs))                                       # Estimated/fixed fishery log_sd (kg)
  #
  # # -- 4.5. Biological reference points
  # NByage0 <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                        # Numbers at age at mean recruitment and F = 0
  # NByageF <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                        # Numbers at age at mean recruitment and F = Flimit
  # DynamicNByage0 <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                 # Numbers at age at F = 0 (accounts for annual recruitment)
  # DynamicNByageF <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                 # Female numbers at age at F = Ftarget (accounts for annual recruitment)
  # DynamicSB0 <- matrix(0, nspp, nyrs)                                         # Estimated dynamic spawning biomass at F = 0 (accounts for S_at_age-R curve)
  # DynamicB0 <- matrix(0, nspp, nyrs)                                          # Estimated dynamic biomass at F = 0 (accounts for S_at_age-R curve)
  # DynamicSBF <- matrix(0, nspp, nyrs)                                         # Estimated dynamic spawning biomass at F = Ftarget (accounts for S_at_age-R curve)
  # NbyageSPR <- array(0, dim = c(4, nspp, max_nages))                            # Estimated numbers at age for spawning biomass per recruit reference points
  # SPRlimit <- rep(0, nspp)                                                     # Estimated Plimit SPR
  # SPRtarget <- rep(0, nspp)                                                    # Estimated Ptarget SPR
  # SPRFinit <- rep(0, nspp)                                                     # Estimated spawning biomass per recruit at Finit
  # SB0 <- matrix(0, nspp, nyrs)                                                # Estimated spawning stock biomass at F = 0 (Accounts for S_at_age-R)
  # SBF <- matrix(0, nspp, nyrs)                                                # Estimated spawning stock biomass at F = target (Accounts for S_at_age-R)
  # B0 <- matrix(0, nspp, nyrs)                                                  # Estimated biomass at F = 0 (Accounts for S_at_age-R)
  # srr_mult <- rep(0, ncol(beta_rec_pars))                                     # Environmental design matrix for rec
  # index_q_mult <- rep(0, ncol(index_q_beta))                                   # Environmental design matrix for q
  # beta_rec_tmp <- rep(0, ncol(beta_rec_pars))                                 # Temporary vector to store beta parameters by species for matrix mult
  # env_rec_tmp <- rep(0, ncol(beta_rec_pars))                                   # Temporary vector to store env data by year for matrix mult
  # beta_q_tmp <- rep(0, ncol(index_q_beta))                                     # Temporary vector to store Q beta parameters by species for matrix mult
  # env_q_tmp <- rep(0, ncol(index_q_beta))                                       # Temporary vector to store Q env data by year for matrix mult
  # proj_F <- matrix(0, nspp, nyrs)                                              # Projected F (Fabc/Ftac/etc) using harvest control rule
  #
  #
  # # -- 4.6. Survey components
  # index_q_sd <- rep(0, n_flt)                                                  # Vector of standard deviation of survey catchability prior
  # index_q_dev_sd <- rep(0, n_flt)                                              # Vector of standard deviation of time-varying survey catchability deviation
  # avgsel_tmp <- 0                                                                # Temporary object for average selectivity across all ages
  # index_hat <- rep(0, nrow(index_obs))                                          # Estimated survey biomass (kg)
  # ln_index_analytical_sd <- rep(0, n_flt)                                      # Temporary vector to save analytical sd follow Ludwig and Walters 1994
  # index_q_analytical <- rep(0, n_flt)                                          # Temporary vector to save analytical sd follow Ludwig and Walters 1994
  # index_n_obs <- rep(0, n_flt)                                                  # Vector to save the number of observations for each survey time series
  #
  #
  # # -- 4.10. Suitability components
  # avail_food <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                     # Available food to predator
  # stom_div_bio <- array(0, dim = c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs)) # Stomach proportion over biomass; U/ (W * N)
  # suit_main <- array(0, dim = c(nspp, max_nsex, max_nages, nspp, max_nsex, max_nages, nyrs))  # Suitability/gamma selectivity of predator age u on prey age a
  # suit_other <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                     # Suitability not accounted for by the included prey
  # suma_suit <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))                      # Sum of suitabilities
  # r_sexes <- matrix(0, nrow = nrow(stom_prop_obs), ncol = 2)
  # k_sexes <- matrix(0, nrow = nrow(stom_prop_obs), ncol = 2)






  # ------------------------------------------------------------------------- #
  # 5. INITIAL CALCULATIONS -----
  # ------------------------------------------------------------------------- #

  # * 5.1. MATURITY AND SEX RATIO ----
  pmature_sexr <- matrix(0, nspp, max_nages)
  for (sp in 1:nspp) {

    # Sex ratio at recruitment
    if(nsex[sp] == 1) {
      R_sexr[sp] <- 1.0
    }

    if(nsex[sp] == 1) {
      # Multiply sex_ratio and pmature for 1 sex models
      pmature_sexr[sp, 1:nages[sp]] <- pmature[sp, 1:nages[sp]] * sex_ratio[sp, 1:nages[sp]]
    }
  }


  # * 5.2. SURVEY CONTROL SWITCHES ----
  flt_ind <- fleet_control[, 2]                    # Temporary fleet index
  flt_type <- fleet_control[flt_ind, 3]            # Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
  flt_spp <- fleet_control[flt_ind, 4]             # Species
  flt_sel_ind <- fleet_control[flt_ind, 5]         # Survey selectivity index
  flt_sel_type <- fleet_control[flt_ind, 6]        # Selectivity type
  flt_nselages <- fleet_control[flt_ind, 7]        # Non-parametric selectivity ages
  flt_varying_sel <- fleet_control[flt_ind, 8]     # Time-varying selectivity type
  flt_sel_age <- fleet_control[flt_ind, 9] - minage[flt_spp] + 1 # First age selected
  flt_sel_maxage <- fleet_control[flt_ind, 10] - minage[flt_spp]  + 1# Age of max selectivity
  comp_ll_type <- fleet_control[flt_ind, 11]       # Index for dirichlet multinomial
  flt_units <- fleet_control[flt_ind, 12]          # Survey units
  flt_wt_index <- fleet_control[flt_ind, 13]       # Dim1 of wt
  flt_age_transition_index <- fleet_control[flt_ind, 14]  # Dim3 of age transition matrix
  flt_q_ind <- fleet_control[flt_ind, 15]          # Index of survey q
  est_index_q <- fleet_control[flt_ind, 16]        # Estimate analytical q?
  index_varying_q <- fleet_control[flt_ind, 17]    # Time varying q type
  est_sigma_index <- fleet_control[flt_ind, 18]    # Whether to estimate standard deviation of survey time series
  est_sigma_fsh <- fleet_control[flt_ind, 19]      # Whether to estimate standard deviation of fishery time series


  # * 5.3. CATCHABILITY -----
  for(flt in 1:n_flt) {
    if(flt_type[flt] == 2){
      for(yr in 1:nyrs_hind) {
        index_q[flt, yr] <- exp(index_ln_q[flt] + index_q_dev[flt, yr])  # Exponentiate

        # Q as a function of environmental index
        if(est_index_q[flt] == 5) {
          beta_q_tmp <- index_q_beta[flt,]
          env_q_tmp <- env_index[yr,]
          index_q_mult <- sum(env_q_tmp * beta_q_tmp)
          index_q[flt, yr] <- exp(index_ln_q[flt] + index_q_mult)
        }

        # QAR1 deviates fit to environmental index (sensu Rogers et al 2024; 10.1093/icesjms/fsae005)
        if(est_index_q[flt] == 6) {
          index_q[flt, yr] <- exp(index_ln_q[flt] + index_q_beta[flt, 1] * index_q_dev[flt, yr])
        }
      }
    }
  }


  # * 5.4. SELECTIVITY -----
  # ** 5.4.1 EMPIRICAL SELECTIVITY ----
  sel <- array(0, dim = c(n_flt, max_nsex, max_nages, nyrs))                           # Estimated selectivity at age

  for(sel_ind in 1:nrow(emp_sel_obs)) {

    # Fishery
    flt <- emp_sel_ctl[sel_ind, 1]            # Temporary index
    sp <- emp_sel_ctl[sel_ind, 2]             # Temporary index of species
    sex <- emp_sel_ctl[sel_ind, 3]            # Temporary index for sex
    yr <- emp_sel_ctl[sel_ind, 4] - styr + 1  # Temporary index for years of data

    # Switch to make sure it doesn't fill out selectivity
    if(flt_sel_type[flt] == 0) {

      # 1 sex model
      if (nsex[sp] == 1) {
        sexes <- c(1)  # Only one sex
      }

      # 2 sex model and emp_sel is for both sex
      if(nsex[sp] == 2 && sex == 0) {
        sexes <- c(1, 2)
      }

      # 2 sex model and emp_sel is for 1 sex
      if(nsex[sp] == 2 && sex > 0) {
        sexes <- sex
      }

      if(yr < nyrs_hind) {
        for(s in sexes) {
          for(age in 1:nages[sp]) {
            if(!is.na(emp_sel_obs[sel_ind, age])) {
              sel[flt, s, age, yr] <- emp_sel_obs[sel_ind, age]
            }
          }
        }
      }
    }
  }



  # ** 5.4.2. ESTIMATED SELECTIVITY ----
  avg_sel <- matrix(0, n_flt, max_nsex)                                             # Average selectivity
  non_par_sel <- array(0, dim = c(n_flt, max_nsex, max_nages))                        # Temporary saved selectivity at age for estimated bits

  for(flt in 1:n_flt) {
    # Temporary indices
    sp <- flt_spp[flt]             # Temporary index of species
    sel_type <- flt_sel_type[flt]
    nselages <- flt_nselages[flt]

    switch(sel_type,
           # Case 1: Logistic selectivity
           {
             for(age in 1:nages[sp]) {
               for(yr in 1:nyrs_hind) {
                 for(sex in 1:nsex[sp]) {
                   # Random walk and block
                   sel[flt, sex, age, yr] <- 1 / (1 + exp(-exp(ln_sel_slp[1, flt, sex] + ln_sel_slp_dev[1, flt, sex, yr]) *
                                                            (age - (sel_inf[1, flt, sex] + sel_inf_dev[1, flt, sex, yr]))))
                 }
               }
             }
           },

           # Case 2: Non-parametric selectivity
           {
             for(sex in 1:nsex[sp]) {
               for(age in 1:nselages) {
                 non_par_sel[flt, sex, age] <- sel_coff[flt, sex, age]
                 avg_sel[flt, sex] <- avg_sel[flt, sex] + exp(sel_coff[flt, sex, age])
               }
               # Average selectivity up to nselages
               avg_sel[flt, sex] <- log(avg_sel[flt, sex] / nselages)

               # Plus group selectivity
               for(age in (nselages + 1):nages[sp]) {
                 non_par_sel[flt, sex, age] <- non_par_sel[flt, sex, nselages]
               }

               # Average selectivity across all ages
               avgsel_tmp <- 0
               for(age in 1:nages[sp]) {
                 avgsel_tmp <- avgsel_tmp + exp(non_par_sel[flt, sex, age])
               }
               avgsel_tmp <- log(avgsel_tmp / nages[sp])

               # Standardize selectivity
               for(age in 1:nages[sp]) {
                 non_par_sel[flt, sex, age] <- non_par_sel[flt, sex, age] - avgsel_tmp
                 non_par_sel[flt, sex, age] <- exp(non_par_sel[flt, sex, age])
               }
             }

             # Move to rest of years
             for(age in 1:nages[sp]) {
               for(sex in 1:nsex[sp]) {
                 for(yr in 1:nyrs_hind) {
                   sel[flt, sex, age, yr] <- non_par_sel[flt, sex, age]
                 }
               }
             }
           },

           # Case 3: Double logistic
           {
             for(age in 1:nages[sp]) {
               for(yr in 1:nyrs_hind) {
                 for(sex in 1:nsex[sp]) {
                   # Random walk and block
                   upper_slope <- 1 / (1 + exp(-exp(ln_sel_slp[1, flt, sex] + ln_sel_slp_dev[1, flt, sex, yr]) *
                                                 (age - (sel_inf[1, flt, sex] + sel_inf_dev[1, flt, sex, yr]))))
                   down_slope <- 1 - 1 / (1 + exp(-exp(ln_sel_slp[2, flt, sex] + ln_sel_slp_dev[2, flt, sex, yr]) *
                                                    (age - (sel_inf[2, flt, sex] + sel_inf_dev[2, flt, sex, yr]))))
                   sel[flt, sex, age, yr] <- upper_slope * down_slope
                 }
               }
             }
           },

           # Case 4: Descending logistic
           {
             for(age in 1:nages[sp]) {
               for(yr in 1:nyrs_hind) {
                 for(sex in 1:nsex[sp]) {
                   # Random walk and block
                   sel[flt, sex, age, yr] <- 1 - 1 / (1 + exp(-exp(ln_sel_slp[2, flt, sex] + ln_sel_slp_dev[2, flt, sex, yr]) *
                                                                (age - (sel_inf[2, flt, sex] + sel_inf_dev[2, flt, sex, yr]))))
                 }
               }
             }
           },

           # Case 5: Non-parametric selectivity (Hake version)
           {
             # Initialize selectivity to zero
             sel[flt, , 1:nages[sp], 1:nyrs_hind] <- 0

             # Sum coefficients
             for(yr in 1:nyrs_hind) {
               for(sex in 1:nsex[sp]) {
                 for(age in flt_sel_age[flt]:nselages) {
                   for(age_tmp in flt_sel_age[flt]:age) {
                     sel[flt, sex, age, yr] <- sel[flt, sex, age, yr] +
                       sel_coff[flt, sex, age_tmp] + sel_coff_dev[flt, sex, age_tmp, yr]
                   }
                 }
               }
             }

             # Normalize by maximum
             for(yr in 1:nyrs_hind) {
               for(sex in 1:nsex[sp]) {
                 max_sel <- 0
                 for(age in flt_sel_age[flt]:nselages) {
                   max_sel <- max2(max_sel, sel[flt, sex, age, yr])
                 }
                 if (max_sel > 0) {
                   sel[flt, sex, flt_sel_age[flt]:nselages, yr] <-
                     exp(sel[flt, sex, flt_sel_age[flt]:nselages, yr] - max_sel)
                 }

                 # Fill in rest of ages
                 sel[flt, sex, (nselages + 1):nages[sp], yr] <- sel[flt, sex, nselages, yr]
               }
             }
           }
    ) # End switch

    # Normalize and account for unselected ages
    if(sel_type > 0) {

      # Zero out ages not selected
      if(!is.na(flt_sel_age[flt])){
        for(age in 1:nages[sp]) {
          if(age < flt_sel_age[flt]) {
            sel[flt, , age, ] <- 0
          }
        }
      }

      # Normalize selectivity
      if(!is.na(flt_sel_maxage[flt])){
        # - Normalize by specific age if specified
        if(flt_sel_maxage[flt] >= 0) {
          for(yr in 1:nyrs_hind) {
            for(sex in 1:nsex[sp]) {
              max_sel <- sel[flt, sex, flt_sel_maxage[flt], yr]
              if (max_sel > 0) {
                sel[flt, sex, , yr] <- sel[flt, sex, , yr] / max_sel
              }
            }
          }
        }

        # - Normalize by max
        if(sel_type < 5 & flt_sel_maxage[flt] < 0) {
          for(yr in 1:nyrs_hind) {
            for(sex in 1:nsex[sp]) {
              max_sel <- 0
              for(age in 1:nages[sp]) {
                max_sel <- max2(max_sel, sel[flt, sex, age, yr])
              }

              sel[flt, sex, , yr] <- sel[flt, sex, , yr] / max_sel
            }
          }
        }
      }
    }

    # Project forward using final year's selectivity
    for(yr in (nyrs_hind + 1):nyrs) {
      sel[flt, , , yr] <- sel[flt, , , nyrs_hind]
    }
  }

  # * 5.5. VULNERABILITY ----
  # Transform predator-prey preference parameters
  # Adopted from https:#github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp (Trijoulet et al 2020)
  # Suitability for other food = 1-sum(predator-prey preference)
  # Criteria for predator-prey preference:
  # 1. predator-prey preference > 0 (hence logs)
  # 2. sum(predator-prey preference) + vuln_other = 1
  # 3. 0 <= sum(predator-prey preference) <= 1 (hence logit transformation)
  #
  # if(suitMode > 0){
  #   sum_phi <- rep(0, nspp)
  #
  #   for(rsp in 1:nspp) {  # Pred loop
  #     for(ksp in 1:nspp) { # Prey loop
  #       sum_phi[rsp] <- sum_phi[rsp] + exp(log_phi[rsp,ksp])
  #     }
  #
  #     for(ksp in 1:nspp) { # Prey loop
  #       vulnerability[rsp,ksp] <- exp(log_phi[rsp,ksp])/(1 + sum_phi[rsp]) # multinomial logistic transformation
  #     }
  #
  #     # vulnerability-other=1-sum-vulnerability but transform so sum_vuln+vuln_other=1
  #     vulnerability_other[rsp] <- sum(vulnerability[rsp,])
  #   }
  # }

  #TEST
  # sum(sel[,1,,] != mod_objects$quantities$sel[,1,,])

  # ------------------------------------------------------------------------- #
  # 6. POPULATION DYNAMICS EQUATIONS -----
  # ------------------------------------------------------------------------- #
  # NOTE: Remember indexing starts at 1
  # Start iterations for multi-species convergence
  for (iter in 1:niter){

    # * 6.1. FISHING MORTALITY ----
    F_results <- calculate_fishing_mortality(n_flt, nspp, nages, max_nages, nsex, max_nsex, nyrs, nyrs_hind,
                                             flt_spp, flt_type, sel, ln_mean_F, F_dev, Ftarget, Fmult,
                                             Flimit, QnormHCR, forecast, proj_F_prop, HCR)

    # * 6.2. TOTAL MORTALITY-AT-AGE ----
    for (sp in 1:nspp) {
      for (age in 1:nages[sp]) {
        for (yr in 1:nyrs) {
          for(sex in 1:nsex[sp]) {
            M_at_age[sp, sex, age, yr] <- M1_at_age[sp, sex, age] + M2_at_age[sp, sex, age, yr]
            Z_at_age[sp, sex, age, yr] <- M1_at_age[sp, sex, age] + F_results$F_spp_age[sp, sex, age, yr] + M2_at_age[sp, sex, age, yr]
            # S_at_age[sp, sex, age, yr] <- exp(-Z_at_age[sp, sex, age, yr])
          }
        }
      }
    }

    #TEST
    # sum(M1_at_age[,1,] - mod_objects$quantities$M1_at_age[,1,], na.rm = TRUE)
    # sum(M_at_age[,1,,] - mod_objects$quantities$M_at_age[,1,,], na.rm = TRUE)
    # sum(F_results$F_spp_age[,1,,] - mod_objects$quantities$F_spp_age[,1,,], na.rm = TRUE)

    # * 6.3. SPR BASED REFERENCE POINTS -----
    spr_results <- calculate_spr_reference_points(nspp, nages, max_nages, nyrs, nyrs_hind,
                                                  initMode, R_sexr, M1_at_age, M2_at_age,
                                                  F_results$Flimit_age_spp, F_results$Ftarget_age_spp,
                                                  wt, ssb_wt_index, pmature_sexr, spawn_month,
                                                  Finit)


    # * 6.4. STOCK-RECRUIT PARAMETERS ----
    # -- For beverton-holt, steepness and R0 are derived from SPR0
    srr_results <- calculate_sr_parameters(nspp, srr_fun, rec_pars, beta_rec_pars,
                                           env_index_srr, spr_results$SPR0, spr_results$SPRFinit)


    # * 6.5. INITIAL NUMBERS AT AGE, BIOMASS, AND SSB (YEAR 1) ----
    biomass_at_age <- array(0, dim = c(nspp, max_nsex, max_nages, nyrs))               # Estimated biomass-at-age (kg)
    biomass <- matrix(0, nrow = nspp, ncol = nyrs)
    ssb_at_age <- array(0, dim = c(nspp, max_nages, nyrs))                             # Spawning biomass at age (kg)
    ssb <- matrix(0, nrow = nspp, ncol = nyrs)

    for (sp in 1:nspp) {

      if(initMode != 3){
        Finit[sp] = 0
      }

      for (age in 1:nages[sp]) {
        for (sex in 1:nsex[sp]) {

          # Handle different estimation dynamics
          # switch(estDynamics[sp] + 1, { # R is 1-based indexing
          # Case 0: Estimated
          # Estimate as free parameters
          if (initMode == 0) {
            R[sp, 1] <- exp(init_dev[sp, 1])
            N_at_age[sp, 1, age, 1] <- exp(init_dev[sp, age]) * R_sexr[sp]
            if(sex == 2){
              N_at_age[sp, 2, age, 1] <- exp(init_dev[sp, age]) * (1 - R_sexr[sp])
            }
          }

          # Equilibrium or non-equilibrium estimated
          if (initMode > 0) {
            # 6.5.1. Amin (recruitment)
            if (age == 1) {
              R[sp, 1] <- srr_results$R_init[sp] * exp(rec_dev[sp, 1])
              N_at_age[sp, 1, 1, 1] <- R[sp, 1] * R_sexr[sp]
              if(sex == 2){
                N_at_age[sp, 2, 1, 1] <- R[sp, 1] * (1 - R_sexr[sp])
              }
            }

            # Sum M1 until age - 1
            mort_sum <- 0
            for (age_tmp in 1:(age-1)) {
              mort_sum <- mort_sum + M1_at_age[sp, sex, age_tmp] + Finit[sp]
            }

            # 6.5.2. Age Amin+1:Amax-1 (initial abundance)
            if (age > 1 && age < nages[sp]) {
              if (sex == 1) {
                N_at_age[sp, 1, age, 1] <- srr_results$R_init[sp] * exp(-mort_sum + init_dev[sp, age - 1]) * R_sexr[sp]
              }
              if (sex == 2) {
                N_at_age[sp, 2, age, 1] <- srr_results$R_init[sp] * exp(-mort_sum + init_dev[sp, age - 1]) * (1 - R_sexr[sp])
              }
            }

            # 6.5.3. Amax
            if (age == nages[sp]) {
              if (sex == 1) {
                N_at_age[sp, 1, age, 1] <- srr_results$R_init[sp] * exp(-mort_sum + init_dev[sp, age - 1]) /
                  (1 - exp(-M1_at_age[sp, sex, nages[sp]])) * R_sexr[sp]
              }
              if (sex == 2) {
                N_at_age[sp, 2, age, 1] <- srr_results$R_init[sp] * exp(-mort_sum + init_dev[sp, age - 1]) /
                  (1 - exp(-M1_at_age[sp, sex, nages[sp]])) * (1 - R_sexr[sp])
              }
            }
          }
          # }, {
          #   # Case 1: Fixed numbers-at-age - fixed scalar
          #   N_at_age[sp, sex, age, 1] <- pop_scalar[sp, 1] * NByageFixed[sp, sex, age, 1]
          # }, {
          #   # Case 2: Fixed numbers-at-age age-independent scalar
          #   N_at_age[sp, sex, age, 1] <- pop_scalar[sp, 1] * NByageFixed[sp, sex, age, 1]
          # }, {
          #   # Case 3: Fixed numbers-at-age age-dependent scalar
          #   N_at_age[sp, sex, age, 1] <- pop_scalar[sp, age] * NByageFixed[sp, sex, age, 1]
          # }, {
          #   stop("Invalid 'estDynamics'")
          # })

          # 6.5.3. Estimate total biomass in year 1
          biomass_at_age[sp, sex, age, 1] <- N_at_age[sp, sex, age, 1] *
            wt[pop_wt_index[sp], sex, age, 1]  # Accessing wt as a data frame
          biomass[sp, 1] <- biomass[sp, 1] + biomass_at_age[sp, sex, age, 1]
        }

        # 6.5.4. Estimated initial female SSB
        ssb_at_age[sp, age, 1] <- N_at_age[sp, 1, age, 1] *
          exp(-Z_at_age[sp, 1, age, 1] * spawn_month[sp] / 12) *
          wt[ssb_wt_index[sp], 1, age, 1] *
          pmature_sexr[sp, age]
        ssb[sp, 1] <- ssb[sp, 1] + ssb_at_age[sp, age, 1]
      }
    }



    # * 6.6. HINDCAST NUMBERS AT AGE, BIOMASS, AND SSB (YEAR 2+) ----
    zero_N_pen <- numeric(nspp)  # Initialize zero_N_pen for each species

    for (sp in 1:nspp) {
      for (yr in 2:nyrs_hind) {

        # Switch for srr (for MSEs if using Ianelli SRR method)
        srr_switch <- srr_fun
        if (yr >= nyrs_srrmean) {
          srr_switch <- srr_pred_fun
        }

        # -- 6.6.1. Recruitment
        R[sp, yr] <- switch(srr_switch + 1,
                            # Random about mean (e.g. Alaska)
                            srr_results$R0[sp] * exp(rec_dev[sp, yr]),

                            # Random about mean with environmental effects
                            {
                              beta_rec_tmp <- beta_rec_pars[sp,]
                              env_rec_tmp <- env_index_srr[yr,]
                              srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
                              srr_results$R0[sp] * exp(rec_dev[sp, yr] + srr_mult)
                            },

                            # Beverton-Holt
                            exp(rec_pars[sp, 1]) * ssb[sp, yr - minage[sp]] * exp(rec_dev[sp, yr]) /
                              (1.0 + exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]]),

                            # Beverton-Holt with environmental impacts on alpha
                            {
                              beta_rec_tmp <- beta_rec_pars[sp,]
                              env_rec_tmp <- env_index_srr[yr,]
                              srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
                              srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
                              srr_alpha * ssb[sp, yr - minage[sp]] * exp(rec_dev[sp, yr]) /
                                (1.0 + exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]])
                            },

                            # Ricker
                            exp(rec_pars[sp, 1]) * ssb[sp, yr - minage[sp]] *
                              exp(-exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]] / 1000000.0) * exp(rec_dev[sp, yr]),

                            # Ricker with environmental impacts on alpha
                            {
                              beta_rec_tmp <- beta_rec_pars[sp,]
                              env_rec_tmp <- env_index_srr[yr,]
                              srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
                              srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
                              srr_alpha * ssb[sp, yr - minage[sp]] *
                                exp(-exp(rec_pars[sp, 2]) * ssb[sp, yr - minage[sp]] / 1000000.0) * exp(rec_dev[sp, yr])
                            },

                            stop("Invalid 'srr_fun'")
        )

        N_at_age[sp, 1, 1, yr] <- R[sp, yr] * R_sexr[sp]

        if(nsex[sp] > 1){
          N_at_age[sp, 2, 1, yr] <- R[sp, yr] * (1.0 - R_sexr[sp])
        }

        # -- 6.6.2. Ages beyond recruitment
        for (age in 1:nages[sp]) {
          for (sex in 1:nsex[sp]) {

            # if (estDynamics[sp] == 0) { # Estimated numbers-at-age
            # -- Where Amin < age < Amax
            if (age < nages[sp]) {
              N_at_age[sp, sex, age + 1, yr] <- N_at_age[sp, sex, age, yr - 1] *
                exp(-Z_at_age[sp, sex, age, yr - 1])
            }

            # -- Plus group where age = Amax
            if (age == nages[sp]) {
              N_at_age[sp, sex, age, yr] <- N_at_age[sp, sex, age - 1, yr - 1] *
                exp(-Z_at_age[sp, sex, age - 1, yr - 1]) +
                N_at_age[sp, sex, age, yr - 1] * exp(-Z_at_age[sp, sex, age, yr - 1])
            }
            # } else if (estDynamics[sp] == 1) { # Fixed numbers-at-age - fixed scalar
            #   N_at_age[sp, sex, age, yr] <- pop_scalar[sp, 1] * NByageFixed[sp, sex, age, yr]
            # } else if (estDynamics[sp] == 2) { # Fixed numbers-at-age age-independent scalar
            #   N_at_age[sp, sex, age, yr] <- pop_scalar[sp, 1] * NByageFixed[sp, sex, age, yr]
            # } else if (estDynamics[sp] == 3) { # Fixed numbers-at-age age-dependent scalar
            #   N_at_age[sp, sex, age, yr] <- pop_scalar[sp, age] * NByageFixed[sp, sex, age, yr]
            # } else {
            #   stop("Invalid 'estDynamics'")
            # }

            # Ensure positive values and calculate penalty
            # pos_tmp <- posfun(N_at_age[sp, sex, age, yr], 0.001)
            # N_at_age[sp, sex, age, yr] <- pos_tmp$ans
            # zero_N_pen[sp] <- zero_N_pen[sp] + pos_tmp$penalty

            # -- 6.6.3. Estimate total biomass
            biomass[sp, yr] <- biomass[sp, yr] +
              N_at_age[sp, sex, age, yr] * wt[pop_wt_index[sp], sex, age, yr]
          }

          # -- 6.6.4. Estimated female ssb
          ssb[sp, yr] <- ssb[sp, yr] +
            N_at_age[sp, 1, age, yr] * exp(-Z_at_age[sp, 1, age, yr] * spawn_month[sp] / 12.0) *
            wt[ssb_wt_index[sp], 1, age, yr] * pmature_sexr[sp, age]
        }
      }
    }

    # ** Calculate mean recruitment from hindcast ----
    avg_R <- rep(0, nspp)
    for(sp in 1:nspp) {
      for(yr in 1:nyrs_srrmean) {
        avg_R[sp] <- avg_R[sp] + R[sp, yr]/nyrs_srrmean # Update mean rec
      }
    }


    # * 6.7. DEPLETION REFERENCE POINTS ----
    rps_results <- calculate_depletion_reference_points(nspp, nyrs, nsex, nages, max_nsex, max_nages, nyrs_hind, proj_mean_rec,
                                                        srr_pred_fun, R, avg_R, R_results$R0, rec_dev,
                                                        beta_rec_pars, env_index_srr, N_at_age, M_at_age,
                                                        F_results$Ftarget_age_spp, pmature, spawn_month,
                                                        wt, ssb_wt_index, R_sexr, minage,
                                                        MSSB0, MSB0, msmMode)


    # * 6.8-6.9. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg) ----
    # Includes Harvest Control Rules
    for (sp in 1:nspp) {
      for (yr in (nyrs_hind+1):nyrs) {

        # ** 6.8. HARVEST CONTROL RULES FOR PROJECTION (i.e. SB0 and dynamic SB0) ----
        # -- Equilibrium Harvest Control Rules
        if(DynamicHCR == 0) {
          # Using switch-like logic with if-else in R
          if(HCR == 0) { # No fishing
            F_results$proj_F[sp, yr] <- 0.0
          }
          if(HCR == 1) { # CMSY
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
          }
          if(HCR == 2) { # Constant F
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
          }
          if(HCR == 3) { # Constant F to achieve X% of SB0
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
          }
          if(HCR == 4) { # Constant Fspr with multiplier
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr] * Fmult[sp]
          }
          if(HCR == 5) { # NPFMC Tier 3 HCR
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
            if(ssb[sp, yr-1] < rps_results$SBF[sp, nyrs-1]) {
              F_results$proj_F[sp, yr] <- Ftarget[sp] * ((ssb[sp, yr-1]/rps_results$SBF[sp, nyrs-1] - Alpha[sp])/(1-Alpha[sp])) # Used Fabc of FtargetSPR%
            }
            if((ssb[sp, yr-1] < rps_results$SB0[sp, nyrs-1] * Plimit[sp]) | (ssb[sp, yr-1] / rps_results$SBF[sp, nyrs-1] < Alpha[sp])) { # If overfished
              F_results$proj_F[sp, yr] <- 0.0
            }
          }
          if(HCR == 6) { # PFMC Category 1 HCR
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
            if(ssb[sp, yr-1] < rps_results$SB0[sp, nyrs-1] * Ptarget[sp]) {
              F_results$proj_F[sp, yr] <- (Flimit[sp] + QnormHCR[sp]) * (rps_results$SB0[sp, nyrs-1] * Ptarget[sp] * (ssb[sp, yr-1] - rps_results$SB0[sp, nyrs-1] * Plimit[sp])) /
                (ssb[sp, yr-1] * (rps_results$SB0[sp, nyrs-1] * (Ptarget[sp] - Plimit[sp])))
            }
            if(ssb[sp, yr-1] < rps_results$SB0[sp, nyrs-1] * Plimit[sp]) { # If overfished
              F_results$proj_F[sp, yr] <- 0.0
            }
          }
          if(HCR == 7) { # SESSF Tier 1 HCR
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
            if(ssb[sp, yr-1] < rps_results$SB0[sp, nyrs-1] * Ptarget[sp]) {
              F_results$proj_F[sp, yr] <- Ftarget[sp] * ((ssb[sp, yr-1]/(rps_results$SB0[sp, nyrs-1] * Plimit[sp]))-1) # Used Fabc of FtargetSPR%
            }
            if(ssb[sp, yr-1] < rps_results$SB0[sp, nyrs-1] * Plimit[sp]) { # If overfished
              F_results$proj_F[sp, yr] <- 0.0
            }
          }
        }

        # Dynamic Harvest Control Rules
        if(DynamicHCR == 1) {
          # Similar structure as above but with Dynamic reference points
          if(HCR == 0) { # No fishing
            F_results$proj_F[sp, yr] <- 0.0
          }
          if(HCR == 1) { # CMSY
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
          }
          if(HCR == 2) { # Constant F
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
          }
          if(HCR == 3) { # Constant F to achieve X% of SB0
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
          }
          if(HCR == 4) { # Constant Fspr with multiplier
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr] * Fmult[sp]
          }
          if(HCR == 5) { # NPFMC Tier 3 HCR
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
            if(ssb[sp, yr-1] < rps_results$DynamicSBF[sp, yr-1]) {
              F_results$proj_F[sp, yr] <- Ftarget[sp] * ((ssb[sp, yr-1]/rps_results$DynamicSBF[sp, yr-1] - Alpha[sp])/(1-Alpha[sp]))
            }
            if((ssb[sp, yr-1] < rps_results$DynamicSB0[sp, yr-1] * Plimit[sp]) | (ssb[sp, yr-1] / rps_results$DynamicSBF[sp, yr-1] < Alpha[sp])) {
              F_results$proj_F[sp, yr] <- 0.0
            }
          }
          if(HCR == 6) { # PFMC Category 1 HCR
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
            if(ssb[sp, yr-1] < rps_results$DynamicSB0[sp, yr-1] * Ptarget[sp]) {
              F_results$proj_F[sp, yr] <- (Flimit[sp] + QnormHCR[sp]) * (rps_results$DynamicSB0[sp, yr-1] * Ptarget[sp] * (ssb[sp, yr-1] - rps_results$DynamicSB0[sp, yr-1] * Plimit[sp])) /
                (ssb[sp, yr-1] * (rps_results$DynamicSB0[sp, yr-1] * (Ptarget[sp] - Plimit[sp])))
            }
            if(ssb[sp, yr-1] < rps_results$DynamicSB0[sp, yr-1] * Plimit[sp]) {
              F_results$proj_F[sp, yr] <- 0.0
            }
          }
          if(HCR == 7) { # SESSF Tier 1 HCR
            F_results$proj_F[sp, yr] <- F_results$proj_F[sp, yr]
            if(ssb[sp, yr-1] < rps_results$DynamicSB0[sp, yr-1] * Ptarget[sp]) {
              F_results$proj_F[sp, yr] <- Ftarget[sp] * ((ssb[sp, yr-1]/(rps_results$DynamicSB0[sp, yr-1] * Plimit[sp]))-1)
            }
            if(ssb[sp, yr-1] < rps_results$DynamicSB0[sp, yr-1] * Plimit[sp]) {
              F_results$proj_F[sp, yr] <- 0.0
            }
          }
        }

        # Set F to 0 if not forecast
        if(forecast[sp] == 0) {
          F_results$proj_F[sp, yr] <- 0.0
        }

        # Adjust F*selex
        # -- 6.8.3. Update F for the projection (account for selectivity and fleets)
        F_results$F_spp_age[sp,,, yr] <- 0.0

        # -- Multiply F from HCR by selectivity and fleet proportion
        F_results$F_spp[sp, yr] <- F_results$proj_F[sp, yr]
        for(flt in 1:n_flt) {
          if(sp == flt_spp[flt]) {
            F_results$F_flt[sp, yr] <- proj_F_prop[flt] * F_results$proj_F[sp, yr]
            for(age in 1:nages[sp]) {
              for(sex in 1:nsex[sp]) {
                F_results$F_flt_age[flt, sex, age, yr] <- sel[flt, sex, age, nyrs_hind] * proj_F_prop[flt] * F_results$proj_F[sp, yr] # FIXME using last year of selectivity
                if(flt_type[flt] == 1) {
                  F_results$F_spp_age[sp, sex, age, yr] <- F_results$F_spp_age[sp, sex, age, yr] + F_results$F_flt_age[flt, sex, age, yr]
                }
              }
            }
          }
        }

        # -- 6.8.4. Update mortality for forecast
        for(age in 1:nages[sp]) {
          for(sex in 1:nsex[sp]) {
            M_at_age[sp, sex, age, yr] <- M1_at_age[sp, sex, age] + M2_at_age[sp, sex, age, yr]
            Z_at_age[sp, sex, age, yr] <- M1_at_age[sp, sex, age] + F_results$F_spp_age[sp, sex, age, yr] + M2_at_age[sp, sex, age, yr]
          }
        }


        # ** 6.9. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg) ----
        # -- 6.9.1. Forecasted recruitment
        # - Option 1: Use mean rec
        if((proj_mean_rec == 1) & (srr_pred_fun > 1)) {
          R[sp, yr] <- exp(log(avg_R[sp]) + rec_dev[sp, yr]) # Projections use mean R given bias in R0
        }

        # - Mean rec and environment
        if((proj_mean_rec == 1) & (srr_pred_fun < 2)) {
          beta_rec_tmp <- beta_rec_pars[sp,]
          env_rec_tmp <- env_index_srr[yr,]
          srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
          R[sp, yr] <- exp(log(avg_R[sp]) + rec_dev[sp, yr]) * exp(srr_mult)
        }

        # - Option 2: Use SRR and rec devs
        if(proj_mean_rec == 0) {
          if(srr_pred_fun == 0) { # Random about mean (e.g. Alaska)
            R[sp, yr] <- R0[sp] * exp(rec_dev[sp, yr])
          }
          if(srr_pred_fun == 1) { # Random about mean with environmental effects
            beta_rec_tmp <- beta_rec_pars[sp,]
            env_rec_tmp <- env_index_srr[yr,]
            srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
            R[sp, yr] <- R0[sp] * exp(rec_dev[sp, yr] + srr_mult)
          }
          if(srr_pred_fun == 2) { # Beverton-Holt
            R[sp, yr] <- exp(rec_pars[sp, 1]) * ssb[sp, yr-minage[sp]] * exp(rec_dev[sp, yr]) /
              (1 + exp(rec_pars[sp, 2]) * ssb[sp, yr-minage[sp]])
          }
          if(srr_pred_fun == 3) { # Beverton-Holt with environmental impacts on alpha
            beta_rec_tmp <- beta_rec_pars[sp,]
            env_rec_tmp <- env_index_srr[yr,]
            srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
            srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
            R[sp, yr] <- srr_alpha * ssb[sp, yr-minage[sp]] * exp(rec_dev[sp, yr]) /
              (1 + exp(rec_pars[sp, 2]) * ssb[sp, yr-minage[sp]])
          }
          if(srr_pred_fun == 4) { # Ricker
            R[sp, yr] <- exp(rec_pars[sp, 1]) * ssb[sp, yr-minage[sp]] *
              exp(-exp(rec_pars[sp, 2]) * ssb[sp, yr-minage[sp]]/1000000.0) * exp(rec_dev[sp, yr])
          }
          if(srr_pred_fun == 5) { # Ricker with environmental impacts on alpha
            beta_rec_tmp <- beta_rec_pars[sp,]
            env_rec_tmp <- env_index_srr[yr,]
            srr_mult <- sum(env_rec_tmp * beta_rec_tmp)
            srr_alpha <- exp(rec_pars[sp, 1] + srr_mult)
            R[sp, yr] <- srr_alpha * ssb[sp, yr-minage[sp]] *
              exp(-exp(rec_pars[sp, 2]) * ssb[sp, yr-minage[sp]]/1000000.0) * exp(rec_dev[sp, yr])
          }
        }

        N_at_age[sp, 1, 1, yr] <- R[sp, yr] * R_sexr[sp]
        if(nsex[sp] > 1){
          N_at_age[sp, 2, 1, yr] <- R[sp, yr] * (1 - R_sexr[sp])
        }

        # -- Ages > recruitment
        for(age in 1:nages[sp]) {
          for(sex in 1:nsex[sp]) {
            # if(estDynamics[sp] == 0) { # Estimated numbers-at-age
            # -- Where age < plus group
            if(age < nages[sp]) {
              N_at_age[sp, sex, age + 1, yr] <- N_at_age[sp, sex, age, yr - 1] * exp(-Z_at_age[sp, sex, age, yr-1])
            }
            # -- Plus group
            if(age == nages[sp]) {
              N_at_age[sp, sex, age, yr] <- N_at_age[sp, sex, age - 1, yr - 1] * exp(-Z_at_age[sp, sex, age-1, yr-1]) +
                N_at_age[sp, sex, age, yr - 1] * exp(-Z_at_age[sp, sex, age, yr-1])
            }
            # } else if(estDynamics[sp] == 1) { # Fixed numbers-at-age - fixed scalar
            #   N_at_age[sp, sex, age, yr] <- pop_scalar[sp, 1] * NByageFixed[sp, sex, age, yr]
            # } else if(estDynamics[sp] == 2) { # Fixed numbers-at-age age-independent scalar
            #   N_at_age[sp, sex, age, yr] <- pop_scalar[sp, 1] * NByageFixed[sp, sex, age, yr]
            # } else if(estDynamics[sp] == 3) { # Fixed numbers-at-age age-dependent scalar
            #   N_at_age[sp, sex, age, yr] <- pop_scalar[sp, age] * NByageFixed[sp, sex, age, yr]
            # } else {
            #   stop("Invalid 'estDynamics'")
            # }

            # Constraint to reduce population collapse
            # pos_tmp <- posfun(N_at_age[sp, sex, age, yr], 0.001)
            # N_at_age[sp, sex, age, yr] <- pos_tmp$ans
            # zero_N_pen[sp] <- zero_N_pen[sp] + pos_tmp$penalty

            # -- 6.9.4. FORECAST ssb BY AGE
            biomass_at_age[sp, sex, age, yr] <- N_at_age[sp, sex, age, yr] * wt[pop_wt_index[sp], sex, age, nyrs_hind] # 6.5.
            biomass[sp, yr] <- biomass[sp, yr] + biomass_at_age[sp, sex, age, yr]
          } # End sex loop

          # -- 6.9.5. FORECAST ssb (SUM ACROSS AGES)
          ssb_at_age[sp, age, yr] <- N_at_age[sp, 1, age, yr] *
            exp(-Z_at_age[sp, 1, age, yr] * (spawn_month[sp]/12.0)) *
            wt[ssb_wt_index[sp], 1, age, nyrs_hind] * pmature_sexr[sp, age] # 6.6.
          ssb[sp, yr] <- ssb[sp, yr] + ssb_at_age[sp, age, yr]
        }
      }
    }


    # * 6.10. ESTIMATE AVERAGE NUMBERS AT AGE ----
    avgN_at_age <- switch(
      avgnMode+1,
      # Case 0: MSVPA approach
      {
        N_at_age *
          (1 - exp(-Z_at_age)) / Z_at_age
      },
      # Case 1: Kinzey and Punt (2009) approximation
      {
        N_at_age *
          exp(-Z_at_age / 2)
      },
      # Case 2: Van Kirk et al (2010) approximation
      {
        N_at_age
      },
      # Default case
      {
        stop("Invalid 'avgnMode'")
      }
    )
    #FIXME -- doesnt work
    # avgN_at_age[is.nan(avgN_at_age)] <- 0 # Divide by 0 for ages not represented

    #TEST
    # sum(N_at_age[1,1,,1:39] - mod_objects$quantities$N_at_age[1,1,,1:39])
    # sum(avgN_at_age[,1,,] - mod_objects$quantities$avgN_at_age[,1,,])


    # ------------------------------------------------------------------------- #
    # 7. RATION AND DIET ----
    # ------------------------------------------------------------------------- #
    # * 7.1. Calculate ration ----
    ration <- calculate_ration(nspp, nyrs, max_nsex, max_nages, nyrs_hind, Ceq, Qc, Tcm, Tco, Tcl,
                               CK1, CK4, CA, CB, fday, env_index, Cindex,
                               wt, pop_wt_index, Pvalue, Pyrs, nsex, nages)
    #TEST
    #sum(ration[1,1,,1:39] - mod_objects$quantities$ration[1,1,,1:39])

    # * 7.2. Reorganize stomach content ----
    diet_prop <- reorganize_stomach_content(stom_prop_obs, stom_prop_ctl, minage, nspp,
                                            nyrs, nsex, nages, max_nsex, max_nages, nyrs_hind, styr)


    #TEST
    #sum(diet_prop[1,1,,1,1,,1:39] - mod_objects$quantities$diet_prop[1,1,,,1:39])

    # * 7.3. Calculate other food stomach content ----
    other_food_diet_prop <- calculate_other_food_diet_prop(nyrs, nspp, nsex, nages, max_nsex, max_nages,
                                                           diet_prop, other_food)


    # ------------------------------------------------------------------------- #
    # 8. START PREDATION ----
    # ------------------------------------------------------------------------- #
    if (msmMode > 0) {

      # ------------------------------------------------------------------------- #
      # * 8.1. SUITABILITY EQUATIONS ----
      # ------------------------------------------------------------------------- #
      # 8.1.1. Holsman and MSVPA based suitability # FIXME - not flexible for interannual variation
      if (suitMode == 0) {
        suit_list <- calculate_MSVPA_suitability(diet_prop, avgN_at_age, wt, pop_wt_index,
                                                 other_food_diet_prop, nspp, nsex, max_nsex, nages, max_nages, nyrs,
                                                 nyrs_hind, suit_styr, suit_endyr, nyrs_suit,
                                                 msmMode)

        #TEST
        #sum(suit_list$suit_main[1,1,,1,1,,1:39] - mod_objects$quantities$suit_main[1,1,,,1:39])
        #TEST
        #sum(suit_list$suit_other[,1,,] - mod_objects$quantities$suit_other[,1,,])
      }

      # 8.1.2. GAMMA suitability
      if(suitMode %in% c(1,2)){
        suit_list <- calculate_gamma_suitability(nspp, nages, nsex, nyrs, nyrs_hind,
                                                 max_nages, max_nsex,
                                                 laa, wt, pop_wt_index, vulnerability,
                                                 vulnerability_other, gam_a, gam_b, suitMode)
      }

      # 8.1.3. Lognormal suitability
      if(suitMode %in% c(3,4)){
        suit_list <- calculate_lognormal_suitability(nspp, nages, nsex, nyrs, nyrs_hind,
                                                     max_nages, max_nsex,
                                                     laa, wt, pop_wt_index, vulnerability,
                                                     vulnerability_other, gam_a, gam_b, suitMode)
      }



      # ------------------------------------------------------------------------- #
      # * 8.2. PREDATION MORTALITY EQUATIONS ----
      # ------------------------------------------------------------------------- #
      # -- 8.2.1. MSVPA PREDATION MORTALITY
      if ((msmMode == 1) | (msmMode == 2)) {
        predation_results <- calculate_predation(nspp, nsex, nages, nyrs, nyrs_hind,
                                                 max_nages, max_nsex,
                                                 avgN_at_age, suit_list$suit_main, suit_list$suit_other, other_food,
                                                 wt, pop_wt_index, ration, msmMode)

        # Update M2 array
        M2_at_age <- predation_results$M2_at_age


        #TEST
        #sum(M2_at_age[,1,,] - mod_objects$quantities$M2_at_age[,1,,])
      }
    } # Predation loop
  } # Population dynamics loop (niter)
  # ------------------------------------------------------------------------- #
  # 8. ENDS PREDATION AND POPULATION DYNAMICS LOOPS ----
  # ------------------------------------------------------------------------- #


  # ------------------------------------------------------------------------- #
  # 9. INDEX EQUATIONS ----
  # ------------------------------------------------------------------------- #
  # * 9.1. Index of abundance/biomass ----
  index_hat <- calculate_abundance_index(index_ctl, index_n, N_at_age, Z_at_age, sel,
                                         wt, flt_units, flt_wt_index, nages, nsex,
                                         nyrs_hind, styr)

  # * 9.2. Analytical survey q following Ludwig and Martell 1994 ----
  index_q_analytical <- calculate_analytical_q(index_ctl, index_n, index_obs, index_hat,
                                               est_sigma_index, est_index_q, n_flt, nyrs_hind, styr)


  for(index in 1:n_flt) {
    # Set index_q to analytical if used
    if(est_index_q[index] == 3 & flt_type[index] == 2) {
      index_q[index, ] <- rep(index_q_analytical[index], nyrs_hind)
    }
  }


  # * 9.3. Survey Biomass - multiply by q ----
  # - Extract necessary columns for vectorized operations
  indices <- index_ctl[, 1]       # Temporary survey indices (1-based)
  flt_yrs <- index_ctl[, 3]       # Temporary index for years of data

  # - Adjust flt_yrs
  flt_yrs[flt_yrs > 0] <- flt_yrs[flt_yrs > 0] - styr + 1
  flt_yrs[flt_yrs < 0] <- -flt_yrs[flt_yrs < 0] - styr + 1

  # - Determine yr_ind based on flt_yrs
  yr_inds <- flt_yrs
  yr_inds[flt_yrs > nyrs_hind] <- nyrs_hind

  # - Ensure indices are within valid bounds
  valid_indices <- indices > 0 & indices <= ncol(index_q) & yr_inds > 0 & yr_inds <= nyrs_hind

  # - Update index_hat only for valid indices
  index_hat[valid_indices] <- index_q[cbind(indices[valid_indices], yr_inds[valid_indices])] * index_hat[valid_indices]

  # - Optionally, handle cases where indices were out of bounds
  if (any(!valid_indices)) {
    warning("Some indices were out of bounds.")
  }

  # * 9.4. Calculate analytical sigma following Ludwig and Walters 1994 ----
  ln_index_analytical_sd <- calculate_analytical_sd(index_ctl, index_obs, index_hat, n_flt, nyrs_hind, styr)

  #TEST
  #sum(abs(index_hat - mod_objects$quantities$index_hat))


  # ------------------------------------------------------------------------- #
  # 10. FISHERY EQUATIONS  ----
  # ------------------------------------------------------------------------- #
  # * 10.1. ESTIMATE CATCH ----

  catch_hat <- estimate_catch(catch_ctl, catch_n, F_results$F_flt_age, Z_at_age, N_at_age,
                              wt, sel, flt_wt_index, proj_F_prop, flt_units,
                              nsex, nages, styr, nyrs_hind)

  #TEST
  #sum(abs(catch_hat - mod_objects$quantities$catch_hat))

  # * 10.2 Exploitable biomass ----
  exploitable_biomass <- calculate_exploitable_biomass(n_flt, flt_spp, flt_type, nyrs, nyrs_hind,
                                                       nages, nsex, N_at_age, wt, sel,
                                                       flt_wt_index, proj_F_prop)

  # ------------------------------------------------------------------------- #
  # 11. COMPOSITION EQUATIONS ----
  # ------------------------------------------------------------------------- #
  comp_hat <- estimate_comp(comp_ctl, comp_n, comp_obs, F_results$F_flt_age, Z_at_age, N_at_age,
                            sel, index_q, age_error, age_trans_matrix,
                            flt_type, nages, nlengths, nsex, styr,
                            nyrs_hind, flt_age_transition_index)

  #TEST
  #sum(abs(comp_hat - mod_objects$quantities$comp_hat))

  # ------------------------------------------------------------------------- #
  # 12. ESTIMATED DIET ----
  # ------------------------------------------------------------------------- #
  #TODO

  # ------------------------------------------------------------------------- #
  # 13. DERIVED QUANTITIES ----
  # ------------------------------------------------------------------------- #

  # 13.1. Depletion
  biomass_depletion <- matrix(0, nrow = nspp, ncol = nyrs)
  ssb_depletion <- matrix(0, nrow = nspp, ncol = nyrs)

  for(sp in 1:nspp) {
    for(yr in 1:nyrs) {

      if(DynamicHCR == 0) {
        biomass_depletion[sp, yr] <- biomass[sp, yr] / rps_results$B0[sp, nyrs]  # Corrected index
        ssb_depletion[sp, yr] <- ssb[sp, yr] / rps_results$SB0[sp, nyrs]        # Corrected index
      }

      if(DynamicHCR == 1) {
        biomass_depletion[sp, yr] <- biomass[sp, yr] / rps_results$DynamicB0[sp, yr]
        ssb_depletion[sp, yr] <- ssb[sp, yr] / rps_results$DynamicSB0[sp, yr]
      }

      # Multi-species and no HCR (MSSB0 is input otherwise)
      if(HCR == 0 && msmMode > 0) {
        biomass_depletion[sp, yr] <- biomass[sp, yr] / biomass[sp, nyrs]  # Corrected index
        ssb_depletion[sp, yr] <- ssb[sp, yr] / ssb[sp, nyrs]              # Corrected index
      }
    }
  }

  RTMB::REPORT(biomass_depletion)
  RTMB::REPORT(ssb_depletion)


  # ------------------------------------------------------------------------- #
  # 14. LIKELIHOOD EQUATIONS ----
  # ------------------------------------------------------------------------- #

  # 14.0. OBJECTIVE FUNCTION
  jnll_comp = matrix(0, 19, n_flt) # matrix of negative log-likelihood components

  # -- Data likelihood components
  # Slot 0 -- Survey biomass
  # Slot 1 -- Total catch (kg)
  # Slot 2 -- Age/length composition
  # Slot 3 -- Sex ratio likelihood TODO
  # Slot 4 -- Selectivity
  # Slot 5 -- Selectivity annual deviates
  # Slot 6 -- Survey selectivity normalization
  # Slot 7 -- Survey catchability prior
  # Slot 8 -- Survey catchability annual deviates
  # -- Priors/penalties
  # Slot 9 -- Stock recruitment parameter (alpha) prior
  # Slot 10 -- Tau -- Annual recruitment deviation
  # Slot 11 -- init_dev -- Initial abundance-at-age
  # Slot 12 -- Epsilon -- Annual fishing mortality deviation
  # Slot 13 -- SPR penalities
  # Slot 14 -- N-at-age < 0 penalty
  # Slot 15 -- M_at_age prior
  # -- M2_at_age likelihood components
  # Slot 16 -- Ration likelihood
  # Slot 17 -- Ration penalties
  # Slot 18 -- Diet proportion by weight likelihood

  # * 14.1. INDEX DATA ----
  jnll_comp <- calculate_index_nll(index_obs, index_ctl, est_sigma_index, index_ln_sd,
                                   ln_index_analytical_sd, jnll_comp,
                                   index_hat, flt_type, endyr)

  # * 14.2. CATCH DATA ----
  jnll_comp <- calculate_catch_nll(catch_obs, catch_ctl, catch_hat, catch_ln_sd,
                                   est_sigma_fsh, jnll_comp,
                                   F_dev, flt_type, styr, endyr)

  # * 14.3. COMPOSITION DATA ----
  jnll_comp <- calculate_comp_nll(comp_obs, comp_hat, comp_ctl, comp_n, nages, nlengths,
                                  DM_pars, flt_type, comp_ll_type, comp_weights, jnll_comp, endyr)

  # #TEST
  # comp_obs - mod_objects$quantities$comp_obs
  # mod_objects$quantities$jnll_comp
  # sum(abs(comp_hat - mod_objects$quantities$comp_hat))

  # * 14.5. SELECTIVITY ----
  jnll_comp <- calculate_selectivity_nll(n_flt, flt_spp, flt_type, flt_sel_type, flt_varying_sel,
                                         non_par_sel, sel_curve_pen, avg_sel, sel_inf_dev,
                                         ln_sel_slp_dev, sel_dev_sd, flt_nselages, nsex, nages,
                                         nyrs_hind, jnll_comp)

  # * 14.5. CATCHABILITY ----
  jnll_comp <- calculate_catchability_nll(n_flt, nyrs_hind, est_index_q, index_varying_q,
                                          flt_type, index_ln_q, index_ln_q_prior,
                                          index_q_sd, index_q_dev, index_q_dev_sd,
                                          env_index, rho_trans, jnll_comp)

  # * 14.6. RECRUITMENT ----
  jnll_comp <- calculate_recruitment_nll(nspp, srr_est_mode, srr_pred_fun, steepness, srr_prior,
                                         srr_prior_sd, rec_pars, Bmsy_lim, initMode, nages,
                                         init_dev, rec_dev, R_sd, R, R_hat, srr_fun,
                                         srr_hat_styr, srr_hat_endyr, nyrs_hind, jnll_comp)

  # * 14.7. MORTALITY ----
  jnll_comp <- calculate_mortality_nll(nspp, nsex, nages, nyrs, M1_model, M1_use_prior, M2_use_prior,
                                       M1_at_age, M_at_age, M_prior, M_prior_sd, jnll_comp)

  jnll = sum(jnll_comp)

  #TEST
  #round(jnll_comp - mod_objects$quantities$jnll_comp,4)

  # ------------------------------------------------------------------------- #
  # 15. REPORT SECTION ----
  # ------------------------------------------------------------------------- #
  RTMB::REPORT(ssb)
  RTMB::REPORT(biomass)
  RTMB::REPORT(R)
  RTMB::REPORT(jnll_comp)
  RTMB::REPORT(comp_hat)
  RTMB::REPORT(jnll)

  RTMB::ADREPORT(ssb)
  RTMB::ADREPORT(biomass)
  RTMB::ADREPORT(R)

  return(jnll)
}

# source("~/Documents/GitHub/Rceattle/R/dev/rtmb dev.R", echo=TRUE)

