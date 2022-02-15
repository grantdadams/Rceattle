
# Not working via TMB so copied from https://rdrr.io/cran/TMB/src/R/TMB.R
#' Test for invalid external pointer
#'
#' @param pointer
#'
#' @return
#' @export
#'
isNullPointer <- function(pointer) {
  .Call("isNullPointer", pointer, PACKAGE="TMB")
}


#' Add external pointer finalizer
#'
#' @param ADFun
#' @param DLL
#'
#' @return
#' @export
#'
registerFinalizer <- function(ADFun, DLL) {
  finalizer <- function(ptr) {
    if ( ! isNullPointer(ptr) ) {
      .Call("FreeADFunObject", ptr, PACKAGE=DLL)
    } else {
      ## Nothing to free
    }
  }
  reg.finalizer(ADFun$ptr, finalizer)
}


#' Model average of derived quantities
#'
#' @param Rceattle list of Rceattle model objects
#' @param weights vector of weights to be used for weighting models
#' @param Uncertainty TRUE/FALSE Sample uncertainty across derived quantities using weighted bootstrap from the asymptotic  distribution of MLEs
#' @param nboot Number of bootstraps taken from asymptotic distribution of MLEs. Default = 10000
#'
#' @return an Rceattle object with derived quantities weighted by the specified weights. The length of the derived quantities spans the years which overlap across all models.
#' @export
#'
model_average <- function(Rceattle, weights = NULL, uncertainty = FALSE, nboot = 10000){
  # --------------------------------------------------------------------------------------------
  # Average derived quantities of models
  # --------------------------------------------------------------------------------------------
  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    stop("Only one model provided")
  }
  if(is.null(weights)){
    weights = rep(1/length(Rceattle), length(Rceattle))
  }
  weights <- weights/sum(weights) # Normalize

  # Extract number of species of each model
  nspp <- sapply(Rceattle, function(x) x$data_list$nspp)
  nspp <- unique(nspp)
  if(length(nspp) > 1){stop("Number of species does not match across models")}

  # Extract timespan of each model
  years <- sapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  endyrs <- sapply(Rceattle, function(x) x$data_list$endyr)
  styrs <- sapply(Rceattle, function(x) x$data_list$styr)
  projyrs <- sapply(Rceattle, function(x) x$data_list$projyr)
  nyrs <- sapply(years, length)

  min_projyr <- min(projyrs, na.rm = TRUE) # Find coverage across all models
  min_hindyr <- min(endyrs, na.rm = TRUE) # Find coverage across all models
  max_styr <- max(styrs, na.rm = TRUE)

  # Initialize model average object
  mod_avg <- Rceattle[[which(nyrs == min(nyrs))[1]]] # Copy a model of the smallest length to fill in
  mod_avg_rel_proj_yrs <- (max_styr - mod_avg$data_list$styr + 1) : (min_projyr - mod_avg$data_list$styr + 1) # Find years of overlap
  mod_avg_rel_hind_yrs <- (max_styr - mod_avg$data_list$styr + 1) : (min_hindyr - mod_avg$data_list$styr + 1) # Find years of overlap
  mod_avg$estimated_params <- mod_avg$initial_params <- mod_avg$opt <- mod_avg$run_time <- mod_avg$obj <- mod_avg$map <- mod_avg$bounds <- NULL
  if(!uncertainty){mod_avg$sdrep <- NULL}

  # -- Set all quantities to zero
    mod_avg$quantities <- sapply(mod_avg$quantities, function(x) replace(x, values = rep(0, length(x))))



  # Average out derived quantities (Note: dimensions may be different)
  for(i in 1:length(Rceattle)){
    sub_rel_proj_yrs <- (max_styr - Rceattle[[i]]$data_list$styr + 1) : (min_projyr - Rceattle[[i]]$data_list$styr + 1) # Find years of overlap
    sub_rel_hind_yrs <- (max_styr - Rceattle[[i]]$data_list$styr + 1) : (min_hindyr - Rceattle[[i]]$data_list$styr + 1) # Find years of overlap

    # Average quantities of interest
    # // -- Biomass components
    mod_avg$quantities$R[, mod_avg_rel_proj_yrs] <- mod_avg$quantities$R[, mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$R[, sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$mn_rec <- mod_avg$quantities$mn_rec + Rceattle[[i]]$quantities$mn_rec * weights[i]

    mod_avg$quantities$NByage[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$NByage[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$NByage[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$AvgN[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$AvgN[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$AvgN[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$biomassByage[,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$biomassByage[,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$biomassByage[,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$biomassSSBByage[,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$biomassSSBByage[,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$biomassSSBByage[,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$biomass[, mod_avg_rel_proj_yrs] <- mod_avg$quantities$biomass[, mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$biomass[, sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$biomassSSB[, mod_avg_rel_proj_yrs] <- mod_avg$quantities$biomassSSB[, mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$biomassSSB[, sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$B_eaten[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$B_eaten[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$B_eaten[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$sel[,,,mod_avg_rel_hind_yrs] <- mod_avg$quantities$sel[,,,mod_avg_rel_hind_yrs] + Rceattle[[i]]$quantities$sel[,,,sub_rel_hind_yrs] * weights[i]

    # // -- Reference points
    mod_avg$quantities$SB0 <- mod_avg$quantities$SB0 + Rceattle[[i]]$quantities$SB0 * weights[i]

    # // -- Survey components
    # mod_avg$quantities$srv_bio_hat );
    # mod_avg$quantities$srv_log_sd_hat );


    # // Fishery components
    # mod_avg$quantities$fsh_bio_hat );
    # mod_avg$quantities$fsh_log_sd_hat );


    # // 12.5. Age/length composition
    # mod_avg$quantities$age_obs_hat );
    # mod_avg$quantities$age_hat );
    # mod_avg$quantities$comp_obs );
    # mod_avg$quantities$comp_hat );
    # mod_avg$quantities$true_age_comp_hat );
    # mod_avg$quantities$n_hat );
    # mod_avg$quantities$comp_n );

    # # -- Ration components
    # mod_avg$quantities$ConsumAge );
    # mod_avg$quantities$LbyAge );
    # mod_avg$quantities$mnWt_obs );
    # mod_avg$quantities$fT );
    # mod_avg$quantities$env_index_hat );
    # mod_avg$quantities$ration2Age );


    # Mortality components
    # mod_avg$quantities$suma_suit );
    # mod_avg$quantities$suit_main );
    # mod_avg$quantities$suit_other );
    # mod_avg$quantities$stom_div_bio2 );
    # mod_avg$quantities$stomKir );
    # mod_avg$quantities$stomKirWt );
    # mod_avg$quantities$avail_food );
    # mod_avg$quantities$othersuit );
    # mod_avg$quantities$of_stomKir );
    mod_avg$quantities$M1 <- mod_avg$quantities$M1 + Rceattle[[i]]$quantities$M1 * weights[i]
    mod_avg$quantities$M2[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$M2[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$M2[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$M[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$M[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$M[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$Zed[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$Zed[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$Zed[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$S[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$S[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$S[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$M2_prop[,,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$M2_prop[,,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$M2_prop[,,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$B_eaten_prop[,,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$B_eaten_prop[,,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$B_eaten_prop[,,,,sub_rel_proj_yrs] * weights[i]
    # mod_avg$quantities$UobsAge_hat );
    # mod_avg$quantities$UobsWtAge_hat );
  }

  # --------------------------------------------------------------------------------------------
  # Sample uncertainty across derived quantities weighted bootstrap from the asymptotic intervals
  # --------------------------------------------------------------------------------------------

  if(uncertainty){

    # # Buckland, S.T., Burnham, K.P., Augustin, N.H., 1997. Model Selection : An Integral Part of Inference. Biometrics 53, 603–618.
    # # - Calculate SD
    # # -- R
    # rec_rows <- which(names(mod_avg$sdrep$value) == "R")
    # mod_avg$sdrep$sd[rec_rows] <- 0
    # mod_avg$sdrep$value[rec_rows] <- 0
    #
    #
    # # -- B
    # biomass_rows <- which(names(mod_avg$sdrep$value) == "biomass")
    # mod_avg$sdrep$sd[biomass_rows] <- 0
    # mod_avg$sdrep$value[biomass_rows] <- 0
    #
    #
    # # -- SSB
    # ssb_rows <- which(names(mod_avg$sdrep$value) == "biomassSSB")
    # mod_avg$sdrep$sd[ssb_rows] <- 0
    # mod_avg$sdrep$value[ssb_rows] <- 0
    #
    # # Loop across models and calculate mean
    # for(i in 1:length(Rceattle)){
    #   # Years
    #   sub_rel_proj_yrs <- (max_styr - Rceattle[[i]]$data_list$styr + 1) : (min_projyr - Rceattle[[i]]$data_list$styr + 1) # Find years of overlap
    #
    #
    #   # Extract Rec, SSB, B
    #   # - R
    #   mu_rec <- array(NA, dim = c(nspp, nyrs[i]))
    #   rec_rows <- which(names(Rceattle[[i]]$sdrep$value) == "R")
    #   mu_rec <- replace(mu_rec, values = Rceattle[[i]]$sdrep$value[rec_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   # - B
    #   mu_biomass <- array(NA, dim = c(nspp, nyrs[i]))
    #   biomass_rows <- which(names(Rceattle[[i]]$sdrep$value) == "biomass")
    #   mu_biomass <- replace(mu_biomass, values = Rceattle[[i]]$sdrep$value[biomass_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   # - SSB
    #   mu_biomassSSB <- array(NA, dim = c(nspp, nyrs[i]))
    #   ssb_rows <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
    #   mu_biomassSSB <- replace(mu_biomassSSB, values = Rceattle[[i]]$sdrep$value[ssb_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #
    #   # - Calculate weighted mean
    #   # -- R
    #   rec_rows <- which(names(mod_avg$sdrep$value) == "R")
    #   mod_avg$sdrep$value[rec_rows] <- mod_avg$sdrep$value[rec_rows] + c(mu_rec * weights[i])
    #
    #
    #   # -- B
    #   biomass_rows <- which(names(mod_avg$sdrep$value) == "biomass")
    #   mod_avg$sdrep$value[biomass_rows] <- mod_avg$sdrep$value[biomass_rows] + c(mu_biomass * weights[i])
    #
    #
    #   # -- SSB
    #   ssb_rows <- which(names(mod_avg$sdrep$value) == "biomassSSB")
    #   mod_avg$sdrep$value[ssb_rows] <- mod_avg$sdrep$value[ssb_rows] + c(mu_biomassSSB * weights[i])
    # }
    #
    #
    # # Loop across models and calculate var
    # for(i in 1:length(Rceattle)){
    #   # - R
    #   sd_rec <- array(NA, dim = c(nspp, nyrs[i]))
    #   mu_rec <- array(NA, dim = c(nspp, nyrs[i]))
    #
    #   rec_rows <- which(names(Rceattle[[i]]$sdrep$value) == "R")
    #   sd_rec <- replace(sd_rec, values = Rceattle[[i]]$sdrep$sd[rec_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #   mu_rec <- replace(mu_rec, values = Rceattle[[i]]$sdrep$value[rec_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   # - B
    #   sd_biomass <- array(NA, dim = c(nspp, nyrs[i]))
    #   mu_biomass <- array(NA, dim = c(nspp, nyrs[i]))
    #
    #   biomass_rows <- which(names(Rceattle[[i]]$sdrep$value) == "biomass")
    #   sd_biomass <- replace(sd_biomass, values = Rceattle[[i]]$sdrep$sd[biomass_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #   mu_biomass <- replace(mu_biomass, values = Rceattle[[i]]$sdrep$value[biomass_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   # - SSB
    #   sd_biomassSSB <- array(NA, dim = c(nspp, nyrs[i]))
    #   mu_biomassSSB <- array(NA, dim = c(nspp, nyrs[i]))
    #
    #   ssb_rows <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
    #   sd_biomassSSB <- replace(sd_biomassSSB, values = Rceattle[[i]]$sdrep$sd[ssb_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #   mu_biomassSSB <- replace(mu_biomassSSB, values = Rceattle[[i]]$sdrep$value[ssb_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #
    #   # - Calculate weighted mean
    #   # -- R
    #   rec_rows <- which(names(mod_avg$sdrep$value) == "R")
    #   mod_avg$sdrep$sd[rec_rows] <- mod_avg$sdrep$sd[rec_rows] + weights[i] * sqrt(sd_rec^2 + (mu_rec - mod_avg$sdrep$value[rec_rows])^2)
    #
    #
    #   # -- B
    #   biomass_rows <- which(names(mod_avg$sdrep$value) == "biomass")
    #   mod_avg$sdrep$sd[biomass_rows] <- mod_avg$sdrep$sd[biomass_rows] + weights[i] * sqrt(sd_biomass^2 + (mu_biomass - mod_avg$sdrep$value[biomass_rows])^2)
    #
    #
    #   # -- SSB
    #   ssb_rows <- which(names(mod_avg$sdrep$value) == "biomassSSB")
    #   mod_avg$sdrep$sd[ssb_rows] <- mod_avg$sdrep$sd[ssb_rows] + weights[i] * sqrt(sd_biomassSSB^2 + (mu_biomassSSB - mod_avg$sdrep$value[ssb_rows])^2)
    # }
    #



    # Using univariate normal asymptotics of ssb
    # weights <- round(weights * nboot)
    # nrows <- 1:sum(weights)
    # rowid <- rep(1:length(weights), weights)
    #
    # # List to save samples of SSB, B, and R
    # samples_rec <- list()
    # samples_biomassSSB <- list()
    # samples_biomass <- list()
    #
    # # Loop across models
    # for(i in 1:length(Rceattle)){
    #   # Years
    #   sub_rel_proj_yrs <- (max_styr - Rceattle[[i]]$data_list$styr + 1) : (min_projyr - Rceattle[[i]]$data_list$styr + 1) # Find years of overlap
    #
    #
    #   # Extract Rec, SSB, B
    #   # - R
    #   sd_rec <- array(NA, dim = c(nspp, nyrs[i]))
    #   mu_rec <- array(NA, dim = c(nspp, nyrs[i]))
    #
    #   rec_rows <- which(names(Rceattle[[i]]$sdrep$value) == "R")
    #   sd_rec <- replace(sd_rec, values = Rceattle[[i]]$sdrep$sd[rec_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #   mu_rec <- replace(mu_rec, values = Rceattle[[i]]$sdrep$value[rec_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   samples_rec[[i]] <- MASS::mvrnorm(weights[i], mu = c(mu_rec), Sigma = diag(c(sd_rec)))
    #
    #
    #   # - B
    #   sd_biomass <- array(NA, dim = c(nspp, nyrs[i]))
    #   mu_biomass <- array(NA, dim = c(nspp, nyrs[i]))
    #
    #   biomass_rows <- which(names(Rceattle[[i]]$sdrep$value) == "biomass")
    #   sd_biomass <- replace(sd_biomass, values = Rceattle[[i]]$sdrep$sd[biomass_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #   mu_biomass <- replace(mu_biomass, values = Rceattle[[i]]$sdrep$value[biomass_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   samples_biomass[[i]] <- MASS::mvrnorm(weights[i], mu = c(mu_biomass), Sigma = diag(c(sd_biomass)))
    #
    #
    #   # - SSB
    #   sd_biomassSSB <- array(NA, dim = c(nspp, nyrs[i]))
    #   mu_biomassSSB <- array(NA, dim = c(nspp, nyrs[i]))
    #
    #   ssb_rows <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
    #   sd_biomassSSB <- replace(sd_biomassSSB, values = Rceattle[[i]]$sdrep$sd[ssb_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #   mu_biomassSSB <- replace(mu_biomassSSB, values = Rceattle[[i]]$sdrep$value[ssb_rows])[,sub_rel_proj_yrs] # Remove years not shared across models
    #
    #   samples_biomassSSB[[i]] <- MASS::mvrnorm(weights[i], mu = c(mu_biomassSSB), Sigma = diag(c(sd_biomassSSB)))
    # }
    #
    #
    # # - Save samples
    # mod_avg$asymptotic_samples <- list(recruitment = do.call(rbind, samples_rec),
    #                                    biomassSSB = do.call(rbind, samples_biomassSSB),
    #                                    biomass = do.call(rbind, samples_biomass))
    #
    # # - Calculate SD
    # # -- R
    # rec_rows <- which(names(mod_avg$sdrep$value) == "R")
    # mod_avg$sdrep$sd[rec_rows] <- sqrt(apply(mod_avg$asymptotic_samples$recruitment, 2, var))
    # check <- data.frame(Sim =  mod_avg$sdrep$sd[rec_rows], Mod3 = Rceattle[[3]]$sdrep$sd[rec_rows])
    # check$diff <- check$Sim - check$Mod3
    # check
    #
    # mod_avg$sdrep$value[rec_rows] <- colMeans(mod_avg$asymptotic_samples$recruitment)
    #
    #
    # # -- B
    # biomass_rows <- which(names(mod_avg$sdrep$value) == "biomass")
    # mod_avg$sdrep$sd[biomass_rows] <- sqrt(apply(mod_avg$asymptotic_samples$biomass, 2, var))
    # mod_avg$sdrep$value[biomass_rows] <- colMeans(mod_avg$asymptotic_samples$biomass)
    #
    #
    # # -- SSB
    # ssb_rows <- which(names(mod_avg$sdrep$value) == "biomassSSB")
    # mod_avg$sdrep$sd[ssb_rows] <- sqrt(apply(mod_avg$asymptotic_samples$biomassSSB, 2, var))
    # mod_avg$sdrep$value[ssb_rows] <- colMeans(mod_avg$asymptotic_samples$biomassSSB)


    # Assuming asymptotic multivariate normal
    # - Must have estimated joint precision matrix
    length_ran <- sapply(Rceattle, function(x) !is.null(x$sdrep$par.random)) # All parameters random effects
    joint_est <- sapply(Rceattle, function(x) is.null(x$sdrep$jointPrecision))
    if(sum(length_ran) > 0){
      if(sum(joint_est) > 0){
        stop(paste("Model(s):", which(joint_est > 0)," do not have joint precision matrices: re-estimate with getJointPrecision = TRUE"))
      }

      if(sum(length_ran) != length(Rceattle)){
        stop("Does not currently support purely fixed and random effects models together")
      }
    }



    # - Pointers
    draws <- round(weights * nboot)
    nrows <- 1:sum(draws)
    rowid <- rep(1:length(draws), draws)

    # List to save samples of SSB, B, and R
    recruitment <-
      array(NA, dim = c(nspp, length(mod_avg_rel_proj_yrs),  sum(draws)))
    biomassSSB <-
      array(NA, dim = c(nspp, length(mod_avg_rel_proj_yrs),  sum(draws)))
    biomass <-
      array(NA, dim = c(nspp, length(mod_avg_rel_proj_yrs),  sum(draws)))


    # Loop across models
    for(i in 1:length(Rceattle)){

      # - Sample parameters from asymptotic normal distribution
      # - Models with random effects
      if(sum(length_ran) > 0){
      mle <- Rceattle[[i]]$obj$env$last.par.best # Includes fixed and random
      vcov <- solve(Rceattle[[i]]$sdrep$jointPrecision) # names(mle) == rownames(vcov)
      vcov[1,1] <- 100
      }


      # - Fixed effects models
      if(sum(length_ran) == 0){
        mle <- Rceattle[[i]]$obj$env$last.par.best # Includes fixed
        vcov <- solve(Rceattle[[i]]$sdrep$cov.fixed) # names(mle) == rownames(vcov)
      }
      samples <- MASS::mvrnorm(draws[i], mu = mle, Sigma = vcov)

      # - Get quantities
      quantities <- apply(samples, 1, function(x) Rceattle[[i]]$obj$report(x)[c("biomass", "biomassSSB", "R")]) # Only want uncertainty in SSB, B, and R

      # - Subset years of interest and assign to objects
      sub_rel_proj_yrs <- (max_styr - Rceattle[[i]]$data_list$styr + 1) : (min_projyr - Rceattle[[i]]$data_list$styr + 1) # Find years of overlap

      # Extract R, B, and SSB and assign to objects
      recruitment[1:nspp, 1:length(mod_avg_rel_proj_yrs), nrows[which(rowid == i)]] <- unlist(lapply(quantities, function(x) x$R[, sub_rel_proj_yrs]))
      biomassSSB[1:nspp, 1:length(mod_avg_rel_proj_yrs), nrows[which(rowid == i)]] <- unlist(lapply(quantities, function(x) x$biomassSSB[, sub_rel_proj_yrs]))
      biomass[1:nspp, 1:length(mod_avg_rel_proj_yrs), nrows[which(rowid == i)]] <- unlist(lapply(quantities, function(x) x$biomass[, sub_rel_proj_yrs]))
    }
    plot_ssb(mod_avg, mod_avg = FALSE, add_ci = TRUE)

    # - Calculate SD
    rec_rows <- which(names(mod_avg$sdrep$value) == "R")
    mod_avg$sdrep$sd[rec_rows] <- sqrt(apply(recruitment, c(1,2), var))

    biomass_rows <- which(names(mod_avg$sdrep$value) == "biomass")
    mod_avg$sdrep$sd[biomass_rows] <- sqrt(apply(biomass, c(1,2), var))

    ssb_rows <- which(names(mod_avg$sdrep$value) == "biomassSSB")
    mod_avg$sdrep$sd[ssb_rows] <- sqrt(apply(biomassSSB, c(1,2), var))

    # - Save samples
    mod_avg$asymptotic_samples <- list(recruitment = recruitment, biomass = biomass, biomassSSB = biomassSSB)
  }

  return(mod_avg)
}
