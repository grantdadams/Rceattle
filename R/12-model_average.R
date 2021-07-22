#' Model average of derived quantities
#'
#' @param Rceattle list of Rceattle model objects
#' @param weights vector of weights to be used for weighting models
#'
#' @return an Rceattle object with derived quantities weighted by the specified weights. The length of the derived quantities spans the years which overlap across all models.
#' @export
#'
model_average <- function(Rceattle, weights = NULL){

  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    stop("Only one model provided")
  }
  if(is.null(weights)){
    weights = rep(1/length(Rceattle), length(Rceattle))
  }
  weights <- weights/sum(weights)

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
  mod_avg$estimated_params <- mod_avg$initial_params <- mod_avg$opt <- mod_avg$run_time <- mod_avg$sdrep <- mod_avg$obj <- mod_avg$map <- mod_avg$bounds <- NULL

  # -- Set all quantities to zero
  for(i in 1:length(mod_avg$quantities)){
    mod_avg$quantities[[i]] <- replace(mod_avg$quantities[[i]], values = rep(0, length(mod_avg$quantities[[i]])))
  }


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
    mod_avg$quantities$NbyageSPR <- mod_avg$quantities$NbyageSPR + mod_avg$quantities$NbyageSPR * weights[i]
    mod_avg$quantities$SB0 <- mod_avg$quantities$SB0 + Rceattle[[i]]$quantities$SB0 * weights[i]
    mod_avg$quantities$SB35 <- mod_avg$quantities$SB35 + Rceattle[[i]]$quantities$SB35 * weights[i]
    mod_avg$quantities$SB40 <- mod_avg$quantities$SB40 + Rceattle[[i]]$quantities$SB40 * weights[i]


    # // -- Survey components
    # mod_avg$quantities$srv_bio_hat );
    # mod_avg$quantities$srv_log_sd_hat );


    # // Fishery components
    mod_avg$quantities$F[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$F[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$F[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$F_tot[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$F_tot[,,,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$F_tot[,,,sub_rel_proj_yrs] * weights[i]
    # mod_avg$quantities$fsh_bio_hat );
    # mod_avg$quantities$fsh_log_sd_hat );
    mod_avg$quantities$proj_FABC[,mod_avg_rel_proj_yrs] <- mod_avg$quantities$proj_FABC[,mod_avg_rel_proj_yrs] + Rceattle[[i]]$quantities$proj_FABC[,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$FSPR <- mod_avg$quantities$FSPR + Rceattle[[i]]$quantities$FSPR * weights[i]
    mod_avg$quantities$F35_tot <- mod_avg$quantities$F35_tot + Rceattle[[i]]$quantities$F35_tot * weights[i]
    mod_avg$quantities$F40_tot <- mod_avg$quantities$F40_tot + Rceattle[[i]]$quantities$F40_tot * weights[i]


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

  return(mod_avg)
}
