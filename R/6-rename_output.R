

#' Function to rename derived quantities from Rceattle
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#' @param quantities list of "report" objects from Rceattle.
#'
#' @export
#'
rename_output = function(data_list = NULL, quantities = NULL){

  # Dimension attributed
  max_age <- max(data_list$nages, na.rm = T)
  max_sex <- max(data_list$nsex, na.rm = T)
  sex_labels <- c("Sex combined or females", "males")
  if(max_sex == 1){
    sex_labels <- "Sex combined"
  }
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)


  # Rename ----
  # Vectors
  # - Biological
  names(quantities$avg_R) <- data_list$spnames
  names(quantities$Flimit) <- data_list$spnames
  names(quantities$Ftarget) <- data_list$spnames
  names(quantities$gam_a) <- data_list$spnames
  names(quantities$gam_b) <- data_list$spnames
  names(quantities$R0) <- data_list$spnames
  names(quantities$R_init) <- data_list$spnames
  names(quantities$SPR0) <- data_list$spnames
  names(quantities$SPRFinit) <- data_list$spnames
  names(quantities$SPRlimit) <- data_list$spnames
  names(quantities$SPRtarget) <- data_list$spnames
  names(quantities$steepness) <- data_list$spnames

  # - Fleets
  names(quantities$ln_catch_sd) <- data_list$catch_data$Fleet_name
  names(quantities$ln_index_sd) <- data_list$index_data$Fleet_name


  # 2D array
  # - Population quantities
  dimnames(quantities$biomass) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$biomass_depletion) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$B0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$DynamicB0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$exploitable_biomass) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$ssb_depletion) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$ssb) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$SB0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$DynamicSB0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$DynamicSBF) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$vulnerability) <- list(paste0("Pred: ", data_list$spnames), paste0("Prey: ", data_list$spnames)) # Pred/prey?

  dimnames(quantities$R) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$F_spp) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$proj_F) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$fT) <- list(data_list$spnames, yrs_proj) # Temperature function of consumption

  dimnames(quantities$pop_scalar) <- list(data_list$spnames, paste0("Age", 1:max_age))


  # - Fleet quantities
  dimnames(quantities$F_flt) <- list(data_list$fleet_control$Fleet_name, yrs_proj) # Sex specific?
  dimnames(quantities$index_q) <- list(data_list$fleet_control$Fleet_name, yrs_hind)

  # 4D array
  # - Biological
  dimnames(quantities$biomass_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$consumption_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$B_eaten_as_prey) <-  list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$F_spp_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M1_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M2_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$Z_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$N_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$NByage0) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$NByageF) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$ration) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)

  # - Fleet
  dimnames(quantities$F_flt_age) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)

  # 5D arrays
  dimnames(quantities$B_eaten) <- list(paste("Pred:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                       paste("Prey:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                       paste0("Pred age", 1:max_age),
                                       paste0("Prey age", 1:max_age),
                                       yrs_proj)
  dimnames(quantities$suitability) <- list(paste("Pred:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                           paste("Prey:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                           paste0("Pred age", 1:max_age),
                                           paste0("Prey age", 1:max_age),
                                           yrs_proj)


  # Rename likelihood components
  quantities$jnll_comp[8,1:data_list$nspp] <- 1:data_list$nspp
  quantities$unweighted_jnll_comp[8,1:data_list$nspp] <- 1:data_list$nspp

  quantities$jnll_comp <- rbind(1:nrow(data_list$fleet_control), quantities$jnll_comp)
  quantities$unweighted_jnll_comp <- rbind(1:nrow(data_list$fleet_control), quantities$unweighted_jnll_comp)

  colnames(quantities$jnll_comp) <- 1:ncol(quantities$jnll_comp)
  colnames(quantities$unweighted_jnll_comp) <- 1:ncol(quantities$unweighted_jnll_comp)

  rownames(quantities$jnll_comp) <- rownames(quantities$unweighted_jnll_comp) <- c(
    "1. Fleet components",
    "Index data",
    "Catch data",
    "Composition data",
    "Non-parametric selectivity",
    "Selectivity deviates",
    "Catchability prior",
    "Catchability deviates",
    "2. Species components", # Empty row
    "Stock-recruit prior",
    "Initial abundance deviates",
    "Recruitment deviates",
    "Stock-recruit penalty",
    "Reference point penalities",
    "Zero n-at-age penalty",
    "M prior",
    "M random effects",
    "Ration",
    "Ration penalties",
    "Stomach content data"
  )

  return(quantities)
}


#' Function to calculate McAllister-Ianelli weights
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#' @param data_list_reorganized reorganized data_list
#' @param quantities list of "report" objects from Rceattle.
#'
#' @export
#'
calc_mcall_ianelli <- function(data_list = NULL, data_list_reorganized = NULL, quantities = NULL){


  # - Calculate Mcallister-Iannelli coefficients
  # Effective sample size for the length data for year y
  eff_n_mcallister <- rowSums(quantities$comp_hat * (1 - quantities$comp_hat), na.rm = TRUE)/rowSums((data_list_reorganized$comp_obs - quantities$comp_hat)^2, na.rm = TRUE) # sum_length (p_hat * (1 - p_hat))/ sum_length ((p - p_hat) ^ 2)


  # Loop fleets and take harmonic mean
  data_list$fleet_control$Est_weights_mcallister <- NA
  for(flt in unique(data_list$comp_data$Fleet_code)){
    comp_sub <- which(data_list$comp_data$Fleet_code == flt & data_list$comp_data$Year > 0)
    data_list$fleet_control$Est_weights_mcallister[which(data_list$fleet_control$Fleet_code == flt)] <- ((1/length(comp_sub))*sum((eff_n_mcallister[comp_sub]/data_list$comp_data$Sample_size[comp_sub])^-1))^-1
  }

  return(data_list)
}

#' Match predicted diet proportions to observed data (Final, Robust Version)
#'
#' @description A helper function that handles any mix of data aggregation
#' within a single diet dataset by processing row by row.
#'
#' @param data_list The Rceattle data_list object.
#' @param quantities The 'quantities' object from a model run.
#'
#' @return The input `diet_data` data frame with a new column, `Est_diet`.
#' @export
match_diet_preds <- function(data_list, quantities) {

  obs_diet <- data_list$diet_data
  pred_diet_array <- quantities$diet_prop_hat

  if (is.null(obs_diet) || is.null(pred_diet_array) || nrow(obs_diet) == 0) {
    return(NULL)
  }

  # Create an empty vector to store the results
  est_diet_vec <- numeric(nrow(obs_diet))

  # Loop through each row of the observed diet data
  for (i in 1:nrow(obs_diet)) {
    obs_row <- obs_diet[i, ]

    # Extract indices for this specific observation
    p <- obs_row$Pred
    k <- obs_row$Prey
    pa <- obs_row$Pred_age
    ka <- obs_row$Prey_age
    yr <- obs_row$Year

    # --- Apply aggregation based on the values in THIS row ---

    if (yr > 0 && pa > 0 && ka > 0) { # Fully disaggregated
      estimated_value <- pred_diet_array[p, k, pa, ka, yr]

    } else if (yr == 0 && pa > 0 && ka > 0) { # Year aggregated only
      estimated_value <- mean(pred_diet_array[p, k, pa, ka, ])

    } else if (yr > 0 && pa > 0 && ka < 0) { # Prey-age aggregated only
      estimated_value <- sum(pred_diet_array[p, k, pa, , yr])

    } else if (yr == 0 && pa > 0 && ka < 0) { # Year AND Prey-age aggregated
      annual_sums <- apply(pred_diet_array[p, k, pa, , ], 2, sum) # Sum over prey_age for each year
      estimated_value <- mean(annual_sums)

    } else {
      # This logic can be expanded to handle the other pred_age < 0 cases if needed
      estimated_value <- NA
    }

    est_diet_vec[i] <- estimated_value
  }

  # Append the final vector of estimated values to the original data frame
  obs_diet$Est_diet <- est_diet_vec

  return(obs_diet)
}


#' Function to calculate McAllister-Ianelli weights for diet data (UPDATED)
#'
#' @description Calculates effective sample size data weights for diet composition
#' data using the McAllister-Ianelli method. This version works with the new
#' multinomial likelihood structure.
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}.
#' @param quantities The 'quantities' object from a model run.
#'
#' @return The original data_list with a new element,
#'  `Diet_weights_mcallister`, containing the calculated weights for each
#'   predator species.
#' @export
calc_mcall_ianelli_diet <- function(data_list = NULL, quantities = NULL){

  # 1. Get a clean data frame with both observed and matched estimated proportions
  diet_comparison <- match_diet_preds(data_list, quantities)

  # Return original list if there's nothing to process
  if (is.null(diet_comparison)) {
    return(data_list)
  }

  # Small constant to avoid division by zero
  epsilon <- 1e-10

  # 2. Calculate effective sample size (eff_n) for each unique stomach sample
  eff_n_per_stomach <- diet_comparison %>%
    # A unique stomach sample is defined by these groups
    dplyr::group_by(Pred, Pred_age, Pred_sex, Year) %>%
    dplyr::summarise(
      Sample_size = dplyr::first(Sample_size),
      # The M-I formula: sum(p_hat*(1-p_hat)) / sum((p_obs - p_hat)^2)
      eff_n = sum(Est_diet * (1 - Est_diet), na.rm = TRUE) /
        (sum((Stomach_proportion_by_weight - Est_diet)^2, na.rm = TRUE) + epsilon),
      .groups = 'drop'
    )

  # 3. Calculate the final weight for each predator species using the harmonic mean
  #    of the ratio of effective N to nominal sample size.
  diet_weights <- eff_n_per_stomach %>%
    dplyr::group_by(Pred) %>%
    dplyr::summarise(
      Est_weights_mcallister = (1/n() * sum((eff_n / Sample_size)^-1, na.rm = TRUE))^-1,
      .groups = 'drop'
    )

  # 4. Add the result to the data_list
  data_list$Diet_weights_mcallister <- diet_weights

  return(data_list)
}
