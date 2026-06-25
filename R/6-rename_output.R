

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

  obs_diet        <- data_list$diet_data
  pred_diet_array <- quantities$diet_prop_hat

  if (is.null(obs_diet) || is.null(pred_diet_array) || nrow(obs_diet) == 0) {
    return(NULL)
  }

  styr <- data_list$styr

  get_suit_yrs <- function(pred_sp) {
    suit_start <- data_list$suit_styr[pred_sp] - styr + 1
    suit_end   <- data_list$suit_endyr[pred_sp] - styr + 1
    suit_start:suit_end
  }

  est_diet_vec <- numeric(nrow(obs_diet))

  for (i in seq_len(nrow(obs_diet))) {

    obs_row <- obs_diet[i, ]
    p    <- obs_row$Pred
    k    <- obs_row$Prey
    psex <- obs_row$Pred_sex
    ksex <- obs_row$Prey_sex
    pa   <- obs_row$Pred_age
    ka   <- obs_row$Prey_age
    yr   <- obs_row$Year

    p <- p + data_list$nspp * max(c(psex - 1, 0))
    k <- k + data_list$nspp * max(c(ksex - 1, 0))

    suit_yrs <- get_suit_yrs(obs_row$Pred)

    estimated_value <- tryCatch({

      if (yr > 0 && pa > 0 && ka > 0) {
        # Fully disaggregated
        pred_diet_array[p, k, pa, ka, yr]

      } else if (yr == 0 && pa > 0 && ka > 0) {
        # Year aggregated only
        mean(sapply(suit_yrs, function(sy)
          pred_diet_array[p, k, pa, ka, sy]))

      } else if (yr > 0 && pa > 0 && ka < 0) {
        # Prey-age aggregated only
        sum(pred_diet_array[p, k, pa, , yr])

      } else if (yr == 0 && pa > 0 && ka < 0) {
        # Year AND prey-age aggregated
        # Loop over suit years explicitly to avoid dimension collapse
        mean(sapply(suit_yrs, function(sy)
          sum(pred_diet_array[p, k, pa, , sy])))

      } else {
        NA_real_
      }

    }, error = function(e) {
      warning(sprintf("Row %d: %s", i, conditionMessage(e)))
      NA_real_
    })

    est_diet_vec[i] <- estimated_value
  }

  obs_diet$Est_diet <- est_diet_vec
  return(obs_diet)
}




#' Function to calculate McAllister-Ianelli weights for diet data
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#' @param quantities list of "report" objects from Rceattle, including diet_hat predictions
#'
#' @export
#'
calc_mcall_ianelli_diet <- function(data_list = NULL, quantities = NULL){

  # - Calculate Mcallister-Iannelli coefficients for diet data
  diet_multiplier <- data_list$Diet_comp_weights

  # Small constant to avoid division by zero
  epsilon <- 1e-10


  # Calculate effective sample size for diet data (predator specific)
  # Using the same formula as for length: sum(p_hat * (1 - p_hat)) / sum((p - p_hat)^2)
  eff_n_mcallister <- data_list$diet_data %>%
    dplyr::mutate(Diet_hat = quantities$diet_hat) %>%
    dplyr::group_by(Pred, Pred_age) %>%
    dplyr::summarise(
      Sample_size = dplyr::first(Sample_size), # Sample size should be the same across predators of the same age
      eff_n_mcallister = sum(Diet_hat * (1 - Diet_hat), na.rm = TRUE) /
        (sum((Stomach_proportion_by_weight - Diet_hat)^2, na.rm = TRUE) + epsilon)
    )

  # Take harmonic mean across predator ages
  data_list$Diet_weights_mcallister <- eff_n_mcallister %>%
    dplyr::group_by(Pred) %>%
    filter(eff_n_mcallister != 0) %>%
    dplyr::summarise(Diet_weights_mcallister = (1/n() * sum((eff_n_mcallister /Sample_size)^-1))^-1 ) %>%
    dplyr::pull(Diet_weights_mcallister)

  return(data_list)
}
