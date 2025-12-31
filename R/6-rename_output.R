

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
  # * Vectors ----
  #  Biological quantities
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

  # * Fleets ----
  names(quantities$ln_catch_sd) <- data_list$catch_data$Fleet_name
  names(quantities$ln_index_sd) <- data_list$index_data$Fleet_name


  # * 2D array ----
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

  # * 4D array ----
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

  dimnames(quantities$weight_hat) <- list(
    c(paste(rep(data_list$spnames, each = 2), rep(c("biomass_weight", "spawn_weight"), data_list$nspp)), data_list$fleet_control$Fleet_name),
    sex_labels, paste0("Age", 1:max_age), yrs_proj)

  # - Fleet
  dimnames(quantities$F_flt_age) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)

  # * 5D arrays ----
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


  # Likelihood components ----
  S = paste0(ss_run$data_list$spnames, " or ")
  X = ss_run$data_list$fleet_control$Fleet_name
  llcolnames <- paste0(append(S, rep("", max(c(0, length(X)-length(S))))), append(X, rep(NA, max(c(0, length(S)-length(X))))))

  colnames(quantities$jnll_comp) <- llcolnames
  colnames(quantities$unweighted_jnll_comp) <- llcolnames

  rownames(quantities$jnll_comp) <- rownames(quantities$unweighted_jnll_comp) <- c(
    "Index data",
    "Catch data",
    "Composition data",
    "CAAL data",
    "Non-parametric selectivity",
    "Selectivity deviates",
    "Catchability prior",
    "Catchability deviates",
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
