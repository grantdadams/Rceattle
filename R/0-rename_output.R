

#' Function to rename derived quantities from Rceattle
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#' @param quantities list of "report" objects from Rceattle.
#'
#' @return
#' @export
#'
rename_output = function(data_list = NULL, quantities = NULL){

  # Dimension attributed
  max_age <- max(data_list$nages, na.rm = T)
  max_sex <- max(data_list$nsex, na.rm = T)
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

  # 3D
  dimnames(quantities$M1_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age))

  # 4D array
  # - Biological
  dimnames(quantities$biomass_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$consumption_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$B_eaten_as_prey) <-  list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$F_spp_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M2_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$Z_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$N_at_age) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$NByage0) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$NByageF) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$ration) <- list(data_list$spnames, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)

  # - Fleet
  dimnames(quantities$F_flt_age) <- list(data_list$fleet_control$Fleet_name, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel) <- list(data_list$fleet_control$Fleet_name, c("Sex combines/females", "males"), paste0("Age", 1:max_age), yrs_proj)

  # 5D arrays
  dimnames(quantities$B_eaten) <- list(paste("Pred:", data_list$spnames, rep(c("Sex combines/females", "males"), each = 3)),
                                       paste("Prey:", data_list$spnames, rep(c("Sex combines/females", "males"), each = 3)),
                                       paste0("Pred age", 1:max_age),
                                       paste0("Prey age", 1:max_age),
                                       yrs_proj)


  # Rename jnll
  colnames(quantities$jnll_comp) <- paste0("Sp/Srv/Fsh_", 1:ncol(quantities$jnll_comp))
  rownames(quantities$jnll_comp) <- c(
    "Index data",
    "Catch data",
    "Composition data",
    "Sex ratio",
    "Non-parametric selectivity",
    "Selectivity deviates",
    "Selectivity normalization",
    "Catchability prior",
    "Catchability deviates",
    "Stock-recruit prior",
    "Recruitment deviates",
    "Initial abundance deviates",
    "Fishing mortality deviates",
    "SPR Calculation",
    "Zero n-at-age penalty",
    "M prior",
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
#' @return
#' @export
#'
#' @examples
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
