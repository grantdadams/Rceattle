#' Function to clean data prior to Rceattle runs
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
clean_data <- function(data_list){

  # Transpose fleet control
  data_list <- transpose_fleet_control(data_list)

  # Update switches if necessary for older files
  data_list <- revert_switches(data_list)

  # --- 1. Filter Data by Year ----
  # Data in likelihood (use absolute Year)
  abs_year_data <- c("index_data", "catch_data", "comp_data", "caal_data")
  for(df_name in abs_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- data_list[[df_name]] %>%
        dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
    }
  }

  # Fixed data (allow Year == 0)
  fixed_year_data <- c("diet_data", "weight", "emp_sel", "NByageFixed", "ration_data")
  for(df_name in fixed_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- as.data.frame(data_list[[df_name]]) %>%
        dplyr::filter((Year >= data_list$styr & Year <= data_list$projyr) | Year == 0)
    }
  }

  # --- 2. Add temp multi-species SB0 ----
  #FIXME:may be redundant now?
  if(is.null(data_list$MSSB0)){
    data_list$MSSB0 <- rep(999, data_list$nspp)
    data_list$MSB0 <- rep(999, data_list$nspp)
  }


  # --- 3. Extend catch data to proj year for projections ----
  if(data_list$projyr > data_list$endyr){
    for(flt in unique(data_list$catch_data$Fleet_code)){
      catch_data_sub <- data_list$catch_data %>% dplyr::filter(Fleet_code == flt)

      yrs_proj <- (data_list$endyr + 1):data_list$projyr
      yrs_proj <- yrs_proj[which(!yrs_proj %in% catch_data_sub$Year)]
      nyrs_proj <- length(yrs_proj)

      if(nyrs_proj > 0) {
        proj_catch_data <- data.frame(
          Fleet_name = rep(catch_data_sub$Fleet_name[1], nyrs_proj),
          Fleet_code = rep(flt, nyrs_proj),
          Species = rep(catch_data_sub$Species[1], nyrs_proj),
          Year = yrs_proj,
          Month = rep(dplyr::last(catch_data_sub$Month), nyrs_proj),
          Selectivity_block = rep(dplyr::last(catch_data_sub$Selectivity_block), nyrs_proj),
          Catch = rep(NA, nyrs_proj),
          Log_sd = rep(dplyr::last(catch_data_sub$Log_sd), nyrs_proj)
        )
        data_list$catch_data <- rbind(data_list$catch_data, proj_catch_data)
      }
    }
  }


  # --- 4. Column names ----
  expected_cols <- c("Species", paste0("Age", 1:max(data_list$nages)))

  if(any(!colnames(data_list$sex_ratio) %in% expected_cols)){
    colnames(data_list$sex_ratio) <- expected_cols
    message("Renaming column names in 'sex_ratio' data to 'Age1', 'Age2', ....")
  }

  if(any(!colnames(data_list$maturity) %in% expected_cols)){
    colnames(data_list$maturity) <- expected_cols
    message("Renaming column names in 'maturity' data to 'Age1', 'Age2', ....")
  }


  # --- 5. Arrange diet data ----
  if(!is.null(data_list$diet_data)){
    data_list$diet_data <- data_list$diet_data %>%
      dplyr::arrange(Pred, Pred_sex, Pred_age, Prey, Prey_sex, Prey_age, Year) %>%
      dplyr::mutate(stratum_id = paste(Pred, Pred_sex, Pred_age, Year, sep = "_"),
                    stomach_id = as.numeric(as.factor(stratum_id)) - 1) %>%
      dplyr::arrange(stomach_id)
  }

  return(data_list)
}



#' Function to check for missing switches for map and parameter functions
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
switch_check <- function(data_list){

  # Helper to set defaults and notify
  set_default <- function(val, default, msg) {
    if(is.null(val)) {
      message(msg)
      return(default)
    }
    return(val)
  }

  if(is.null(data_list$srr_fun)){
    data_list$srr_fun <- 0
    data_list$srr_pred_fun <- 0
    data_list$srr_est_mode <- 0
    message("'srr_fun' are not included in data, assuming 0")
  }

  # Model and multi-species switches
  data_list$estDynamics <- set_default(data_list$estDynamics, rep(0, data_list$nspp), "'estDynamics' are not included in data, assuming 0")
  data_list$Diet_comp_weights <- set_default(data_list$Diet_comp_weights, rep(1, data_list$nspp), "'Diet_comp_weights' are not included in data, assuming 1")
  data_list$Diet_loglike <- set_default(data_list$Diet_loglike, rep(0, data_list$nspp), "'Diet_loglike' are not included in data, assuming multinomial")
  data_list$alpha_wt_len <- set_default(data_list$alpha_wt_len, 1e-6, "'alpha_wt_len' not specified in data, assuming 1e-6")
  data_list$beta_wt_len <- set_default(data_list$beta_wt_len, 3, "'beta_wt_len' not specified in data, assuming 3")
  data_list$M1_model <- set_default(data_list$M1_model, rep(0, data_list$nspp), "'M1_model' is not included in data, assuming 0")
  data_list$msmMode <- set_default(data_list$msmMode, 0, "'msmMode' is not included in data, assuming single-species (0)")
  data_list$M1_re <- set_default(data_list$M1_re, rep(0, data_list$nspp), "'M1_re' is not in data, assuming 0 for all species")
  data_list$initMode <- set_default(data_list$initMode, 2, "'initMode' is not in the data, setting to 2 (default)")

  # Fleet Control defaults
  data_list$fleet_control$Sel_norm_bin1 <- set_default(data_list$fleet_control$Sel_norm_bin1, NA, "'Sel_norm_bin1' not specified in 'fleet_control', assuming 'NA'")
  data_list$fleet_control$Sel_norm_bin2 <- set_default(data_list$fleet_control$Sel_norm_bin2, NA, "'Sel_norm_bin2' not specified in 'fleet_control', assuming 'NA'")
  data_list$fleet_control$Sel_curve_pen1 <- set_default(data_list$fleet_control$Sel_curve_pen1, 0, "'Sel_curve_pen1' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Sel_curve_pen2 <- set_default(data_list$fleet_control$Sel_curve_pen2, 0, "'Sel_curve_pen2' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Selectivity_dimension <- set_default(data_list$fleet_control$Selectivity_dimension, "Age", "'Selectivity_dimension' not specified in 'fleet_control', assuming 'Age'")
  data_list$fleet_control$Comp_loglike <- set_default(data_list$fleet_control$Comp_loglike, -1, "'Comp_loglike' not specified in 'fleet_control', assuming multinomial")
  data_list$fleet_control$CAAL_loglike <- set_default(data_list$fleet_control$CAAL_loglike, 0, "'CAAL_loglike' not specified in 'fleet_control', assuming multinomial")
  data_list$fleet_control$CAAL_weights <- set_default(data_list$fleet_control$CAAL_weights, 1, "'CAAL_weights' not specified in 'fleet_control', assuming 1")
  data_list$fleet_control$Month <- set_default(data_list$fleet_control$Month, 0, "'Month' not specified in 'fleet_control', assuming 0")

  # Format adjustment for NonParametric
  np_idx <- data_list$fleet_control$Selectivity %in% c(2, "NonParametric", "Non-parametric")
  if(any(np_idx & !is.na(data_list$fleet_control$Time_varying_sel) & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1)))){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::mutate(
        Sel_curve_pen1 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1)), Time_varying_sel, Sel_curve_pen1),
        Sel_curve_pen2 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1)), Time_varying_sel_sd_prior, Sel_curve_pen2),
        Time_varying_sel = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1)), 0, Time_varying_sel),
        Time_varying_sel_sd_prior = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1)), 0, Time_varying_sel_sd_prior)
      )
    message("Updating format where 'Selectivity == Non-parametric'. Moving non-parametric penalties to 'Sel_curve_pen1' and 'Sel_curve_pen2'.")
  }

  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen1))) stop("'Sel_curve_pen1' is NA in 'fleet_control' for fleet with non-parametric selectivity")
  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen2))) stop("'Sel_curve_pen2' is NA in 'fleet_control' for fleet with non-parametric selectivity")

  return(data_list)
}


#' Convert integer switches to intuitive text strings. Maintains backwards compatability.
#'
#' @param data_list Rceattle data list
#'
revert_switches <- function(data_list) {

  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(
      Selectivity = ifelse(as.character(Selectivity) %in% names(sel_rev_map),
                           sel_rev_map[as.character(Selectivity)],
                           Selectivity)
      # ,Catchability = ifelse(as.character(Catchability) %in% names(q_rev_map),
      #                       q_rev_map[as.character(Catchability)],
      #                       Catchability)
    )

  return(data_list)
}

#' Function to transpose fleet_control if long format
#'
#' @param data_list Rceattle data list
#'
#' @export
transpose_fleet_control <- function(data_list){

  if(sum(colnames(data_list$fleet_control)[1:2] == c("Fleet_name", "Fleet_code")) != 2){ #, "Fleet_type", "Species", "Selectivity_index", "Selectivity")) != 6){
    data_list$fleet_control <- as.data.frame(t(data_list$fleet_control))
    colnames(data_list$fleet_control) <- data_list$fleet_control[1,]
    data_list$fleet_control <- data_list$fleet_control[-1,]
    data_list$fleet_control <- cbind(data.frame(Fleet_name = rownames(data_list$fleet_control)),
                                     data_list$fleet_control)
    rownames(data_list$fleet_control) = NULL

    data_list$fleet_control[,-which(colnames(data_list$fleet_control) %in% c("Fleet_name", "Time_varying_q"))] <- apply(
      data_list$fleet_control[,-which(colnames(data_list$fleet_control) %in% c("Fleet_name", "Time_varying_q"))], 2, as.numeric)
  }
  return(data_list)
}

