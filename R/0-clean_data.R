#' Function to clean data prior to Rceattle runs
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
clean_data <- function(data_list){

  # --- 0. Default-fill optional data.frame fields ----
  # These fields can be NULL when the user is not using the corresponding feature
  # (e.g., comp_data in a model without composition likelihoods, ration_data /
  # diet_data in single-species mode, NByageFixed when estDynamics == 0).
  # Filling with empty data.frames that carry the metadata columns the
  # downstream code expects lets rearrange_data / build_params use uniform
  # `nrow > 0` checks without separate NULL guards. Whether the missing data is
  # actually a problem is enforced by data_check() based on the model
  # configuration (msmMode, growth_model, estDynamics, Selectivity, etc.).
  default_dfs <- list(
    comp_data   = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Age0_Length1 = integer(),
                             Year = integer(), Month = integer(),
                             Sample_size = numeric()),
    caal_data   = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Year = integer(),
                             Length = numeric(), Sample_size = numeric()),
    emp_sel     = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Year = integer()),
    NByageFixed = data.frame(Species = integer(), Sex = integer(),
                             Year = integer()),
    ration_data = data.frame(Species = integer(), Sex = integer(),
                             Year = integer()),
    diet_data   = data.frame(Pred = integer(), Prey = integer(),
                             Pred_sex = integer(), Prey_sex = integer(),
                             Pred_age = integer(), Prey_age = integer(),
                             Year = integer(),
                             Sample_size = numeric(),
                             Stomach_proportion_by_weight = numeric())
  )
  for(df_name in names(default_dfs)){
    if(is.null(data_list[[df_name]])){
      data_list[[df_name]] <- default_dfs[[df_name]]
    }
  }

  # env_data: default to a Year-only data.frame so downstream code that does
  # `ncol(env_data) - 1` (number of indices) gets 0 when no environmental data
  # are supplied, and `merge(env_data, ...)` in rearrange_dat() works.
  # data_check() still errors when a model feature (env-q catchability,
  # temperature-dependent consumption, env linkages, srr/M1 indices) needs
  # actual indices.
  if(is.null(data_list$env_data)){
    yrs <- if(!is.null(data_list$styr) && !is.null(data_list$projyr))
             data_list$styr:data_list$projyr else integer(0)
    data_list$env_data <- data.frame(Year = yrs)
  }

  # --- 1. Filter Data by Year ----
  # Data in likelihood (use absolute Year)
  abs_year_data <- c("index_data", "catch_data", "comp_data", "caal_data")
  for(df_name in abs_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- data_list[[df_name]] |>
        dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
    }
  }

  # Fixed data (allow Year == 0)
  fixed_year_data <- c("diet_data", "weight", "emp_sel", "NByageFixed", "ration_data")
  for(df_name in fixed_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- as.data.frame(data_list[[df_name]]) |>
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
      catch_data_sub <- data_list$catch_data |> dplyr::filter(Fleet_code == flt)

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
    data_list$diet_data <- data_list$diet_data |>
      dplyr::arrange(Pred, Pred_sex, Pred_age, Prey, Prey_sex, Prey_age, Year) |>
      dplyr::mutate(stratum_id = paste(Pred, Pred_sex, Pred_age, Year, sep = "_"),
                    stomach_id = as.numeric(as.factor(stratum_id)) - 1) |>
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
  data_list$Diet_loglike <- set_default(data_list$Diet_loglike, rep(0, data_list$nspp), "'Diet_loglike' are not included in data, assuming 'Multinomial'")
  data_list$alpha_wt_len <- set_default(data_list$alpha_wt_len, 1e-6, "'alpha_wt_len' not specified in data, assuming 1e-6")
  data_list$beta_wt_len <- set_default(data_list$beta_wt_len, 3, "'beta_wt_len' not specified in data, assuming 3")
  data_list$M1_model <- set_default(data_list$M1_model, rep(0, data_list$nspp), "'M1_model' is not included in data, assuming 0")
  data_list$msmMode <- set_default(data_list$msmMode, 0, "'msmMode' is not included in data, assuming single-species (0)")
  data_list$M1_re <- set_default(data_list$M1_re, rep(0, data_list$nspp), "'M1_re' is not in data, assuming 0 for all species")
  data_list$initMode <- set_default(data_list$initMode, 2, "'initMode' is not in the data, setting to 2 (default)")

  # Bioenergetics scalars: TMB declares them as DATA_VECTOR length-nspp, so
  # they must exist even when not used. In single-species mode (msmMode == 0)
  # they are never read by the consumption code, so silently fill any missing
  # entries with safe sentinels. When msmMode > 0 we leave them untouched and
  # let data_check() report which ones are missing or wrong-length.
  if(data_list$msmMode == 0){
    bioenergetics_defaults <- list(
      Ceq    = rep(1L,  data_list$nspp),
      Cindex = rep(0L,  data_list$nspp),
      Pvalue = rep(1,   data_list$nspp),
      fday   = rep(365, data_list$nspp),
      CA     = rep(0,   data_list$nspp),
      CB     = rep(0,   data_list$nspp),
      Qc     = rep(0,   data_list$nspp),
      Tco    = rep(0,   data_list$nspp),
      Tcm    = rep(0,   data_list$nspp),
      Tcl    = rep(0,   data_list$nspp),
      CK1    = rep(0,   data_list$nspp),
      CK4    = rep(0,   data_list$nspp)
    )
    for(nm in names(bioenergetics_defaults)){
      if(is.null(data_list[[nm]])){
        data_list[[nm]] <- bioenergetics_defaults[[nm]]
      }
    }
  }

  # 1. Fleet Control defaults ----
  data_list$fleet_control$Sel_norm_bin1 <- set_default(data_list$fleet_control$Sel_norm_bin1, NA, "'Sel_norm_bin1' not specified in 'fleet_control', assuming 'NA'")
  data_list$fleet_control$Sel_norm_bin2 <- set_default(data_list$fleet_control$Sel_norm_bin2, NA, "'Sel_norm_bin2' not specified in 'fleet_control', assuming 'NA'")
  data_list$fleet_control$Sel_curve_pen1 <- set_default(data_list$fleet_control$Sel_curve_pen1, 0, "'Sel_curve_pen1' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Sel_curve_pen2 <- set_default(data_list$fleet_control$Sel_curve_pen2, 0, "'Sel_curve_pen2' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Selectivity_dimension <- set_default(data_list$fleet_control$Selectivity_dimension, "Age", "'Selectivity_dimension' not specified in 'fleet_control', assuming 'Age'")
  data_list$fleet_control$Comp_loglike <- set_default(data_list$fleet_control$Comp_loglike, "MultinomialAFSC", "'Comp_loglike' not specified in 'fleet_control', assuming 'MultinomialAFSC'")
  data_list$fleet_control$CAAL_loglike <- set_default(data_list$fleet_control$CAAL_loglike, "Multinomial", "'CAAL_loglike' not specified in 'fleet_control', assuming 'Multinomial'")
  data_list$fleet_control$CAAL_weights <- set_default(data_list$fleet_control$CAAL_weights, 1, "'CAAL_weights' not specified in 'fleet_control', assuming 1")
  data_list$fleet_control$Month <- set_default(data_list$fleet_control$Month, 0, "'Month' not specified in 'fleet_control', assuming 0")

  # Format adjustment for NonParametric
  np_idx <- data_list$fleet_control$Selectivity %in% c(2, "NonParametric", "Non-parametric")
  if(any(np_idx & !is.na(data_list$fleet_control$Time_varying_sel) & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1, "Off", "IID")))){
    data_list$fleet_control <- data_list$fleet_control |>
      dplyr::mutate(
        Sel_curve_pen1 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID")), Time_varying_sel, Sel_curve_pen1),
        Sel_curve_pen2 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID")), Time_varying_sel_sd_prior, Sel_curve_pen2),
        Time_varying_sel = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID")), 0, Time_varying_sel),
        Time_varying_sel_sd_prior = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID")), 0, Time_varying_sel_sd_prior)
      )
    message("Updating format where 'Selectivity == Non-parametric'. Moving non-parametric penalties to 'Sel_curve_pen1' and 'Sel_curve_pen2'.")
  }

  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen1))) stop("'Sel_curve_pen1' is NA in 'fleet_control' for fleet with non-parametric selectivity")
  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen2))) stop("'Sel_curve_pen2' is NA in 'fleet_control' for fleet with non-parametric selectivity")


  # 2. Sel bins ----
  for(flt in 1:nrow(data_list$fleet_control)){

    sp_idx <- data_list$fleet_control$Species[flt]
    age_selex = data_list$fleet_control$Selectivity_dimension[flt] == "Age"
    selex_text <- ifelse(age_selex, "nages", "nlengths")
    max_bin <- ifelse(age_selex,
                      data_list$nages[sp_idx],
                      data_list$nlengths[sp_idx])

    # - Sel normalization bin
    if(any(data_list$fleet_control$Sel_norm_bin1[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$Sel_norm_bin1[flt] <- max_bin
      message(paste0("'Sel_norm_bin1' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }

    # - Upper sel normalization bin
    if(any(data_list$fleet_control$Sel_norm_bin2[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$Sel_norm_bin2[flt] <- max_bin
      message(paste0("'Sel_norm_bin2' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }

    # - N bins
    if(any(data_list$fleet_control$N_sel_bins[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$N_sel_bins[flt] <- max_bin
      message(paste0("'N_sel_bins' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }
  }

  # Update switches to text-based if necessary for older files
  data_list <- revert_switches(data_list)

  # Auto-Off fleets with no observations contributing to the likelihood.
  # catch_data, index_data, comp_data,or caal_data should have
  # at least one row with Year > 0 referencing its
  # Fleet_code.
  active_codes <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(integer(0))
    if (!all(c("Fleet_code", "Year") %in% colnames(df))) return(integer(0))
    df <- df[!is.na(df$Year) & df$Year > 0 & !is.na(df$Fleet_code), , drop = FALSE]
    if (nrow(df) == 0) return(integer(0))
    unique(as.integer(df$Fleet_code))
  }
  active <- unique(c(
    active_codes(data_list$catch_data),
    active_codes(data_list$index_data),
    active_codes(data_list$comp_data),
    active_codes(data_list$caal_data)
  ))
  fc <- data_list$fleet_control
  inactive_idx <- which(!fc$Fleet_code %in% active & fc$Fleet_type != "Off")
  if (length(inactive_idx) > 0) {
    message(paste0(
      "Auto-Off fleet(s) with no observations in the likelihood (all Year < 0 or missing): ",
      paste(fc$Fleet_name[inactive_idx], collapse = ", "),
      ". Setting Fleet_type = 'Off'."
    ))
    data_list$fleet_control$Fleet_type[inactive_idx]   <- "Off"
    data_list$fleet_control$proj_F_prop[inactive_idx]  <- NA_real_
    data_list$fleet_control$Catchability[inactive_idx] <- NA
  }

  return(data_list)
}


#' Convert integer switches to intuitive text strings. Maintains backwards compatability.
#'
#' @param data_list Rceattle data list
#' @importFrom rlang .data
revert_switches <- function(data_list) {

  # Helper: convert integer codes to map names. Numeric-looking strings ("0",
  # " 0", "0.000000000") are canonicalized via as.numeric so they match the
  # integer keys in `map`. Strings that are already a map name (or anything
  # else) pass through unchanged.
  .conv <- function(x, map) {
    x_char <- gsub(" ", "", as.character(x))
    x_num <- suppressWarnings(as.numeric(x_char))
    is_num <- !is.na(x_num)
    x_char[is_num] <- as.character(x_num[is_num])
    idx <- match(x_char, as.character(map))
    matched <- !is.na(idx)
    if (any(matched)) {
      x_char[matched] <- names(map)[idx[matched]]
    }
    return(x_char)
  }

  # - Fleet switches
  data_list$fleet_control <- data_list$fleet_control |>
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Catchability = .conv(.data$Catchability, q_map),
      Comp_loglike = .conv(.data$Comp_loglike, comp_loglike_map),
      CAAL_loglike = .conv(.data$CAAL_loglike, comp_loglike_map)
    )

  # Time_varying_q doubles as an environmental-index column when Catchability
  # is "AR1" or "Environmental", so only convert the rows that hold a switch.
  non_env_idx <- !data_list$fleet_control$Catchability %in% c("AR1", "Environmental")
  if (any(non_env_idx)) {
    data_list$fleet_control$Time_varying_q[non_env_idx] <-
      .conv(data_list$fleet_control$Time_varying_q[non_env_idx], tv_q_map)
  }

  # - Population dynamics switches
  data_list$initMode <- .conv(data_list$initMode, initMode_map)
  data_list$HCR <- .conv(data_list$HCR, hcr_map)

  return(data_list)
}

