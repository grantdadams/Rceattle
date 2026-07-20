#' Function to check data for errors. Does not update the data set!
#'
#' @param data_list Rceattle data list
#'
#' @keywords internal
data_check <- function(data_list) {
  errors <- character(0)

  # ---- Helpers ----
  has_data <- function(df) !is.null(df) && nrow(df) > 0
  fc_num <- function(fc, col, flt){
    if(!col %in% colnames(fc)) return(NA_real_)
    val <- suppressWarnings(as.numeric(fc[[col]][flt]))
    if(length(val) == 0) NA_real_ else val
  }

  # =======================================================================
  # 1. Top-level scalars and run-time switches ----
  # =======================================================================

  # msmMode: Kinzey & Punt (2009) functional responses are deprecated.
  if(!is.null(data_list$msmMode) && data_list$msmMode %in% 3:9){
    errors <- c(errors, paste0(
      "msmMode = ", data_list$msmMode, " selects a Kinzey & Punt (2009) ",
      "predation formulation (Holling I/II/III, predator interference, ",
      "preemption, Hassell-Varley, or Ecosim). These are deprecated and ",
      "have not been validated in the current code base. Use msmMode = 0 ",
      "(single-species), 1 (Holsman et al. 2015 MSVPA), or 2 (Holling ",
      "Type III MSVPA)."
    ))
  }

  # suitMode: length-based modes are not yet implemented.
  # 1 = GammaLength, 3 = LognormalLength, 5 = NormalLength.
  if(any(data_list$suitMode %in% c(1, 3, 5))){
    errors <- c(errors, "Length based suitability not yet implemented")
  }

  # Catchability = "PowerEquation" is not yet implemented: the power coefficient
  # (index_q_pow) is not built as a parameter and the template does not apply it,
  # so the fleet would silently get a plain estimated q instead.
  if(!is.null(data_list$fleet_control$Catchability) &&
     any(data_list$fleet_control$Catchability %in% c("PowerEquation", 4), na.rm = TRUE)){
    errors <- c(errors, "'PowerEquation' catchability not yet implemented")
  }

  # Catchability = "Environmental" (Estimate_q = 5) is superseded by a q
  # linkage. The old form still fits -- it keeps its own C++ path and exact
  # numerics -- but names covariates by position in Time_varying_q rather than
  # by a formula, and cannot carry priors, bounds or a phase. Redirect the
  # user without changing their result.
  if(!is.null(data_list$fleet_control$Catchability) &&
     any(data_list$fleet_control$Catchability %in% c("Environmental", 5),
         na.rm = TRUE)){
    env_flts <- data_list$fleet_control$Fleet_name[
      data_list$fleet_control$Catchability %in% c("Environmental", 5)]
    warning(paste0(
      "Catchability = 'Environmental' is deprecated. Fleet(s) ",
      paste(env_flts, collapse = ", "), " still fit with their current ",
      "numerics, but the environmental effect on q is better expressed as a ",
      "linkage:\n  build_catchability(linkages = list(q = ",
      "linkage_spec(~ your_covariate, by = ~ fleet)))\n",
      "which names the covariate by a formula and can carry priors, bounds ",
      "and an estimation phase."), call. = FALSE)
  }

  # minage: < 0 error
  if(any(data_list$minage < 0)){
    errors <- c(errors, "Minimum age is < 0. Check 'minage'.")
  }

  # Year ordering
  if(!is.null(data_list$styr) && !is.null(data_list$endyr) && data_list$styr > data_list$endyr){
    errors <- c(errors, paste0("styr (", data_list$styr, ") must be <= endyr (", data_list$endyr, ")"))
  }
  if(!is.null(data_list$projyr) && !is.null(data_list$endyr) && data_list$projyr < data_list$endyr){
    errors <- c(errors, paste0("projyr (", data_list$projyr, ") must be >= endyr (", data_list$endyr, ")"))
  }
  if(!is.null(data_list$suit_styr) && !is.null(data_list$suit_endyr) && data_list$suit_styr > data_list$suit_endyr){
    errors <- c(errors, paste0("suit_styr (", data_list$suit_styr, ") must be <= suit_endyr (", data_list$suit_endyr, ")"))
  }

  # =======================================================================
  # 2. Per-species dimensions ----
  # =======================================================================

  if(data_list$nspp != max(data_list$weight$Species)){
    errors <- c(errors, "`nspp` does not match the number of species in the weight data. Check `nspp` or `weight`")
  }
  if(length(data_list$spnames)     != data_list$nspp) errors <- c(errors, "species names not included for all species")
  if(length(data_list$spawn_month) != data_list$nspp) errors <- c(errors, "'spawn_month' not included for all species")
  if(length(data_list$nages)       != data_list$nspp) errors <- c(errors, "'nages' not included for all species")
  if(length(data_list$nlengths)    != data_list$nspp) errors <- c(errors, "'nlengths' not included for all species")
  if(length(data_list$other_food)  != data_list$nspp) errors <- c(errors, "'other_food' not included for all species")

  # spawn_month is used as a fraction of year (spawn_month/12 in TMB), so 0 is
  # a valid sentinel meaning "spawn at start of year".
  if(!is.null(data_list$spawn_month) &&
     any(data_list$spawn_month < 0 | data_list$spawn_month > 12, na.rm = TRUE)){
    errors <- c(errors, "spawn_month values must be in 0:12")
  }

  # =======================================================================
  # 3. Biology: weight, maturity, sex_ratio, M1, ration, ALK, age error, ----
  #             NByageFixed, bioenergetics
  # =======================================================================

  if(length(data_list$M1_base) == 1)    errors <- c(errors, "M1 is a single value, please make it age/species specific")
  if(sum(data_list$other_food < 0) > 0) errors <- c(errors, "Other food for one species is negative")

  # Weight: year coverage (only when time-varying)
  wt_yr <- data_list$weight |>
    dplyr::group_by(Wt_index, Sex) |>
    dplyr::distinct(Year) |>
    dplyr::mutate(Tmp_ind = paste0("index = ", Wt_index, " & sex = ", Sex))
  for(ind in unique(wt_yr$Tmp_ind)){
    tmp_wt <- wt_yr |> dplyr::filter(Tmp_ind == ind) |> dplyr::distinct(Year) |> dplyr::pull(Year)
    if(length(tmp_wt) > 1 && any(!(data_list$styr:data_list$endyr) %in% tmp_wt)){
      errors <- c(errors, paste0("Weight data for ", ind, " does not span all hindcast years"))
    }
  }

  # Weight: pop_wt_index / ssb_wt_index reference valid Wt_index values
  wt_index <- data_list$weight |> dplyr::distinct(Wt_index, Species, Sex)
  if(any(!data_list$pop_wt_index %in% wt_index$Wt_index)){
    errors <- c(errors, "Check population weight index, not in weight file")
  }
  if(any(!data_list$ssb_wt_index %in% wt_index$Wt_index)){
    errors <- c(errors, "Check SSB weight index, not in weight file")
  }

  # Weight: Wt_index must be unique per species
  for(sp in 1:data_list$nspp){
    wt_no  <- data_list$weight |> dplyr::filter(Species != sp) |> dplyr::distinct(Wt_index) |> dplyr::pull(Wt_index)
    wt_yes <- data_list$weight |> dplyr::filter(Species == sp) |> dplyr::distinct(Wt_index) |> dplyr::pull(Wt_index)
    if(any(wt_no %in% wt_yes)){
      errors <- c(errors, "Check weight indices (Wt_index), the same weight index was used for multiple species")
    }
  }

  # Weight / maturity / sex_ratio age coverage
  if(any(data_list$weight |>
         dplyr::select(-c(Wt_name, Wt_index, Species, Sex, Year)) |>
         ncol() < data_list$nages)){
    errors <- c(errors, "Weight data does not span range of ages")
  }
  if(ncol(data_list$maturity)  <= max(data_list$nages)) errors <- c(errors, "Maturity-at-age (maturity) does not span all ages")
  if(ncol(data_list$sex_ratio) <= max(data_list$nages)) errors <- c(errors, "Sex ratio does not span all ages")

  # Maturity / sex_ratio value ranges
  mat_vals <- data_list$maturity[, grep("^Age", colnames(data_list$maturity)), drop = FALSE]
  if(any(mat_vals < 0, na.rm = TRUE)){ # | mat_vals > 1
    errors <- c(errors, "maturity values must be > 0")
  }
  sr_vals <- data_list$sex_ratio[, grep("^Age", colnames(data_list$sex_ratio)), drop = FALSE]
  if(any(sr_vals < 0 | sr_vals > 1, na.rm = TRUE)){
    errors <- c(errors, "sex_ratio values must be in [0, 1]")
  }

  # Sex consistency: max Sex in M1_base / weight / ration_data must be <= nsex
  m1_sex <- data_list$M1_base |> dplyr::group_by(Species) |>
    dplyr::summarise(max_sex = max(Sex)) |> dplyr::arrange(Species)
  if(any(m1_sex$max_sex > data_list$nsex)){
    errors <- c(errors, "'M1_base' has more sexes than specified in 'nsex'")
  }
  wt_sex <- data_list$weight |> dplyr::group_by(Species) |>
    dplyr::summarise(max_sex = max(Sex)) |> dplyr::arrange(Species)
  if(any(wt_sex$max_sex > data_list$nsex)){
    errors <- c(errors, "'weight' has more sexes than specified in 'nsex'")
  }

  # ration_data
  if(has_data(data_list$ration_data)){
    if(any(data_list$ration_data |>
           dplyr::select(-c(Species, Sex, Year)) |>
           ncol() < data_list$nages)){
      errors <- c(errors, "'ration_data' data does not span range of ages")
    }
    ration_sex <- data_list$ration_data |> dplyr::group_by(Species) |>
      dplyr::summarise(max_sex = max(Sex)) |> dplyr::arrange(Species)
    if(any(ration_sex$max_sex > data_list$nsex)){
      errors <- c(errors, "'ration_data' has more sexes than specified in 'nsex'")
    }
  } else if(data_list$msmMode > 0){
    message("No ration data")
  }

  # Age transition matrix: length coverage when growth is fixed and length data are used
  if(any(data_list$growth_model == 0) &
     has_data(data_list$comp_data |> dplyr::filter(Age0_Length1 == 1)) &
     any(data_list$age_trans_matrix |>
         dplyr::select(-c(Age_transition_name, Age_transition_index, Species, Sex, Age)) |>
         ncol() < data_list$nlengths)){
    errors <- c(errors, "`age_trans_matrix` data does not span range of lengths")
  }

  # Age error matrix: observed-age column count
  if(any(data_list$age_error |>
         as.data.frame() |>
         dplyr::select(-c(Species, True_age)) |>
         ncol() < data_list$nages)){
    errors <- c(errors, "`age_error` observed ages do not span range of ages")
  }
  # ALK & age_error: per-species age coverage (fillable with 0s downstream — message-level)
  for(sp in 1:data_list$nspp){
    expected_ages <- data_list$minage[sp]:(data_list$minage[sp] + data_list$nages[sp] - 1)

    atm_ages <- data_list$age_trans_matrix |> as.data.frame() |>
      dplyr::filter(Species == sp) |> dplyr::pull(Age)
    if(!all(expected_ages %in% atm_ages)){
      message(paste("`age_trans_matrix` data does not span range of age for species", sp, "will fill with 0s"))
    }

    ae_ages <- data_list$age_error |> as.data.frame() |>
      dplyr::filter(Species == sp) |> dplyr::pull(True_age)
    if(!all(expected_ages %in% ae_ages)){
      message(paste("`age_error` data does not span range of true ages for species", sp, "will fill with 0s"))
    }
  }
  # ALK / age_error: row sums (warning — rearrange may renormalize)
  if(has_data(data_list$age_error)){
    ae_cols <- setdiff(colnames(data_list$age_error), c("Species", "True_age"))
    if(length(ae_cols) > 0){
      ae_sums <- rowSums(data_list$age_error[, ae_cols, drop = FALSE], na.rm = TRUE)
      if(any(ae_sums > 0 & abs(ae_sums - 1) > 1e-3)){
        warning("Some `age_error` rows do not sum to 1")
      }
    }
  }
  if(has_data(data_list$age_trans_matrix) &
     has_data(data_list$comp_data |> dplyr::filter(Age0_Length1 == 1))
  ){
    atm_cols <- setdiff(colnames(data_list$age_trans_matrix),
                        c("Age_transition_name", "Age_transition_index", "Species", "Sex", "Age"))
    if(length(atm_cols) > 0){
      atm_sums <- rowSums(data_list$age_trans_matrix[, atm_cols, drop = FALSE], na.rm = TRUE)
      if(any(atm_sums > 0 & abs(atm_sums - 1) > 1e-3)){
        warning("Some `age_trans_matrix` rows do not sum to 1")
      }
    }
  }

  # NByageFixed: required when estDynamics > 0
  if(any(data_list$estDynamics > 0)){
    if(!has_data(data_list$NByageFixed)){
      errors <- c(errors, "estDynamics > 0 requires NByageFixed data to be provided")
    } else {
      expected_cols <- 4 + max(data_list$nages)
      if(ncol(data_list$NByageFixed) != expected_cols){
        errors <- c(errors, paste0("NByageFixed should have ", expected_cols,
                                   " columns (Species_name, Species, Sex, Year, Age1...Age",
                                   max(data_list$nages), "), but has ", ncol(data_list$NByageFixed)))
      }
    }
  }

  # Bioenergetics: temperature-dependent consumption requires environmental driver
  if(!is.null(data_list$Ceq)){
    for(sp in 1:data_list$nspp){
      if(data_list$Ceq[sp] > 1){
        if(is.null(data_list$env_data) || ncol(data_list$env_data) < (data_list$Cindex[sp] + 1)){
          errors <- c(errors, paste0("Species ", sp, " uses temperature-dependent consumption (Ceq = ",
                                     data_list$Ceq[sp], ") but Cindex (", data_list$Cindex[sp],
                                     ") exceeds available environmental indices"))
        }
      }
    }
  }

  # =======================================================================
  # 4. Fleet control (structure, references, per-fleet bin/q checks) ----
  # =======================================================================
  if(!has_data(data_list$fleet_control)) {
    stop("Missing 'fleet_control' from data.")
  } else{
    fc <- data_list$fleet_control

    # Fleet_code uniqueness; Species in 1:nspp
    fcodes <- if("Fleet_code" %in% colnames(fc)) suppressWarnings(as.numeric(fc$Fleet_code)) else integer(0)
    fsp    <- if("Species"    %in% colnames(fc)) suppressWarnings(as.numeric(fc$Species))    else integer(0)
    if(any(duplicated(fcodes[!is.na(fcodes)]))){
      errors <- c(errors, "fleet_control$Fleet_code values must be unique")
    }
    # Fleet_code is used directly as the fleet slot of the per-fleet parameter
    # and map arrays, which are built in fleet_control row order. The two must
    # therefore agree, or parameters are silently attached to the wrong fleet.
    if(length(fcodes) == nrow(fc) && !identical(as.integer(fcodes), seq_len(nrow(fc)))){
      errors <- c(errors, paste0(
        "fleet_control$Fleet_code must equal the row number (1, 2, ... ", nrow(fc),
        "); got ", paste(fcodes, collapse = ", "),
        ". Fleet_code indexes the per-fleet parameter/map arrays, which are built in row order."))
    }
    if(any(!is.na(fsp) & (fsp < 1 | fsp > data_list$nspp))){
      errors <- c(errors, paste0("fleet_control$Species values must be in 1:", data_list$nspp))
    }

    # Month range (0 = unspecified sentinel)
    if(!is.null(fc$Month)){
      fm <- suppressWarnings(as.numeric(fc$Month))
      if(any(!is.na(fm) & (fm < 0 | fm > 12))){
        errors <- c(errors, "fleet_control$Month values must be in 0:12")
      }
    }

    # proj_F_prop: must be 0 (NoFishing) or sum to 1 across fleets
    flt_spp <- unique(fc$Species)
    for(sp in flt_spp){
      fc_spp <- fc |> dplyr::filter(Species == sp &
                                      Fleet_type %in% c(1, "Fishery"))
      total_proj_f <- sum(fc_spp$proj_F_prop, na.rm = TRUE)

      if(total_proj_f == 0 && !(data_list$HCR %in% c(0, "NoFishing"))){
        errors <- c(errors, "HCR is > 0 and 'proj_F_prop' is 0")
      }

      if(total_proj_f > 0 && abs(total_proj_f - 1) > 1e-6){
        errors <- c(errors, paste0("proj_F_prop values should sum to 1 for species ", sp, ", but sum to ", total_proj_f))
      }
    }

    # Per-fleet checks
    for(flt in 1:nrow(fc)){
      sp_idx <- fc_num(fc, "Species", flt)
      if(is.na(sp_idx) || sp_idx < 1 || sp_idx > data_list$nspp) next

      flt_name   <- fc$Fleet_name[flt]
      dim_is_age <- !isTRUE(fc$Selectivity_dimension[flt] == "Length")
      max_bin    <- if(dim_is_age) data_list$nages[sp_idx] else data_list$nlengths[sp_idx]

      # Bin_first_selected
      bfs <- fc_num(fc, "Bin_first_selected", flt)
      if(!is.na(bfs) && (bfs < 1 || bfs > max_bin)){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Bin_first_selected (", bfs, ") must be in 1:", max_bin))
      }

      # N_sel_bins
      nsb <- fc_num(fc, "N_sel_bins", flt)
      if(!is.na(nsb) && (nsb < 1 || nsb > max_bin)){
        errors <- c(errors, paste0("Fleet '", flt_name, "': N_sel_bins (", nsb, ") must be in 1:", max_bin))
      }

      # Non-parametric shape-penalty range and cap, given on the fleet's own
      # selectivity dimension: an age (from minage) for age-based fleets, a
      # 1-based length-bin ordinal for length-based. Out-of-range values would
      # index past the selectivity array in the template.
      bin_lo <- if(dim_is_age) data_list$minage[sp_idx] else 1L
      bin_hi <- bin_lo + max_bin - 1L
      for(col in c("Sel_pen_first_bin", "Sel_pen_last_bin", "Sel_cap_bin")){
        val <- fc_num(fc, col, flt)
        if(!is.na(val) && (val < bin_lo || val > bin_hi)){
          errors <- c(errors, paste0("Fleet '", flt_name, "': ", col, " (", val, ") must be in ",
                                     bin_lo, ":", bin_hi,
                                     if(dim_is_age) " (age-based selectivity)" else " (length-based selectivity)"))
        }
      }
      pf <- fc_num(fc, "Sel_pen_first_bin", flt); pl <- fc_num(fc, "Sel_pen_last_bin", flt)
      if(!is.na(pf) && !is.na(pl) && pf > pl){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Sel_pen_first_bin (", pf,
                                   ") must be <= Sel_pen_last_bin (", pl, ")"))
      }

      # Accumulation ages
      lo <- fc_num(fc, "Accumulation_age_lower", flt)
      hi <- fc_num(fc, "Accumulation_age_upper", flt)
      if(!is.na(lo) && (lo < 1 || lo > max_bin)){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Accumulation_age_lower (", lo, ") must be in 1:", max_bin))
      }
      if(!is.na(hi) && (hi < 1 || hi > max_bin)){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Accumulation_age_upper (", hi, ") must be in 1:", max_bin))
      }
      if(!is.na(lo) && !is.na(hi) && lo > hi){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Accumulation_age_lower (", lo,
                                   ") must be <= Accumulation_age_upper (", hi, ")"))
      }

      # AR1 catchability: Time_varying_q must be a valid env_data column index (1..ncol-1)
      if(!is.na(fc$Catchability[flt]) && fc$Catchability[flt] == "AR1"){
        if((fc$Time_varying_q[flt] > (ncol(data_list$env_data) - 1) ||
            fc$Time_varying_q[flt] < 1)){
          errors <- c(errors, "For 'AR1' catchability, environmental index specified in 'Time_varying_q' is greater than number of indices in 'env_data'")
        }
      }

      # Time-varying form is selectivity-type specific:
      #  - NonParametric (Ianelli, type 2) / NonParametricPM (type 9): the cpp
      #    penalizes year-to-year log selectivity-at-age, i.e. a RANDOM WALK ->
      #    allow only "Off"/"RandomWalk".
      #  - Hake (Taylor, type 5): IID coefficient deviates -> allow only "Off"/"IID".
      if(fc$Selectivity[flt] %in% c("NonParametric", "NonParametricPM") &&
         !fc$Time_varying_sel[flt] %in% c("Off", "RandomWalk")){
        errors <- c(errors, "For 'NonParametric'/'NonParametricPM' selectivity, 'Time_varying_sel' must be 'Off' or 'RandomWalk'")
      }
      if(fc$Selectivity[flt] == "Hake" &&
         !fc$Time_varying_sel[flt] %in% c("Off", "IID")){
        errors <- c(errors, "For 'Hake' selectivity, 'Time_varying_sel' must be 'Off' or 'IID'")
      }
      #  - LogisticPM (ADMB AMAK "pm" BTS, type 11): random-walk deviates on
      #    slope/inflection/age-1 -> allow only "Off"/"RandomWalk".
      if(fc$Selectivity[flt] == "LogisticPM" &&
         !fc$Time_varying_sel[flt] %in% c("Off", "RandomWalk")){
        errors <- c(errors, "For 'LogisticPM' selectivity, 'Time_varying_sel' must be 'Off' or 'RandomWalk'")
      }

      # Non-parametric (Ianelli) selectivity penalties (Sel_curve_pen1 =
      # decreasing penalty, Sel_curve_pen2 = curvature) must be present and
      # numeric to identify the free coefficients. Catch the case where they are
      # missing / non-numeric (e.g. a Time_varying_sel mode string accidentally
      # written into Sel_curve_pen) before it surfaces as a cryptic
      # "inits not within bounds" error in build_bounds.
      if(fc$Selectivity[flt] %in% c("NonParametric", "NonParametricPM")){
        cp1 <- suppressWarnings(as.numeric(fc$Sel_curve_pen1[flt]))
        cp2 <- suppressWarnings(as.numeric(fc$Sel_curve_pen2[flt]))
        if(is.na(cp1) || is.na(cp2)){
          errors <- c(errors, paste0(
            "Fleet '", fc$Fleet_name[flt], "' has Selectivity = '", fc$Selectivity[flt], "' but ",
            "'Sel_curve_pen1'/'Sel_curve_pen2' are missing or non-numeric (got '",
            fc$Sel_curve_pen1[flt], "' / '", fc$Sel_curve_pen2[flt],
            "'). Non-parametric selectivity requires numeric curvature/smoothness penalties."))
        }
      }

      # LogisticPM: Sel_curve_pen1/2/3 are the random-walk weights on the
      # slope / inflection / age-1 deviates (ADMB 50 / 50 / 8). Require numeric
      # when time-varying so a stray mode string is caught early.
      if(fc$Selectivity[flt] == "LogisticPM" && fc$Time_varying_sel[flt] == "RandomWalk"){
        cps <- suppressWarnings(as.numeric(c(fc$Sel_curve_pen1[flt], fc$Sel_curve_pen2[flt], fc$Sel_curve_pen3[flt])))
        if(any(is.na(cps))){
          errors <- c(errors, paste0(
            "Fleet '", fc$Fleet_name[flt], "' has Selectivity = 'LogisticPM' with time-varying ",
            "selectivity but 'Sel_curve_pen1'/'Sel_curve_pen2'/'Sel_curve_pen3' (slope/inflection/age-1 ",
            "random-walk weights) are missing or non-numeric."))
        }
      }

      # 2DAR1/3DAR1: 'Sel_curve_pen1'/'Sel_curve_pen2' are reused as logit-scale AR1
      # correlation parameters via rho_trans(x) = 2/(1+exp(-2x)) - 1. |x| > ~10 saturates
      # rho at +-1 and the AR1 log-density evaluates to NaN.
      if(fc$Selectivity[flt] %in% c("2DAR1", "3DAR1")){
        for(col in c("Sel_curve_pen1", "Sel_curve_pen2")){
          val <- suppressWarnings(as.numeric(fc[[col]][flt]))
          if(!is.na(val) && abs(val) > 10){
            errors <- c(errors, sprintf(
              "Fleet '%s' has Selectivity = '%s'. For 2DAR1/3DAR1, '%s' is reused as a logit-scale AR1 correlation parameter (rho_trans maps to (-1, 1)); abs(value) > 10 saturates rho at +-1 and produces NaN likelihoods. Got %s = %s. Suggested replacement: %s = 0 (rho = 0). Use a value in roughly [-5, 5] to set non-trivial correlation.",
              flt_name, fc$Selectivity[flt], col, col, val, col
            ))
          }
        }
      }
    }

    # emp_sel required when any fleet has Selectivity = "Fixed"
    fixed_flts <- fc[!is.na(fc$Selectivity) & fc$Selectivity == "Fixed", , drop = FALSE]
    if(nrow(fixed_flts) > 0 && !has_data(data_list$emp_sel)){
      errors <- c(errors, paste0(
        "Fleet(s) with Selectivity = 'Fixed' (",
        paste(fixed_flts$Fleet_name, collapse = ", "),
        ") require emp_sel data"
      ))
    }

    # Estimated selectivity (Selectivity != "Fixed" and Fleet_type != "Off")
    # requires comp or CAAL data with Year > 0 to be identifiable. Otherwise
    # the selectivity parameters are unconstrained and the optimizer wanders.
    # EXCEPTION: a fleet whose Selectivity_index is shared (mirrored) with
    # another fleet that DOES have active comp/CAAL data is identifiable through
    # the master fleet's data, so it is not flagged.
    has_active_age_data <- function(flt_code, df) {
      if (!has_data(df) || !all(c("Fleet_code", "Year") %in% colnames(df))) return(FALSE)
      any(df$Fleet_code == flt_code & !is.na(df$Year) & df$Year > 0 & df$Sample_size > 0)
    }
    # Fleets sharing a Selectivity_index share one deviation block, so a
    # differing Sel_start_year within the group is resolved to the group minimum
    # (the earliest year any sharing fleet has data). Surface it so the
    # resolution is deliberate rather than silent.
    if (all(c("Selectivity_index", "Sel_start_year") %in% colnames(fc))) {
      for (si in unique(fc$Selectivity_index[!is.na(fc$Selectivity_index)])) {
        rows <- which(!is.na(fc$Selectivity_index) & fc$Selectivity_index == si)
        if (length(rows) > 1) {
          ys <- suppressWarnings(as.integer(fc$Sel_start_year[rows]))
          ys[is.na(ys)] <- as.integer(data_list$styr)
          if (length(unique(ys)) > 1) {
            warning(paste0("Fleets sharing Selectivity_index ", si, " (",
                           paste(fc$Fleet_name[rows], collapse = ", "),
                           ") have different Sel_start_year (", paste(ys, collapse = ", "),
                           "); the shared selectivity deviations use the earliest (", min(ys), ")."))
          }
        }
      }
    }

    est_sel_flts <- fc[!is.na(fc$Selectivity) &
                         fc$Selectivity != "Fixed" &
                         (!"Fleet_type" %in% colnames(fc) | fc$Fleet_type != "Off"),
                       , drop = FALSE]
    # Selectivity_index values that have active age data in ANY sharing fleet
    sel_idx_has_data <- if ("Selectivity_index" %in% colnames(fc)) {
      vapply(unique(fc$Selectivity_index), function(si) {
        codes <- fc$Fleet_code[!is.na(fc$Selectivity_index) & fc$Selectivity_index == si]
        any(vapply(codes, function(cc)
          has_active_age_data(cc, data_list$comp_data) ||
            has_active_age_data(cc, data_list$caal_data), logical(1)))
      }, logical(1))
    } else NULL
    if (!is.null(sel_idx_has_data)) names(sel_idx_has_data) <- as.character(unique(fc$Selectivity_index))
    if (nrow(est_sel_flts) > 0) {
      missing_age_data <- vapply(seq_len(nrow(est_sel_flts)), function(i) {
        fc_code <- est_sel_flts$Fleet_code[i]
        own_data <- has_active_age_data(fc_code, data_list$comp_data) ||
          has_active_age_data(fc_code, data_list$caal_data)
        # mirrored identifiability: another fleet sharing this Selectivity_index
        # supplies the comp/CAAL data
        si <- est_sel_flts$Selectivity_index[i]
        mirror_data <- !is.null(sel_idx_has_data) && !is.na(si) &&
          isTRUE(sel_idx_has_data[[as.character(si)]])
        !own_data && !mirror_data
      }, logical(1))
      if (any(missing_age_data)) {
        errors <- c(errors, paste0(
          "Fleet(s) with estimated Selectivity but no comp_data or caal_data ",
          "rows in the likelihood (all Year/Sample_size == 0 or missing): ",
          paste(est_sel_flts$Fleet_name[missing_age_data], collapse = ", "),
          ". Either provide composition / CAAL data, mark Selectivity = 'Fixed' ",
          "with emp_sel, or set Fleet_type = 'Off'."
        ))
      }
    }


    # Mirroring (informational)
    mirror_sel <- fc |> dplyr::group_by(Selectivity_index) |>
      dplyr::filter(dplyr::n() > 1) |> dplyr::ungroup()
    if(nrow(mirror_sel) > 0){
      message(paste0("Selectivity for ", paste(mirror_sel$Fleet_name, collapse = ", "),
                     " is mirrored with another fleet"))
    }
    mirror_q <- fc |> dplyr::filter(!is.na(Catchability)) |>
      dplyr::group_by(Q_index) |> dplyr::filter(dplyr::n() > 1) |> dplyr::ungroup()
    if(nrow(mirror_q) > 0){
      message(paste0("Catchability for ", paste(mirror_q$Fleet_name, collapse = ", "),
                     " is mirrored with another fleet"))
    }
  }

  # =======================================================================
  # 5. Observation tables: catch_data, index_data, comp_data, caal_data ----
  # =======================================================================

  # Required columns
  required_cols_by_table <- list(
    index_data = c("Fleet_code", "Year", "Observation"),
    catch_data = c("Fleet_code", "Year", "Catch"),
    comp_data  = c("Fleet_code", "Species", "Sex", "Year")
  )
  for(df_name in names(required_cols_by_table)){
    df <- data_list[[df_name]]
    if(has_data(df)){
      missing_cols <- setdiff(required_cols_by_table[[df_name]], colnames(df))
      if(length(missing_cols) > 0){
        errors <- c(errors, paste(df_name, "is missing required columns:",
                                  paste(missing_cols, collapse = ", ")))
      }
    }
  }

  # Fleet_code referential integrity
  if(has_data(data_list$fleet_control)){
    fcodes <- suppressWarnings(as.numeric(data_list$fleet_control$Fleet_code))
    valid_codes <- fcodes[!is.na(fcodes)]
    for(df_name in c("catch_data", "index_data", "comp_data", "caal_data")){
      df <- data_list[[df_name]]
      if(has_data(df) && "Fleet_code" %in% colnames(df)){
        bad <- setdiff(unique(df$Fleet_code), valid_codes)
        if(length(bad) > 0){
          errors <- c(errors, paste0(df_name, " references Fleet_code(s) not in fleet_control: ",
                                     paste(bad, collapse = ", ")))
        }
      }
    }
  }

  # Lognormal SDs must be > 0
  for(df_name in c("index_data", "catch_data")){
    df <- data_list[[df_name]]
    if(has_data(df) && "Log_sd" %in% colnames(df) && any(!is.na(df$Log_sd) & df$Log_sd <= 0 & df$Year > 0)){
      errors <- c(errors, paste0(df_name, " has 'Log_sd' <= 0; lognormal likelihood requires positive SD"))
    }
  }

  # Sample_size >= 0
  for(df_name in c("comp_data", "caal_data", "diet_data")){
    df <- data_list[[df_name]]
    if(has_data(df) && "Sample_size" %in% colnames(df) && any(df$Sample_size < 0 & df$Year > 0, na.rm = TRUE)){
      errors <- c(errors, paste0(df_name, " contains 'Sample_size' values < 0"))
    }
  }

  # Observation values
  if(has_data(data_list$index_data) && "Observation" %in% colnames(data_list$index_data) &&
     any(!is.na(data_list$index_data$Observation) & data_list$index_data$Observation <= 0 & data_list$index_data$Year > 0)){
    errors <- c(errors, "index_data$Observation must be > 0 (lognormal likelihood breaks at 0)")
  }
  # NA Catch is allowed (clean_data sets NA for projection rows)
  catch_df <- data_list$catch_data |>
    dplyr::filter(Year <= data_list$endyr)
  if(has_data(catch_df) && "Catch" %in% colnames(catch_df) &&
     any(is.na(catch_df$Catch) | catch_df$Catch < 0)){
    errors <- c(errors, "catch_data$Catch must be >= 0")
  }

  # MVN survey covariance requirement ----
  # A fleet using Index_loglike == "MVN" or "MVNORM" must supply a square,
  # symmetric variance-covariance matrix in data_list$index_cov (keyed by
  # Fleet_name or Fleet_code) whose dimension equals the number of fitted survey
  # observations for that fleet. Validated here so the requirement surfaces with a
  # clear message at the flag rather than as a cryptic error in rearrange_data().
  fc <- data_list$fleet_control
  if(has_data(fc) && "Index_loglike" %in% colnames(fc)){
    mvn_flts <- which(fc$Index_loglike %in% c("MVN", "MVNORM", 1, 2, "1", "2"))
    for(flt in mvn_flts){
      flt_name <- fc$Fleet_name[flt]
      flt_code <- fc$Fleet_code[flt]
      Sigma <- NULL
      if(!is.null(data_list$index_cov)){
        Sigma <- data_list$index_cov[[as.character(flt_name)]]
        if(is.null(Sigma)) Sigma <- data_list$index_cov[[as.character(flt_code)]]
      }
      if(is.null(Sigma)){
        errors <- c(errors, paste0("Fleet '", flt_name, "' has Index_loglike == 'MVN' but no covariance matrix ",
                                   "was found in data_list$index_cov (expected an element named '", flt_name,
                                   "' or '", flt_code, "')."))
        next
      }
      Sigma <- as.matrix(Sigma)
      n_fit <- 0
      if(has_data(data_list$index_data)){
        n_fit <- sum(data_list$index_data$Fleet_code == flt_code &
                       data_list$index_data$Year > 0 &
                       data_list$index_data$Year <= data_list$endyr &
                       data_list$index_data$Observation > 0, na.rm = TRUE)
      }
      if(nrow(Sigma) != ncol(Sigma)){
        errors <- c(errors, paste0("index_cov matrix for fleet '", flt_name, "' must be square (got ",
                                   nrow(Sigma), " x ", ncol(Sigma), ")."))
      } else if(nrow(Sigma) != n_fit){
        errors <- c(errors, paste0("index_cov matrix for fleet '", flt_name, "' is ", nrow(Sigma), " x ",
                                   ncol(Sigma), " but the fleet has ", n_fit, " fitted survey observations ",
                                   "(Year in [styr, endyr], Observation > 0). Sigma must match those rows in index_data order."))
      } else if(!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-4, check.attributes = FALSE))){
        errors <- c(errors, paste0("index_cov matrix for fleet '", flt_name, "' is not symmetric."))
      }
    }
  }

  # Duplicates per (Fleet_code, Year, Month)
  for(df_name in c("catch_data", "index_data")){
    df <- data_list[[df_name]]
    if(has_data(df) && all(c("Fleet_code", "Year", "Month") %in% colnames(df))){
      n_dup <- sum(duplicated(df[, c("Fleet_code", "Year", "Month")]))
      if(n_dup > 0){
        errors <- c(errors, paste0(df_name, " has ", n_dup,
                                   " duplicated row(s) for (Fleet_code, Year, Month)"))
      }
    }
  }

  # catch_data must span hindcast years (use 0 where no catch occurred);
  # index_data gaps are normal (biennial / triennial surveys, missed years).
  if(!is.null(data_list$styr) && !is.null(data_list$endyr) && has_data(data_list$catch_data)){
    missing_years <- setdiff(data_list$styr:data_list$endyr, unique(data_list$catch_data$Year))
    if(length(missing_years) > 0){
      errors <- c(errors, paste0("catch_data is missing data for years: ",
                                 paste(missing_years, collapse = ", ")))
    }
  }

  # Empirical growth (growth_model == 0) does NOT populate growth_matrix from
  # age_trans_matrix; the C++ pred_CAAL = N * sel * growth_matrix then sees
  # growth_matrix = 0 and pred_CAAL collapses to ~0. The multinomial NLL
  # evaluates to a constant offset with no useful gradient, and the optimizer
  # finds spurious minima. Check per-species since growth_model is per-species.
  if(has_data(data_list$caal_data) && "Species" %in% colnames(data_list$caal_data)){
    bad_sp <- integer(0)
    for(sp in seq_along(data_list$growth_model)){
      if(data_list$growth_model[sp] == 0 &&
         any(data_list$caal_data$Species == sp & data_list$caal_data$Year > 0  & data_list$caal_data$Sample_size > 0)){
        bad_sp <- c(bad_sp, sp)
      }
    }
    if(length(bad_sp) > 0){
      sp_names <- if(!is.null(data_list$spnames)) data_list$spnames[bad_sp] else as.character(bad_sp)
      errors <- c(errors, paste0(
        "Empirical growth (growth_model == 0) is incompatible with CAAL data ",
        "for species: ", paste(sp_names, collapse = ", "), ". The C++ ",
        "growth_matrix is not populated from age_trans_matrix in the empirical ",
        "branch, so pred_CAAL = 0 and the CAAL likelihood gradient is ",
        "uninformative. Either (a) switch to parametric growth for these ",
        "species via build_growth(fun = 'vonBertalanffy'), or (b) drop ",
        "caal_data rows for these species (set to empty, Sample_size = 0, or Year < 0)."
      ))
    }
  }

  # CAAL: required when growth is being estimated
  if(any(data_list$growth_model > 0)){
    if(!has_data(data_list$caal_data)){
      errors <- c(errors, "Growth estimation (growth_model > 0) requires caal_data to be provided")
    } else {
      missing_cols <- setdiff(c("Fleet_code", "Species", "Year", "Length", "Sample_size"),
                              colnames(data_list$caal_data))
      if(length(missing_cols) > 0){
        errors <- c(errors, paste("caal_data is missing required columns:",
                                  paste(missing_cols, collapse = ", ")))
      }
      for(sp in 1:data_list$nspp){
        sp_lengths <- unique(data_list$caal_data$Length[data_list$caal_data$Species == sp])
        if(length(sp_lengths) != data_list$nlengths[sp]){
          errors <- c(errors, paste0("Species ", sp, " has ", length(sp_lengths),
                                     " unique lengths in caal_data, but nlengths[", sp, "] = ",
                                     data_list$nlengths[sp]))
        }
      }
      caal_cols <- grep("^CAAL_", colnames(data_list$caal_data), value = TRUE)
      if(length(caal_cols) == 0){
        errors <- c(errors, "caal_data is missing CAAL_ columns (CAAL_1, CAAL_2, etc.)")
      } else {
        missing_caal <- setdiff(paste0("CAAL_", 1:max(data_list$nages)), caal_cols)
        if(length(missing_caal) > 0){
          errors <- c(errors, paste("caal_data is missing CAAL columns:",
                                    paste(missing_caal, collapse = ", ")))
        }
      }
    }
  }

  # =======================================================================
  # 6. Diet & predation ----
  # =======================================================================

  if(has_data(data_list$diet_data)){
    dd <- data_list$diet_data

    # Pred age bounds
    pred_max <- dd |> dplyr::group_by(Pred) |>
      dplyr::summarise(Max_age = max(Pred_age)) |> dplyr::arrange(Pred)
    if(any(pred_max$Max_age > data_list$nages)){
      errors <- c(errors, "Pred ages in 'diet_data' > 'nages'")
    }
    # Prey age bounds
    prey_max <- dd |> dplyr::group_by(Prey) |>
      dplyr::summarise(Max_age = max(Prey_age)) |> dplyr::arrange(Prey)
    if(any(prey_max$Max_age > data_list$nages)){
      errors <- c(errors, "Prey ages in 'diet_data' > 'nages'")
    }
    # Duplicates
    if(sum(duplicated(dd)) > 0){
      errors <- c(errors, "Diet data includes duplicated rows")
    }
    # Stomach proportion ranges
    if("Stomach_proportion_by_weight" %in% colnames(dd)){
      if(any(dd$Stomach_proportion_by_weight < 0, na.rm = TRUE)){
        errors <- c(errors, "diet_data contains negative Stomach_proportion_by_weight values")
      }
      diet_sum <- dd |> dplyr::group_by(Pred, Pred_age, Pred_sex, Year) |>
        dplyr::summarise(diet_sum = sum(Stomach_proportion_by_weight))
      if(any(diet_sum$diet_sum > 1)){
        errors <- c(errors, "Stomach proportion in `diet_data` for some predators-at-age/sex/year is > 1")
      }
    }
  } else if(data_list$msmMode > 0){
    errors <- c(errors, "No diet data included")
  }

  # Diet composition likelihood weights & types (only enforced when suitMode > 0)
  if(any(is.na(data_list$Diet_comp_weights) & data_list$suitMode > 0)){
    errors <- c(errors, "Diet composition likelihood weight for a species with estimated suitability is NA")
  }
  if(any(!(data_list$Diet_loglike %in% c(0, 1)) & data_list$suitMode > 0)){
    errors <- c(errors, "Diet composition likelihood for a species with estimated suitability is not 0 or 1")
  }

  # other_food = 0 causes divide-by-zero in suitability
  if(!is.null(data_list$msmMode) && data_list$msmMode > 0 &&
     !is.null(data_list$other_food) && any(data_list$other_food == 0, na.rm = TRUE)){
    errors <- c(errors, "msmMode > 0 requires other_food > 0 for all species; zero values cause divide-by-zero in suitability")
  }

  # Bioenergetics scalars: required when msmMode > 0. switch_check fills these
  # with safe sentinels in single-species mode, so the only way they're missing
  # or wrong-length here is in multispecies mode.
  if(!is.null(data_list$msmMode) && data_list$msmMode > 0){
    bio_scalars <- c("Ceq", "Cindex", "Pvalue", "fday", "CA", "CB",
                     "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4")
    bad <- bio_scalars[vapply(bio_scalars, function(nm){
      is.null(data_list[[nm]]) || length(data_list[[nm]]) != data_list$nspp
    }, logical(1))]
    if(length(bad) > 0){
      errors <- c(errors, paste0(
        "msmMode > 0 requires bioenergetics scalars of length nspp = ",
        data_list$nspp, "; missing or wrong length: ",
        paste(bad, collapse = ", ")
      ))
    }
  }

  # =======================================================================
  # 7. Environmental data ----
  # =======================================================================

  # env_data may have just a Year column (no indices) — clean_data fills that
  # default. Downstream checks (Cindex when Ceq > 1, env-q catchability,
  # srr_indices, M1_indices) report when an index is actually needed but
  # missing.
  if(has_data(data_list$env_data)){
    if(!"Year" %in% colnames(data_list$env_data)){
      errors <- c(errors, "env_data is missing required 'Year' column")
    }
  }
  if(any(data_list$M1_indices > ncol(data_list$env_index))){
    errors <- c(errors, "'M1_indices' greater than the number of indices included")
  }

  # =======================================================================
  # 8. Switches (delegated to validate_switches) ----
  # =======================================================================

  errors <- c(errors, validate_switches(data_list))

  if(length(errors) > 0){
    stop(paste(errors, collapse = "\n"))
  }
}
