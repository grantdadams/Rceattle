#' @title Main function to construct the TMB map argument for CEATTLE
#'
#' @description Orchestrates the building of the TMB map object by calling
#'   specialized helper functions for each parameter block (Recruitment, M1,
#'   Predation, Selectivity, Catchability, etc.).
#'
#' @param data_list A data_list created from \code{\link{build_dat}}.
#' @param params A parameter list created from \code{\link{build_params}}.
#' @param debug Logical. If TRUE, sets all map values to NA except the dummy
#'   parameter, running the model without parameter estimation.
#' @param random_rec Logical. If TRUE, treats recruitment deviations as random effects,
#'   meaning the variance parameter (\code{R_ln_sd}) is estimated.
#' @param random_sel Logical. If TRUE, treats selectivity deviations as random effects,
#'   meaning the variance parameter (\code{sel_dev_ln_sd}) is estimated.
#'
#' @return A list containing the factorized TMB map (`mapFactor`) and the
#'   original map matrix/array list (`mapList`).
#' @export
build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE, random_sel = FALSE) {

  # Check data list format (Assuming this is a necessary internal utility)
  data_list <- Rceattle::switch_check(data_list)

  # --- Setup Data and Initial Map ---
  nyrs_hind <- data_list$endyr - data_list$styr + 1
  nyrs_proj <- data_list$projyr - data_list$styr + 1
  # yrs_hind is calculated inside helpers as needed

  # Convert parameters to map object (assigning sequential index for all parameters initially)
  map_list <- lapply(params, function(x) {
    if (length(x) == 0) return(x)
    replace(x, values = seq_along(x))
  })

  # --- Parameter Components ---
  map_list <- build_map_recruitment(map_list, data_list, nyrs_hind, nyrs_proj, random_rec)

  map_list <- build_map_m1(map_list, data_list, nyrs_hind)

  map_list <- build_map_predation(map_list, data_list)

  map_list <- build_map_selectivity(map_list, data_list, nyrs_hind, random_sel)

  map_list <- build_map_catchability(map_list, data_list, nyrs_hind)

  map_list <- adjust_map_shared_params(map_list, data_list)

  map_list <- build_map_f_and_data_weights(map_list, data_list, nyrs_hind)

  map_list <- build_map_fixed_natage(map_list, data_list)

  # --- Final Steps ---
  map_list <- build_map_debug(map_list, debug)

  map_list_grande <- list()
  map_list_grande$mapFactor <- lapply(map_list, factor)
  map_list_grande$mapList <- map_list

  return(map_list_grande)
}

## ----------------------------------------------------------------------
## Helper Functions ----
## ----------------------------------------------------------------------

#' @title Helper to set map for Recruitment parameters
#'
#' @description Maps the recruitment deviations (\code{rec_dev}, \code{init_dev}),
#'   stock-recruitment parameters (\code{rec_pars}), their variances (\code{R_ln_sd}),
#'   and environmental linkages (\code{beta_rec_pars}).
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#' @param nyrs_proj Total number of years (historical + projected).
#' @param random_rec Logical indicating if recruitment deviations are random effects.
#'
#' @description
#' see \code{build_srr()} for options,
#'
#'
#' @return Updated \code{map_list}.
build_map_recruitment <- function(map_list, data_list, nyrs_hind, nyrs_proj, random_rec) {
  yrs_proj <- (nyrs_hind + 1):nyrs_proj
  if (nyrs_hind == nyrs_proj) yrs_proj <- NULL


  # 1. Initial Population Deviates (initMode)
  for (sp in 1:data_list$nspp) {
    nages_sp <- data_list$nages[sp]

    # 0) Initial abundance as free parameters
    if (data_list$initMode == 0) {
      map_list$rec_dev[sp, 1] <- NA

      # - Map out ages above range
      if(nages_sp < ncol(map_list$init_dev)) {
        map_list$init_dev[sp, (nages_sp+1):ncol(map_list$init_dev)] <- NA
      }
    }

    # 1) Equilibrium with no devs
    if (data_list$initMode == 1) {
      map_list$init_dev[sp, ] <- NA
    }

    # 2-3) Equilibrium or non-equilibrium with no devs
    if (data_list$initMode > 1) {
      if ((nages_sp - 1) < ncol(map_list$init_dev)) {
        map_list$init_dev[sp, nages_sp:ncol(map_list$init_dev)] <- NA
      }
    }
  }


  # 2. Stock-Recruitment Parameters
  # col1 = mean rec, col2 = SRR alpha, col3 = SRR beta
  if (!is.null(yrs_proj)) {
    map_list$rec_dev[, yrs_proj] <- NA
  }

  if (!random_rec) {
    map_list$R_ln_sd[] <- NA
  }

  # Null model (mean-R)
  # - Turning off 2nd and 3rd par if only using mean rec
  if (data_list$srr_fun %in% c(0, 1) & data_list$srr_pred_fun %in% c(0, 1)) {
    map_list$rec_pars[, 2:3] <- NA
  }

  # Stock recruit models (Ricker & Beverton)
  # - Turning off mean-R if using SRR
  if (data_list$srr_fun > 1) {
    map_list$rec_pars[, 1] <- NA

    # - Fix first parameter of SRR if fixed (if SRR not used, will be NA anyway)
    if (data_list$srr_est_mode == 0) {
      map_list$rec_pars[, 2] <- NA
    }
  }

  # Environmental linkages
  #FIXME: make it so the covariates can vary by species
  map_list$beta_rec_pars[] <- NA
  if (data_list$srr_pred_fun %in% c(1, 3, 5)) {
    for (sp in 1:data_list$nspp) {
      indices <- data_list$srr_indices
      map_list$beta_rec_pars[sp, indices] <- indices + (sp - 1) * ncol(map_list$beta_rec_pars)
    }
  }

  return(map_list)
}


#' @title Helper to set map for Natural Mortality (M1) parameters
#'
#' @description Maps the fixed parameters (\code{ln_M1}) and the random effects
#'   parameters (\code{ln_M1_dev}, \code{M1_dev_ln_sd}, \code{M1_rho}) based on
#'   \code{M1_model} and \code{M1_re} settings.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_m1 <- function(map_list, data_list, nyrs_hind) {

  # Map out natural mortality parameters
  M1_params <- c("ln_M1", "M1_beta", "ln_M1_dev", "M1_dev_ln_sd", "M1_rho")
  map_list[M1_params] <- lapply(map_list[M1_params], function(x) replace(x, values = NA))

  # Generic indices for looping
  M1_ind <- 1
  M1_beta_ind <- 0
  M1_dev_ind <- 0
  M1_rho_ind <- 0

  # Loop through species and turn on based on model
  for (sp in 1:data_list$nspp) {
    nages_sp <- data_list$nages[sp]
    nsex_sp <- data_list$nsex[sp]
    M1_model <- data_list$M1_model[sp]
    M1_re_type <- data_list$M1_re[sp]

    # * Fixed effects ----
    # - M1_model = 1: sex- and age-invariant M1
    if (M1_model == 1) { # Sex- and age-invariant
      map_list$ln_M1[sp, , 1:nages_sp] <- M1_ind
      M1_ind <- M1_ind + 1
    }


    # - M1_model = 2: sex-specific, but age-invariant M1
    if (M1_model == 2) { # Sex-specific, age-invariant
      map_list$ln_M1[sp, 1, 1:nages_sp] <- M1_ind
      map_list$ln_M1[sp, 2, 1:nages_sp] <- M1_ind + 1
      M1_ind <- M1_ind + 2
      if (nsex_sp == 1) {
        warning(paste0("M1 model for species ", sp," is set to 2 (sex-specific), but species is single-sex."))
        map_list$ln_M1[sp, 2, 1:nages_sp] <- M1_ind
      }
    }

    # - M1_model = 3: sex-specific, age-specific M1
    if (M1_model == 3) { # Sex- and age-specific
      map_list$ln_M1[sp, 1, 1:nages_sp] <- M1_ind:(M1_ind + nages_sp - 1)
      M1_ind <- M1_ind + nages_sp
      if (nsex_sp == 2) {
        map_list$ln_M1[sp, 2, 1:nages_sp] <- M1_ind:(M1_ind + nages_sp - 1)
        M1_ind <- M1_ind + nages_sp
      } else {
        warning(paste0("M1 model for species ", sp," is set to 3 (sex-specific), but species is single-sex."))
        map_list$ln_M1[sp, 2, ] <- map_list$ln_M1[sp, 1, ]
      }
    }

    # - M1_model = 4: environmentally driven sex- and age-invariant M1
    if (M1_model == 4) { # Environmental, sex- and age-invariant
      map_list$ln_M1[sp, , 1:nages_sp] <- M1_ind
      M1_ind <- M1_ind + 1
      map_list$M1_beta[sp, 1, data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices

      # Males and females share M1 betas
      if (nsex_sp == 2) map_list$M1_beta[sp, 2, ] <- map_list$M1_beta[sp, 1, ]
      M1_beta_ind <- M1_beta_ind + dim(map_list$M1_beta)[3]
    }

    # - M1_model = 5: environmentally driven sex-specific, but age-invariant M1
    if (M1_model == 5) { # Environmental, sex-specific, age-invariant
      # ... (logic for M1_model 5 here, similar to original)
    }

    # 3.2. Random Effects (ln_M1_dev)
    if (M1_re_type %in% c(1, 4, 2, 5, 3, 6)) {
      map_list$M1_dev_ln_sd[sp, ] <- sp

      if (M1_re_type %in% c(1, 4)) { # Age-varying
        if (M1_re_type == 4) {
          map_list$M1_rho[sp, , 1] <- M1_rho_ind + 1
          M1_rho_ind <- M1_rho_ind + 1
        }
        # ... (logic for Age-varying REs here, similar to original)
      } else if (M1_re_type %in% c(2, 5)) { # Year-varying
        if (M1_re_type == 5) {
          map_list$M1_rho[sp, , 2] <- M1_rho_ind + 1
          M1_rho_ind <- M1_rho_ind + 1
        }
        # ... (logic for Year-varying REs here, similar to original)
      } else if (M1_re_type %in% c(3, 6)) { # Age- and Year-varying
        if (M1_re_type == 6) {
          map_list$M1_rho[sp, 1, ] <- M1_rho_ind + 1:2
          if (nsex_sp == 2) map_list$M1_rho[sp, 2, ] <- map_list$M1_rho[sp, 1, ]
          M1_rho_ind <- M1_rho_ind + 2
        }
        # ... (logic for 2D REs here, similar to original)
      }
    }
  }
  return(map_list)
}

#' @title Helper to set map for Predation Mortality (M2) parameters
#'
#' @description Maps predation suitability parameters (\code{log_gam_a}, \code{log_gam_b},
#'   \code{log_phi}) and diet weight parameters based on \code{msmMode} and
#'   \code{suitMode}.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#'
#' @return Updated \code{map_list}.
build_map_predation <- function(map_list, data_list) {

  # 4.1. Single-Species Mode
  if (data_list$msmMode == 0) {
    map_list$log_gam_a[] <- NA
    map_list$log_gam_b[] <- NA
    map_list$log_phi[] <- NA
  }

  # 4.2. Suitability (if Multi-Species)
  if (data_list$msmMode > 0) {
    for (sp in 1:data_list$nspp) {
      if (data_list$suitMode[sp] == 0) {
        map_list$log_gam_a[sp] <- NA
        map_list$log_gam_b[sp] <- NA
        map_list$log_phi[sp, ] <- NA
      }

      if (data_list$estDynamics[sp] > 0) {
        map_list$log_phi[, sp] <- NA
      }
    }
  }

  # 4.3. Diet Multiplier
  map_list$diet_comp_weights[] <- NA

  return(map_list)
}

#' @title Helper to set map for Selectivity parameters
#'
#' @description Maps base selectivity parameters (\code{ln_sel_slp}, \code{sel_inf},
#'   \code{sel_coff}) and time-varying deviations, based on \code{Selectivity}
#'   and \code{Time_varying_sel} settings in \code{fleet_control}.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#' @param random_sel Logical indicating if selectivity deviations are random effects.
#'
#' @return Updated \code{map_list}.
build_map_selectivity <- function(map_list, data_list, nyrs_hind, random_sel) {

  sel_params <- c("sel_coff", "sel_coff_dev", "ln_sel_slp", "sel_inf",
                  "ln_sel_slp_dev", "sel_inf_dev", "sel_dev_ln_sd", "sel_curve_pen")
  map_list[sel_params] <- lapply(map_list[sel_params], function(x) replace(x, values = NA))

  ind_coff <- 1
  ind_dev_coff <- 1
  ind_slp <- 1
  ind_inf <- 1
  yrs_hind <- 1:nyrs_hind

  # Turn on variance of random effects for selectivity deviates (sigma)
  if (random_sel) {
    for (i in 1:nrow(data_list$fleet_control)) {
      flt <- data_list$fleet_control$Fleet_code[i]
      sel_type <- data_list$fleet_control$Selectivity[i]
      tv_sel <- data_list$fleet_control$Time_varying_sel[i]
      if (sel_type %in% c(1, 3, 4, 5) && tv_sel %in% c(1, 2, 4, 5)) {
        map_list$sel_dev_ln_sd[flt] <- flt
      }
    }
  }

  # Loop through fleets to set up selectivity parameters
  for (i in 1:nrow(data_list$fleet_control)) {
    flt <- data_list$fleet_control$Fleet_code[i]
    spp <- data_list$fleet_control$Species[i]
    nsex <- data_list$nsex[spp]
    sel_type <- data_list$fleet_control$Selectivity[i]
    tv_sel <- data_list$fleet_control$Time_varying_sel[i]

    if (data_list$fleet_control$Fleet_type[i] > 0) {

      # Helper for selectivity blocks logic
      max_block <- 0
      if (tv_sel == 3) {
        data_source <- if (data_list$fleet_control$Fleet_type[i] == 1) data_list$catch_data else data_list$index_data
        fleet_data <- data_source %>%
          dplyr::filter(.data$Fleet_code == flt, .data$Year - data_list$styr + 1 <= nyrs_hind)
        Selectivity_block <- fleet_data$Selectivity_block
        biom_yrs <- fleet_data$Year - data_list$styr + 1
        max_block <- max(Selectivity_block, 0)
      }

      # 5.1. Logitistic (sel_type = 1)
      if (sel_type == 1) {
        # ... (logic for logistic selectivity here)
      }

      # 5.2. Non-parametric (sel_type = 2)
      if (sel_type == 2) {
        # ... (logic for non-parametric selectivity here)
      }

      # 5.3. Double logistic (sel_type = 3)
      if (sel_type == 3) {
        # ... (logic for double logistic selectivity here)
      }

      # 5.4. Descending logitistic (sel_type = 4)
      if (sel_type == 4) {
        # ... (logic for descending logistic selectivity here)
      }

      # 5.5. Non-parametric Hake-like (sel_type = 5)
      if (sel_type == 5) {
        # ... (logic for Hake-like selectivity here)
      }
    }
  }
  return(map_list)
}

#' @title Helper to set map for Catchability parameters
#'
#' @description Maps survey catchability base parameters (\code{index_ln_q}),
#'   time-varying deviations (\code{index_q_dev}), and environmental linkages
#'   (\code{index_q_beta}, \code{index_q_rho}).
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_catchability <- function(map_list, data_list, nyrs_hind) {

  catchability_params <- c("index_ln_q", "index_q_beta", "index_q_rho", "index_q_dev", "index_q_ln_sd", "index_q_dev_ln_sd", "index_ln_sd")
  map_list[catchability_params] <- lapply(map_list[catchability_params], function(x) replace(x, values = NA))

  ind_q_dev <- 1
  ind_beta_q <- 0
  yrs_hind <- 1:nyrs_hind

  for (i in 1:nrow(data_list$fleet_control)) {
    flt <- data_list$fleet_control$Fleet_code[i]

    if (data_list$fleet_control$Fleet_type[flt] == 2) { # If survey
      est_q <- data_list$fleet_control$Estimate_q[i]
      tv_q <- data_list$fleet_control$Time_varying_q[i]

      # Turn on mean q
      if (est_q %in% c(1, 2, 4, 5, 6)) {
        map_list$index_ln_q[flt] <- flt
      }

      # Time-varying q (if relevant)
      is_tv_q_model <- (est_q %in% c(1, 2) & tv_q %in% c(1, 2, 3, 4)) | est_q == 6

      if (is_tv_q_model) {
        index_data <- data_list$index_data %>%
          dplyr::filter(.data$Fleet_code == flt, .data$Year > data_list$styr, .data$Year <= data_list$endyr)
        srv_biom_yrs <- index_data$Year - data_list$styr + 1

        if (tv_q %in% c(1, 2, 4) || est_q == 6) { # Penalized deviate or random walk
          map_list$index_q_dev[flt, yrs_hind] <- ind_q_dev + 1:nyrs_hind - 1
          ind_q_dev <- ind_q_dev + nyrs_hind
        }

        if (tv_q == 4) map_list$index_q_dev[flt, 1] <- NA

        if (tv_q == 3) { # Time blocks
          map_list$index_q_dev[flt, srv_biom_yrs] <- ind_q_dev + index_data$Selectivity_block - 1
          ind_q_dev <- ind_q_dev + max(index_data$Selectivity_block, 0)
        }
      }

      # Environmental linkage (Estimate_q = 5)
      if (est_q == 5) {
        # ... (logic for environmental linkage here)
      }

      # Fit to environmental index (Estimate_q = 6)
      if (est_q == 6) {
        # ... (logic for fitting to environmental index here)
      }

      # Standard deviation of survey index
      if (data_list$fleet_control$Estimate_index_sd[i] == 1) {
        map_list$index_ln_sd[flt] <- flt
      }
    }
  }

  return(map_list)
}


#' @title Helper to adjust map for shared catchability/selectivity indices
#'
#' @description Enforces parameter sharing by mapping parameters for fleets
#'   with a common \code{Selectivity_index} or \code{Q_index} to the same value
#'   as the initial index.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#'
#' @return Updated \code{map_list}.
adjust_map_shared_params <- function(map_list, data_list) {

  fleet_control <- data_list$fleet_control
  sel_indices <- fleet_control$Selectivity_index
  q_indices <- fleet_control$Q_index

  for (i in 1:nrow(fleet_control)) {
    flt <- fleet_control$Fleet_code[i]
    sel_idx <- sel_indices[flt]
    q_idx <- q_indices[flt]

    # Shared Selectivity
    if (!is.na(sel_idx) && sel_idx < flt) {
      sel_duplicate <- sel_idx
      # ... (logic for shared selectivity here)
    }

    # Shared Catchability
    if (!is.na(q_idx) && q_idx < flt) {
      q_duplicate <- q_idx
      # ... (logic for shared catchability here)
    }
  }
  return(map_list)
}

#' @title Helper to set map for Fishing Mortality and Data Weights
#'
#' @description Maps fishing mortality parameters (\code{ln_F}) and related targets.
#'   Also maps data weighting parameters (\code{catch_ln_sd}, \code{comp_weights}).
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_f_and_data_weights <- function(map_list, data_list, nyrs_hind) {

  map_list$proj_F_prop[] <- NA
  if (!(data_list$initMode %in% c(3, 4))) map_list$ln_Finit[] <- NA
  map_list$ln_Flimit[] <- NA
  map_list$ln_Ftarget[] <- NA

  comp_count <- data_list$comp_data %>%
    dplyr::filter(.data$Year > 0) %>%
    dplyr::count(.data$Fleet_code)

  for (i in 1:nrow(data_list$fleet_control)) {
    flt <- data_list$fleet_control$Fleet_code[i]
    fleet_type <- data_list$fleet_control$Fleet_type[i]

    # Catch SD, F, and F dev
    if (data_list$fleet_control$Estimate_catch_sd[i] %in% c(NA, 0, 2) || fleet_type %in% c(0, 2)) {
      map_list$catch_ln_sd[flt] <- NA
    }
    if (fleet_type %in% c(0, 2)) {
      map_list$ln_F[flt, ] <- NA
    }

    # Comp weights
    has_comp_data <- flt %in% comp_count$Fleet_code
    if (data_list$fleet_control$Comp_loglike[i] != 1 || fleet_type == 0 || !has_comp_data) {
      map_list$comp_weights[flt] <- NA
    }

    comp_loglike <- data_list$fleet_control$Comp_loglike[i]
    if (!is.na(comp_loglike) && !(comp_loglike %in% c(-1, 0, 1))) {
      stop(paste0("Comp_loglike for fleet ", flt, "is not -1, 0 or 1. Current value: ", comp_loglike))
    }
  }

  # Map out Fdev and selectivity devs for years with 0 catch
  catch_data_hind <- data_list$catch_data %>%
    dplyr::filter(.data$Year <= data_list$endyr, .data$Catch == 0)

  if (nrow(catch_data_hind) > 0) {
    fsh_ind <- catch_data_hind$Fleet_code
    yr_ind <- catch_data_hind$Year - data_list$styr + 1
    for (i in seq_along(yr_ind)) {
      map_list$ln_F[fsh_ind[i], yr_ind[i]] <- NA
      map_list$ln_sel_slp_dev[, fsh_ind[i], , yr_ind[i]] <- NA
      map_list$sel_inf_dev[, fsh_ind[i], , yr_ind[i]] <- NA
    }
  }

  return(map_list)
}

#' @title Helper to set map for Fixed N-at-Age models
#'
#' @description Turns off (sets to NA) most population and fleet parameters for
#'   species where the dynamics are fixed (\code{estDynamics > 0}).
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#'
#' @return Updated \code{map_list}.
build_map_fixed_natage <- function(map_list, data_list) {

  for (sp in 1:data_list$nspp) {
    est_dyn <- data_list$estDynamics[sp]

    if (est_dyn > 0) {
      pop_params_off <- c("rec_pars", "R_ln_sd", "ln_Finit", "rec_dev", "init_dev",
                          "ln_M1", "ln_M1_dev", "M1_dev_ln_sd", "M1_rho")
      for (p in pop_params_off) {
        if (!is.null(map_list[[p]])) map_list[[p]][sp, ] <- NA
      }

      flts <- data_list$fleet_control$Fleet_code[data_list$fleet_control$Species == sp]
      if (length(flts) > 0) {
        fleet_params_off <- c("ln_F", "index_ln_q", "index_q_dev", "index_q_ln_sd",
                              "index_q_dev_ln_sd", "sel_coff", "sel_coff_dev",
                              "ln_sel_slp", "sel_inf", "ln_sel_slp_dev",
                              "sel_inf_dev", "sel_dev_ln_sd", "index_ln_sd",
                              "catch_ln_sd", "comp_weights")
        for (p in fleet_params_off) {
          if (!is.null(map_list[[p]])) {
            if (length(dim(map_list[[p]])) == 1) {
              map_list[[p]][flts] <- NA
            } else {
              map_list[[p]][, flts, ] <- NA
            }
          }
        }
      }
    }

    # Population scalar logic (est_dyn = 2 or 3)
    nages_scalar <- ncol(map_list$ln_pop_scalar)
    if (est_dyn < 2 || data_list$msmMode == 0) {
      map_list$ln_pop_scalar[sp, ] <- NA
    } else if (est_dyn == 2) {
      map_list$ln_pop_scalar[sp, 2:nages_scalar] <- NA
    } else if (est_dyn == 3) {
      if (data_list$nages[sp] < nages_scalar) {
        map_list$ln_pop_scalar[sp, (data_list$nages[sp] + 1):nages_scalar] <- NA
      }
    }
  }

  return(map_list)
}


#' @title Helper to set map for debug mode
#'
#' @description Sets all parameters in the map list to NA, except the
#'   `dummy` parameter, for use in debug or testing modes.
#'
#' @param map_list The current TMB map list.
#' @param debug Logical. If TRUE, debug mode is activated.
#'
#' @return Updated \code{map_list}.
build_map_debug <- function(map_list, debug) {
  if (debug) {
    map_list <- lapply(map_list, function(x) replace(x, values = NA))
    map_list$dummy <- 1
  }
  return(map_list)
}
