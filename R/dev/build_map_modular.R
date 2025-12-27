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


  # * 1. Initial Population Deviates (initMode) ----
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


  # * 2. Stock-Recruitment Parameters ----
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

  # * 3. Environmental linkages ----
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

    # * Switches ----
    nages_sp <- data_list$nages[sp]
    nsex_sp <- data_list$nsex[sp]
    M1_model <- data_list$M1_model[sp] # Fixed effects model
    M1_re_model <- data_list$M1_re[sp] # Random effects model

    # * 1. Fixed effects ----
    # ** M1_model = 1: sex- and age-invariant M1
    if (M1_model == 1) {
      map_list$ln_M1[sp, , 1:nages_sp] <- M1_ind
      M1_ind <- M1_ind + 1
    }


    # ** M1_model = 2: sex-specific, but age-invariant M1
    if (M1_model == 2) {
      map_list$ln_M1[sp, 1, 1:nages_sp] <- M1_ind
      map_list$ln_M1[sp, 2, 1:nages_sp] <- M1_ind + 1
      M1_ind <- M1_ind + 2
      if (nsex_sp == 1) {
        warning(paste0("M1 model for species ", sp," is set to 2 (sex-specific), but species is single-sex."))
        map_list$ln_M1[sp, 2, 1:nages_sp] <- M1_ind
      }
    }

    # ** M1_model = 3: sex-specific, age-specific M1
    if (M1_model == 3) {
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

    # ** M1_model = 4: environmentally driven sex- and age-invariant M1
    if (M1_model == 4| (M1_model == 5 & nsex_sp == 1)) {

      # - Mean M
      map_list$ln_M1[sp, , 1:nages_sp] <- M1_ind
      M1_ind <- M1_ind + 1
      map_list$M1_beta[sp, 1, data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices

      # Males and females share M1 betas
      if (nsex_sp == 2){
        map_list$M1_beta[sp, 2, ] <- map_list$M1_beta[sp, 1, ]
      }
      M1_beta_ind <- M1_beta_ind + dim(map_list$M1_beta)[3]

      if(M1_model == 5){
        warning(paste0("M1 model for species ", sp," is set to 5 (sex-specific), but species is single-sex."))
      }
    }

    # ** M1_model = 5: environmentally driven sex-specific, but age-invariant M1
    if (M1_model == 5 & nsex_sp == 2) {

      # - Mean M
      map_list$ln_M1[sp, 1, 1:nages_sp] <- M1_ind # Females
      map_list$ln_M1[sp, 2, 1:nages_sp] <- M1_ind + 1 # Males
      M1_ind = M1_ind + 2

      # - Betas
      # -- Females
      map_list$M1_beta[sp,1,data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices;
      M1_beta_ind = M1_beta_ind + dim(map_list$M1_beta)[3]

      # -- Males
      map_list$M1_beta[sp,2,data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices;
      M1_beta_ind = M1_beta_ind + dim(map_list$M1_beta)[3]
    }

    # * 2. Random Effects ----
    # - M1_re = 0: No random effects (default).
    # - M1_re = 1: Random effects varies by age, but uncorrelated (IID) and constant over years.
    # - M1_re = 2: Random effects varies by year, but uncorrelated (IID) and constant over ages.
    # - M1_re = 3: Random effects varies by year and age, but uncorrelated (IID).
    # - M1_re = 4: Correlated AR1 random effects varies by age, but constant over years.
    # - M1_re = 5: Correlated AR1 random effects varies by year, but constant over ages.
    # - M1_re = 6: Correlated 2D-AR1 random effects varies by year and age.
    # "ln_M1_dev"
    # - M1_re = 1/4: Random effects varies by age (IID or AR1) and constant over years.
    if(M1_re_model %in% c(1, 4)){
      if(M1_model == 1){ # Sex-invariant
        # - Random effects
        map_list$ln_M1_dev[sp,1, 1:nages_sp,] <- M1_dev_ind + 1:nages_sp

        # Males mapped the same, if present
        if(nsex_sp == 2){
          map_list$ln_M1_dev[sp,2,,] <- map_list$ln_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + nages_sp
      }

      if(M1_model == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$ln_M1_dev[sp, 1, 1:nages_sp,] <- M1_dev_ind + 1:nages_sp
        M1_dev_ind = M1_dev_ind + nages_sp

        # -- Males
        map_list$ln_M1_dev[sp, 2, 1:nages_sp,] <- M1_dev_ind + 1:nages_sp
        M1_dev_ind = M1_dev_ind + nages_sp
      }


      # - Standard deviation (shared across sexes)
      map_list$M1_dev_ln_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(M1_re_model == 4){
        map_list$M1_rho[sp,,1] =  sp
      }
    }

    # - M1_re = 2/5: Random effects varies by year (IID or AR1) and constant over ages
    if(M1_re_model %in% c(2, 5)){
      if(M1_model == 1){ # Sex-invariant
        # - Random effects
        map_list$ln_M1_dev[sp,1,1:nages_sp, 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = nages_sp)

        # Males mapped the same, if present
        if(nsex_sp == 2){
          map_list$ln_M1_dev[sp,2,,] <- map_list$ln_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + nyrs_hind
      }

      if(nsex_sp == 2 & M1_model == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$ln_M1_dev[sp,1, 1:nages_sp, 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = nages_sp)
        M1_dev_ind = M1_dev_ind + nyrs_hind

        # -- Males
        map_list$ln_M1_dev[sp,2, 1:nages_sp, 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = nages_sp)
        M1_dev_ind = M1_dev_ind + nyrs_hind
      }

      # - Standard deviation (shared across sexes)
      map_list$M1_dev_ln_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(M1_re_model == 5){
        map_list$M1_rho[sp,,2] =  sp #FIXME: may want sex-varying?? Hard to estimate
      }
    }

    # - M1_re = 3/6: Random effects varies by age and year (IID or 2D-AR1)
    if(M1_re_model %in% c(3, 6)){
      if(M1_model == 1){ # Sex-invariant
        # - Random effects
        map_list$ln_M1_dev[sp,1,1:nages_sp, 1:nyrs_hind] <- M1_dev_ind + (1:nyrs_hind * nages_sp)

        # Males mapped the same, if present
        if(nsex_sp == 2){
          map_list$ln_M1_dev[sp,2,,] <- map_list$ln_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + (nyrs_hind * nages_sp)
      }

      if(nsex_sp == 2 & M1_model == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$ln_M1_dev[sp,1, 1:nages_sp, 1:nyrs_hind] <- M1_dev_ind + (1:nyrs_hind * nages_sp)
        M1_dev_ind = M1_dev_ind + (nyrs_hind * nages_sp)

        # -- Males
        map_list$ln_M1_dev[sp,2, 1:nages_sp, 1:nyrs_hind] <- M1_dev_ind + (1:nyrs_hind * nages_sp)
        M1_dev_ind = M1_dev_ind + (nyrs_hind * nages_sp)
      }

      # - Standard deviation (shared across sexes)
      map_list$M1_dev_ln_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(M1_re_model == 5){
        map_list$M1_rho[sp,1,] = M1_dev_ln_sd_ind + 1:2
        map_list$M1_rho[sp,2,] = map_list$M1_rho[sp,1,]
        M1_dev_ln_sd_ind = M1_dev_ln_sd_ind + 2  #FIXME: may want sex-varying?? Hard to estimate
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

  # * 1. Single-Species Mode ----
  if (data_list$msmMode == 0) {
    map_list$log_gam_a[] <- NA
    map_list$log_gam_b[] <- NA
    map_list$log_phi[] <- NA

    # # Multispecies kinzey parameters
    # map_list$logH_1 <- map_list$logH_1 * NA
    # map_list$logH_1a <- map_list$logH_1a * NA
    # map_list$logH_1b <- map_list$logH_1b * NA
    #
    # map_list$logH_2 <- map_list$logH_2 * NA
    # map_list$logH_3 <- map_list$logH_3 * NA
    # map_list$H_4 <- map_list$H_4 * NA
  }

  # # * 2. MSVPA form ----
  # # Turn off all functional form parameters
  # if (data_list$msmMode %in% c(1,2)) {
  #
  #   # Multispecies kinzey parameters
  #   map_list$logH_1 <- map_list$logH_1 * NA
  #   map_list$logH_1a <- map_list$logH_1a * NA
  #   map_list$logH_1b <- map_list$logH_1b * NA
  #
  #   map_list$logH_2 <- map_list$logH_2 * NA
  #   map_list$logH_3 <- map_list$logH_3 * NA
  #   map_list$H_4 <- map_list$H_4 * NA
  #
  # }
  #
  # # * 3. Kinzey and Punt predation equations ----
  # if (data_list$msmMode > 2) {
  #   # Holling Type 1
  #   if (data_list$msmMode == 3) {
  #     map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
  #     map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
  #     map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
  #   }
  #
  #   # Holling Type 2
  #   if (data_list$msmMode == 4) {
  #     map_list$logH_3 <- map_list$logH_3 * NA
  #     map_list$H_4 <- map_list$H_4 * NA
  #   }
  #
  #   # Holling Type 3
  #   if (data_list$msmMode == 5) {
  #     map_list$logH_3 <- map_list$logH_3 * NA
  #   }
  #
  #   # Predator interference
  #   if (data_list$msmMode == 6) {
  #     map_list$H_4 <- map_list$H_4 * NA
  #   }
  #
  #   # Predator preemption
  #   if (data_list$msmMode == 7) {
  #     map_list$H_4 <- map_list$H_4 * NA
  #   }
  #
  #   # Hassell-Varley
  #   if (data_list$msmMode == 8) {
  #     map_list$logH_3 <- map_list$logH_3 * NA
  #   }
  #
  #   # Ecosim
  #   if (data_list$msmMode == 9) {
  #     map_list$logH_2 <- map_list$logH_2 * NA
  #     map_list$H_4 <- map_list$H_4 * NA
  #   }
  # }

  # * 4. Suitability (if Multi-Species) ----
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

  # * 5. Diet Multiplier ----
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
#' @description
#'   \code{Selectivity} in \code{fleet_control} of the data determines shape of selectivity curve:
#' 0 = empirical selectivity provided in \code{emp_sel} in the data
#' 1 = logistic selectivity
#' 2 = non-parametric selecitivty sensu Ianelli et al 2018
#' 3 = double logistic
#' 4 = descending logistic
#' 5 = non-parametric selectivity sensu Taylor et al 2014 (Hake)
#'
#' \code{N_sel_bins}	Number of age/length bins to estimate non-parametric selectivity when Selectivity = 2 & 5. Not used otherwise
#'
#' \code{Time_varying_sel}	determines if time-varying selectivity should be estimated for logistic, double logistic selectivity,  descending logistic , or non-parametric (\code{Selectivity = 1, 3, 4, or 5}).
#' 0 = no
#' 1 = penalized deviates given \code{sel_sd_prior}
#' 3 = time blocks with no penality
#' 4 = random walk following Dorn
#' 5 = random walk on ascending portion of double logistic only.
#' NOTE: If selectivity is set to type = 2 (non-parametric) \code{sel_sd_prior} will be the 1st penalty on selectivity. \code{random_sel} treats random deviates and random walk parameters as random effects, estimating the variance.
#'
#'
#' @return Updated \code{map_list}.
build_map_selectivity <- function(map_list, data_list, nyrs_hind, random_sel) {

  # -- Map out parameters (then turned on)
  sel_params <- c("sel_coff", "sel_coff_dev", "ln_sel_slp", "sel_inf",
                  "ln_sel_slp_dev", "sel_inf_dev", "sel_dev_ln_sd", "sel_curve_pen")
  map_list[sel_params] <- lapply(map_list[sel_params], function(x) replace(x, values = NA))

  # -- Selectivity  indices
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
      if (sel_type > 0 && tv_sel %in% c(1, 2, 4, 5)) {
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
          dplyr::filter(Fleet_code == flt, Year - data_list$styr + 1 <= nyrs_hind)
        Selectivity_block <- fleet_data$Selectivity_block
        biom_yrs <- fleet_data$Year - data_list$styr + 1
        max_block <- max(Selectivity_block, 0)
      }

      # * 1. Logitistic (sel_type = 1) ----
      if (sel_type == 1) {

        # Turn on slp and asymptote for each sex
        for (sex in 1:nsex) {
          map_list$ln_sel_slp[1, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
          map_list$sel_inf[1, flt, sex] <- ind_inf; ind_inf <- ind_inf + 1
        }

        # Time-varying parameters
        if (tv_sel %in% c(1, 2, 4)) { # Random walk or deviate
          for (sex in 1:nsex) {
            map_list$ln_sel_slp_dev[1, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
            map_list$sel_inf_dev[1, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1

            ind_slp <- ind_slp + nyrs_hind
            ind_inf <- ind_inf + nyrs_hind
          }
          if (tv_sel == 4) { # Random walk: fix first deviate (start at mean)
            map_list$ln_sel_slp_dev[1, flt, , 1] <- NA
            map_list$sel_inf_dev[1, flt, , 1] <- NA
          }
        } else if (tv_sel == 3 && max_block > 0) { # Selectivity blocks
          for (sex in 1:nsex) {
            map_list$ln_sel_slp_dev[1, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
            map_list$sel_inf_dev[1, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

            ind_slp <- ind_slp + max_block
            ind_inf <- ind_inf + max_block

            # Turn off main parameters
            #FIXME will fail if random_sel = TRUE?
            map_list$ln_sel_slp[1, flt, sex] <- NA
            map_list$sel_inf[1, flt, sex] <- NA
          }
        }
      }


      # * 2. Non-parametric (sel_type = 2) ----
      if (sel_type == 2) {
        if (tv_sel %in% c(2, 4, 5)) { # Error check
          stop(paste0("'Time_varying_sel' for fleet ", flt, " with non-parametric selectivity is not 0 or 1. Current value: ", tv_sel))
        }

        age_first_selected <- data_list$fleet_control$Age_first_selected[i]
        minage_spp <- data_list$minage[spp]
        N_sel_bins <- data_list$fleet_control$N_sel_bins[i]

        if (is.na(age_first_selected)) age_first_selected <- minage_spp
        ages_on <- (age_first_selected - minage_spp + 1):N_sel_bins
        max_age_on <- max(ages_on)

        for (sex in 1:nsex) {
          map_list$sel_coff[flt, sex, ages_on] <- ind_coff + ages_on
          ind_coff <- ind_coff + max_age_on

          if (tv_sel == 1) { # Time-varying deviates
            map_list$sel_coff[flt, , ] <- NA # Must turn off mean parameter
            dev_indices <- ind_dev_coff + 1:(length(ages_on) * nyrs_hind)
            map_list$sel_coff_dev[flt, sex, ages_on, yrs_hind] <- dev_indices
            ind_dev_coff <- ind_dev_coff + length(dev_indices)
          }
        }
      }


      # * 5.3. Double logistic (sel_type = 3)
      if (sel_type == 3) {

        # Base parameters (j=1 ascending, j=2 descending)
        for (j in 1:2) {
          for (sex in 1:nsex) {
            map_list$ln_sel_slp[j, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
            map_list$sel_inf[j, flt, sex] <- ind_inf; ind_inf <- ind_inf + 1
          }
        }

        # Time varying parameters
        if (tv_sel %in% c(1, 2, 4, 5)) { # Random walk or deviate

          j_range <- if (tv_sel == 5) 1 else 1:2 # 5 only does ascending portion (j=1)

          for (j in j_range) {
            for (sex in 1:nsex) {
              map_list$ln_sel_slp_dev[j, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
              map_list$sel_inf_dev[j, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1

              ind_slp <- ind_slp + nyrs_hind
              ind_inf <- ind_inf + nyrs_hind
            }
          }

          # Random walk: fix first deviate
          if (tv_sel %in% c(4, 5)) {
            for (j in j_range) {
              map_list$ln_sel_slp_dev[j, flt, , 1] <- NA
              map_list$sel_inf_dev[j, flt, , 1] <- NA
            }
          }
        } else if (tv_sel == 3 && max_block > 0) { # Selectivity blocks
          for (j in 1:2) {
            for (sex in 1:nsex) {
              map_list$ln_sel_slp_dev[j, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
              map_list$sel_inf_dev[j, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

              ind_slp <- ind_slp + max_block
              ind_inf <- ind_inf + max_block
            }
          }
        }
      }


      # * 5.4. Descending logitistic (sel_type = 4)
      if (sel_type == 4) { # Descending portion only (j=2)

        # Base parameters
        for (sex in 1:nsex) {
          map_list$ln_sel_slp[2, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
          map_list$sel_inf[2, flt, sex] <- ind_inf; ind_inf <- ind_inf + 1
        }

        # Time varying parameters
        if (tv_sel %in% c(1, 2, 4)) { # Random walk or deviate
          for (sex in 1:nsex) {
            map_list$ln_sel_slp_dev[2, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
            map_list$sel_inf_dev[2, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1

            ind_slp <- ind_slp + nyrs_hind
            ind_inf <- ind_inf + nyrs_hind
          }

          if (tv_sel == 4) { # Random walk: fix first deviate
            map_list$ln_sel_slp_dev[2, flt, , 1] <- NA
            map_list$sel_inf_dev[2, flt, , 1] <- NA
          }
        } else if (tv_sel == 3 && max_block > 0) { # Selectivity blocks
          for (sex in 1:nsex) {
            map_list$ln_sel_slp_dev[2, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
            map_list$sel_inf_dev[2, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

            ind_slp <- ind_slp + max_block
            ind_inf <- ind_inf + max_block
          }
        }
      }


      # * 5.5. Non-parametric Hake-like (sel_type = 5)
      if (sel_type == 5) {
        if (tv_sel > 1) {
          warning(paste("Time_varying_sel for fleet", flt, "is not compatible (select NA, 0, or 1). Current value:", tv_sel))
        }

        age_first_selected <- data_list$fleet_control$Age_first_selected[i]
        minage_spp <- data_list$minage[spp]
        N_sel_bins <- data_list$fleet_control$N_sel_bins[i]

        if (is.na(age_first_selected)) age_first_selected <- minage_spp
        # +2 because first parameter is not-identifiable and is not estimated
        ages_on <- (age_first_selected - minage_spp + 2):N_sel_bins
        max_age_on <- max(ages_on, 0)

        for (sex in 1:nsex) {
          map_list$sel_coff[flt, sex, ages_on] <- ind_coff + ages_on
          ind_coff <- ind_coff + max_age_on

          if (tv_sel == 1) { # Time-varying deviates
            dev_indices <- ind_dev_coff + 1:(length(ages_on) * nyrs_hind)
            map_list$sel_coff_dev[flt, sex, ages_on, yrs_hind] <- dev_indices
            ind_dev_coff <- ind_dev_coff + length(dev_indices)
          }
        }
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
