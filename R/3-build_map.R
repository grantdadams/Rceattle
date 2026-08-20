#' @title Main function to construct the TMB map argument for CEATTLE
#'
#' @description Builds the TMB map, which tells TMB which parameters to estimate
#'   and which to hold fixed (a fixed parameter is mapped to `NA`), and which
#'   parameters share a single estimated value. One helper handles each process
#'   block (recruitment, M1, predation, selectivity, catchability, ...).
#'
#' @param data_list an Rceattle data_list
#' @param params A parameter list created from \code{\link{build_params}}.
#' @param debug Logical. If TRUE, fixes every parameter except the dummy
#'   (maps all to `NA`), so the model runs with no parameters estimated.
#' @param random_rec Logical. If TRUE, treats recruitment deviations as random effects,
#'   meaning the variance parameter (\code{R_log_sd}) is estimated.
#' @param random_sel Logical. If TRUE, treats selectivity deviations as random effects,
#'   meaning the variance parameter (\code{sel_dev_log_sd}) is estimated.
#'
#' @return A list containing the factorized TMB map (`mapFactor`) and the
#'   original map matrix/array list (`mapList`).
#' @export
build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE, random_sel = FALSE) {

  # Fill in defaulted switches and upgrade any deprecated column names
  data_list <- Rceattle::switch_check(data_list)

  # --- Setup Data and Initial Map ---
  nyrs_hind <- data_list$endyr - data_list$styr + 1
  nyrs_proj <- data_list$projyr - data_list$styr + 1
  # yrs_hind is calculated inside helpers as needed

  # Start every parameter estimated: give each element its own index. The
  # per-process helpers below then fix (NA) or share indices as needed.
  map_list <- lapply(params, function(x) {
    if (length(x) == 0) return(x)
    replace(x, values = seq_along(x))
  })

  # --- Parameter Components ---
  map_list <- build_map_recruitment(map_list, data_list, nyrs_hind, nyrs_proj, random_rec)

  map_list <- build_map_m1(map_list, data_list, nyrs_hind)

  map_list <- build_map_growth(map_list, data_list, nyrs_hind)

  map_list <- build_map_linkages(map_list, data_list)

  map_list <- build_map_predation(map_list, data_list)

  map_list <- build_map_selectivity(map_list, data_list, nyrs_hind, random_sel)

  map_list <- build_map_catchability(map_list, data_list, nyrs_hind)

  # After the per-process builders, not before: a linkage supplies the level
  # of the parameter it targets, so it must have the last word on whether
  # that base parameter stays estimable. Ordered ahead of them, a linkage on
  # catchability or selectivity would be silently overwritten by
  # build_map_catchability() / build_map_selectivity().
  map_list <- map_linkage_adjuster(map_list, data_list)

  map_list <- adjust_map_shared_params(map_list, data_list)

  map_list <- build_map_f_and_data_weights(map_list, data_list, nyrs_hind)

  map_list <- build_map_fixed_natage(map_list, data_list)

  # --- Debug Mode ---
  map_list <- build_map_debug(map_list, debug)

  # --- Final Steps ---
  map_list_grande <- list()
  map_list_grande$mapFactor <- lapply(map_list, factor)
  map_list_grande$mapList <- map_list

  return(map_list_grande)
}

## Helper Functions ----

#' @title Helper to set map for Recruitment parameters
#'
#' @description Maps the recruitment deviations (\code{rec_dev}, \code{init_dev}),
#'   stock-recruitment parameters (\code{rec_pars}) and Rec dev variances (\code{R_log_sd})
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

    # 1) Equilibrium with no devs. OffsetEquilibrium (5) likewise turns init_dev
    # off entirely -- the initial age-structure is the deterministic equilibrium
    # seeded by first-year recruitment, so every element the cpp reads (columns
    # 1:(nages-1)) is fixed at its build_params() value of 0.
    if (data_list$initMode %in% c("Equilibrium", "OffsetEquilibrium")) {
      map_list$init_dev[sp, ] <- NA
    }

    # 0, 2-3) Equilibrium or non-equilibrium with no devs
    if (!data_list$initMode %in% c("Equilibrium", "OffsetEquilibrium")) {
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
    map_list$R_log_sd[] <- NA
  }

  # Null model (mean-R)
  # - Turning off 2nd and 3rd par if only using mean rec
  if (data_list$srr_fun %in% c(0, 1) & data_list$srr_pred_fun %in% c(0, 1)) {
    map_list$rec_pars[, 2:3] <- NA
  }

  # Stock recruit models (Ricker & Beverton)
  # - Turning off mean-R if the hindcast itself uses the SRR, where R0 is derived
  #   from alpha and beta rather than estimated.
  if (data_list$srr_fun > 1) {
    map_list$rec_pars[, 1] <- NA
  }

  # - Fix alpha at the prior mean. Keyed on srr_pred_fun, the curve alpha belongs
  #   to, so it also covers the Ianelli configuration (srr_fun = 0,
  #   srr_pred_fun > 1), where the curve is estimated as a recruitment penalty
  #   and alpha is a free parameter just the same.
  if (data_list$srr_pred_fun > 1 && data_list$srr_est_mode == 0) {
    map_list$rec_pars[, 2] <- NA
  }

  return(map_list)
}


#' @title Helper to set map for Natural Mortality (M1) parameters
#'
#' @description Maps the fixed parameters (\code{log_M1}) and the random effects
#'   parameters (\code{log_M1_dev}, \code{M1_dev_log_sd}, \code{M1_rho}) based on
#'   \code{M1_model} and \code{M1_re} settings.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_m1 <- function(map_list, data_list, nyrs_hind) {

  # Map out natural mortality parameters
  M1_params <- c("log_M1", "M1_beta", "log_M1_dev", "M1_dev_log_sd", "M1_rho")
  map_list[M1_params] <- lapply(map_list[M1_params], function(x) replace(x, values = NA))

  # Generic indices for looping
  M1_ind <- 1
  M1_beta_ind <- 0
  M1_dev_ind <- 0
  M1_rho_ind <- 0

  # Loop through species and turn on based on model
  for (sp in 1:data_list$nspp) {

    # * Switches
    nages_sp <- data_list$nages[sp]
    nsex_sp <- data_list$nsex[sp]
    M1_model <- data_list$M1_model[sp] # Fixed effects model
    M1_re_model <- data_list$M1_re[sp] # Random effects model

    # * 1. Fixed effects ----
    # ** M1_model = 1: sex- and age-invariant M1
    if (M1_model == 1) {
      map_list$log_M1[sp, , 1:nages_sp] <- M1_ind
      M1_ind <- M1_ind + 1
    }


    # ** M1_model = 2: sex-specific, but age-invariant M1
    if (M1_model == 2) {
      map_list$log_M1[sp, 1, 1:nages_sp] <- M1_ind
      map_list$log_M1[sp, 2, 1:nages_sp] <- M1_ind + 1
      M1_ind <- M1_ind + 2
      if (nsex_sp == 1) {
        warning(paste0("M1 model for species ", sp," is set to 2 (sex-specific), but species is single-sex."))
        map_list$log_M1[sp, 2, 1:nages_sp] <- M1_ind
      }
    }

    # ** M1_model = 3: sex-specific, age-specific M1
    if (M1_model == 3) {
      map_list$log_M1[sp, 1, 1:nages_sp] <- M1_ind:(M1_ind + nages_sp - 1)
      M1_ind <- M1_ind + nages_sp
      if (nsex_sp == 2) {
        map_list$log_M1[sp, 2, 1:nages_sp] <- M1_ind:(M1_ind + nages_sp - 1)
        M1_ind <- M1_ind + nages_sp
      } else {
        warning(paste0("M1 model for species ", sp," is set to 3 (sex-specific), but species is single-sex."))
        map_list$log_M1[sp, 2, ] <- map_list$log_M1[sp, 1, ]
      }
    }

    # ** M1_model = 4: environmentally driven sex- and age-invariant M1
    if (M1_model == 4| (M1_model == 5 & nsex_sp == 1)) {

      # - Mean M
      map_list$log_M1[sp, , 1:nages_sp] <- M1_ind
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
      map_list$log_M1[sp, 1, 1:nages_sp] <- M1_ind # Females
      map_list$log_M1[sp, 2, 1:nages_sp] <- M1_ind + 1 # Males
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
    # "log_M1_dev"
    # - M1_re = 1/4: Random effects varies by age (IID or AR1) and constant over years.
    if(M1_re_model %in% c(1, 4)){
      if(M1_model == 1){ # Sex-invariant
        # - Random effects
        map_list$log_M1_dev[sp,1, 1:nages_sp,] <- M1_dev_ind + 1:nages_sp

        # Males mapped the same, if present
        if(nsex_sp == 2){
          map_list$log_M1_dev[sp,2,,] <- map_list$log_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + nages_sp
      }

      if(M1_model == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$log_M1_dev[sp, 1, 1:nages_sp,] <- M1_dev_ind + 1:nages_sp
        M1_dev_ind = M1_dev_ind + nages_sp

        # -- Males
        map_list$log_M1_dev[sp, 2, 1:nages_sp,] <- M1_dev_ind + 1:nages_sp
        M1_dev_ind = M1_dev_ind + nages_sp
      }


      # - Standard deviation (shared across sexes)
      map_list$M1_dev_log_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(M1_re_model == 4){
        map_list$M1_rho[sp,,1] =  sp
      }
    }

    # - M1_re = 2/5: Random effects varies by year (IID or AR1) and constant over ages
    if(M1_re_model %in% c(2, 5)){
      if(M1_model == 1){ # Sex-invariant
        # - Random effects
        map_list$log_M1_dev[sp,1,1:nages_sp, 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = nages_sp)

        # Males mapped the same, if present
        if(nsex_sp == 2){
          map_list$log_M1_dev[sp,2,,] <- map_list$log_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + nyrs_hind
      }

      if(nsex_sp == 2 & M1_model == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$log_M1_dev[sp,1, 1:nages_sp, 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = nages_sp)
        M1_dev_ind = M1_dev_ind + nyrs_hind

        # -- Males
        map_list$log_M1_dev[sp,2, 1:nages_sp, 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = nages_sp)
        M1_dev_ind = M1_dev_ind + nyrs_hind
      }

      # - Standard deviation (shared across sexes)
      map_list$M1_dev_log_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(M1_re_model == 5){
        map_list$M1_rho[sp,,2] =  sp #FIXME: may want sex-varying?? Hard to estimate
      }
    }

    # - M1_re = 3/6: Random effects varies by age and year (IID or 2D-AR1)
    # One free deviation per age-year cell: the index vector must be
    # 1:(nyrs_hind * nages_sp), not 1:nyrs_hind. The shorter form recycled into
    # the age dimension without warning (an exact multiple), leaving the field
    # with only nyrs_hind distinct parameters on a stride pattern while the
    # density and the simulator both treated it as a full age x year grid.
    if(M1_re_model %in% c(3, 6)){
      if(M1_model == 1){ # Sex-invariant
        # - Random effects
        map_list$log_M1_dev[sp,1,1:nages_sp, 1:nyrs_hind] <- M1_dev_ind + 1:(nyrs_hind * nages_sp)

        # Males mapped the same, if present
        if(nsex_sp == 2){
          map_list$log_M1_dev[sp,2,,] <- map_list$log_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + (nyrs_hind * nages_sp)
      }

      if(nsex_sp == 2 & M1_model == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$log_M1_dev[sp,1, 1:nages_sp, 1:nyrs_hind] <- M1_dev_ind + 1:(nyrs_hind * nages_sp)
        M1_dev_ind = M1_dev_ind + (nyrs_hind * nages_sp)

        # -- Males
        map_list$log_M1_dev[sp,2, 1:nages_sp, 1:nyrs_hind] <- M1_dev_ind + 1:(nyrs_hind * nages_sp)
        M1_dev_ind = M1_dev_ind + (nyrs_hind * nages_sp)
      }

      # - Standard deviation (shared across sexes)
      map_list$M1_dev_log_sd[sp,] = sp

      # AR1 correlations (age AND year, shared across sexes) for the separable
      # 2D-AR1 (mode 6). The gate was `== 5` -- impossible inside this
      # `%in% c(3, 6)` block, so mode 6's rho parameters were never estimated and
      # the separable AR1 silently collapsed to IID. Estimate both here.
      # `M1_rho_ind` is offset past nspp so these two map values never collide
      # with the single `sp`-valued rho of a mode-4/5 species.
      if(M1_re_model == 6){
        # dim 3: 1 = age correlation, 2 = year correlation. `sp,,k` sets every
        # sex at once (shared across sexes), matching the mode-4/5 style and
        # working for 1- or 2-sex models. Two distinct map values so age and
        # year rho are estimated separately.
        map_list$M1_rho[sp,,1] = data_list$nspp + M1_rho_ind + 1  # age rho
        map_list$M1_rho[sp,,2] = data_list$nspp + M1_rho_ind + 2  # year rho
        M1_rho_ind = M1_rho_ind + 2  #FIXME: may want sex-varying?? Hard to estimate
      }
    }
  }
  return(map_list)
}

#' @title Helper to set map for growth parameters
#'
#' @description Maps the growth parameters (\code{log_growth_pars},
#'   \code{growth_log_sd}, \code{weight_length_pars}) according to
#'   \code{growth_model}. Everything starts mapped off and each growth function
#'   turns on the parameters it uses.
#'
#'   Time-varying growth comes from the linkage grammar
#'   (\code{build_growth(linkages = )}), whose random effects carry their own
#'   density and map.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_growth <- function(map_list, data_list, nyrs_hind) {

  # Map out growth parameters
  growth_params <- c("log_growth_pars", "weight_length_pars", "growth_log_sd")
  map_list[growth_params] <- lapply(map_list[growth_params], function(x) replace(x, values = NA))

  # Loop through species and turn on based on model
  for (sp in 1:data_list$nspp) {

    # * Switches ----
    nages_sp <- data_list$nages[sp]
    nsex_sp <- data_list$nsex[sp]
    growth_model <- data_list$growth_model[sp]

    # * 1. Fixed effects ----
    if(growth_model %in% c(1, "vonBertalanffy")){ # Von Bertalanffy
      # k, l1, linf
      map_list$log_growth_pars[sp, 1, 1:3] <- (sp - 1) * 4 + 1:3 # Females/sex combined
      map_list$growth_log_sd[sp, 1, 1:2] <- (sp - 1) * 2 + 1:2 # Growth SD
      if(nsex_sp == 2){
        map_list$log_growth_pars[sp, 2, 1:3] <- (sp - 1) * 4 + 5:7 # Males
        map_list$growth_log_sd[sp, 2, 1:2] <- (sp - 1) * 2 + 1:2 # Growth SD same as females
      }
    }

    if(growth_model %in% c(2, "Richards")){ # Richards
      # k, l1, linf, m
      map_list$log_growth_pars[sp, 1, 1:4] <- (sp - 1) * 4 + 1:4 # Females/sex combined
      map_list$growth_log_sd[sp, 1, 1:2] <- (sp - 1) * 2 + 1:2 # Growth SD
      if(nsex_sp == 2){
        map_list$log_growth_pars[sp, 2, 1:4] <- (sp - 1) * 4 + 5:8 # Males
        map_list$growth_log_sd[sp, 2, 1:2] <- (sp - 1) * 2 + 1:2 # Growth SD same as females
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
  map_list$diet_comp_weights[] <- NA
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

      # Dirichlet-multinomial overdispersion for the diet composition. Only
      # estimated where the diet likelihood is actually fit: empirical
      # suitability (suitMode = 0) takes the diet proportions as given and the
      # cpp skips that predator, leaving nothing to inform the weight.
      if (data_list$suitMode[sp] > 0 && data_list$Diet_distribution[sp] == 1) {
        map_list$diet_comp_weights[sp] <- sp
      }
    }
  }

  return(map_list)
}

#' @title Helper to set map for Selectivity parameters
#'
#' @description Maps base selectivity parameters (\code{log_sel_slp}, \code{sel_inf},
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
#' `0` = "Fixed" empirical selectivity provided in \code{emp_sel} in the data
#' `1` = "Logistic"
#' `2` = "NonParametric" selecitivty sensu Ianelli et al 2018
#' `3` = "DoubleLogistic"
#' `4` = "DescendingLogistic"
#' `5` = "Hake" non-parametric selectivity sensu Taylor et al 2014 (Hake)
#' `6` = "2DAR1" across age x year
#' `7` = "3DAR1" across age x cohort x year (Cheng et al 2024)
#' `8` = "DoubleNormal" Gaussian ascending/descending limbs blended at peak (analogous to SS3 pattern 24).
#'     Parameters: `sel_inf[1]` = peak; `sel_inf[2]` = logit(right_floor) (right-tail floor, SS3 P6/end_logit);
#'     `log_sel_slp[1]` = log(sigma_asc); `log_sel_slp[2]` = log(sigma_desc).
#'     right_floor->0: dome-shaped; right_floor->1: logistic ascending only.
#'
#' \code{N_sel_bins}	Number of age/length bins to estimate non-parametric selectivity when Selectivity = 2 or 5. Not used otherwise
#'
#' \code{Time_varying_sel}	determines if time-varying selectivity should be estimated for logistic, double logistic selectivity,  descending logistic , non-parametric, or hake (\code{Selectivity = 1, 2, 3, 4, or 5}).
#' `0` = 'Off'
#' `1` = 'IID' penalized deviates given \code{sel_sd_prior}
#' `3` = 'Block' time blocks with no penalty
#' `4` = 'RandomWalk' random walk
#' `5` = 'RandomWalkAscending' random walk on ascending portion of double logistic only.
#' \code{random_sel} in \code{fit_mod} treats random deviates and random walk parameters as random effects, estimating the variance.
#'
#'
#' @return Updated \code{map_list}.
build_map_selectivity <- function(map_list, data_list, nyrs_hind, random_sel) {

  # -- Map out parameters (then turned on)
  sel_params <- c("sel_coff", "sel_coff_dev", "log_sel_slp", "sel_inf",
                  "log_sel_slp_dev", "sel_inf_dev", "sel_dev_log_sd", "sel_curve_pen")
  map_list[sel_params] <- lapply(map_list[sel_params], function(x) replace(x, values = NA))

  # -- Selectivity  indices
  ind_coff <- 1
  ind_dev_coff <- 1
  ind_slp <- 1
  ind_inf <- 1
  yrs_hind <- 1:nyrs_hind
  n_flt <- nrow(data_list$fleet_control)

  # Turn on variance of random effects for selectivity deviates (sigma)
  if (random_sel) {
    for (i in 1:n_flt) {
      flt <- data_list$fleet_control$Fleet_code[i]
      sel_type <- data_list$fleet_control$Selectivity[i]
      tv_sel <- data_list$fleet_control$Time_varying_sel[i]
      if (sel_type != "Fixed" && tv_sel %in% c("IID", "AR1", "RandomWalk", "RandomWalkAscending")) {
        map_list$sel_dev_log_sd[flt] <- flt
      }

      if (sel_type %in% c("2DAR1", "3DAR1")) {
        map_list$sel_dev_log_sd[flt] <- flt
      }
    }
  }


  # Selectivity start year per fleet, resolved to the minimum across each
  # Selectivity_index group: mirrored fleets share one deviation block, so the
  # mask must start at the earliest year any of them has data (otherwise the
  # result depends on fleet_control row order).
  sel_start_yr_grp <- .sel_start_year_by_group(data_list$fleet_control, data_list$styr)

  # Loop through fleets to set up selectivity parameters
  for (i in 1:n_flt) {
    flt <- data_list$fleet_control$Fleet_code[i]
    spp <- data_list$fleet_control$Species[i]
    nsex <- data_list$nsex[spp]
    sel_type <- data_list$fleet_control$Selectivity[i]
    tv_sel <- data_list$fleet_control$Time_varying_sel[i]

    if (data_list$fleet_control$Fleet_type[i] != "Off") {

      # Helper for selectivity blocks logic
      max_block <- 0
      if (tv_sel == "Block") {
        data_source <- if (data_list$fleet_control$Fleet_type[i] == "Fishery") data_list$catch_data else data_list$index_data
        fleet_data <- data_source |>
          dplyr::filter(Fleet_code == flt, Year - data_list$styr + 1 <= nyrs_hind)
        Selectivity_block <- fleet_data$Selectivity_block
        block_yrs <- fleet_data$Year - data_list$styr + 1
        max_block <- max(Selectivity_block, 0)
      }

      # * Logitistic ----
      # - sel_type = 1 (age-based), 8 (length-based)
      if (sel_type == "Logistic") {

        # Turn on slp and asymptote for each sex
        for (sex in 1:nsex) {
          map_list$log_sel_slp[1, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
          map_list$sel_inf[1, flt, sex] <- ind_inf; ind_inf <- ind_inf + 1
        }

        # Time-varying parameters
        if (tv_sel %in% c('IID', 'AR1', 'RandomWalk')) { # Random walk or deviate
          for (sex in 1:nsex) {
            map_list$log_sel_slp_dev[1, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
            map_list$sel_inf_dev[1, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1

            ind_slp <- ind_slp + nyrs_hind
            ind_inf <- ind_inf + nyrs_hind
          }
          if (tv_sel == "RandomWalk") { # Random walk: fix first deviate (start at mean)
            map_list$log_sel_slp_dev[1, flt, , 1] <- NA
            map_list$sel_inf_dev[1, flt, , 1] <- NA
          }
        } else if (tv_sel == "Block" && max_block > 0) { # Selectivity blocks
          for (sex in 1:nsex) {
            map_list$log_sel_slp_dev[1, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_slp
            map_list$sel_inf_dev[1, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_inf

            ind_slp <- ind_slp + max_block
            ind_inf <- ind_inf + max_block

            # Turn off main parameters
            #FIXME will fail if random_sel = TRUE?
            map_list$log_sel_slp[1, flt, sex] <- NA
            map_list$sel_inf[1, flt, sex] <- NA
          }
        }
      }


      # * LogisticPM ----
      # ---- sel_type = 11 (ADMB AMAK "pm" bottom-trawl survey).
      #      Logistic (slope = log_sel_slp[1], inflection = sel_inf[1]) with
      #      MULTIPLICATIVE time-varying deviates, PLUS a free first-bin (age-1)
      #      log-selectivity carried in sel_inf[2] (base) and sel_inf_dev[2] (deviates).
      #      Time_varying_sel must be "Off" or "RandomWalk".
      if (sel_type == "LogisticPM") {

        # Base parameters: logistic slope + inflection, and the free age-1 log-selectivity
        for (sex in 1:nsex) {
          map_list$log_sel_slp[1, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
          map_list$sel_inf[1, flt, sex]     <- ind_inf; ind_inf <- ind_inf + 1
          map_list$sel_inf[2, flt, sex]     <- ind_inf; ind_inf <- ind_inf + 1  # age-1 log-sel base
        }

        if (tv_sel %in% c("IID", "AR1", "RandomWalk")) {
          for (sex in 1:nsex) {
            map_list$log_sel_slp_dev[1, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1  # slope devs
            map_list$sel_inf_dev[1, flt, sex, yrs_hind]     <- ind_inf + yrs_hind - 1  # inflection devs
            ind_slp <- ind_slp + nyrs_hind
            ind_inf <- ind_inf + nyrs_hind
            map_list$sel_inf_dev[2, flt, sex, yrs_hind]     <- ind_inf + yrs_hind - 1  # age-1 devs
            ind_inf <- ind_inf + nyrs_hind
          }
          if (tv_sel == "RandomWalk") {
            # Fix every deviate through the fleet's start year, so the first
            # estimated deviate is Sel_start_year + 1. This pins the random walk's
            # level (otherwise the base parameter and deviates are confounded) and
            # drops the pre-start deviates, which have neither data nor a penalty
            # (the cpp penalties also begin at start_yr + 1) and would otherwise
            # leave flat directions that stall the optimizer.
            sel_start_yr <- sel_start_yr_grp[i]  # group-resolved (mirrored fleets share one block)
            start_idx <- if (is.null(sel_start_yr) || is.na(sel_start_yr)) 1L else
              max(1L, min(nyrs_hind, as.integer(sel_start_yr) - data_list$styr + 1L))
            map_list$log_sel_slp_dev[1, flt, , 1:start_idx] <- NA
            map_list$sel_inf_dev[1, flt, , 1:start_idx]     <- NA
            map_list$sel_inf_dev[2, flt, , 1:start_idx]     <- NA
          }
        }
      }


      # * Non-parametric ----
      # ---- sel_type = 2 (Ianelli penalty), 9 (NonParametricPM, ADMB AMAK penalty).
      #      Both share identical parameters / mapping; they differ only in the
      #      selectivity penalty form (see ceattle.cpp).
      if (sel_type %in% c("NonParametric", "NonParametricPM")) {
        if (tv_sel %in% c("AR1", "IID", "RandomWalkAscending")) { # Error check
          stop(paste0("'Time_varying_sel' for fleet ", flt, " with non-parametric selectivity is not 'Off'(0) or 'RandomWalk'(4). Current value: ", tv_sel))
        }

        bin_first_selected <- data_list$fleet_control$Bin_first_selected[i]
        N_sel_bins <- data_list$fleet_control$N_sel_bins[i]
        if (is.na(N_sel_bins)) stop(paste0("N_sel_bins is NA for fleet ", i))

        if (is.na(bin_first_selected)) bin_first_selected <- 1
        bins_on <- bin_first_selected:N_sel_bins
        max_bin_on <- max(bins_on)

        for (sex in 1:nsex) {
          map_list$sel_coff[flt, sex, bins_on] <- ind_coff + bins_on
          ind_coff <- ind_coff + max_bin_on

          if (tv_sel == "RandomWalk") { # Time-varying deviates
            map_list$sel_coff[flt, , ] <- NA # Must turn off mean parameter
            dev_indices <- ind_dev_coff + 1:(length(bins_on) * nyrs_hind)
            map_list$sel_coff_dev[flt, sex, bins_on, yrs_hind] <- dev_indices
            ind_dev_coff <- ind_dev_coff + length(dev_indices)

            # Fix the deviates BEFORE the start year. The mean parameter
            # (sel_coff) is mapped off for a random walk, so the deviate AT the
            # start year carries the base shape and stays estimated; only the
            # earlier ones are dropped, having neither data nor a penalty (all
            # NonParametricPM penalties begin at start_yr).
            sel_start_yr <- sel_start_yr_grp[i]  # group-resolved (mirrored fleets share one block)
            start_idx <- if (is.null(sel_start_yr) || is.na(sel_start_yr)) 1L else
              max(1L, min(nyrs_hind, as.integer(sel_start_yr) - data_list$styr + 1L))
            if (start_idx > 1L) {
              map_list$sel_coff_dev[flt, sex, bins_on, 1:(start_idx - 1L)] <- NA
            }
          }
        }
      }


      # * Double logistic ----
      # ---- sel_type = 3 (age-based), 10 (length-based)
      if (sel_type == "DoubleLogistic") {

        # Base parameters (j=1 ascending, j=2 descending)
        for (j in 1:2) {
          for (sex in 1:nsex) {
            map_list$log_sel_slp[j, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
            map_list$sel_inf[j, flt, sex] <- ind_inf; ind_inf <- ind_inf + 1
          }
        }

        # Time varying parameters
        if (tv_sel %in% c("IID", "AR1", "RandomWalk", "RandomWalkAscending")) { # Random walk or deviate

          j_range <- if (tv_sel == "RandomWalkAscending") 1 else 1:2 # RandomWalkAscending only does ascending portion (j=1)

          for (j in j_range) {
            for (sex in 1:nsex) {
              map_list$log_sel_slp_dev[j, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
              map_list$sel_inf_dev[j, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1

              ind_slp <- ind_slp + nyrs_hind
              ind_inf <- ind_inf + nyrs_hind
            }
          }

          # Random walk: fix first deviate
          if (tv_sel %in% c("RandomWalk", "RandomWalkAscending")) {
            for (j in j_range) {
              map_list$log_sel_slp_dev[j, flt, , 1] <- NA
              map_list$sel_inf_dev[j, flt, , 1] <- NA
            }
          }
        } else if (tv_sel == "Block" && max_block > 0) { # Selectivity blocks
          for (j in 1:2) {
            for (sex in 1:nsex) {
              map_list$log_sel_slp_dev[j, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_slp
              map_list$sel_inf_dev[j, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_inf

              ind_slp <- ind_slp + max_block
              ind_inf <- ind_inf + max_block
            }
          }
        }
      }


      # * Double Normal ----
      # Parameters reuse existing arrays (same dimensions as DoubleLogistic):
      #   sel_inf[1]     = peak bin/length (mode)
      #   sel_inf[2]     = logit(right_floor) -- right-tail floor (analogous to SS3 P6 / end_logit)
      #   log_sel_slp[1] = log(sigma_ascending)
      #   log_sel_slp[2] = log(sigma_descending)
      if (sel_type == "DoubleNormal") {

        # Base parameters
        for (sex in 1:nsex) {
          map_list$sel_inf[1, flt, sex]     <- ind_inf; ind_inf <- ind_inf + 1
          map_list$sel_inf[2, flt, sex]     <- ind_inf; ind_inf <- ind_inf + 1
          map_list$log_sel_slp[1, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
          map_list$log_sel_slp[2, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
        }

        # Time-varying parameters
        if (tv_sel %in% c("IID", "AR1", "RandomWalk")) {
          for (sex in 1:nsex) {
            # Peak (sel_inf[1]) and right-floor (sel_inf[2]) deviates
            map_list$sel_inf_dev[1, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1
            ind_inf <- ind_inf + nyrs_hind
            map_list$sel_inf_dev[2, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1
            ind_inf <- ind_inf + nyrs_hind

            # Ascending-SD (log_sel_slp[1]) and descending-SD (log_sel_slp[2]) deviates
            map_list$log_sel_slp_dev[1, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
            ind_slp <- ind_slp + nyrs_hind
            map_list$log_sel_slp_dev[2, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
            ind_slp <- ind_slp + nyrs_hind
          }

          # Random walk: fix first deviate
          if (tv_sel == "RandomWalk") {
            map_list$sel_inf_dev[1, flt, , 1]     <- NA
            map_list$sel_inf_dev[2, flt, , 1]     <- NA
            map_list$log_sel_slp_dev[1, flt, , 1] <- NA
            map_list$log_sel_slp_dev[2, flt, , 1] <- NA
          }

        } else if (tv_sel == "Block" && max_block > 0) {
          for (sex in 1:nsex) {
            # Peak and right-floor block deviates
            map_list$sel_inf_dev[1, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_inf
            ind_inf <- ind_inf + max_block
            map_list$sel_inf_dev[2, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_inf
            ind_inf <- ind_inf + max_block

            # Ascending-SD and descending-SD block deviates
            map_list$log_sel_slp_dev[1, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_slp
            ind_slp <- ind_slp + max_block
            map_list$log_sel_slp_dev[2, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_slp
            ind_slp <- ind_slp + max_block

            # Fix base parameters -- blocks fully replace them (deviates are absolute values)
            map_list$sel_inf[1, flt, sex]     <- NA
            map_list$sel_inf[2, flt, sex]     <- NA
            map_list$log_sel_slp[1, flt, sex] <- NA
            map_list$log_sel_slp[2, flt, sex] <- NA
          }
        }
      }


      # * Descending logitistic ----
      # ---- sel_type = 4 (age-based), 11 (length-based)
      if (sel_type == "DescendingLogistic") {

        # Base parameters
        for (sex in 1:nsex) {
          map_list$log_sel_slp[2, flt, sex] <- ind_slp; ind_slp <- ind_slp + 1
          map_list$sel_inf[2, flt, sex] <- ind_inf; ind_inf <- ind_inf + 1
        }

        # Time varying parameters
        if (tv_sel %in% c("IID", "AR1", "RandomWalk")) { # Random walk or deviate
          for (sex in 1:nsex) {
            map_list$log_sel_slp_dev[2, flt, sex, yrs_hind] <- ind_slp + yrs_hind - 1
            map_list$sel_inf_dev[2, flt, sex, yrs_hind] <- ind_inf + yrs_hind - 1

            ind_slp <- ind_slp + nyrs_hind
            ind_inf <- ind_inf + nyrs_hind
          }

          if (tv_sel == "RandomWalk") { # Random walk: fix first deviate
            map_list$log_sel_slp_dev[2, flt, , 1] <- NA
            map_list$sel_inf_dev[2, flt, , 1] <- NA
          }
        } else if (tv_sel == "Block" && max_block > 0) { # Selectivity blocks
          for (sex in 1:nsex) {
            map_list$log_sel_slp_dev[2, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_slp
            map_list$sel_inf_dev[2, flt, sex, block_yrs] <- Selectivity_block - 1 + ind_inf

            ind_slp <- ind_slp + max_block
            ind_inf <- ind_inf + max_block
          }
        }
      }


      # * Non-parametric Hake-like ----
      # ---- sel_type = 5 (age-based), 12 (length-based)
      if (sel_type == "Hake") {
        # "Off" (0) is legal here -- the bin coefficients are estimated with no
        # annual deviates -- and data_check() enforces exactly {"Off", "IID"}
        # for this form. Warn only on a mode the Hake form cannot represent
        # ("Block", "AR1", "RandomWalk", ...); the old `tv_sel != "IID"` guard
        # fired on "Off" too, which is why a valid Time_varying_sel = 0 fleet
        # warned about itself. Guard NA explicitly -- `NA != "IID"` is NA, which
        # errors inside `if()` rather than falling through.
        if (!is.na(tv_sel) && !tv_sel %in% c("Off", "IID")) {
          warning(paste("Time_varying_sel for fleet", flt, "is not compatible (select 'Off' (0) or 'IID' (1)). Current value:", tv_sel))
        }

        bin_first_selected <- data_list$fleet_control$Bin_first_selected[i]
        N_sel_bins <- data_list$fleet_control$N_sel_bins[i]
        if (is.na(N_sel_bins)) stop(paste0("N_sel_bins is NA for fleet ", i))

        if (is.na(bin_first_selected)) bin_first_selected <- 1
        # +1 because first parameter is not-identifiable and is not estimated
        bins_on <- (bin_first_selected + 1):N_sel_bins
        max_bin_on <- max(bins_on, 0)

        for (sex in 1:nsex) {
          map_list$sel_coff[flt, sex, bins_on] <- ind_coff + bins_on
          ind_coff <- ind_coff + max_bin_on

          if (tv_sel == 'IID') { # Time-varying deviates
            dev_indices <- ind_dev_coff + 1:(length(bins_on) * nyrs_hind)
            map_list$sel_coff_dev[flt, sex, bins_on, yrs_hind] <- dev_indices
            ind_dev_coff <- ind_dev_coff + length(dev_indices)
          }
        }
      }

      # * 2DAR1 ----
      # ---- sel_type = 6 (age-based), 13 (length-based)
      if (sel_type == "2DAR1") {
        if (!is.na(tv_sel)) {
          warning(paste("Time_varying_sel for fleet", flt, "is ignored for 2DAR1 selectivity."))
        }
        if(!random_sel){
          warning("Laplace approximation for selectivity is turned off via `random_sel` in `fit_mod`.  Variance will not be estimated.")
        }

        bin_first_selected <- data_list$fleet_control$Bin_first_selected[i]
        N_sel_bins <- data_list$fleet_control$N_sel_bins[i]
        if (is.na(N_sel_bins)) stop(paste0("N_sel_bins is NA for fleet ", i))

        if (is.na(bin_first_selected)) bin_first_selected <- 1
        bins_on <- (bin_first_selected):N_sel_bins
        max_bin_on <- max(bins_on, 0)

        for (sex in 1:nsex) {
          map_list$sel_coff[flt, sex, bins_on] <- ind_coff + bins_on
          ind_coff <- ind_coff + max_bin_on

          # Time-varying deviates
          dev_indices <- ind_dev_coff + 1:(length(bins_on) * nyrs_hind)
          map_list$sel_coff_dev[flt, sex, bins_on, yrs_hind] <- dev_indices
          ind_dev_coff <- ind_dev_coff + length(dev_indices)
        }

        # Rho. Slot 1 is the correlation across selectivity BINS, slot 2
        # across YEARS -- the array passed to SEPARABLE() is (bin, year) and
        # SEPARABLE(f, g) puts g on the fastest-running dimension. These two
        # comments said the opposite until 5.9.0, as did the density itself.
        map_list$sel_curve_pen[flt,1] <- flt # bin
        map_list$sel_curve_pen[flt,2] <- flt + n_flt # year
      }


      # * 3DAR1 ----
      # ---- sel_type = 7 (age-based), 14 (length-based)
      if (sel_type == "3DAR1") {
        if (!is.na(tv_sel)) {
          warning(paste("Time_varying_sel for fleet", flt, "is ignored for 3DAR1 selectivity."))
        }
        if(!random_sel){
          warning("Laplace approximation for selectivity is turned off via `random_sel` in `fit_mod`. Variance will not be estimated.")
        }

        bin_first_selected <- data_list$fleet_control$Bin_first_selected[i]
        N_sel_bins <- data_list$fleet_control$N_sel_bins[i]
        if (is.na(N_sel_bins)) stop(paste0("N_sel_bins is NA for fleet ", i))

        if (is.na(bin_first_selected)) bin_first_selected <- 1
        bins_on <- (bin_first_selected):N_sel_bins
        max_bin_on <- max(bins_on, 0)

        for (sex in 1:nsex) {
          map_list$sel_coff[flt, sex, bins_on] <- ind_coff + bins_on
          ind_coff <- ind_coff + max_bin_on

          # Time-varying deviates
          dev_indices <- ind_dev_coff + 1:(length(bins_on) * nyrs_hind)
          map_list$sel_coff_dev[flt, sex, bins_on, yrs_hind] <- dev_indices
          ind_dev_coff <- ind_dev_coff + length(dev_indices)
        }

        # Rho, as for 2DAR1 above: slot 1 across BINS, slot 2 across YEARS,
        # slot 3 across COHORTS. construct_Q()'s own argument names are
        # inverted relative to what they multiply, and Rceattle's ay_index fill
        # is the mirror of WHAM's, so the two cancel -- the mapping here has
        # always been bin/year/cohort, whatever these comments used to say.
        map_list$sel_curve_pen[flt,1] <- flt # bin
        map_list$sel_curve_pen[flt,2] <- flt + n_flt # year
        map_list$sel_curve_pen[flt,3] <- flt + n_flt * 2 # cohort
      }
    }
  }
  return(map_list)
}


#' @title Helper to set map for Catchability parameters
#'
#' @description Maps catchability base parameters (\code{index_log_q}),
#'   time-varying deviations (\code{index_q_dev}), and environmental linkages
#'   (\code{index_q_beta}, \code{index_q_rho}) for every fleet that carries
#'   fitted \code{index_data} -- a fishery with a CPUE series as much as a
#'   survey. A fleet with no index rows gets none of them, whatever its
#'   \code{Catchability} says, since a q with no index to inform it is a flat
#'   direction. Sharing overrides this: \code{adjust_map_shared_params()} then
#'   copies the lead fleet's slice across each \code{Catchability_index} group.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_catchability <- function(map_list, data_list, nyrs_hind) {

  # -- Catchability indices
  ind_q_dev <- 1
  ind_beta_q <- 0
  yrs_hind <- 1:nyrs_hind


  catchability_params <- c("index_log_q", "index_q_beta", "index_q_rho", "index_q_dev", "index_q_log_sd", "index_q_dev_log_sd", "index_log_sd") # "index_q_pow"
  map_list[catchability_params] <- lapply(map_list[catchability_params], function(x) replace(x, values = rep(NA, length(x))))

  # Fleets whose catchability block is estimable: those carrying fitted index
  # observations, whatever their Fleet_type. The template fits an index row for
  # any non-Off fleet, so a fishery CPUE series is scored like a survey's index
  # and needs its q the same way; keying on Fleet_type == "Survey" instead left a
  # fishery's q, time-varying q and index sd mapped out, so Catchability =
  # "Estimated" did nothing.
  #
  # The converse holds too, and is why this is the data and not the fleet type: a
  # q with no index to inform it is a flat direction in the likelihood. A survey
  # with no index rows no longer gets one either.
  #
  # A fleet with no index of its own can still end up estimated, by sharing a
  # Catchability_index group whose LEAD estimates -- adjust_map_shared_params()
  # copies the lead's slice over the group afterwards, which is intended.
  q_fleets <- .fleets_with_index(data_list)

  # Loop through fleets
  for( i in 1: nrow(data_list$fleet_control)){
    flt = data_list$fleet_control$Fleet_code[i]
    if(flt %in% q_fleets){
      # Q
      # - 0 = fixed at prior
      # - 1 = Estimate single parameter
      # - 2 = Estimate single parameter with prior
      # - 3 = Estimate analytical q
      # - 4 = Estimate power equation
      # - 5 = Use env index ln(q_y) = q_mu + beta * index_y
      # - 6 = Fit to env index dnorm(d_y, env_index, sigma) [Rogers et al 2024]


      # - Turn on mean q for:
      # - 1 = Estimate single parameter
      # - 2 = Estimate single parameter with prior
      # - 4 = Estimate power equation
      # - 5 = Use env index ln(q_y) = q_mu + beta * index_y
      # - 6 = Fit to env index
      # "Analytical" (geometric) and "AnalyticalArith" (arithmetic) both SOLVE q
      # from the data rather than estimating it, so index_log_q is unused and must
      # be mapped out -- otherwise it is a free parameter that never enters the
      # objective, leaving a flat direction that makes the Hessian singular.
      if(!data_list$fleet_control$Catchability[i] %in% c("Fixed", "Analytical", "AnalyticalArith")){
        map_list$index_log_q[flt] <- flt
      }

      # - Turn on power param for:
      # - 4 = Estimate power equation
      if (data_list$fleet_control$Catchability[i] == "PowerEquation") {
        # map_list$index_q_pow[flt] <- flt
      }

      # Time- varying q parameters "Time_varying_q"
      # - 0 = "Off",
      # - 1 = "IID" penalized deviate or random effect
      # - 2 = "AR1"
      # - 3 = "Block" time blocks with no penalty
      # - 4 = "RandomWalk" random walk from mean following Dorn 2018 (dnorm(q_y - q_y-1, 0, sigma)
      # - If estimate_q == 5 or 6; "Time_varying_q" determines the environmental indices to be used in the equation log(q_y) = q_mu + beta * index_y or to fit to.
      # - Catchability = 6 turns on time-varying deviates

      # -- Set up time varying catchability if used (account for missing years)
      if((data_list$fleet_control$Catchability[i] %in% c("Estimated", "Estimated-with-prior") &
          data_list$fleet_control$Time_varying_q[i] %in% c("IID", "Block", "AR1", "RandomWalk")) |
         data_list$fleet_control$Catchability[i] == "AR1"){

        # Extract survey years where data is provided
        index_data <- data_list$index_data[which(data_list$index_data$Fleet_code == flt & data_list$index_data$Year > data_list$styr & data_list$index_data$Year <= data_list$endyr),]
        block_yrs <- index_data$Year - data_list$styr + 1

        # Penalized deviate or random walk
        if(data_list$fleet_control$Time_varying_q[i] %in% c("IID", "AR1", "RandomWalk")){
          map_list$index_q_dev[flt, yrs_hind] <- ind_q_dev + (1:nyrs_hind) - 1
          ind_q_dev <- ind_q_dev + nyrs_hind
        }

        # Turn on first deviate for random walk
        if(data_list$fleet_control$Time_varying_q[i] == "RandomWalk"){
          map_list$index_q_dev[flt, 1] <- NA
        }

        # Time blocks
        if(data_list$fleet_control$Time_varying_q[i] == "Block"){
          map_list$index_q_dev[flt, block_yrs] <- ind_q_dev + index_data$Selectivity_block - 1
          ind_q_dev <- ind_q_dev + max(index_data$Selectivity_block)
        }
      }

      # - Turn on regression coefficients for:
      # - 5 = Estimate environmental linkage
      # FIXME: use formula
      if (data_list$fleet_control$Catchability[i] == "Environmental") {
        if(nchar(data_list$fleet_control$Time_varying_q[i]) == 1){
          turn_on <- as.numeric(data_list$fleet_control$Time_varying_q[i])
        }else{
          turn_on <- as.numeric(unlist(strsplit(data_list$fleet_control$Time_varying_q[i],","))) # Parameters to turn on
        }

        if(any(turn_on > ncol(data_list$env_data))){
          stop("For 'Environmental' catchability, index specified in 'Time_varying_q' is greater than number of indices in 'env_data'")
        }

        map_list$index_q_beta[flt, turn_on] <- turn_on + ind_beta_q
        ind_beta_q <- ind_beta_q + max(turn_on)
      }

      # - 6 = Fit to environmental index
      if (data_list$fleet_control$Catchability[i] == "AR1") {
        if(!nchar(data_list$fleet_control$Time_varying_q[i]) == 1){
          warning("Cant fit catchability deviates to multiple indices")
        }
        map_list$index_q_beta[flt, 1] <- 1 + ind_beta_q # The effect size
        ind_beta_q <- ind_beta_q + 1

        map_list$index_q_rho[flt] <- flt # Correlation coeff

        # Turn on standard deviations
        map_list$index_q_log_sd[flt] <- flt # Obseration error
        map_list$index_q_dev_log_sd[flt] <- flt # AR1 process error
      }

      # Standard deviation of surveys index
      # - 0 = use CV from index_data
      # - 1 = estimate a free parameter
      # - 2 = analytically estimate following (Ludwig and Walters 1994)
      if (data_list$fleet_control$Estimate_index_sd[i] == 1) {
        map_list$index_log_sd[flt] <- flt
      }
    }
  }
  return(map_list)
}


#' @title Helper to adjust map for shared catchability/selectivity indices
#'
#' @description Enforces parameter sharing by mapping parameters for fleets
#'   with a common \code{Selectivity_index} or \code{Catchability_index} to the
#'   same value as the initial index.
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#'
#' @return Updated \code{map_list}.
adjust_map_shared_params <- function(map_list, data_list) {

  # Based on `Selectivity_index` or `Catchability_index` in `fleet_control`
  sel_index <- data_list$fleet_control$Selectivity_index
  sel_index_tested <- c()

  q_index <- data_list$fleet_control$Catchability_index
  q_index_tested <- c()
  rows_tests <- c()

  # A fleet with Fleet_type "Off" has no estimated parameters, so its (all-NA)
  # map slice must not be copied onto the fleets sharing its index -- that would
  # silently stop estimating their selectivity/catchability. Prefer the first
  # estimated fleet in the group as the donor.
  flt_off <- data_list$fleet_control$Fleet_type == "Off"
  first_est <- function(rows) {
    est <- rows[!flt_off[rows]]
    if (length(est)) est[1] else NA_integer_
  }

  for(i in 1: nrow(data_list$fleet_control)){
    flt = data_list$fleet_control$Fleet_code[i]
    sel_test <- sel_index[i] %in% sel_index_tested
    if(!is.na(q_index[i])){ # Make sure not using fishery data
      q_test <- q_index[i] %in% q_index_tested
    } else {
      q_test <- FALSE
    }

    # If selectivity is the same as a previous index
    if(sel_test){
      sel_duplicate_vec <- c(which(sel_index_tested == sel_index[i]), i)
      sel_duplicate <- first_est(which(sel_index_tested == sel_index[i]))

      # Per-fleet settings a shared Selectivity_index does not reconcile
      # (Selectivity, Selectivity_dimension, Bin_first_selected, N_sel_bins,
      # Sel_norm_bin*, Time_varying_sel) are checked in data_check(): fit_mod()
      # wraps this call in suppressWarnings(), so a warning raised here is never
      # seen.

      # FIXME add checks for surveys sel sigma

      # Make selectivity maps the same if selectivity is the same
      if(!is.na(sel_duplicate)){
        map_list$log_sel_slp[1:2, flt,] <- map_list$log_sel_slp[1:2, sel_duplicate,]
        map_list$sel_inf[1:2, flt,] <- map_list$sel_inf[1:2, sel_duplicate,]
        map_list$sel_coff[flt,,] <- map_list$sel_coff[sel_duplicate,,]
        map_list$sel_coff_dev[flt,,,] <- map_list$sel_coff_dev[sel_duplicate,,,]
        map_list$log_sel_slp_dev[1:2, flt,,] <- map_list$log_sel_slp_dev[1:2, sel_duplicate,,]
        map_list$sel_inf_dev[1:2, flt,,] <- map_list$sel_inf_dev[1:2, sel_duplicate,,]
        map_list$sel_dev_log_sd[flt] <- map_list$sel_dev_log_sd[sel_duplicate]
        map_list$sel_curve_pen[flt,] <- map_list$sel_curve_pen[sel_duplicate,]
      }
    }


    # If catchability is the same as a previous index
    if(q_test){
      q_duplicate_vec <- c(which(q_index_tested == q_index[i]), i)
      q_duplicate <- first_est(which(q_index_tested == q_index[i]))

      # Error check selectivity type
      if(length(unique(data_list$fleet_control$Catchability[q_duplicate_vec])) > 1){
        warning("Survey catchability of surveys with same Catchability_index is not the same")
        warning(paste0("Double check Catchability in fleet_control of surveys:", paste(data_list$fleet_control$Fleet_name[q_duplicate_vec])))
      }


      # Error check time-varying selectivity type
      if(length(unique(data_list$fleet_control$Time_varying_q[q_duplicate_vec])) > 1){
        warning("Time varying survey catchability of surveys with same Catchability_index is not the same")
        warning(paste0("Double check Time_varying_q in fleet_control of surveys:", paste(data_list$fleet_control$Fleet_name[q_duplicate_vec])))
      }

      # FIXME add checks for surveys q sigma

      # Make catchability maps the same.
      #
      # The group shares ONE q parameter, so it can carry only one answer to
      # "is q estimated?", and the LEAD fleet -- first_est(), the group's first
      # non-Off fleet -- decides it for everyone regardless of fleet type. That
      # is intended: a fishery sharing a group whose lead estimates q follows the
      # lead and is estimated too, even with no index_data of its own, and a
      # fishery whose lead is Fixed stays fixed whatever its own Catchability
      # says.
      #
      # TODO: what is NOT intended is that the disagreement is invisible. The
      # two warnings above fire when the group's Catchability / Time_varying_q
      # settings differ, but fit_mod() wraps build_map() in suppressWarnings()
      # (see the TODO there), so a fleet whose own setting was overridden by the
      # lead never hears about it.
      if(!is.na(q_duplicate)){
        map_list$index_log_q[flt] <- map_list$index_log_q[q_duplicate]
        # map_list$index_q_pow[flt] <- map_list$index_q_pow[q_duplicate]
        map_list$index_q_rho[flt] <- map_list$index_q_rho[q_duplicate]
        map_list$index_q_beta[flt,] <- map_list$index_q_beta[q_duplicate,]
        map_list$index_q_dev[flt,] <- map_list$index_q_dev[q_duplicate,]
        map_list$index_q_log_sd[flt] <- map_list$index_q_log_sd[q_duplicate]
        map_list$index_q_dev_log_sd[flt] <- map_list$index_q_dev_log_sd[q_duplicate]
      }
    }


    # Add index
    sel_index_tested <- c(sel_index_tested, sel_index[i])
    q_index_tested <- c(q_index_tested, q_index[i])
  }

  return(map_list)
}

#' @title Helper to set map for Fishing Mortality and Data Weights
#'
#' @description Maps fishing mortality parameters (\code{log_F}) and related targets.
#'   Also maps data weighting parameters (\code{catch_log_sd}, \code{comp_weights}, \code{caal_weights}).
#'
#' @param map_list The current TMB map list.
#' @param data_list The data list containing model settings.
#' @param nyrs_hind Number of historical years.
#'
#' @return Updated \code{map_list}.
build_map_f_and_data_weights <- function(map_list, data_list, nyrs_hind) {

  # -- Map out future fishing mortality
  map_list$proj_F_prop <- map_list$proj_F_prop * NA

  # -- Map out initial F if starting at equilibrium
  if(!(data_list$initMode %in% c("FishedNonEquilibrium", "FishedNonEquilibriumScaled"))){
    map_list$log_Finit <- rep(NA, data_list$nspp)
  }

  # -- FSPR mapped out
  map_list$log_Flimit <- rep(NA, data_list$nspp)
  map_list$log_Ftarget <- rep(NA, data_list$nspp)


  comp_count <- data_list$comp_data |> # Count comp obs by fleet
    dplyr::filter(Year > 0) |>
    dplyr::count(Fleet_code)

  caal_count <- data_list$caal_data |> # Count CAAL obs by fleet
    dplyr::filter(Year > 0) |>
    dplyr::count(Fleet_code)

  for (i in 1:nrow(data_list$fleet_control)) {
    flt = data_list$fleet_control$Fleet_code[i]
    # Fishery catch time-series SD: fix it unless it is being estimated
    if (data_list$fleet_control$Estimate_catch_sd[i] %in% c(NA, 0, 2)) {
      map_list$catch_log_sd[flt] <- NA
    }

    # Turn off F and F dev if not estimating or it is a Survey
    if (data_list$fleet_control$Fleet_type[i] != "Fishery") {
      map_list$catch_log_sd[flt] <- NA
      map_list$log_F[flt, ] <- NA
    }

    # Map out comp weights if using multinomial
    if(data_list$fleet_control$Comp_distribution[i] !=  "DirichletMultinomial") {
      map_list$comp_weights[i] <- NA
    }

    # Map out CAAL weights if using multinomial
    if(data_list$fleet_control$CAAL_distribution[i] != "DirichletMultinomial") {
      map_list$caal_weights[i] <- NA
    }

    # Map out comp weights and sigma if fleet is turned off or there are no comp data
    if(data_list$fleet_control$Fleet_type[i] == "Off") {
      map_list$comp_weights[i] <- NA
      map_list$caal_weights[i] <- NA
      map_list$sel_dev_log_sd[i] <- NA
    }
    if(!data_list$fleet_control$Fleet_code[i] %in% comp_count$Fleet_code){
      map_list$comp_weights[i] <- NA
    }
    if(!data_list$fleet_control$Fleet_code[i] %in% caal_count$Fleet_code){
      map_list$caal_weights[i] <- NA
    }
  }


  # - Map out Fdev and sel-devs for years with 0 catch to very low number
  catch_data <- data_list$catch_data[which(data_list$catch_data$Year <= data_list$endyr),]
  fsh_ind <- catch_data$Fleet_code[which(catch_data$Catch == 0)]
  yr_ind <- catch_data$Year[which(catch_data$Catch == 0)] - data_list$styr + 1

  for(i in 1:length(yr_ind)){
    map_list$log_F[fsh_ind[i], yr_ind[i]] <- NA
    map_list$log_sel_slp_dev[1:2, fsh_ind[i], , yr_ind[i]] <- NA
    map_list$sel_inf_dev[1:2, fsh_ind[i], , yr_ind[i]] <- NA
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

  # For any species whose numbers-at-age are fixed (estDynamics > 0), fix its
  # population and fleet parameters -- there is nothing to estimate for it.
  for(sp in 1:data_list$nspp){

    # Fixed numbers-at-age: fix (map out) most parameters for this species
    if(data_list$estDynamics[sp] > 0){

      # Population parameters
      map_list$rec_pars[sp,] <- NA
      map_list$R_log_sd[sp] <- NA
      map_list$log_Finit[sp] <- NA
      # map_list$sex_ratio_log_sd[sp] <- NA
      map_list$rec_dev[sp,] <- NA
      map_list$init_dev[sp,] <- NA
      map_list$log_M1[sp,,] <- NA
      map_list$log_M1_dev[sp,,,] <- NA
      map_list$M1_dev_log_sd[sp,] <- NA
      map_list$M1_rho[sp,,] <- NA
      map_list$log_Finit[sp] <- NA

      # Survey and fishery fleet parameters
      flts <- data_list$fleet_control$Fleet_code[which(data_list$fleet_control$Species == sp)]


      map_list$log_F[flts,] <- NA
      map_list$index_log_q[flts] <- NA
      # map_list$index_q_pow[flts] <- NA
      map_list$index_q_dev[flts,] <- NA
      map_list$index_q_log_sd[flts] <- NA
      map_list$index_q_dev_log_sd[flts] <- NA
      map_list$sel_coff[flts,,] <- NA
      map_list$sel_coff_dev[flts,,,] <- NA
      map_list$log_sel_slp[, flts, ] <- NA
      map_list$sel_inf[, flts, ] <- NA
      map_list$log_sel_slp_dev[, flts, ,] <- NA
      map_list$sel_inf_dev[, flts, ,] <- NA
      map_list$sel_dev_log_sd[flts] <- NA
      map_list$index_log_sd[flts] <- NA
      map_list$catch_log_sd[flts] <- NA
      map_list$comp_weights[flts] <- NA
      map_list$caal_weights[flts] <- NA
    }

    # Don't estimate the scalar
    if(data_list$estDynamics[sp] < 2 | data_list$msmMode == 0){
      map_list$log_pop_scalar[sp,] <- NA
    }

    # Age-independent scalar
    if(data_list$estDynamics[sp] == 2 | data_list$msmMode != 0){
      map_list$log_pop_scalar[sp,2:ncol(map_list$log_pop_scalar)] <- NA # Only estimate first parameter
    }

    # Age-dependent scalar
    if(data_list$estDynamics[sp] == 3 | data_list$msmMode != 0){
      if(data_list$nages[sp] < ncol(map_list$log_pop_scalar)){ # Map out ages beyond maxage of the species
        map_list$log_pop_scalar[sp,(data_list$nages[sp]+1):ncol(map_list$log_pop_scalar)] <- NA # Only estimate parameters for each age of species
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
  } else{
    map_list$dummy <- NA
  }
  return(map_list)
}


#' @title Helper to set map for linkage-table parameters
#'
#' @description Maps `beta_linkage` (one entry per row of
#'   `data_list$linkage_table`). Rows whose `est_phase == 0` are fixed
#'   at their initial values via `NA`; `(Intercept)` rows are also
#'   fixed (their value stays at 0 -- the base parameter carries the
#'   level). Everything else is estimated.  Phased estimation honoring
#'   nonzero phase ordinals can layer on later via the `phase` argument
#'   to [fit_control()].
#'
#' @param map_list The current TMB map list.
#' @param data_list an Rceattle data_list (with the pooled
#'   `linkage_table` from `pool_linkages()`).
#'
#' @return Updated \code{map_list}.
#' @keywords internal
build_map_linkages <- function(map_list, data_list) {
  if (is.null(data_list$linkage_table) ||
      nrow(data_list$linkage_table) == 0L) {
    return(map_list)
  }
  tbl <- data_list$linkage_table
  est_phase   <- as.integer(tbl$est_phase)
  is_intercept <- tbl$design_col == "(Intercept)"
  # Random-effect indicator rows carry their deviation in beta_linkage_re (which
  # holds the density), so the fixed beta_linkage entry is pinned at 0. Key on
  # re_index -- the registry marker -- so a fixed row that merely inherited a
  # spec-level re_group stays estimable.
  is_re_row <- !is.na(tbl$re_index)
  m <- map_list$beta_linkage
  m[est_phase == 0L | is_intercept | is_re_row] <- NA
  map_list$beta_linkage <- m
  # beta_linkage_re keeps the blanket "all estimable" map (the density damps
  # it), except a random walk fixes its FIRST deviate for identifiability: the
  # walk's mean level is carried by the base parameter the intercept re-targets,
  # exactly as the legacy RandomWalk fixes index_q_dev[flt, 1]. Fixed means held
  # at its `inits` value, which is 0 by default but need not be -- supplying a
  # non-zero first deviate is how a caller reproduces a reference model that
  # estimates it freely, without shifting the level into the base (which would
  # move the point any prior on that base is evaluated at).
  rw_rows <- which(!is.na(tbl$re_struct) & tbl$re_struct == "rw")
  if (length(rw_rows) > 0L) {
    # first slot of each rw group = smallest re_index (earliest time, since
    # re_index is assigned in (group, elapsed-time) order). The slot is a GLOBAL
    # index, so translate it to the vector that actually holds it -- pinning
    # position n of beta_linkage_re when the walk is penalized would fix an
    # unrelated deviation and leave this walk's level unidentified.
    first_slot <- as.integer(tapply(tbl$re_index[rw_rows],
                                    tbl$sigma_index[rw_rows], min))
    rt <- .re_slot_routing(tbl)
    for (s in first_slot) {
      r <- rt[rt$slot == s, , drop = FALSE]
      nm <- if (r$integrate) "beta_linkage_re" else "beta_linkage_re_pen"
      m_re <- map_list[[nm]]
      m_re[r$pos + 1L] <- NA
      map_list[[nm]] <- m_re
    }
  }

  # A group's log_sigma_linkage is estimable unless the spec supplied a fixed
  # input SD (`init = list(sigma = )` with no sigma prior) -- then it is held at
  # that value, reproducing the reference Time_varying_*_sd_prior.
  gt <- .re_group_table(tbl)
  if (!is.null(gt) && any(gt$sigma_fixed)) {
    m_sig <- map_list$log_sigma_linkage
    m_sig[gt$sigma_index[gt$sigma_fixed] + 1L] <- NA   # sigma_index 0-based
    map_list$log_sigma_linkage <- m_sig
  }
  # ar1 correlation: estimable unless the spec supplied a fixed input rho
  # (`init = list(rho = )` with no rho prior). trans_rho_linkage is ordered by
  # ar1 group (gt order among ar1 rows), matching build_params / linkage_re_rho.
  if (!is.null(gt)) {
    ar1 <- gt[gt$re_struct == "ar1", , drop = FALSE]
    if (nrow(ar1) > 0L && any(ar1$rho_fixed)) {
      m_rho <- map_list$trans_rho_linkage
      m_rho[which(ar1$rho_fixed)] <- NA
      map_list$trans_rho_linkage <- m_rho
    }
  }
  # Rogers QAR1 observation SD: one log-SD per observed group. Held FIXED at the
  # input `obs_sd` by default (mapped to NA); estimated when the spec sets
  # `obs_sd_est = TRUE` (a distinct free level). Estimation is opt-in because a
  # freely-estimated obs_sd collapses toward 0 on a smooth covariate (the
  # beta/obs_sd identifiability degeneracy) -- keep it fixed unless the covariate
  # is informative.
  if (!is.null(gt) && any(gt$observed)) {
    obs_est <- gt$obs_sd_est[gt$observed]
    m <- rep(NA_integer_, length(obs_est))
    if (any(obs_est)) m[obs_est] <- seq_len(sum(obs_est))
    map_list$log_obs_sd_linkage <- m
  }
  map_list
}


#' @title Helper to turn off base linked parameters
#'
#' @description Maps the base parameter (`rec_pars`, `log_M1`,
#'   `log_growth_pars`) out of estimation only for stratum groups
#'   whose linkage formula carries *no* intercept. With an intercept
#'   in the formula (`~ 1`, `~ temp`, ...) and a nonzero `est_phase` the
#'   base parameter holds the level and stays estimable; the linkage
#'   `(Intercept)` row is fixed at 0 instead. When an intercept row has
#'   `est_phase == 0` (the documented "fix at init" contract), the base
#'   parameter is mapped off as well, so a fixed intercept truly fixes the
#'   parameter at its `build_params()` initial value. Slope-only formulas
#'   (`~ 0 + temp`) emit no intercept row, so we mask the base parameter to
#'   keep it at its `build_params()` default and let the slope rows define
#'   the year-by-year offset.
#'
#' @param map_list The current TMB map list.
#' @param data_list an Rceattle data_list (with the pooled
#'   `linkage_table` from `pool_linkages()`).
#'
#' @return Updated \code{map_list}.
#' @keywords internal
map_linkage_adjuster <- function(map_list, data_list) {

  tbl <- data_list$linkage_table
  if (is.null(tbl) || nrow(tbl) == 0L) return(map_list)

  group_key <- function(row) {
    paste(row$process, row$param,
          ifelse(is.na(row$species), "*", row$species),
          ifelse(is.na(row$sex),     "*", row$sex),
          ifelse(is.na(row$age_bin), "*", row$age_bin),
          ifelse(is.na(row$fleet),   "*", row$fleet),
          sep = "|")
  }
  keys <- vapply(seq_len(nrow(tbl)),
                 function(i) group_key(tbl[i, , drop = FALSE]),
                 character(1))
  is_intercept <- tbl$design_col == "(Intercept)"
  groups_with_intercept <- unique(keys[is_intercept])

  for (i in seq_len(nrow(tbl))) {
    row <- tbl[i, , drop = FALSE]
    if (is_intercept[i]) {
      # The intercept row carries the base parameter's level. est_phase == 0
      # means "fix at init" (the documented prior/fix/init contract), so the
      # base parameter must be mapped out of estimation as well. With a
      # nonzero phase the base stays estimable and holds the level.
      if (as.integer(row$est_phase) != 0L) next
    } else {
      # Slope rows only mask the base in slope-only groups (no intercept);
      # in intercept-bearing groups the intercept holds the level.
      if (keys[i] %in% groups_with_intercept) next
    }
    idx <- .linkage_row_indices(row, data_list)
    switch(row$process,
      growth = {
        # Mean-growth params live on log_growth_pars[sp, sex, k];
        # SD endpoints live on growth_log_sd[sp, sex, k']. Same
        # slope-only-mask logic applies to both: mask the base so the
        # slope rows in beta_linkage define the offset alone.
        # (SD specs are pre-validated to be intercept-bearing, so this
        # SD branch is only reachable if a future caller bypasses
        # `.validate_growth_linkages`.)
        mean_idx <- .GROWTH_PARAM_TO_INDEX[row$param]
        sd_idx   <- .GROWTH_SD_PARAM_TO_INDEX[row$param]
        if (!is.na(mean_idx)) {
          for (s in idx$species) {
            map_list$log_growth_pars[s, idx$per_sp[[as.character(s)]]$sex, mean_idx] <- NA
          }
        } else if (!is.na(sd_idx)) {
          for (s in idx$species) {
            map_list$growth_log_sd[s, idx$per_sp[[as.character(s)]]$sex, sd_idx] <- NA
          }
        }
      },
      M = {
        for (s in idx$species) {
          sx <- idx$per_sp[[as.character(s)]]$sex
          ag <- idx$per_sp[[as.character(s)]]$age
          map_list$log_M1[s, sx, ag] <- NA
        }
      },
      recruitment = {
        par_idx <- .REC_PARAM_TO_INDEX[row$param]
        if (is.na(par_idx)) next
        map_list$rec_pars[idx$species, par_idx] <- NA
      },
      q = {
        # Catchability is indexed by fleet, so the stratum to mask is
        # idx$fleet rather than the species/sex/age triple.
        map_list$index_log_q[idx$fleet] <- NA
      },
      sel = {
        # Selectivity is indexed by (slot, fleet, sex). A slope-only formula
        # masks the base slot across all sexes of the linked fleet; the coff
        # form masks every bin.
        m <- .SEL_PARAM_TO_SLOT[[row$param]]
        if (!is.null(m)) {
          if (m$arr == "log_sel_slp") {
            map_list$log_sel_slp[m$slot, idx$fleet, ] <- NA
          } else if (m$arr == "sel_inf") {
            map_list$sel_inf[m$slot, idx$fleet, ] <- NA
          } else if (m$arr == "sel_coff") {
            map_list$sel_coff[idx$fleet, , ] <- NA
          }
        }
      }
    )
  }

  map_list
}
