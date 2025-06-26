#' Function to construct the TMB map argument for CEATTLE
#'
#' @description Reads a parameter list to construct map
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param params a parameter list created from \code{\link{build_params}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA except dummy
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE, which sets the map for R_ln_sd to NA
#' @param random_sel logical. If TRUE, treats selectivity deviations as random effects.The default is FALSE, which sets the map for sel_dev_ln_sd to NA. Only viable for logisitc, Double Logistic, Descending Logistic, and Hake Non-parametric with Random walk or deviates.
#'
#' @description
#' TODO: turn on selectivity and catchability deviance variance parameters
#'
#'
#' @return a list of map arguments for each parameter
#' @export
build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE, random_sel = FALSE) {

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Setup ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  data_list <- Rceattle::switch_check(data_list)

  # Get year objects
  nyrs_hind <- data_list$endyr - data_list$styr + 1
  nyrs_proj <- data_list$projyr - data_list$styr + 1
  yrs_proj <- (nyrs_hind + 1):nyrs_proj
  yrs_hind <- 1:nyrs_hind
  if(nyrs_hind == nyrs_proj){
    yrs_proj = NULL
  }

  # Convert parameters to map object and
  # - Set each item in map_list to seperate value
  map_list <- sapply(params, function(x) replace(x, values = c(1:length(x))))


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 1. Recruitment ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # * 1.1. Initial population deviates ----

  # -- Map out first year rec devs if estimating initial abundance as free parameters
  if(data_list$initMode == 0){
    map_list$rec_dev[, 1] <- NA
  }

  # -- Map out initial devs if starting at equilibrium with no devs
  if(data_list$initMode == 1){
    map_list$init_dev[] <- NA
  }

  # -- Map out initial population deviations not to be estimated - map out last age and ages not seen
  for(sp in 1:data_list$nspp) {
    if(data_list$initMode > 1){ # Unfinished or fished equilibrium
      if((data_list$nages[sp] - 1) < ncol(map_list$init_dev)) {
        map_list$init_dev[sp, (data_list$nages[sp]):ncol(map_list$init_dev)] <- NA
      }
    }else{ # Free parameters
      if((data_list$nages[sp]) < ncol(map_list$init_dev)) {
        map_list$init_dev[sp, (data_list$nages[sp]+1):ncol(map_list$init_dev)] <- NA
      }
    }
  }


  # * 1.2. Stock-recruitment parameters ----
  # --  Map out future recruitment deviations
  map_list$rec_dev[, yrs_proj] <- as.numeric(replace(map_list$rec_dev[, yrs_proj],
                                                     values = rep(NA, length(map_list$rec_dev[, yrs_proj]))))

  # -- Recruitment deviation sigmas - turn off if not estimating
  if(random_rec == FALSE){
    map_list$R_ln_sd <- map_list$R_ln_sd * NA
  }

  # -- Stock recruit relationship (SRR) parameters:
  # col1 = mean rec, col2 = SRR alpha, col3 = SRR beta
  # - Turning off 2nd and 3rd par if only using mean rec
  if(data_list$srr_fun %in% c(0, 1) & data_list$srr_pred_fun  %in% c(0, 1)){
    map_list$rec_pars[, 2:3] <- NA
  }

  # - Turning off mean rec par if using SRR
  if(data_list$srr_fun > 1){
    map_list$rec_pars[, 1] <- NA
  }

  # - Fix first parameter in SRR (if SRR not used, will be NA anyway)
  if(data_list$srr_est_mode == 0){
    map_list$rec_pars[, 2] <- NA
  }

  # - Environmental linkages
  #FIXME: make it so the covariates can vary by species
  map_list$beta_rec_pars[] <- NA
  if(data_list$srr_pred_fun %in% c(1, 3, 5)){
    for(sp in 1:data_list$nspp){
      map_list$beta_rec_pars[sp, data_list$srr_indices] <- data_list$srr_indices + (sp-1) * ncol(map_list$beta_rec_pars)
    }
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 2. Natural mortality (M1) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  M1_ind = 1 # Generic indices for looping
  M1_dev_ind = 0
  M1_beta_ind = 0
  M1_dev_ln_sd_ind = 0

  # -- Map out natural mortality parameters
  map_list$ln_M1 <- replace(map_list$ln_M1,
                            values = rep(NA, length(map_list$ln_M1)))
  map_list$M1_beta <- replace(map_list$M1_beta,
                              values = rep(NA, length(map_list$M1_beta)))
  map_list$ln_M1_dev <- replace(map_list$ln_M1_dev,
                                values = rep(NA, length(map_list$ln_M1_dev)))
  map_list$M1_dev_ln_sd <- replace(map_list$M1_dev_ln_sd,
                                   values = rep(NA, length(map_list$M1_dev_ln_sd)))
  map_list$M1_rho <- replace(map_list$M1_rho,
                             values = rep(NA, length(map_list$M1_rho)))

  # -- Loop through and turn on based on model
  for(sp in 1:data_list$nspp){

    # * Fixed effects ----
    # - M1_model = 1: sex- and age-invariant M1
    if(data_list$M1_model[sp] == 1){
      map_list$ln_M1[sp,,1:data_list$nages[sp]] <- M1_ind
      M1_ind = M1_ind + 1
    }

    # - M1_model = 2: sex-specific, but age-invariant M1
    if(data_list$M1_model[sp] == 2){
      map_list$ln_M1[sp,1,1:data_list$nages[sp]] <- M1_ind # Females
      map_list$ln_M1[sp,2,1:data_list$nages[sp]] <- M1_ind + 1 # Males
      M1_ind = M1_ind + 2
    }

    # - M1_model = 3: sex-specific, age-specific M1
    if(data_list$M1_model[sp] == 3){
      if(data_list$nsex[sp] == 1){ # One sex population
        map_list$ln_M1[sp,1, 1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1) # Females (one-sex though)
        map_list$ln_M1[sp,2,] <- map_list$ln_M1[sp,1,]
        M1_ind = M1_ind + data_list$nages[sp]
      }
      if(data_list$nsex[sp] == 2){ # Two sex population
        # Females
        map_list$ln_M1[sp,1, 1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1)
        M1_ind = M1_ind + data_list$nages[sp]

        # Males
        map_list$ln_M1[sp,2, 1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1)
        M1_ind = M1_ind + data_list$nages[sp]
      }
    }

    # - M1_model = 4: environmentally driven sex- and age-invariant M1
    if(data_list$M1_model[sp] == 4){
      # Mean M
      map_list$ln_M1[sp,,1:data_list$nages[sp]] <- M1_ind # Females and Males
      M1_ind = M1_ind + 1

      # - Betas
      map_list$M1_beta[sp,1,data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices
      if(data_list$nsex[sp] == 2){ # Males share M1 betas
        map_list$M1_beta[sp,2,] <- map_list$M1_beta[sp,1,]
      }
      M1_beta_ind = M1_beta_ind + dim(map_list$M1_beta)[3]
    }

    # - M1_model = 5: environmentally driven sex-specific, but age-invariant M1
    if(data_list$M1_model[sp] == 5){
      if(data_list$nsex[sp] == 1){ # One sex population
        # - Mean M
        map_list$ln_M1[sp,,1:data_list$nages[sp]] <- M1_ind
        M1_ind = M1_ind + 1

        # - Betas
        map_list$M1_beta[sp,1,data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices
        M1_beta_ind = M1_beta_ind + dim(map_list$M1_beta)[3]
      }
      if(data_list$nsex[sp] == 2){ # Two sex population
        # - Mean M
        map_list$ln_M1[sp, 1, 1:data_list$nages[sp]] <- M1_ind # Females
        map_list$ln_M1[sp, 2, 1:data_list$nages[sp]] <- M1_ind + 1 # Males
        M1_ind = M1_ind + 2

        # - Betas
        # -- Females
        map_list$M1_beta[sp,1,data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices;
        M1_beta_ind = M1_beta_ind + dim(map_list$M1_beta)[3]

        # -- Males
        map_list$M1_beta[sp,2,data_list$M1_indices] <- M1_beta_ind + data_list$M1_indices;
        M1_beta_ind = M1_beta_ind + dim(map_list$M1_beta)[3]
      }
    }


    # * Random effects ----
    # - M1_re = 0: No random effects (default).
    # - M1_re = 1: Random effects varies by age, but uncorrelated (IID) and constant over years.
    # - M1_re = 2: Random effects varies by year, but uncorrelated (IID) and constant over ages.
    # - M1_re = 3: Random effects varies by year and age, but uncorrelated (IID).
    # - M1_re = 4: Correlated AR1 random effects varies by age, but constant over years.
    # - M1_re = 5: Correlated AR1 random effects varies by year, but constant over ages.
    # - M1_re = 6: Correlated 2D-AR1 random effects varies by year and age.

    # - M1_re = 1/4: Random effects varies by age (IID or AR1) and constant over years.
    if(data_list$M1_re[sp] %in% c(1, 4)){
      if(data_list$M1_model[sp] == 1){ # Sex-invariant
        # - Random effects
        map_list$ln_M1_dev[sp,1, 1:data_list$nages[sp],] <- M1_dev_ind + 1:data_list$nages[sp]

        # Males mapped the same, if present
        if(data_list$nsex[sp] == 2){
          map_list$ln_M1_dev[sp,2,,] <- map_list$ln_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + data_list$nages[sp]
      }

      if(data_list$M1_model[sp] == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$ln_M1_dev[sp, 1, 1:data_list$nages[sp],] <- M1_dev_ind + 1:data_list$nages[sp]
        M1_dev_ind = M1_dev_ind + data_list$nages[sp]

        # -- Males
        map_list$ln_M1_dev[sp, 2, 1:data_list$nages[sp],] <- M1_dev_ind + 1:data_list$nages[sp]
        M1_dev_ind = M1_dev_ind + data_list$nages[sp]
      }


      # - Standard deviation (shared across sexes)
      map_list$M1_dev_ln_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(data_list$M1_re[sp] == 4){
        map_list$M1_rho[sp,,1] =  sp
      }
    }

    # - M1_re = 2/5: Random effects varies by year (IID or AR1) and constant over ages
    if(data_list$M1_re[sp] %in% c(2, 5)){
      if(data_list$M1_model[sp] == 1){ # Sex-invariant
        # - Random effects
        map_list$ln_M1_dev[sp,1,1:data_list$nages[sp], 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = data_list$nages[sp])

        # Males mapped the same, if present
        if(data_list$nsex[sp] == 2){
          map_list$ln_M1_dev[sp,2,,] <- map_list$ln_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + nyrs_hind
      }

      if(data_list$nsex[sp] == 2 & data_list$M1_model[sp] == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$ln_M1_dev[sp,1, 1:data_list$nages[sp], 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = data_list$nages[sp])
        M1_dev_ind = M1_dev_ind + nyrs_hind

        # -- Males
        map_list$ln_M1_dev[sp,2, 1:data_list$nages[sp], 1:nyrs_hind] <- rep(M1_dev_ind + 1:nyrs_hind, each = data_list$nages[sp])
        M1_dev_ind = M1_dev_ind + nyrs_hind
      }

      # - Standard deviation (shared across sexes)
      map_list$M1_dev_ln_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(data_list$M1_re[sp] == 5){
        map_list$M1_rho[sp,,2] =  sp #FIXME: may want sex-varying?? Hard to estimate
      }
    }

    # - M1_re = 3/6: Random effects varies by age and year (IID or 2D-AR1)
    if(data_list$M1_re[sp] %in% c(3, 6)){
      if(data_list$M1_model[sp] == 1){ # Sex-invariant
        # - Random effects
        map_list$ln_M1_dev[sp,1,1:data_list$nages[sp], 1:nyrs_hind] <- M1_dev_ind + (1:nyrs_hind * data_list$nages[sp])

        # Males mapped the same, if present
        if(data_list$nsex[sp] == 2){
          map_list$ln_M1_dev[sp,2,,] <- map_list$ln_M1_dev[sp,1,,]
        }

        M1_dev_ind = M1_dev_ind + (nyrs_hind * data_list$nages[sp])
      }

      if(data_list$nsex[sp] == 2 & data_list$M1_model[sp] == 2){ # Two sex population and sex-specific
        # - Random effects
        # -- Females
        map_list$ln_M1_dev[sp,1, 1:data_list$nages[sp], 1:nyrs_hind] <- M1_dev_ind + (1:nyrs_hind * data_list$nages[sp])
        M1_dev_ind = M1_dev_ind + (nyrs_hind * data_list$nages[sp])

        # -- Males
        map_list$ln_M1_dev[sp,2, 1:data_list$nages[sp], 1:nyrs_hind] <- M1_dev_ind + (1:nyrs_hind * data_list$nages[sp])
        M1_dev_ind = M1_dev_ind + (nyrs_hind * data_list$nages[sp])
      }

      # - Standard deviation (shared across sexes)
      map_list$M1_dev_ln_sd[sp,] = sp

      # AR1 correlation (shared across sexes)
      if(data_list$M1_re[sp] == 5){
        map_list$M1_rho[sp,1,] = M1_dev_ln_sd_ind + 1:2
        map_list$M1_rho[sp,2,] = map_list$M1_rho[sp,1,]
        M1_dev_ln_sd_ind = M1_dev_ln_sd_ind + 2  #FIXME: may want sex-varying?? Hard to estimate
      }
    }
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 3. Predation mortality (M2) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # -- Turn off all predation parameters for single species
  if (data_list$msmMode == 0) { # Single-species

    # Suitability parameters
    map_list$log_gam_a <- map_list$log_gam_a * NA
    map_list$log_gam_b <- map_list$log_gam_b * NA
    map_list$log_phi <- map_list$log_phi * NA

    # # Multispecies kinzey parameters
    # map_list$logH_1 <- map_list$logH_1 * NA
    # map_list$logH_1a <- map_list$logH_1a * NA
    # map_list$logH_1b <- map_list$logH_1b * NA
    #
    # map_list$logH_2 <- map_list$logH_2 * NA
    # map_list$logH_3 <- map_list$logH_3 * NA
    # map_list$H_4 <- map_list$H_4 * NA
  }

  # # * 3.1. Functional form ----
  # # ** MSVPA based predation ----
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
  # # ** Kinzey and Punt predation equations ----
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


  # * 3.2. Suitability ----
  if (data_list$msmMode > 0) {
    # -- Empirical suitability

    for(sp in 1:data_list$nspp){
      if (data_list$suitMode[sp] == 0) {
        # Turn off suitability parameters
        map_list$log_gam_a[sp] <- NA
        map_list$log_gam_b[sp] <- NA
        map_list$log_phi[sp,] <- NA
      }

      # Turn off predation for fixed-prey species
      if(data_list$estDynamics[sp] > 0) {
        map_list$log_phi[,sp] <- NA
      }
    }
  }

  # * 3.4. Diet multiplier ----
  # TODO add in dirichlet multinomial switch
  map_list$diet_comp_weights[] <- NA


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 4. Selectivity ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # SETTINGS
  # "Selectivity" determines shape of selectivity curve:
  # - 0 = empirical selectivity provided in "emp_sel"
  # - 1 = logistic selectivity
  # - 2 = non-parametric selecitivty sensu Ianelli et al 2018
  # - 3 = double logistic
  # - 4 = descending logistic
  # - 5 = non-parametric selectivity sensu Taylor et al 2014 (Hake)

  # "Nselages"	Number of ages to estimate non-parametric selectivity when Selectivity = 2 & 5. Not used otherwise

  # "Time_varying_sel"	determines if time-varying selectivity should be estimated for logistic, double logistic selectivity,  descending logistic , or non-parametric (Selectivity = 1, 3, 4, or 5).
  # - 0 = no
  # - 1 = penalized deviates given "sel_sd_prior"
  # - 3 = time blocks with no penality
  # - 4 = random walk following Dorn
  # - 5 = random walk on ascending portion of double logistic only.
  # NOTE: If selectivity is set to type = 2 (non-parametric) "Sel_sd_prior" will be the 1st penalty on selectivity. "random_sel" treats random deviates and random walk parameters as random effects.

  # -- Selectivity  indices
  ind_coff <- 1
  ind_dev_coff <- 1
  ind_slp <- 1
  ind_inf <- 1
  ind_inf_re <- 1
  ind_slp_re <- 1

  # -- Map out parameters (then turned on)
  # - non-parametric
  map_list$sel_coff <- replace(map_list$sel_coff, values = rep(NA, length(map_list$sel_coff)))
  map_list$sel_coff_dev <- replace(map_list$sel_coff_dev, values = rep(NA, length(map_list$sel_coff_dev)))

  # - logistic and double logistic
  map_list$ln_sel_slp <- replace(map_list$ln_sel_slp, values = rep(NA, length(map_list$ln_sel_slp)))
  map_list$sel_inf <- replace(map_list$sel_inf, values = rep(NA, length(map_list$sel_inf)))

  map_list$ln_sel_slp_dev <- replace(map_list$ln_sel_slp_dev, values = rep(NA, length(map_list$ln_sel_slp_dev)))
  map_list$sel_inf_dev <- replace(map_list$sel_inf_dev, values = rep(NA, length(map_list$sel_inf_dev)))

  # - time-varying selectivity variance
  map_list$sel_dev_ln_sd <- map_list$sel_dev_ln_sd * NA

  # - non-parametric selectivity penalties. Leaving as parameters in case we want to estimate down the line
  map_list$sel_curve_pen <- map_list$sel_curve_pen * NA


  # -- Turn on parameters
  # --- Variance of random effects for selectivity deviates (turn on sigma)
  if(random_sel){
    for (i in 1:nrow(data_list$fleet_control)) {
      flt = data_list$fleet_control$Fleet_code[i]
      # - Logisitc, Double Logistic, Descending Logistic, and Hake Non-parametric
      if (data_list$fleet_control$Selectivity[i] %in% c(1,2,3,4,5) & data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4,5)) {
        map_list$sel_dev_ln_sd[flt] <- flt
      }
    }
  }


  # Loop through fleets
  for (i in 1:nrow(data_list$fleet_control)) {
    flt = data_list$fleet_control$Fleet_code[i]

    # -- Turn off sex-specific parameters if 1 sex model
    spp <- data_list$fleet_control$Species[i]
    nsex <- data_list$nsex[spp]

    if(data_list$fleet_control$Fleet_type[i] > 0){ # If estimating observation model for fleet i

      # * 4.1. Logitistic ----
      # - sel_type = 1
      if (data_list$fleet_control$Selectivity[i] == 1) {

        # Turn on slp and asymptote for each sex
        for(sex in 1:nsex){
          map_list$ln_sel_slp[1, flt, sex] <- ind_slp; ind_slp = ind_slp + 1
          map_list$sel_inf[1, flt, sex] <- ind_inf; ind_inf = ind_inf + 1
        }


        # Turn on time-varying parameters
        # ** Random walk or deviate ----
        if(data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4)){
          for(sex in 1:nsex){
            map_list$ln_sel_slp_dev[1, flt, sex, yrs_hind] <- ind_slp + 1:length(yrs_hind) - 1
            map_list$sel_inf_dev[1, flt, sex, yrs_hind] <- ind_inf + 1:length(yrs_hind) - 1

            ind_slp <- ind_slp + length(yrs_hind)
            ind_inf <- ind_inf + length(yrs_hind)
          }
        }

        # Turn off first deviate parameters for random walk (start at mean)
        # - Ascending
        if(data_list$fleet_control$Time_varying_sel[i] == 4){
          # map_list$ln_sel_slp[, flt,] <- NA
          # map_list$sel_inf[, flt,] <- NA

          map_list$ln_sel_slp_dev[1, flt, sex, 1] <- NA
          map_list$sel_inf_dev[1, flt, sex, 1] <- NA
        }

        # ** Selectivity blocks ----
        if(data_list$fleet_control$Time_varying_sel[i] == 3){

          # If a fishery use the years from the fishery
          if(data_list$fleet_control$Fleet_type[i] == 1){
            catch_data <- data_list$catch_data[which(data_list$catch_data$Fleet_code == flt),]
            Selectivity_block <- catch_data$Selectivity_block
            biom_yrs <- catch_data$Year - data_list$styr + 1

            Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
            biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
          }

          # if a survey use the survey years
          if(data_list$fleet_control$Fleet_type[i] == 2){
            index_data <- data_list$index_data[which(data_list$index_data$Fleet_code == flt),]
            Selectivity_block <- index_data$Selectivity_block
            biom_yrs <- index_data$Year - data_list$styr + 1

            Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
            biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
          }

          # Turn on selectivity blocks
          for(sex in 1:nsex){
            map_list$ln_sel_slp_dev[1, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
            map_list$sel_inf_dev[1, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

            ind_slp <- ind_slp + max(Selectivity_block)
            ind_inf <- ind_inf + max(Selectivity_block)
          }
        }
      }


      # * 4.2. Non-parametric ----
      # - sel_type = 2 (Ianelli et al 20??)
      if(data_list$fleet_control$Selectivity[i] == 2){ # Non-parametric at age

        # Ages to turn on
        # Age_first_selected until (age_first_selected + nselages)
        if(is.na(data_list$fleet_control$Age_first_selected[i])){data_list$fleet_control$Age_first_selected[i] = data_list$minage[spp]}
        ages_on <- (data_list$fleet_control$Age_first_selected[i] - data_list$minage[spp] + 1):data_list$fleet_control$Nselages[i]

        # Turn on parameters for each sex
        for(sex in 1:nsex){
          map_list$sel_coff[flt, sex, ages_on] <- ind_coff + ages_on; ind_coff = ind_coff + max(ages_on)

          # -- Time-varying deviates
          if(data_list$fleet_control$Time_varying_sel[i] == 1){
            map_list$sel_coff[flt, , ] <- NA
            map_list$sel_coff_dev[flt,sex, ages_on, yrs_hind] <- ind_dev_coff + 1:(length(ages_on) * length(yrs_hind))
            ind_dev_coff = ind_dev_coff + (length(ages_on) * length(yrs_hind))
          }

          # Switch check
          if(!data_list$fleet_control$Time_varying_sel[i] %in% c(NA, 0, 1)){
            stop(paste0("'Time_varying_sel' for fleet ", i, " with non-parametric selectivity is not 0 or 1"))
          }
        }
      }



      # * 4.3. Double logistic ----
      # - sel_type = 3
      if(data_list$fleet_control$Selectivity[i] == 3){ # Double logistic

        # Turn on slp and asymptote for each sex
        for(j in 1:2){
          for(sex in 1:nsex){
            map_list$ln_sel_slp[j, flt, sex] <- ind_slp; ind_slp = ind_slp + 1
            map_list$sel_inf[j, flt, sex] <- ind_inf; ind_inf = ind_inf + 1
          }
        }

        # -- Time varying parameters
        # ** Random walk or deviate ----
        if(data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4,5)){
          for(j in 1:2){
            for(sex in 1:nsex){
              map_list$ln_sel_slp_dev[j, flt, sex,yrs_hind] <- ind_slp + 1:length(yrs_hind) - 1
              map_list$sel_inf_dev[j, flt, sex,yrs_hind] <- ind_inf + 1:length(yrs_hind) - 1

              ind_slp <- ind_slp + length(yrs_hind)
              ind_inf <- ind_inf + length(yrs_hind)
            }
          }

          # If only doing the ascending portion
          if(data_list$fleet_control$Time_varying_sel[i] %in% c(5)){
            map_list$ln_sel_slp_dev[2, flt, ,] <- NA
            map_list$sel_inf_dev[2, flt, ,] <- NA
          }
        }

        # Turn off first deviate for random walk
        # - Ascending and descending
        if(data_list$fleet_control$Time_varying_sel[i] == 4){
          map_list$ln_sel_slp_dev[j, flt, sex, 1] <- NA
          map_list$sel_inf_dev[j, flt, sex, 1] <- NA
        }

        # - Ascending only
        if(data_list$fleet_control$Time_varying_sel[i] == 5){
          map_list$ln_sel_slp_dev[1, flt, sex, 1] <- NA
          map_list$sel_inf_dev[1, flt, sex, 1] <- NA
        }


        # ** Selectivity blocks ----
        if(data_list$fleet_control$Time_varying_sel[i] == 3){

          # If a fishery use the years from the fishery
          if(data_list$fleet_control$Fleet_type[i] == 1){
            catch_data <- data_list$catch_data[which(data_list$catch_data$Fleet_code == flt),]
            Selectivity_block <- catch_data$Selectivity_block
            biom_yrs <- catch_data$Year - data_list$styr + 1

            Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
            biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
          }

          # if a survey use the survey years
          if(data_list$fleet_control$Fleet_type[i] == 2){
            index_data <- data_list$index_data[which(data_list$index_data$Fleet_code == flt),]
            Selectivity_block <- index_data$Selectivity_block
            biom_yrs <- index_data$Year - data_list$styr + 1

            Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
            biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
          }

          # Loop through upper and lower
          for(j in 1:2){
            for(sex in 1:nsex){
              map_list$ln_sel_slp_dev[j, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
              map_list$sel_inf_dev[j, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

              ind_slp <- ind_slp + max(Selectivity_block)
              ind_inf <- ind_inf + max(Selectivity_block)
            }
          }
        }
      }


      # * 4.4. Descending logitistic ----
      # - sel_type = 4
      if (data_list$fleet_control$Selectivity[i] == 4) {

        # Turn on descending slp and asymptote for each sex
        for(sex in 1:nsex){
          map_list$ln_sel_slp[2, flt, sex] <- ind_slp; ind_slp = ind_slp + 1
          map_list$sel_inf[2, flt, sex] <- ind_inf; ind_inf = ind_inf + 1
        }

        # ** Random walk or deviate ----
        if(data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4)){
          for(sex in 1:nsex){
            map_list$ln_sel_slp_dev[2, flt, sex, yrs_hind] <- ind_slp + 1:length(yrs_hind) - 1
            map_list$sel_inf_dev[2, flt, sex, yrs_hind] <- ind_inf + 1:length(yrs_hind) - 1

            ind_slp <- ind_slp + length(yrs_hind)
            ind_inf <- ind_inf + length(yrs_hind)
          }
        }

        # Turn off first deviate for random walk
        # - Descending
        if(data_list$fleet_control$Time_varying_sel[i] == 4){
          map_list$ln_sel_slp_dev[2, flt, sex, 1] <- NA
          map_list$sel_inf_dev[2, flt, sex, 1] <- NA
        }

        # ** Selectivity blocks ----
        if(data_list$fleet_control$Time_varying_sel[i] == 3){

          # If a fishery use the years from the fishery
          if(data_list$fleet_control$Fleet_type[i] == 1){
            catch_data <- data_list$catch_data[which(data_list$catch_data$Fleet_code == flt),]
            Selectivity_block <- catch_data$Selectivity_block
            biom_yrs <- catch_data$Year - data_list$styr + 1

            Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
            biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
          }

          # if a survey use the survey years
          if(data_list$fleet_control$Fleet_type[i] == 2){
            index_data <- data_list$index_data[which(data_list$index_data$Fleet_code == flt),]
            Selectivity_block <- index_data$Selectivity_block
            biom_yrs <- index_data$Year - data_list$styr + 1

            Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
            biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
          }

          # Map selectivity blocks
          for(sex in 1:nsex){
            map_list$ln_sel_slp_dev[2, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
            map_list$sel_inf_dev[2, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

            ind_slp <- ind_slp + max(Selectivity_block)
            ind_inf <- ind_inf + max(Selectivity_block)
          }
        }
      }


      # * 4.5. Non-parametric similar to Hake ----
      # (Taylor et al 2014) - sel_type = 5
      if(data_list$fleet_control$Selectivity[i] == 5){ # Non-parametric at age
        # Ages to turn on
        # Age_first_selected until (age_first_selected + nselages)
        if(is.na(data_list$fleet_control$Age_first_selected[i])){data_list$fleet_control$Age_first_selected[i] = data_list$minage[spp]}
        ages_on <- (data_list$fleet_control$Age_first_selected[i] - data_list$minage[spp] + 2):data_list$fleet_control$Nselages[i] # + 2 because first parameter is not-identifiable and is not estimated

        # Turn on parameters for each sex
        for(sex in 1:nsex){
          map_list$sel_coff[flt, sex, ages_on] <- ind_coff + ages_on; ind_coff = ind_coff + max(ages_on)

          # -- time-varying deviates
          if(data_list$fleet_control$Time_varying_sel[i] == 1){
            map_list$sel_coff_dev[flt, sex, ages_on, yrs_hind] <- ind_dev_coff + 1:(length(ages_on) * length(yrs_hind))
            ind_dev_coff = ind_dev_coff + (length(ages_on) * length(yrs_hind))
          }

          if(!(data_list$fleet_control$Time_varying_sel[i] %in% c(NA, 0, 1))){
            warning(paste("Time_varying_sel for fleet", i, "is not compatible (select NA, 0, or 1)"))
          }
        }
      }
    }
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 5. Catchability ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # -- Catchability indices
  ind_q_dev <- 1
  ind_beta_q <- 0


  catchability_params <- c("index_ln_q", "index_q_beta", "index_q_rho", "index_q_dev", "index_q_ln_sd", "index_q_dev_ln_sd", "index_ln_sd") # "index_q_pow"
  map_list[catchability_params] <- lapply(map_list[catchability_params], function(x) replace(x, values = rep(NA, length(x))))

  # Loop through fleets
  for( i in 1: nrow(data_list$fleet_control)){
    flt = data_list$fleet_control$Fleet_code[i]

    if(data_list$fleet_control$Fleet_type[flt] == 2){ # If survey
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
      if(data_list$fleet_control$Estimate_q[i] %in% c(1, 2, 4, 5, 6)){
        map_list$index_ln_q[flt] <- flt
      }

      # - Turn on power param for:
      # - 4 = Estimate power equation
      if (data_list$fleet_control$Estimate_q[i] %in% c(4)) {
        # map_list$index_q_pow[flt] <- flt
      }

      # Time- varying q parameters "Time_varying_q"
      # - 0 = no,
      # - 1 = penalized deviate
      # - 2 = random effect
      # - 3 = time blocks with no penalty
      # - 4 = random walk from mean following Dorn 2018 (dnorm(q_y - q_y-1, 0, sigma)
      # - If estimate_q == 5 or 6; "Time_varying_q" determines the environmental indices to be used in the equation log(q_y) = q_mu + beta * index_y or to fit to.
      # - Estimate_q = 6 turns on time-varying deviates

      # -- Set up time varying catchability if used (account for missing years)
      if((data_list$fleet_control$Estimate_q[i] %in% c(1, 2) &
          as.numeric(data_list$fleet_control$Time_varying_q[i]) %in% c(1, 2, 3, 4)) |
         data_list$fleet_control$Estimate_q[i] == 6){

        # Extract survey years where data is provided
        index_data <- data_list$index_data[which(data_list$index_data$Fleet_code == flt & data_list$index_data$Year > data_list$styr & data_list$index_data$Year <= data_list$endyr),]
        srv_biom_yrs <- index_data$Year - data_list$styr + 1

        # Penalized deviate or random walk
        if(data_list$fleet_control$Time_varying_q[i] %in% c(1,2,4)){
          map_list$index_q_dev[flt, yrs_hind] <- ind_q_dev + (1:nyrs_hind) - 1
          ind_q_dev <- ind_q_dev + nyrs_hind
        }

        # Turn on first deviate for random walk
        if(data_list$fleet_control$Time_varying_q[i] == 4){
          map_list$index_q_dev[flt, 1] <- NA
        }

        # Time blocks
        if(data_list$fleet_control$Time_varying_q[i] == 3){
          map_list$index_q_dev[flt, srv_biom_yrs] <- ind_q_dev + index_data$Selectivity_block - 1
          ind_q_dev <- ind_q_dev + max(index_data$Selectivity_block)
        }
      }

      # - Turn on regression coefficients for:
      # - 5 = Estimate environmental linkage
      if (data_list$fleet_control$Estimate_q[i] == 5) {
        if(nchar(data_list$fleet_control$Time_varying_q[i]) == 1){
          turn_on <- as.numeric(data_list$fleet_control$Time_varying_q[i])
        }else{
          turn_on <- as.numeric(unlist(strsplit(data_list$fleet_control$Time_varying_q[i],","))) # Parameters to turn on
        }
        map_list$index_q_beta[flt, turn_on] <- turn_on + ind_beta_q
        ind_beta_q <- ind_beta_q + max(turn_on)
      }

      # - 6 = Fit to environmental index
      if (data_list$fleet_control$Estimate_q[i] == 6) {
        if(!nchar(data_list$fleet_control$Time_varying_q[i]) == 1){
          warning("Cant fit catchability deviates to multiple indices")
        }
        map_list$index_q_beta[flt, 1] <- 1 + ind_beta_q # The effect size
        ind_beta_q <- ind_beta_q + 1

        map_list$index_q_rho[flt] <- flt # Correlation coeff

        # Turn on standard deviations
        map_list$index_q_ln_sd[flt] <- flt # Obseration error
        map_list$index_q_dev_ln_sd[flt] <- flt # AR1 process error
      }

      # Standard deviation of surveys index
      # - 0 = use CV from index_data
      # - 1 = estimate a free parameter
      # - 2 = analytically estimate following (Ludwig and Walters 1994)
      if (data_list$fleet_control$Estimate_index_sd[i] == 1) {
        map_list$index_ln_sd[flt] <- flt
      }
    }
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 6. Adjust map  for shared q/selectivity ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Based on `Selectivity_index` or `Q_index` in `fleet_control`
  sel_index <- data_list$fleet_control$Selectivity_index
  sel_index_tested <- c()

  q_index <- data_list$fleet_control$Q_index
  q_index_tested <- c()
  rows_tests <- c()

  for(i in 1: nrow(data_list$fleet_control)){
    flt = data_list$fleet_control$Fleet_code[i]
    sel_test <- sel_index[flt] %in% sel_index_tested
    if(!is.na(q_index[flt])){ # Make sure not using fishery data
      q_test <- q_index[flt] %in% q_index_tested
    } else {
      q_test <- FALSE
    }

    # If selectivity is the same as a previous index
    if(sel_test){
      sel_duplicate <- which(sel_index_tested == sel_index[flt])[1]
      sel_duplicate_vec <- c(which(sel_index_tested == sel_index[flt]), flt)

      # Error check selectivity type
      if(length(unique(data_list$fleet_control$Selectivity[sel_duplicate_vec])) > 1){
        warning("Survey selectivity of surveys with same Selectivity_index is not the same")
        warning(paste0("Double check Selectivity in fleet_control of surveys:", paste(data_list$fleet_control$Fleet_name[sel_duplicate_vec])))
      }


      # Error check time-varying selectivity type
      if(length(unique(data_list$fleet_control$Time_varying_sel[sel_duplicate_vec])) > 1){
        warning("Time varying survey selectivity of surveys with same Selectivity_index is not the same")
        warning(paste0("Double check Time_varying_sel in fleet_control of surveys:", paste(data_list$fleet_control$Fleet_name[sel_duplicate_vec])))
      }

      # FIXME add checks for surveys sel sigma

      # Make selectivity maps the same if selectivity is the same
      map_list$ln_sel_slp[1:2, flt,] <- map_list$ln_sel_slp[1:2, sel_duplicate,]
      map_list$sel_inf[1:2, flt,] <- map_list$sel_inf[1:2, sel_duplicate,]
      map_list$sel_coff[flt,,] <- map_list$sel_coff[sel_duplicate,,]
      map_list$sel_coff_dev[flt,,,] <- map_list$sel_coff_dev[sel_duplicate,,,]
      map_list$ln_sel_slp_dev[1:2, flt,,] <- map_list$ln_sel_slp_dev[1:2, sel_duplicate,,]
      map_list$sel_inf_dev[1:2, flt,,] <- map_list$sel_inf_dev[1:2, sel_duplicate,,]
      map_list$sel_dev_ln_sd[flt] <- map_list$sel_dev_ln_sd[sel_duplicate]
      map_list$sel_curve_pen[flt,] <- map_list$sel_curve_pen[sel_duplicate,]
    }


    # If catchability is the same as a previous index
    if(q_test){
      q_duplicate <- which(q_index_tested == q_index[flt])[1]
      q_duplicate_vec <- c(which(q_index_tested == q_index[flt]), flt)

      # Error check selectivity type
      if(length(unique(data_list$fleet_control$Estimate_q[q_duplicate_vec])) > 1){
        warning("Survey catchability of surveys with same Q_index is not the same")
        warning(paste0("Double check Estimate_q in fleet_control of surveys:", paste(data_list$fleet_control$Fleet_name[q_duplicate_vec])))
      }


      # Error check time-varying selectivity type
      if(length(unique(data_list$fleet_control$Time_varying_q[q_duplicate_vec])) > 1){
        warning("Time varying survey catchability of surveys with same Q_index is not the same")
        warning(paste0("Double check Time_varying_q in fleet_control of surveys:", paste(data_list$fleet_control$Fleet_name[q_duplicate_vec])))
      }

      # FIXME add checks for surveys q sigma

      # Make catchability maps the same
      map_list$index_ln_q[flt] <- map_list$index_ln_q[q_duplicate]
      map_list$index_ln_q[flt] <- map_list$index_ln_q[q_duplicate]
      # map_list$index_q_pow[flt] <- map_list$index_q_pow[q_duplicate]
      map_list$index_q_rho[flt] <- map_list$index_q_rho[q_duplicate]
      map_list$index_q_beta[flt,] <- map_list$index_q_beta[q_duplicate,]
      map_list$index_q_dev[flt,] <- map_list$index_q_dev[q_duplicate,]
      map_list$index_q_ln_sd[flt] <- map_list$index_q_ln_sd[q_duplicate]
      map_list$index_q_dev_ln_sd[flt] <- map_list$index_q_dev_ln_sd[q_duplicate]
    }


    # Add index
    sel_index_tested <- c(sel_index_tested, sel_index[flt])
    q_index_tested <- c(q_index_tested, q_index[flt])
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 7. Fishing mortality and data weights ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # -- Map out future fishing mortality
  map_list$proj_F_prop <- map_list$proj_F_prop * NA

  # -- Map out initial F if starting at equilibrium
  if(!(data_list$initMode %in% c(3,4))){
    map_list$ln_Finit <- rep(NA, data_list$nspp)
  }

  # -- FSPR mapped out
  map_list$ln_Flimit <- rep(NA, data_list$nspp)
  map_list$ln_Ftarget <- rep(NA, data_list$nspp)


  comp_count <- data_list$comp_data %>% # Count comp obs by fleet
    dplyr::filter(Year > 0) %>%
    dplyr::count(Fleet_code)

  for (i in 1:nrow(data_list$fleet_control)) {
    flt = data_list$fleet_control$Fleet_code[i]
    # Standard deviation of fishery time series If not estimating turn of
    if (data_list$fleet_control$Estimate_catch_sd[i] %in% c(NA, 0, 2)) {
      map_list$catch_ln_sd[flt] <- NA
    }

    # Turn of F and F dev if not estimating of it is a Survey
    if (data_list$fleet_control$Fleet_type[i] %in% c(0, 2)) {
      map_list$catch_ln_sd[flt] <- NA
      map_list$ln_F[flt, ] <- NA
    }

    # Map out comp weights if using multinomial
    if(data_list$fleet_control$Comp_loglike[i] != 1) {
      map_list$comp_weights[i] <- NA
    }

    # Map out comp weights if fleet is turned off or there are no comp data
    if(data_list$fleet_control$Fleet_type[i] == 0) {
      map_list$comp_weights[i] <- NA
    }
    if(!data_list$fleet_control$Fleet_code[i] %in% comp_count$Fleet_code){
      map_list$comp_weights[i] <- NA
    }

    if(!data_list$fleet_control$Comp_loglike[i] %in% c(-1, 0, 1)){
      if(!is.na(data_list$fleet_control$Comp_loglike[i])){
        stop(paste0("Comp_loglike for fleet", i, "is not 0 or 1"))
      }
    }
  }


  # - Map out Fdev for years with 0 catch to very low number
  catch_data <- data_list$catch_data[which(data_list$catch_data$Year <= data_list$endyr),]
  fsh_ind <- catch_data$Fleet_code[which(catch_data$Catch == 0)]
  yr_ind <- catch_data$Year[which(catch_data$Catch == 0)] - data_list$styr + 1

  for(i in 1:length(yr_ind)){
    map_list$ln_F[fsh_ind[i], yr_ind[i]] <- NA
    map_list$ln_sel_slp_dev[1:2, fsh_ind[i], , yr_ind[i]] <- NA
    map_list$sel_inf_dev[1:2, fsh_ind[i], , yr_ind[i]] <- NA
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 8. Set up fixed n-at-age ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # - I.E. turn off all parameters besides for species
  for(sp in 1:data_list$nspp){

    # Fixed n-at-age: Turn off most parameters
    if(data_list$estDynamics[sp] > 0){

      # Population parameters
      map_list$rec_pars[sp,] <- NA
      map_list$R_ln_sd[sp] <- NA
      map_list$ln_Finit[sp] <- NA
      # map_list$sex_ratio_ln_sd[sp] <- NA
      map_list$rec_dev[sp,] <- NA
      map_list$init_dev[sp,] <- NA
      map_list$ln_M1[sp,,] <- NA
      map_list$ln_M1_dev[sp,,,] <- NA
      map_list$M1_dev_ln_sd[sp,] <- NA
      map_list$M1_rho[sp,,] <- NA
      map_list$ln_Finit[sp] <- NA

      # Survey and fishery fleet parameters
      flts <- data_list$fleet_control$Fleet_code[which(data_list$fleet_control$Species == sp)]


      map_list$ln_F[flts,] <- NA
      map_list$index_ln_q[flts] <- NA
      # map_list$index_q_pow[flts] <- NA
      map_list$index_q_dev[flts,] <- NA
      map_list$index_q_ln_sd[flts] <- NA
      map_list$index_q_dev_ln_sd[flts] <- NA
      map_list$sel_coff[flts,,] <- NA
      map_list$sel_coff_dev[flts,,,] <- NA
      map_list$ln_sel_slp[, flts, ] <- NA
      map_list$sel_inf[, flts, ] <- NA
      map_list$ln_sel_slp_dev[, flts, ,] <- NA
      map_list$sel_inf_dev[, flts, ,] <- NA
      map_list$sel_dev_ln_sd[flts] <- NA
      map_list$index_ln_sd[flts] <- NA
      map_list$catch_ln_sd[flts] <- NA
      map_list$comp_weights[flts] <- NA
    }

    # Don't estimate the scalar
    if(data_list$estDynamics[sp] < 2 | data_list$msmMode == 0){
      map_list$ln_pop_scalar[sp,] <- NA
    }

    # Age-independent scalar
    if(data_list$estDynamics[sp] == 2 | data_list$msmMode != 0){
      map_list$ln_pop_scalar[sp,2:ncol(map_list$ln_pop_scalar)] <- NA # Only estimate first parameter
    }

    # Age-dependent scalar
    if(data_list$estDynamics[sp] == 3 | data_list$msmMode != 0){
      if(data_list$nages[sp] < ncol(map_list$ln_pop_scalar)){ # Map out ages beyond maxage of the species
        map_list$ln_pop_scalar[sp,(data_list$nages[sp]+1):ncol(map_list$ln_pop_scalar)] <- NA # Only estimate parameters for each age of species
      }
    }
  }






  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 9. Debug ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # - I.E. turn off all parameters besides dummy
  map_list$dummy <- NA
  if(debug){
    map_list <- sapply(map_list, function(x) replace(x, values = rep(NA, length(x))))
    map_list$dummy = 1
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 10. Convert to factor ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  map_list_grande <- list()
  map_list_grande$mapFactor <- sapply(map_list, factor)
  map_list_grande$mapList <- map_list

  return(map_list_grande)
}
