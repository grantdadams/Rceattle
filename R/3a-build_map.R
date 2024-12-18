#' Function to construct the TMB map argument for CEATTLE
#'
#' @description Reads a parameter list to construct map
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param params a parameter list created from \code{\link{build_params}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA except dummy
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects.The default is FALSE, which sets the map for ln_rec_sigma to NA
#' @param random_sel logical. If TRUE, treats selectivity deviations as random effects.The default is FALSE, which sets the map for ln_sigma_sel to NA. Only viable for logisitc, Double Logistic, Descending Logistic, and Hake Non-parametric with Random walk or deviates.
#'
#' @return a list of map arguments for each parameter
#' @export
build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE, random_sel = FALSE) {
  # functions
  '%!in%' <- function(x,y)!('%in%'(x,y))

  # -----------------------------------------------------------
  # -- Setup
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


  # -----------------------------------------------------------
  # STEP 1: Base population dynamics ----
  # -----------------------------------------------------------
  # -- 1.1. Map out future fishing mortality and sex ratio variance
  map_list$proj_F_prop <- map_list$proj_F_prop * NA
  map_list$ln_sex_ratio_sigma <- map_list$ln_sex_ratio_sigma * NA

  # -- 1.2. Map out comp weights
  map_list$comp_weights <- map_list$comp_weights * NA

  # -- 1.3. Map out future recruitment deviations
  map_list$rec_dev[, yrs_proj] <- as.numeric(replace(map_list$rec_dev[, yrs_proj],
                                                     values = rep(NA, length(map_list$rec_dev[, yrs_proj]))))

  # -- Map out first year rec devs if estimating initial abundance as free parameters
  if(data_list$initMode == 0){
    map_list$rec_dev[, 1] <- NA
  }

  if(data_list$initMode != 2){
    map_list$ln_Finit <- rep(NA, data_list$nspp)
  }

  # -- FSPR mapped out
  map_list$ln_Flimit <- rep(NA, data_list$nspp)
  map_list$ln_Ftarget <- rep(NA, data_list$nspp)


  # -- 1.4. Map out initial population deviations not to be estimated - map out last age and ages not seen
  for(sp in 1:data_list$nspp) {
    if(data_list$initMode > 0){ # Unfinished or fished equilibrium
      if((data_list$nages[sp] - 1) < ncol(map_list$init_dev)) {
        map_list$init_dev[sp, (data_list$nages[sp]):ncol(map_list$init_dev)] <- NA
      }
    }else{ # Free parameters
      if((data_list$nages[sp]) < ncol(map_list$init_dev)) {
        map_list$init_dev[sp, (data_list$nages[sp]+1):ncol(map_list$init_dev)] <- NA
      }
    }
  }


  # -- 1.5. Map out natural mortality
  M1_ind = 1 # Generic ind for looping
  for(sp in 1:data_list$nspp){
    # Dim = nspp, nsex (2), nages

    # Turn off all
    map_list$ln_M1[sp,,] <- NA

    # Turn on
    # - M1_model = 1: sex- and age-invariant M1
    if(data_list$M1_model[sp] == 1){
      map_list$ln_M1[sp,,1:data_list$nages[sp]] <- M1_ind
      M1_ind = M1_ind + 1
    }

    # - M1_model = 2: sex-specific, age-invariant M1
    if(data_list$M1_model[sp] == 2){
      if(data_list$nsex[sp] == 1){ # One sex population
        map_list$ln_M1[sp,,1:data_list$nages[sp]] <- M1_ind
        M1_ind = M1_ind + 1
      }
      if(data_list$nsex[sp] == 2){ # Two sex population
        map_list$ln_M1[sp,1,1:data_list$nages[sp]] <- M1_ind # Females
        map_list$ln_M1[sp,2,1:data_list$nages[sp]] <- M1_ind + 1 # Males
        M1_ind = M1_ind + 2
      }
    }

    # - M1_model = 3: sex-specific, age-specific M1
    if(data_list$M1_model[sp] == 3){
      if(data_list$nsex[sp] == 1){ # One sex population
        map_list$ln_M1[sp,1,1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1) # Females (one-sex though)
        map_list$ln_M1[sp,2,1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1) # Males (one-sex though)
        M1_ind = M1_ind + data_list$nages[sp]
      }
      if(data_list$nsex[sp] == 2){ # Two sex population
        # Females
        map_list$ln_M1[sp,1,1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1)
        M1_ind = M1_ind + data_list$nages[sp]

        # Males
        map_list$ln_M1[sp,2,1:data_list$nages[sp]] <- M1_ind: (M1_ind + data_list$nages[sp] - 1)
        M1_ind = M1_ind + data_list$nages[sp]
      }
    }
  }


  # -----------------------------------------------------------
  # Selectivity ----
  # -----------------------------------------------------------
  # -- Selectivity  inds
  ind_slp <- 1
  ind_inf <- 1
  ind_inf_re <- 1
  ind_slp_re <- 1

  # -- Random effects for selectivity
  # -- Map out time-varying selectivity variance
  map_list$ln_sigma_sel <- map_list$ln_sigma_sel * NA

  if(random_sel){
    for (i in 1:nrow(data_list$fleet_control)) {
      flt = data_list$fleet_control$Fleet_code[i]
      # - Logisitc, Double Logistic, Descending Logistic, and Hake Non-parametric
      # - Random walk or deviate
      if (data_list$fleet_control$Selectivity[i] %in% c(1,3,4,5) & data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4)) {
        map_list$ln_sigma_sel[flt] <- flt
      }
    }
  }


  # -- Map out non-parametric selectivity penalties. Leaving as parameters in case we want to estimate down the line
  map_list$sel_curve_pen <- map_list$sel_curve_pen * NA

  # Loop through fleets
  for (i in 1:nrow(data_list$fleet_control)) {
    flt = data_list$fleet_control$Fleet_code[i]

    # -- Turn off sex-specific parameters if 1 sex model
    spp <- data_list$fleet_control$Species[i]
    nsex <- data_list$nsex[spp]

    # -- Turn of selectivity deviates for years not in composition data
    # comp_sub <- data_list$comp_data[which(data_list$comp_data$Fleet_code == flt),]
    # comp_sub <- comp_sub[which(comp_sub$Year <= data_list$endyr),]
    # comp_sub_yrs <- unique(comp_sub$Year) - data_list$styr + 1
    # comp_sub_yrs <- comp_sub_yrs[which(comp_sub_yrs > 0)]
    # comp_miss_yrs <- yrs_hind[which(yrs_hind %!in% comp_sub_yrs)]

    # -- 4.0.  Empirical or not Fit - sel_type = 0
    if (data_list$fleet_control$Selectivity[i] == 0 | data_list$fleet_control$Fleet_type[i] == 0) {

      # Map out non-parametric
      map_list$sel_coff[flt,, ] <- NA
      map_list$sel_coff_dev[flt,,,] <- NA

      # Map out logistic and double logistic
      map_list$ln_sel_slp[1:2, flt, ] <- NA
      map_list$sel_inf[1:2, flt, ] <- NA

      # Map out logistic and double logistic deviates
      map_list$ln_sel_slp_dev[1:2, flt, ,] <- NA
      map_list$sel_inf_dev[1:2, flt, ,] <- NA

      # Map out selectivity var
      map_list$ln_sigma_sel[flt] <- NA
    }


    # -- 4.1. Logitistic - sel_type = 1
    if (data_list$fleet_control$Selectivity[i] == 1) {

      # Map out non-parametric
      map_list$sel_coff[flt,, ] <- NA
      map_list$sel_coff_dev[flt,,,] <- NA

      # Map out double logistic
      map_list$ln_sel_slp[2, flt, ] <- NA
      map_list$sel_inf[2, flt, ] <- NA

      # Map out double logistic deviates
      map_list$ln_sel_slp_dev[2, flt, ,] <- NA
      map_list$sel_inf_dev[2, flt, ,] <- NA

      # Map out time varying parameters if not used
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(NA, 0)){
        map_list$ln_sel_slp_dev[1, flt, ,] <- NA
        map_list$sel_inf_dev[1, flt, ,] <- NA
      }

      # Penalized deviate or random walk
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4)){
        for(sex in 1:nsex){
          map_list$ln_sel_slp_dev[1, flt, sex, ] <- ind_slp + 1:length(map_list$ln_sel_slp_dev[1, flt, sex, ]) - 1
          map_list$sel_inf_dev[1, flt, sex, ] <- ind_inf + 1:length(map_list$sel_inf_dev[1, flt, sex, ]) - 1

          ind_slp <- ind_slp + length(map_list$ln_sel_slp_dev[1, flt, sex, ])
          ind_inf <- ind_inf + length(map_list$sel_inf_dev[1, flt, sex, ])
        }
      }

      # Map selectivity blocks together
      if(data_list$fleet_control$Time_varying_sel[i] == 3){

        # If a fishery use the years from the fishery
        if(data_list$fleet_control$Fleet_type[i] == 1){
          fsh_biom <- data_list$fsh_biom[which(data_list$fsh_biom$Fleet_code == flt),]
          Selectivity_block <- fsh_biom$Selectivity_block
          biom_yrs <- fsh_biom$Year - data_list$styr + 1

          Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
          biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
        }

        # if a survey use the survey years
        if(data_list$fleet_control$Fleet_type[i] == 2){
          srv_biom <- data_list$srv_biom[which(data_list$srv_biom$Fleet_code == flt),]
          Selectivity_block <- srv_biom$Selectivity_block
          biom_yrs <- srv_biom$Year - data_list$styr + 1

          Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
          biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
        }

        # Map out double logistic
        map_list$ln_sel_slp[1:2, flt ,] <- NA
        map_list$sel_inf[1:2, flt ,] <- NA

        # Map selectivity blocks
        for(sex in 1:nsex){
          map_list$ln_sel_slp_dev[1, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
          map_list$sel_inf_dev[1, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

          ind_slp <- ind_slp + max(Selectivity_block)
          ind_inf <- ind_inf + max(Selectivity_block)
        }
      }
    }


    # -- 4.2. Non-parametric - sel_type = 2 (Ianelli et al 20??)
    if(data_list$fleet_control$Selectivity[i] == 2){ # Non-parametric at age
      # -- Map out time-varying deviates
      map_list$sel_coff_dev[flt,,,] <- NA

      # Map out nselages  < max(nselages)
      if(data_list$fleet_control$Nselages[i] < max(data_list$fleet_control$Nselages, na.rm = TRUE)){
        mapped_ages <- (data_list$fleet_control$Nselages[i] + 1):max(data_list$fleet_control$Nselages, na.rm = T)
        map_list$sel_coff[flt, , mapped_ages]  <- NA
      }

      # -- Map out ages < Amin
      if(!is.na(data_list$fleet_control$Age_first_selected[i])){
        if(data_list$minage[spp] < data_list$fleet_control$Age_first_selected[i]){
          mapped_ages <- (data_list$minage[spp]:(data_list$fleet_control$Age_first_selected[i]-1))- data_list$minage[spp] + 1
          map_list$sel_coff[flt, , mapped_ages]  <- NA
        }
      }

      # 1-sex model
      if(nsex == 1){
        map_list$sel_coff[flt, 2, ]  <- NA
      }

      # Map out logistic and double logistic
      map_list$ln_sel_slp[1:2, flt ,] <- NA
      map_list$sel_inf[1:2, flt ,] <- NA

      # Map out logistic and double logistic deviates
      map_list$ln_sel_slp_dev[1:2, flt, ,] <- NA
      map_list$sel_inf_dev[1:2, flt, ,] <- NA

      # Map out selectivity var
      map_list$ln_sigma_sel[flt] <- NA
    }

    # -- 4.3. Double logistic - sel_type = 3
    if(data_list$fleet_control$Selectivity[i] == 3){ # Double logistic

      # Map out non-parametric
      map_list$sel_coff[flt,, ] <- NA
      map_list$sel_coff_dev[flt,,,] <- NA

      # Map out time varying parameters if not used
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(NA, 0)){
        map_list$ln_sel_slp_dev[1:2, flt, ,] <- NA
        map_list$sel_inf_dev[1:2, flt, ,] <- NA
      }

      # Penalized likelihood or random walk
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4,5)){
        for(j in 1:2){
          for(sex in 1:nsex){
            map_list$ln_sel_slp_dev[j, flt, sex,] <- ind_slp + 1:length(map_list$ln_sel_slp_dev[j, flt, sex,]) - 1
            map_list$sel_inf_dev[j, flt, sex,] <- ind_inf + 1:length(map_list$sel_inf_dev[j, flt, sex,]) - 1

            ind_slp <- ind_slp + length(map_list$ln_sel_slp_dev[j, flt, sex,])
            ind_inf <- ind_inf + length(map_list$sel_inf_dev[j, flt, sex,])
          }
        }

        # If only doing the ascending portion
        if(data_list$fleet_control$Time_varying_sel[i] %in% c(5)){
          map_list$ln_sel_slp_dev[2, flt, ,] <- NA
          map_list$sel_inf_dev[2, flt, ,] <- NA
        }
      }


      # Selectivity blocks
      if(data_list$fleet_control$Time_varying_sel[i] == 3){

        # If a fishery use the years from the fishery
        if(data_list$fleet_control$Fleet_type[i] == 1){
          fsh_biom <- data_list$fsh_biom[which(data_list$fsh_biom$Fleet_code == flt),]
          Selectivity_block <- fsh_biom$Selectivity_block
          biom_yrs <- fsh_biom$Year - data_list$styr + 1

          Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
          biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
        }

        # if a survey use the survey years
        if(data_list$fleet_control$Fleet_type[i] == 2){
          srv_biom <- data_list$srv_biom[which(data_list$srv_biom$Fleet_code == flt),]
          Selectivity_block <- srv_biom$Selectivity_block
          biom_yrs <- srv_biom$Year - data_list$styr + 1

          Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
          biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
        }

        # Map out logistic
        map_list$ln_sel_slp[1:2, flt ,] <- NA
        map_list$sel_inf[1:2, flt, ] <- NA

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


    # -- 4.4. Descending logitistic - sel_type = 4
    if (data_list$fleet_control$Selectivity[i] == 4) {

      # Map out non-parametric
      map_list$sel_coff[flt,, ] <- NA
      map_list$sel_coff_dev[flt,,,] <- NA

      # Map out ascending logistic
      map_list$ln_sel_slp[1, flt, ] <- NA
      map_list$sel_inf[1, flt, ] <- NA

      # Map out ascending logistic deviates
      map_list$ln_sel_slp_dev[1, flt, ,] <- NA
      map_list$sel_inf_dev[1, flt, ,] <- NA

      # Map out time varying parameters if not used
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(NA,0)){
        map_list$ln_sel_slp_dev[2, flt, ,] <- NA
        map_list$sel_inf_dev[2, flt, ,] <- NA
      }

      # Penalized/random deviate or random walk
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(1,2,4)){
        for(sex in 1:nsex){
          map_list$ln_sel_slp_dev[2, flt, sex, ] <- ind_slp + 1:length(map_list$ln_sel_slp_dev[2, flt, sex, ]) - 1
          map_list$sel_inf_dev[2, flt, sex, ] <- ind_inf + 1:length(map_list$sel_inf_dev[2, flt, sex, ]) - 1

          ind_slp <- ind_slp + length(map_list$ln_sel_slp_dev[2, flt, sex, ])
          ind_inf <- ind_inf + length(map_list$sel_inf_dev[2, flt, sex, ])
        }
      }

      # Selectivity block
      if(data_list$fleet_control$Time_varying_sel[i] == 3){

        # If a fishery use the years from the fishery
        if(data_list$fleet_control$Fleet_type[i] == 1){
          fsh_biom <- data_list$fsh_biom[which(data_list$fsh_biom$Fleet_code == flt),]
          Selectivity_block <- fsh_biom$Selectivity_block
          biom_yrs <- fsh_biom$Year - data_list$styr + 1

          Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
          biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
        }

        # if a survey use the survey years
        if(data_list$fleet_control$Fleet_type[i] == 2){
          srv_biom <- data_list$srv_biom[which(data_list$srv_biom$Fleet_code == flt),]
          Selectivity_block <- srv_biom$Selectivity_block
          biom_yrs <- srv_biom$Year - data_list$styr + 1

          Selectivity_block <- Selectivity_block[which(biom_yrs <= nyrs_hind)]
          biom_yrs <- biom_yrs[which(biom_yrs <= nyrs_hind)]
        }

        # Map out double logistic
        map_list$ln_sel_slp[1:2, flt ,] <- NA
        map_list$sel_inf[1:2, flt ,] <- NA

        # Map selectivity blocks
        for(sex in 1:nsex){
          map_list$ln_sel_slp_dev[2, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_slp
          map_list$sel_inf_dev[2, flt, sex, biom_yrs] <- Selectivity_block - 1 + ind_inf

          ind_slp <- ind_slp + max(Selectivity_block)
          ind_inf <- ind_inf + max(Selectivity_block)
        }
      }
    }

    # -- 4.5. Non-parametric similar to Hake (Taylor et al 2014) - sel_type = 5
    if(data_list$fleet_control$Selectivity[i] == 5){ # Non-parametric at age

      # -- Turn off time-varying coefs if not used
      if(data_list$fleet_control$Time_varying_sel[i] %in% c(0, NA, 3)){ # 0 = off, 1 = penalized like, 2 = random effects (NOT Implemented),...
        map_list$sel_coff_dev[flt,,,] <- NA
      }

      # -- If nselages is  < max(nselages) map out
      if(data_list$fleet_control$Nselages[i] < max(data_list$fleet_control$Nselages, na.rm = TRUE)){
        mapped_ages <- (data_list$fleet_control$Nselages[i] + 1):max(data_list$fleet_control$Nselages, na.rm = T)

        map_list$sel_coff[flt, , mapped_ages]  <- NA
        map_list$sel_coff_dev[flt, , mapped_ages,]  <- NA
      }

      # -- Map out ages <= Age first selected
      if(!is.na(data_list$fleet_control$Age_first_selected[i])){
        mapped_ages <- (data_list$minage[spp]:data_list$fleet_control$Age_first_selected[i])- data_list$minage[spp] + 1
        map_list$sel_coff[flt, , mapped_ages]  <- NA
        map_list$sel_coff_dev[flt, , mapped_ages,]  <- NA
      }

      # -- Map out years where comp data do not exist (set to average selectivity)
      map_list$sel_coff_dev[flt, , ,]  <- NA

      # Map out logistic and double logistic
      map_list$ln_sel_slp[1:2, flt ,] <- NA
      map_list$sel_inf[1:2, flt ,] <- NA

      # Map out logistic and double logistic deviates
      map_list$ln_sel_slp_dev[1:2, flt, ,] <- NA
      map_list$sel_inf_dev[1:2, flt, ,] <- NA
    }


    # Map out 2nd sex parameters if 1 sex model
    if(nsex == 1){

      # Map out logistic
      map_list$ln_sel_slp[1:2, flt, 2] <- NA
      map_list$sel_inf[1:2, flt, 2] <- NA

      map_list$ln_sel_slp_dev[1:2, flt, 2,] <- NA
      map_list$sel_inf_dev[1:2, flt, 2,] <- NA

      # -- Non-parametric
      map_list$sel_coff[flt, 2, ]  <- NA
      map_list$sel_coff_dev[flt,2,,] <- NA
    }
  }

  # -----------------------------------------------------------
  # Catchability ----
  # -----------------------------------------------------------
  # -- Catchability indices
  ind_q_re <- 1
  ind_q <- 1

  # Loop through fleets
  for( i in 1: nrow(data_list$fleet_control)){
    flt = data_list$fleet_control$Fleet_code[i]
    # Q
    # - 0 = fixed at prior
    # - 1 = Estimate single parameter
    # - 2 = Estimate single parameter with prior
    # - 3 = Estimate analytical q
    # - 4 = Estimate power equation
    # - 5 = Use env index ln(q_y) = q_mu + beta * index_y

    # Turn off catchability variance
    map_list$ln_sigma_srv_q[flt] <- NA

    # Catchability of surveys If not estimating turn of
    if (data_list$fleet_control$Estimate_q[i] %in% c(NA, 0, 3) | data_list$fleet_control$Fleet_type[i] == 0) {
      map_list$ln_srv_q[flt] <- NA
      # map_list$srv_q_pow[flt] <- NA
      map_list$ln_srv_q_dev[flt,] <- NA
      map_list$ln_sigma_srv_q[flt] <- NA
      map_list$ln_sigma_time_varying_srv_q[flt] <- NA
    }

    # Just turn off catchability power coeficient for single parameter.
    if (data_list$fleet_control$Estimate_q[i] %in% c(1,2, 3) | data_list$fleet_control$Fleet_type[i] == 0) {
      # map_list$srv_q_pow[flt] <- NA
    }

    # Time-varying catchability of surveys
    # If not estimating turn of. Also turns off if using environmental index Est_q = 5.
    if(data_list$fleet_control$Time_varying_q[i] %in% c(0, NA) | data_list$fleet_control$Estimate_q[i] %in% c(5)){
      map_list$ln_srv_q_dev[flt,] <- NA
      map_list$ln_sigma_srv_q[flt] <- NA
      map_list$ln_sigma_time_varying_srv_q[flt] <- NA
    }

    # Random walk or penalized likelihood - map out q_sd
    if(data_list$fleet_control$Time_varying_q[flt] %in% c(1, 4)){
      map_list$ln_sigma_srv_q[flt] <- NA
      map_list$ln_sigma_time_varying_srv_q[flt] <- NA
    }

    # Time block - map out q_sd
    if(data_list$fleet_control$Time_varying_q[i] %in% c(3)){
      map_list$ln_sigma_srv_q[flt] <- NA
      map_list$ln_sigma_time_varying_srv_q[flt] <- NA
      map_list$ln_srv_q[flt] <- NA
    }

    # Set up time varying catchability if used (account for missing years)
    if(data_list$fleet_control$Time_varying_q[i] %!in% c(NA, 0)){

      # Set all to NA
      map_list$ln_srv_q_dev[flt,] <- NA

      # Extract survey years where data is provided
      srv_biom <- data_list$srv_biom[which(data_list$srv_biom$Fleet_code == flt & data_list$srv_biom$Year > data_list$styr & data_list$srv_biom$Year <= data_list$endyr),]
      srv_biom_yrs <- srv_biom$Year - data_list$styr + 1
      srv_biom_yrs_miss <- yrs_hind[which(yrs_hind %!in% srv_biom_yrs)]

      # Penalized deviate or random walk
      if(data_list$fleet_control$Time_varying_q[i] %in% c(1,2,4)){
        map_list$ln_srv_q_dev[flt, srv_biom_yrs] <- ind_q + (1:length(srv_biom_yrs)) - 1
        ind_q <- ind_q + length(srv_biom_yrs)
      }

      # Time blocks
      if(data_list$fleet_control$Time_varying_q[i] == 3){
        map_list$ln_srv_q_dev[flt, srv_biom_yrs] <- ind_q + srv_biom$Selectivity_block - 1
        ind_q <- ind_q + max(srv_biom$Selectivity_block)
      }
    }

    # Standard deviation of surveys index If not estimating turn of
    if (data_list$fleet_control$Estimate_survey_sd[i] %in% c(NA, 0, 2) | data_list$fleet_control$Fleet_type[flt] == 0) {
      map_list$ln_sigma_srv_index[flt] <- NA
    }
  }


  # -- 7. Share survey q and selectivity
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
      map_list$ln_sigma_sel[flt] <- map_list$ln_sigma_sel[sel_duplicate]
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

      # Make catchability maps the same if selectivity is the same
      map_list$ln_srv_q[flt] <- map_list$ln_srv_q[q_duplicate]
      # map_list$srv_q_pow[flt] <- map_list$srv_q_pow[q_duplicate]
      map_list$ln_srv_q_dev[flt,] <- map_list$ln_srv_q_dev[q_duplicate,]
      map_list$ln_sigma_srv_q[flt] <- map_list$ln_sigma_srv_q[q_duplicate]
      map_list$ln_sigma_time_varying_srv_q[flt] <- map_list$ln_sigma_time_varying_srv_q[q_duplicate]
    }


    # Add index
    sel_index_tested <- c(sel_index_tested, sel_index[flt])
    q_index_tested <- c(q_index_tested, q_index[flt])
  }


  # -----------------------------------------------------------
  # - Fishery control ----
  # -----------------------------------------------------------
  for (i in 1:nrow(data_list$fleet_control)) {
    flt = data_list$fleet_control$Fleet_code[i]
    # Standard deviation of fishery time series If not estimating turn of
    if (data_list$fleet_control$Estimate_catch_sd[i] %in% c(NA, 0, 2)) {
      map_list$ln_sigma_fsh_catch[flt] <- NA
    }

    # Turn of F and F dev if not estimating of it is a Survey
    if (data_list$fleet_control$Fleet_type[i] %in% c(0, 2)) {
      map_list$ln_sigma_fsh_catch[flt] <- NA
      map_list$F_dev[flt, ] <- NA
      map_list$ln_mean_F[flt] <- NA
    }
  }


  # - Map out Fdev for years with 0 catch to very low number
  fsh_biom <- data_list$fsh_biom[which(data_list$fsh_biom$Year <= data_list$endyr),]
  fsh_ind <- fsh_biom$Fleet_code[which(fsh_biom$Catch == 0)]
  yr_ind <- fsh_biom$Year[which(fsh_biom$Catch == 0)] - data_list$styr + 1

  for(i in 1:length(yr_ind)){
    map_list$F_dev[fsh_ind[i], yr_ind[i]] <- NA
    map_list$ln_sel_slp_dev[1:2, fsh_ind[i], 1:2, yr_ind[i]] <- NA
    map_list$sel_inf_dev[1:2, fsh_ind[i], 1:2, yr_ind[i]] <- NA
  }
  #map_list$ln_sel_slp_dev_re[1:2, fsh_ind, 1:2, yr_ind] <- NA
  #map_list$sel_inf_dev_re[1:2, fsh_ind, 1:2, yr_ind] <- NA


  # -----------------------------------------------------------
  # Recruitment ----
  # -----------------------------------------------------------
  # -- Recruitment deviation sigmas - turn off if not estimating
  if(random_rec == FALSE){
    map_list$ln_rec_sigma <- map_list$ln_rec_sigma * NA
  }

  # -- Stock recruit relationship (SRR) parameters:
  # col1 = mean rec, col2 = SRR alpha, col3 = SRR beta
  # - Turning off 2nd and 3rd par if only using mean rec
  if(data_list$srr_fun == 0 & data_list$srr_pred_fun == 0){
    map_list$rec_pars[, 2:3] <- NA
  }

  # - Turning off mean rec par if using SRR
  if(data_list$srr_fun > 0){
    map_list$rec_pars[, 1] <- NA
  }

  # - Fix first parameter in SRR (if SRR not used, will be NA anyway)
  if(data_list$srr_est_mode == 0){
    map_list$rec_pars[, 2] <- NA
  }


  # -----------------------------------------------------------
  # Predation bits ----
  # (e.g. Turn off all predation parameters for single species)
  # -----------------------------------------------------------
  if (data_list$msmMode == 0) { # Single-species

    # Suitability parameters
    map_list$log_gam_a <- map_list$log_gam_a * NA
    map_list$log_gam_b <- map_list$log_gam_b * NA
    map_list$log_phi <- map_list$log_phi * NA

    # # Multispecies kinzey parameters
    map_list$logH_1 <- map_list$logH_1 * NA
    map_list$logH_1a <- map_list$logH_1a * NA
    map_list$logH_1b <- map_list$logH_1b * NA

    map_list$logH_2 <- map_list$logH_2 * NA
    map_list$logH_3 <- map_list$logH_3 * NA
    map_list$H_4 <- map_list$H_4 * NA

  }

  # 2. MSVPA based predation
  if (data_list$msmMode %in% c(1,2)) {

    # Multispecies kinzey parameters
    map_list$logH_1 <- map_list$logH_1 * NA
    map_list$logH_1a <- map_list$logH_1a * NA
    map_list$logH_1b <- map_list$logH_1b * NA

    map_list$logH_2 <- map_list$logH_2 * NA
    map_list$logH_3 <- map_list$logH_3 * NA
    map_list$H_4 <- map_list$H_4 * NA

  }

  # 3. Kinzey and Punt predation equations
  if (data_list$msmMode > 2) {
    # Holling Type 1
    if (data_list$msmMode == 3) {
      map_list$logH_2 <- replace(map_list$logH_2, values = rep(NA, length(map_list$logH_2)))
      map_list$logH_3 <- replace(map_list$logH_3, values = rep(NA, length(map_list$logH_3)))
      map_list$H_4 <- replace(map_list$H_4, values = rep(NA, length(map_list$H_4)))
    }

    # Holling Type 2
    if (data_list$msmMode == 4) {
      map_list$logH_3 <- map_list$logH_3 * NA
      map_list$H_4 <- map_list$H_4 * NA
    }

    # Holling Type 3
    if (data_list$msmMode == 5) {
      map_list$logH_3 <- map_list$logH_3 * NA
    }

    # Predator interference
    if (data_list$msmMode == 6) {
      map_list$H_4 <- map_list$H_4 * NA
    }

    # Predator preemption
    if (data_list$msmMode == 7) {
      map_list$H_4 <- map_list$H_4 * NA
    }

    # Hassell-Varley
    if (data_list$msmMode == 8) {
      map_list$logH_3 <- map_list$logH_3 * NA
    }

    # Ecosim
    if (data_list$msmMode == 9) {
      map_list$logH_2 <- map_list$logH_2 * NA
      map_list$H_4 <- map_list$H_4 * NA
    }
  }


  # -----------------------------------------------------------
  # Suitability bits ----
  # -----------------------------------------------------------
  if (data_list$msmMode > 0) {

    # 2.1. Empirical suitability
    if (data_list$suitMode == 0) {
      # Suitability parameters
      map_list$log_gam_a <- map_list$log_gam_a * NA
      map_list$log_gam_b <- map_list$log_gam_b * NA
      map_list$log_phi <- map_list$log_phi * NA
    }

    # 2.2. GAMMA or lognormal suitability
    if (data_list$suitMode %in% c(1:4)) {
      # Use all the parameters
    }
  }

  # STEP 3 - set up debug - I.E. turn off all parameters besides dummy
  map_list$dummy <- NA
  if(debug){
    map_list <- sapply(map_list, function(x) replace(x, values = rep(NA, length(x))))
    map_list$dummy = 1
  }

  # -----------------------------------------------------------
  # Set up fixed n-at-age ----
  # -----------------------------------------------------------
  # - I.E. turn off all parameters besides for species
  for(sp in 1:data_list$nspp){

    # Fixed n-at-age: Turn off most parameters
    if(data_list$estDynamics[sp] > 0){

      # Population parameters
      map_list$rec_pars[sp,] <- NA
      map_list$ln_rec_sigma[sp] <- NA
      map_list$ln_sex_ratio_sigma[sp] <- NA
      map_list$rec_dev[sp,] <- NA
      map_list$init_dev[sp,] <- NA
      map_list$ln_M1[sp,,] <- NA

      # Survey and fishery fleet parameters
      flts <- data_list$fleet_control$Fleet_code[which(data_list$fleet_control$Species == sp)]


      map_list$ln_mean_F[flts] <- NA
      map_list$F_dev[flts,] <- NA
      map_list$ln_srv_q[flts] <- NA
      # map_list$srv_q_pow[flts] <- NA
      map_list$ln_srv_q_dev[flts,] <- NA
      map_list$ln_sigma_srv_q[flts] <- NA
      map_list$ln_sigma_time_varying_srv_q[flts] <- NA
      map_list$sel_coff[flts,,] <- NA
      map_list$sel_coff_dev[flts,,,] <- NA
      map_list$ln_sel_slp[, flts, ] <- NA
      map_list$sel_inf[, flts, ] <- NA
      map_list$ln_sel_slp_dev[, flts, ,] <- NA
      map_list$sel_inf_dev[, flts, ,] <- NA
      map_list$ln_sigma_sel[flts] <- NA
      map_list$ln_sigma_srv_index[flts] <- NA
      map_list$ln_sigma_fsh_catch[flts] <- NA
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


  # -----------------------------------------------------------
  # STEP 5 -- Convert to factor ----
  # -----------------------------------------------------------
  map_list_grande <- list()
  map_list_grande$mapFactor <- sapply(map_list, factor)
  map_list_grande$mapList <- map_list

  return(map_list_grande)
}
