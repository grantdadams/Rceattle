# -----------------------------------------------------------
# Catchability ----
# -----------------------------------------------------------
# -- Catchability indices
ind_q_re <- 1
ind_q <- 1
ind_beta_q <- 0

map_list$ln_srv_q <- replace(map_list$ln_srv_q, values = rep(NA, length(x)))
# map_list$srv_q_pow <- replace(map_list$srv_q_beta, values = rep(NA, length(x)))
map_list$srv_q_beta <- replace(map_list$srv_q_beta, values = rep(NA, length(x)))
map_list$ln_srv_q_dev <- replace(map_list$ln_srv_q_dev, values = rep(NA, length(x)))
map_list$ln_sigma_srv_q <- replace(map_list$ln_sigma_srv_q, values = rep(NA, length(x)))
map_list$ln_sigma_time_varying_srv_q <- replace(map_list$ln_sigma_time_varying_srv_q, values = rep(NA, length(x)))
map_list$ln_sigma_srv_index <- replace(map_list$ln_sigma_srv_index, values = rep(NA, length(x)))

# Loop through fleets
for( i in 1: nrow(data_list$fleet_control)){
  if(data_list$fleet_control$Fleet_type[flt] == 1){ # If survey

    flt = data_list$fleet_control$Fleet_code[i]
    # Q
    # - 0 = fixed at prior
    # - 1 = Estimate single parameter
    # - 2 = Estimate single parameter with prior
    # - 3 = Estimate analytical q
    # - 4 = Estimate power equation
    # - 5 = Use env index ln(q_y) = q_mu + beta * index_y


    # - Turn on q for:
    # - 1 = Estimate single parameter
    # - 2 = Estimate single parameter with prior
    # - 4 = Estimate power equation
    # - 5 = Use env index ln(q_y) = q_mu + beta * index_y
    if(data_list$fleet_control$Estimate_q[i] %in% c(1, 2, 4, 5)){
      map_list$ln_srv_q[flt] <- flt
    }

    # - Turn on power param for:
    # - 4 = Estimate power equation
    if (data_list$fleet_control$Estimate_q[i] %in% c(4)) {
      # map_list$srv_q_pow[flt] <- flt
    }

    # Time- varying q parameters
    # - 0 = no,
    # - 1 = penalized deviate
    # - 2 = random effect
    # - 3 = time blocks with no penalty
    # - 4 = random walk from mean following Dorn 2018 (dnorm(q_y - q_y-1, 0, sigma)
    # - If estimate_q == 5; this determines the environmental indices to be used in the equation log(q_y) = q_mu + beta * index_y

    # -- Set up time varying catchability if used (account for missing years)
    if(data_list$fleet_control$Estimate_q[i] %in% c(1, 2) &
       as.numeric(data_list$fleet_control$Time_varying_q[i]) %!in% c(1, 2, 3, 4)){

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

    # - Turn on regression coefficients for:
    # - 5 = Estimate environmental linkage
    if (data_list$fleet_control$Estimate_q[i] == 5) {
      turn_on <- as.numeric(unlist(strsplit(data_list$fleet_control$Time_varying_q[i],","))) # Parameters to turn on
      map_list$srv_q_beta[flt, turn_on] <- turn_on + ind_beta_q
      ind_beta_q <- ind_beta_q + max(turn_on)
    }

    # Standard deviation of surveys index
    # - 0 = use CV from srv_biom
    # - 1 = estimate a free parameter
    # - 2 = analytically estimate following (Ludwig and Walters 1994)
    if (data_list$fleet_control$Estimate_survey_sd[i] == 1) {
      map_list$ln_sigma_srv_index[flt] <- flt
    }
  }
}
