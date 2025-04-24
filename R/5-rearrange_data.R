#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB CEATTLE
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @export
rearrange_dat <- function(data_list){

  # Data dimensions
  max_sex <- max(data_list$nsex, na.rm = T)
  max_age <- max(data_list$nages, na.rm = T)
  max_length <- max(data_list$nlengths, na.rm = T)
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)

  # 1 - remove non-integer objects from control ----
  data_list$index_ln_q_prior <- log(data_list$fleet_control$Q_prior)

  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::select(Fleet_name,
                  Fleet_code,           # 1) Temporary survey index
                  Fleet_type,           # 2) Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
                  Species,              # 3) Species
                  Selectivity_index,    # 4) Survey selectivity index
                  Selectivity,          # 5) Selectivity type
                  Nselages,             # 6) Non-parametric selectivity ages
                  Time_varying_sel,     # 7) Time-varying selectivity type.
                  Age_first_selected,   # 8) First age selected
                  Age_max_selected,     # 9) Age of max selectivity (used for normalization). If NA, does not normalize
                  Comp_loglike,         # 10) Index indicating wether to do dirichlet multinomial for a multinomial)
                  Weight1_Numbers2,     # 11) Survey units
                  Weight_index,         # 12) Dim1 of weight (what weight-at-age data set)
                  Age_transition_index, # 13) Dim3 of age transition matrix (what ALK to use)
                  Q_index,              # 14) Index of survey q
                  Estimate_q,           # 15) Parametric form of q
                  Time_varying_q,       # 16) Time varying q type
                  Estimate_index_sd,    # 17) Wether to estimate standard deviation of survey time series
                  Estimate_catch_sd     # 18) Wether to estimate standard deviation of fishery time series
    )
  # Don't want: "Sel_sd_prior", "Q_prior", "Q_sd_prior", "Time_varying_q_sd_prior", "Survey_sd_prior", "proj_F", "Catch_sd_prior", "Comp_weights", "proj_F_prop"

  data_list$fleet_control$Time_varying_sel <- round(data_list$fleet_control$Time_varying_sel)
  data_list$fleet_control$Fleet_name <- suppressWarnings(as.numeric(as.character(data_list$fleet_control$Fleet_name)))
  data_list$fleet_control$Time_varying_q <- suppressWarnings(as.numeric(as.character(data_list$fleet_control$Time_varying_q)))

  # Species names
  data_list$spnames <- NULL

  # * F reference points ----
  data_list$Ftarget_percent <- data_list$Ftarget
  data_list$Flimit_percent <- data_list$Flimit

  # - Input missing nselages and age first selected (use age range)
  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(Nselages = ifelse(is.na(Nselages), -999, Nselages),
                  Age_max_selected = ifelse(is.na(Age_max_selected), -999, Age_max_selected),
                  Age_first_selected = ifelse(is.na(Age_first_selected), data_list$minage[Species], Age_first_selected)
    )


  # 2 -  Seperate survey biomass info from observation ----
  data_list$index_ctl <- data_list$index_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$index_n <- as.matrix(data_list$index_data[,c("Month")])

  data_list$index_obs <- data_list$index_data %>%
    dplyr::select(Observation, Log_sd) %>%
    dplyr::mutate_all(as.numeric)


  # 3 -  Seperate catch biomass info from observation ----
  data_list$catch_ctl <- data_list$catch_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$catch_n <- as.matrix(data_list$catch_data[,c("Month")])

  data_list$catch_obs <- data_list$catch_data %>%
    dplyr::select(Catch, Log_sd) %>%
    dplyr::mutate_all(as.numeric)


  # 4 -  Seperate comp data from observation ----
  data_list$comp_ctl <- data_list$comp_data %>%
    dplyr::select(Fleet_code, Species, Sex, Age0_Length1, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$comp_n <- data_list$comp_data %>%
    dplyr::select(Month, Sample_size) %>%
    dplyr::mutate_all(as.numeric)

  data_list$comp_obs <- data_list$comp_data %>%
    dplyr::select(contains("Comp_")) %>%
    dplyr::mutate_all(as.numeric)

  data_list <- check_composition_data(data_list)


  # 5 -  Seperate diet info from observation ----
  data_list$diet_ctl <- data_list$diet_data %>%
    dplyr::select(Pred, Prey, Pred_sex, Prey_sex, Pred_age, Prey_age, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$diet_obs <- data_list$diet_data %>%
    dplyr::select(Sample_size, Stomach_proportion_by_weight) %>%
    dplyr::mutate_all(as.numeric)


  # 6 -  Seperate survey empirical selectivity info from observation ----
  yrs <- data_list$styr:data_list$endyr
  if(nrow(data_list$emp_sel) > 0 ){
    for(i in 1:nrow(data_list$emp_sel)){
      # Fill in years
      if(data_list$emp_sel$Year[i] == 0){

        # Change first year
        data_list$emp_sel$Year[i] <- yrs[1]

        # Change the rest
        emp_sel_tmp <- data_list$emp_sel[i,]
        emp_sel_tmp$Year <- yrs[2]
        emp_sel <- emp_sel_tmp

        for(yr in 3:length(yrs)){
          emp_sel_tmp$Year <- yrs[yr]
          emp_sel <- rbind(emp_sel, emp_sel_tmp)
        }

        data_list$emp_sel <- rbind(data_list$emp_sel, emp_sel)
      }
    }
  }

  data_list$emp_sel_ctl <- data_list$emp_sel %>%
    dplyr::select(Fleet_code, Species, Sex, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$emp_sel_obs <- data_list$emp_sel %>%
    dplyr::select(contains("Comp_")) %>%
    dplyr::mutate_all(as.numeric)


  # 7 - Rearrange age-transition matrix ----
  age_trans_matrix <- data_list$age_trans_matrix
  unique_age_transition <- unique(as.character(age_trans_matrix$Age_transition_index))
  if(sum(!data_list$pop_age_transition_index %in% unique_age_transition) > 0){
    stop("Check population age_transition index, not in age_transition file")
  }
  age_transition <- array(0, dim = c(length(unique_age_transition), max_sex, max_age, max_length))


  for (i in 1:nrow(age_trans_matrix)) {
    age_transition_ind <- as.numeric(as.character(age_trans_matrix$Age_transition_index[i]))
    sex <- as.numeric(as.character(age_trans_matrix$Sex[i]))
    sp <- as.numeric(as.character(age_trans_matrix$Species[i]))
    age <- as.numeric(as.character(age_trans_matrix$Age[i])) - data_list$minage[sp] + 1

    if (age > data_list$nages[sp]) {
      message(paste0("Error: number of ages in age_trans_matrix for species: ", sp, " and index: ", age_transition_ind))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    # Handle sex == 0 case for 2-sex species
    sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

    for(j in 1:length(sex_values)){
      age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]] <- as.numeric(as.character(age_trans_matrix[i, (1:data_list$nlengths[sp]) + 5]))

      # Normalize
      age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]] <- age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]] / sum(age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]], na.rm = TRUE)
    }
  }
  data_list$age_trans_matrix <- age_transition


  # 8 - Rearrange age_error matrices ----
  arm <- array(0, dim = c(data_list$nspp, max_age, max_age))

  data_list$age_error <- as.data.frame(data_list$age_error) # FIXME: somethin is up with data.frames
  for (i in 1:nrow(data_list$age_error)) {
    sp <- as.numeric(as.character(data_list$age_error$Species[i]))
    true_age <- as.numeric(as.character(data_list$age_error$True_age[i])) - data_list$minage[sp] + 1

    if (true_age > data_list$nages[sp]) {
      message(paste0("Error: number of true ages specified in age_error for species: ", sp))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    arm[sp, true_age, 1:data_list$nages[sp]] <- as.numeric(as.character(data_list$age_error[i, (1:data_list$nages[sp]) + 2]))

    # Normalize
    arm[sp, true_age, 1:data_list$nages[sp]] <- arm[sp, true_age, 1:data_list$nages[sp]] / sum(arm[sp, true_age, 1:data_list$nages[sp]], na.rm = TRUE)
  }
  data_list$age_error <- arm


  # 10 - Normalize comp data ----
  data_list$comp_obs <- t(apply(data_list$comp_obs, 1, function(x) as.numeric(x) / sum(as.numeric(x), na.rm = TRUE)))


  # 11 - Set up weight-at-age ----
  # - Convert to array
  data_list$weight <- data_list$weight %>%
    dplyr::mutate(
      Wt_index = as.numeric(as.character(Wt_index)),
      Species = as.numeric(as.character(Species)),
      Sex = as.numeric(as.character(Sex)),
      Year = as.numeric(as.character(Year)),
      Year = ifelse(Year == 0,
                    0,
                    Year - data_list$styr + 1)
    )

  # Pre-allocate the array
  unique_wt <- unique(as.numeric(data_list$weight$Wt_index))
  weight <- array(0, dim = c(length(unique_wt), max_length, max_age, nyrs_hind))

  # Convert weight data to numeric once, outside the loop
  weight_matrix <- data_list$weight[, (1:max_age) + 5]
  weight_matrix <- apply(weight_matrix, 2, function(x) as.numeric(as.character(x)))
  weight_matrix <- matrix(weight_matrix, ncol = max_age)

  # Check for spaces in the weight data
  if (any(grepl("[[:space:]]", weight_matrix))) {
    stop("Space found in weight data")
  }

  # Main processing
  for (i in seq_len(nrow(data_list$weight))) {
    wt_ind <- data_list$weight$Wt_index[i]
    sp <- data_list$weight$Species[i]
    sex <- data_list$weight$Sex[i]
    yr <- data_list$weight$Year[i]

    # Check if the year is within the valid range
    if (yr <= data_list$endyr - data_list$styr + 1) {
      # Handle year == 0 case
      if(yr == 0) {
        yr <- seq_len(nyrs_hind)
      }

      # Handle sex == 0 case for 2-sex species
      sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

      # Assign weights to the array
      weight[wt_ind, sex_values, 1:data_list$nages[sp], yr] <- rep(weight_matrix[i, 1:data_list$nages[sp]], each = length(sex_values))
    }
  }

  # Update the data_list with the new weight array
  data_list$weight <- weight


  # 12 - Set up NByageFixed ----
  # - Convert to array
  data_list$NByageFixed <- data_list$NByageFixed %>%
    mutate(Species = as.numeric(as.character(Species)),
           Sex = as.numeric(as.character(Sex)),
           Year = as.numeric(as.character(Year)) - data_list$styr + 1)

  NByageFixed <- array(0, dim = c(data_list$nspp, max_sex, max_age, nyrs_proj))

  if(nrow(data_list$NByageFixed) > 0){
    for (i in 1:nrow(data_list$NByageFixed)) {

      sp <- data_list$NByageFixed$Species[i]
      sex <- data_list$NByageFixed$Sex[i]
      yr <- data_list$NByageFixed$Year[i]

      # Handle sex == 0 case for 2-sex species
      sex_values <- if(sex == 0) 1:data_list$nsex[sp] else sex

      for(j in 1:length(sex_values)){
        if(yr > 0){
          NByageFixed[sp, sex_values[j], 1:max_age, yr] <- as.numeric(data_list$NByageFixed[i,c(1:max_age)+4])
        }
      }
    }
  }

  data_list$NByageFixed <- NByageFixed


  # 13 - Set up environmental indices ----
  # - Fill in missing years with column mean
  data_list$env_index <- merge(data_list$env_data, data.frame(Year = data_list$styr:data_list$projyr), all = TRUE)
  data_list$env_index <-  data_list$env_index %>%
    dplyr::select(-Year) %>%
    dplyr::mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%
    as.matrix()

  # - Create matrix for srr curve
  if(is.null(data_list$srr_indices)){
    data_list$srr_indices <- 1:(ncol(data_list$env_index)-1)
  }
  if(any(is.na(data_list$srr_indices))){
    data_list$srr_indices <- 1:(ncol(data_list$env_index)-1)
  }
  if(sum(sapply(data_list$srr_indices, function(x) x>ncol(data_list$env_index))) > 0){stop("'srr_indices' greater than the number of indices included")}
  data_list$env_index_srr <-  data_list$env_index[, data_list$srr_indices] %>%
    as.matrix()


  # - Create matrix for M1
  if(is.null(data_list$M1_indices)){
    data_list$M1_indices <- 1:(ncol(data_list$env_index)-1)
  }
  if(any(is.na(data_list$M1_indices))){
    data_list$M1_indices <- 1:(ncol(data_list$env_index)-1)
  }
  if(sum(sapply(data_list$M1_indices, function(x) x>ncol(data_list$env_index))) > 0){stop("'M1_indices' greater than the number of indices included")}
  data_list$env_index_M1 <-  data_list$env_index[, data_list$M1_indices] %>%
    as.matrix()


  # 14 - Set up pyrs ----
  # - Convert to array
  data_list$Pyrs <- data_list$Pyrs %>%
    dplyr::mutate(Species = as.numeric(as.character(Species)),
                  Sex = as.numeric(as.character(Sex)),
                  Year = as.numeric(as.character(Year)),
                  Year = ifelse(Year == 0,
                                0,
                                Year - data_list$styr + 1)
    )

  Pyrs <- array(0, dim = c(data_list$nspp, max_sex, max_age, nyrs_hind)) # FIXME: Change for forecast

  if( nrow(data_list$Pyrs) > 0){
    for (i in 1:nrow(data_list$Pyrs)) {
      sp <- data_list$Pyrs$Species[i]
      sex <- data_list$Pyrs$Sex[i]
      yr <- data_list$Pyrs$Year[i]

      # Handle sex == 0 case for 2-sex species
      sex_values <- if(sex == 0) 1:data_list$nsex[sp] else sex

      if(yr <= (data_list$endyr - data_list$styr + 1)){

        # Handle year == 0 case
        if(yr == 0) {
          yr <- seq_len(nyrs_hind)
        }

        # Fill in years
        for(j in 1:length(sex_values)){
          Pyrs[sp, sex_values[j], 1:data_list$nages[sp], yr] <- as.numeric(as.character(data_list$Pyrs[i, (1:data_list$nages[sp]) + 3]))
        }
      }
    }
  }
  data_list$Pyrs <- Pyrs


  # 15 - Remove species column from alw, maturity, sex_ratio ----
  data_list$sex_ratio <- data_list$sex_ratio[,-1]
  data_list$maturity <- data_list$maturity[,-1]
  # data_list$aLW <- data_list$aLW[,-1]


  # 16 - Make data.frames into matrices ----
  df_to_mat <- which(sapply(data_list, function(x) class(x)[1]) == "data.frame")
  data_list[df_to_mat] <- lapply(data_list[df_to_mat], as.matrix)

  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "catch_data",    "index_data", "comp_data", "env_data", "spnames",
                       "aLW", "diet_data", # "NByageFixed", "estDynamics", "Ceq",
                       "avgnMode", "minNByage")
  data_list[items_to_remove] <- NULL

  return(data_list)
}



#' Check and Clean Composition Data
#'
#' This function checks the composition data for zero sum rows, handles NA values,
#' and verifies that the composition data spans the required range of ages/lengths.
#'
#' @param data_list A list containing the following components:
#'   - comp_obs: A matrix or data frame of composition observations.
#'   - comp_data: A data frame with metadata including species, sex, and age/length information.
#'   - nsex: A vector of sex counts by species.
#'   - nages: A vector of age counts by species.
#'   - nlengths: A vector of length counts by species.
#' throws An error if any rows of `comp_obs` sum to 0 or if the composition data does not span the range of ages/lengths.
#'
#' @return The modified `data_list` with NA values in `comp_obs` converted to 0.
#' @examples
#' # Example usage:
#' data_list <- list(
#'   comp_obs = matrix(c(1, 2, 3, 0, 4, 5), nrow = 2),
#'   comp_data = data.frame(Species = c(1, 2), Sex = c(1, 3), Age0_Length1 = c(0, 1)),
#'   nsex = c(1, 2),
#'   nages = c(5, 6),
#'   nlengths = c(10, 12)
#' )
#' cleaned_data_list <- check_composition_data(data_list)
#'
#' @export
check_composition_data <- function(data_list) {
  # Check for zero sum rows in composition data
  if (any(rowSums(data_list$comp_obs, na.rm = TRUE) == 0)) {
    stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")
  }

  # Calculate adjustments for sex and age/length
  joint_adjust <- ifelse(data_list$nsex[data_list$comp_data$Species] == 2 &
                           data_list$comp_data$Sex == 3, 2, 1)

  col_adjust <- ifelse(data_list$comp_data$Age0_Length1 == 0,
                       data_list$nages[data_list$comp_data$Species],
                       data_list$nlengths[data_list$comp_data$Species])

  col_adjust <- col_adjust * joint_adjust

  # Check for NAs in composition data and convert them to 0
  na_check <- apply(cbind(col_adjust, data_list$comp_obs), 1, function(x) sum(x[2:(x[1] + 1)]))

  if (any(is.na(na_check))) {
    na_rows <- which(is.na(na_check))
    warning(sprintf("Composition data have NAs in row(s): %s. Converting to 0s.",
                    paste(na_rows, collapse = ", ")))
    data_list$comp_obs[is.na(data_list$comp_obs)] <- 0
  }

  # Verify composition data size
  if (ncol(data_list$comp_obs) < max(col_adjust)) {
    stop("Comp data does not span range of ages/lengths")
  }

  return(data_list)
}

