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

  # 1 - Fleet control ----
  # - 0) Vector to save  species
  data_list$flt_spp <- data_list$fleet_control %>%
    dplyr::pull(Species) %>% as.integer() - 1

  # - 1) Fleet pointer
  data_list$flt_sel_ind <- data_list$fleet_control %>%
    dplyr::pull(Fleet_code) %>% as.integer() - 1

  # - 2) Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
  data_list$flt_type <- data_list$fleet_control %>%
    dplyr::pull(Fleet_type) %>% as.integer()

  # - 3) Month of observation
  data_list$flt_month <- data_list$fleet_control %>%
    dplyr::pull(Month)

  # - 4) Selectivity type
  data_list$flt_sel_type <- data_list$fleet_control %>%
    dplyr::pull(Selectivity) %>% as.integer()

  # - 5) Number of ages/lengths for non-parametric selectivity
  data_list$flt_n_sel_bins <- data_list$fleet_control %>%
    dplyr::mutate(N_sel_bins = ifelse(is.na(N_sel_bins), -999, N_sel_bins)) %>%
    dplyr::pull(N_sel_bins) %>% as.integer()

  # - 6) Time-varying selectivity type.
  data_list$flt_varying_sel <- data_list$fleet_control %>%
    dplyr::pull(Time_varying_sel) %>% as.integer()

  # - 7) First age selected
  data_list$age_first_selected <- data_list$fleet_control %>%
    dplyr::mutate(Age_first_selected = Age_first_selected - data_list$minage[Species],
                  Age_first_selected = ifelse(is.na(Age_first_selected), 0, Age_first_selected)) %>%
    dplyr::pull(Age_first_selected) %>% as.integer()

  # - 8) Age of max selectivity (used for normalization). If NA, does not normalize
  data_list$sel_norm_bin1 <- data_list$fleet_control %>%
    dplyr::mutate(
      Sel_norm_bin1 = Sel_norm_bin1 - data_list$minage[Species],
      Sel_norm_bin1 = ifelse(Sel_norm_bin1 < 0, -99, Sel_norm_bin1),         # Less than zero, normalize by max
      Sel_norm_bin1 = ifelse(is.na(Sel_norm_bin1), -999, Sel_norm_bin1)) %>% # NA, do not normalize (unless type = 2)
    dplyr::pull(Sel_norm_bin1) %>% as.integer()

  # - 9) upper age of max selectivity (used for normalization). If NA, does not normalize
  data_list$sel_norm_bin2 <- data_list$fleet_control %>%
    dplyr::mutate(Sel_norm_bin2 = Sel_norm_bin2 - data_list$minage[Species],
                  Sel_norm_bin2 = ifelse(is.na(Sel_norm_bin2), -999, Sel_norm_bin2)) %>%
    dplyr::pull(Sel_norm_bin2) %>% as.integer()

  # - 10) Index indicating wether to do dirichlet multinomial or a multinomial
  data_list$comp_ll_type <- data_list$fleet_control %>%
    dplyr::pull(Comp_loglike) %>% as.integer()
  data_list$caal_ll_type <- data_list$fleet_control %>%
    dplyr::pull(CAAL_loglike) %>% as.integer()

  # - 11) Index units (1 = weight, 2 = numbers)
  data_list$flt_units <- data_list$fleet_control %>%
    dplyr::pull(Weight1_Numbers2) %>% as.integer()

  # - 12) Dim1 of weight (what weight-at-age data set)
  data_list$flt_wt_index <- data_list$fleet_control %>%
    dplyr::pull(Weight_index) %>% as.integer() - 1

  # - 13) Dim3 of age transition matrix (what ALK to use)
  data_list$flt_age_transition_index <- data_list$fleet_control %>%
    dplyr::pull(Age_transition_index) %>% as.integer() - 1

  # - 14) Parametric form of q
  data_list$est_index_q <- data_list$fleet_control %>%
    dplyr::pull(Estimate_q) %>% as.integer()

  # - 15) Time varying q type
  data_list$index_varying_q <- data_list$fleet_control %>%
    dplyr::pull(Time_varying_q) %>% as.integer()

  # - 16) Wether to estimate standard deviation of index time series
  data_list$est_sigma_index <- data_list$fleet_control %>%
    dplyr::pull(Estimate_index_sd) %>% as.integer()

  # - 17) Wether to estimate standard deviation of fishery time series
  data_list$est_sigma_fsh <- data_list$fleet_control %>%
    dplyr::pull(Estimate_catch_sd) %>% as.integer()

  data_list$index_ln_q_prior <- log(data_list$fleet_control$Q_prior)

  # Species names
  data_list$spnames <- NULL

  # * F reference points ----
  data_list$Ftarget_percent <- data_list$Ftarget
  data_list$Flimit_percent <- data_list$Flimit


  # 2 -  Index data ----
  # - Seperate index metadata from observation
  data_list$index_ctl <- data_list$index_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$index_n <- as.matrix(data_list$index_data[,c("Month")])

  data_list$index_obs <- data_list$index_data %>%
    dplyr::select(Observation, Log_sd) %>%
    dplyr::mutate_all(as.numeric)


  # 3 -  Catch data ----
  # - Seperate catch metadata from observation
  data_list$catch_ctl <- data_list$catch_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$catch_n <- as.matrix(data_list$catch_data[,c("Month")])

  data_list$catch_obs <- data_list$catch_data %>%
    dplyr::select(Catch, Log_sd) %>%
    dplyr::mutate_all(as.numeric)


  # 4 -  Comp data ----
  # - Seperate comp metadata from observation
  data_list$comp_ctl <- data_list$comp_data %>%
    dplyr::select(Fleet_code, Species, Sex, Age0_Length1, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$comp_n <- data_list$comp_data %>%
    dplyr::select(Month, Sample_size) %>%
    dplyr::mutate_all(as.numeric)

  data_list$comp_obs <- data_list$comp_data %>%
    dplyr::select(contains("Comp_")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()
  data_list$comp_obs <- t(apply(data_list$comp_obs, 1, function(x) as.numeric(x) / sum(as.numeric(x), na.rm = TRUE))) # Normalize

  data_list <- check_composition_data(data_list)

  # 5 - CAAL data ----
  data_list$caal_ctl <- data_list$caal_data %>%
    dplyr::select(Fleet_code, Species, Sex, Year, Length) %>%
    dplyr::mutate_all(as.integer)

  data_list$caal_n <- data_list$caal_data %>%
    dplyr::select(Sample_size) %>%
    dplyr::mutate_all(as.numeric)

  data_list$caal_obs <- data_list$caal_data %>%
    dplyr::select(contains("CAAL_")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()
  data_list$caal_obs <- (apply(data_list$caal_obs, 1, function(x) as.numeric(x) / sum(as.numeric(x), na.rm = TRUE))) # Normalize

  data_list <- check_caal_data(data_list)

  # * Get lengths from CAAL ----
  data_list$lengths <- matrix(rep(1:max(data_list$nlengths), data_list$nspp), nrow = data_list$nspp, byrow = TRUE)

  if(nrow(data_list$caal_data) > 0){
    caal_lengths <- data_list$caal_data %>%
      dplyr::distinct(Species, Length) %>%
      dplyr::arrange(Species, Length) %>%
      dplyr::group_by(Species) %>%
      dplyr::mutate(Bin = paste0("Bin", 1:n())) %>%
      tidyr::pivot_wider(names_from = Bin, values_from = Length)
    data_list$lengths[caal_lengths$species, 1:ncol(caal_lengths[,-1])] <- as.matrix(caal_lengths[,-1])
  }


  # 6 -  Diet data ----
  # - Seperate diet metadata from observation
  if(!is.null(data_list$diet_data)){ # Add a check in case there's no diet data

    # Create a temporary diet data frame to work with
    diet_dat <- data_list$diet_data %>%
      # - Sort by stomach_id so that rows for the same predator are contiguous
      dplyr::arrange(Pred, Pred_sex, Pred_age, Prey, Prey_sex, Prey_age, Year) %>%
      # - Create the unique stratum identifier and the zero-indexed stomach_id
      dplyr::mutate(stratum_id = paste(Pred, Pred_sex, Pred_age, Year, sep = "_"),
                    stomach_id = as.numeric(as.factor(stratum_id)) - 1)

    # Add n_stomach_obs and the stomach_id vector to the main data_list
    data_list$n_stomach_obs <- max(diet_dat$stomach_id) + 1
    data_list$stomach_id <- diet_dat$stomach_id

    # Create diet_ctl (without the stomach_id column)
    data_list$diet_ctl <- diet_dat %>%
      dplyr::select(Pred, Prey, Pred_sex, Prey_sex, Pred_age, Prey_age, Year) %>%
      dplyr::mutate_all(as.integer)

    # Create diet_obs as before
    data_list$diet_obs <- diet_dat %>%
      dplyr::select(Sample_size, Stomach_proportion_by_weight) %>%
      dplyr::mutate_all(as.numeric)
  }


  # 7 -  Seperate empirical selectivity info from observation ----
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


  # 8 - Rearrange age-transition matrix ----
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


  # 9 - Rearrange age_error matrices ----
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
  weight <- array(0, dim = c(max(unique_wt), max_sex, max_age, nyrs_hind))

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
  data_list$weight_obs <- weight


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


  # 14 - Set up ration_data ----
  # - Convert to array
  data_list$ration_data <- data_list$ration_data %>%
    dplyr::mutate(Species = as.numeric(as.character(Species)),
                  Sex = as.numeric(as.character(Sex)),
                  Year = as.numeric(as.character(Year)),
                  Year = ifelse(Year == 0,
                                0,
                                Year - data_list$styr + 1)
    )

  ration_data <- array(0, dim = c(data_list$nspp, max_sex, max_age, nyrs_hind)) # FIXME: Change for forecast

  if( nrow(data_list$ration_data) > 0){
    for (i in 1:nrow(data_list$ration_data)) {
      sp <- data_list$ration_data$Species[i]
      sex <- data_list$ration_data$Sex[i]
      yr <- data_list$ration_data$Year[i]

      # Handle sex == 0 case for 2-sex species
      sex_values <- if(sex == 0) 1:data_list$nsex[sp] else sex

      if(yr <= (data_list$endyr - data_list$styr + 1)){

        # Handle year == 0 case
        if(yr == 0) {
          yr <- seq_len(nyrs_hind)
        }

        # Fill in years
        for(j in 1:length(sex_values)){
          ration_data[sp, sex_values[j], 1:data_list$nages[sp], yr] <- as.numeric(
            as.character(
              as.matrix(data_list$ration_data[i, (1:data_list$nages[sp]) + 3])
            )
          )
        }
      }
    }
  }
  data_list$ration_data <- ration_data


  # 15 - Remove species column from alw, maturity, sex_ratio ----
  data_list$sex_ratio <- data_list$sex_ratio[,-1]
  data_list$maturity <- data_list$maturity[,-1]
  # data_list$aLW <- data_list$aLW[,-1]


  # 16 - Make data.frames into matrices ----
  df_to_mat <- which(sapply(data_list, function(x) class(x)[1]) == "data.frame")
  data_list[df_to_mat] <- lapply(data_list[df_to_mat], as.matrix)

  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "catch_data",    "index_data", "comp_data", "env_data", "spnames",
                       "aLW", "diet_data", # "NByageFixed", "estDynamics", "Ceq",
                       "avgnMode", "minNByage", "weight", "fleet_control")
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

  # If no data, convert to empty matrix
  if(is.null(dim(data_list$comp_obs))){
    data_list$comp_obs <- matrix(NA, ncol = 10, nrow = 0)
  } else{


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
  }

  return(data_list)
}


#' Check and Clean CAAL Data
#'
#' This function checks the CAAL data for zero sum rows, handles NA values,
#' and verifies that the CAAL data spans the required range of ages
#'
#' @param data_list A list containing the following components:
#'   - caal_obs: A matrix or data frame of composition observations.
#'   - caal_data: A data frame with metadata including species, sex, and age/length information.
#'   - nsex: A vector of sex counts by species.
#'   - nages: A vector of age counts by species.
#'   - nlengths: A vector of length counts by species.
#' throws An error if any rows of `caal_obs` sum to 0 or if the composition data does not span the range of ages
#'
#' @return The modified `data_list` with NA values in `caal_obs` converted to 0.
#' cleaned_data_list <- check_caal_data(data_list)
#'
#' @export
check_caal_data <- function(data_list) {

  # If no data, convert to empty matrix
  if(is.null(dim(data_list$caal_obs))){
    data_list$caal_obs <- matrix(NA, ncol = 10, nrow = 0)
  } else{

    # Check for zero sum rows in composition data
    if (any(rowSums(data_list$caal_obs, na.rm = TRUE) == 0)) {
      stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")
    }

    # Calculate adjustments for joint sex and CAAL data
    joint_adjust <- ifelse(data_list$nsex[data_list$caal_data$Species] == 2 &
                             data_list$caal_data$Sex == 3, 2, 1)

    col_adjust <- data_list$nages[data_list$comp_data$Species]
    col_adjust <- col_adjust * joint_adjust

    # Check for NAs in composition data and convert them to 0
    na_check <- apply(cbind(col_adjust, data_list$caal_obs), 1, function(x) sum(x[2:(x[1] + 1)]))

    if (any(is.na(na_check))) {
      na_rows <- which(is.na(na_check))
      warning(sprintf("Composition data have NAs in row(s): %s. Converting to 0s.",
                      paste(na_rows, collapse = ", ")))
      data_list$caal_obs[is.na(data_list$caal_obs)] <- 0
    }

    # Verify CAAL data size
    if (ncol(data_list$caal_obs) < max(col_adjust)) {
      stop("CAAL data does not span range of ages")
    }
  }

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

