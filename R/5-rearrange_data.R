#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB
#'
#' @param data_list an Rceattle data_list
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr n
#' @importFrom tidyselect contains
rearrange_dat <- function(data_list){

  # Convert text to integer for switches used in TMB
  data_list <- convert_switches(data_list)

  # Data dimensions
  max_sex <- max(data_list$nsex, na.rm = TRUE)
  max_age <- max(data_list$nages, na.rm = TRUE)
  max_length <- max(data_list$nlengths, na.rm = TRUE)
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)

  # 1 - Fleet control ----
  # - 0) Vector to save  species
  data_list$flt_spp <- data_list$fleet_control %>%
    dplyr::pull(.data$Species) %>% as.integer() - 1

  # - 1) Fleet pointer
  data_list$flt_sel_ind <- data_list$fleet_control %>%
    dplyr::pull(.data$Fleet_code) %>% as.integer() - 1

  # - 2) Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
  data_list$flt_type <- data_list$fleet_control %>%
    dplyr::pull(.data$Fleet_type) %>% as.integer()

  # - 3) Month of observation
  data_list$flt_month <- data_list$fleet_control %>%
    dplyr::pull(.data$Month)

  # - 4) Selectivity type
  data_list$flt_sel_type <- data_list$fleet_control %>%
    dplyr::pull(.data$Selectivity) %>% as.integer()

  # - 4b) Selectivity penalty "lead": 0 if an EARLIER fleet shares this fleet's
  #       Selectivity_index AND selectivity type (i.e. this fleet mirrors an
  #       earlier one's estimated selectivity), else 1. The selectivity penalty is
  #       only accumulated for lead fleets so a mirrored selectivity is penalized
  #       once (matching ADMB). Requiring the same type means a fleet with empirical
  #       selectivity (type 0) sharing an index with a parametric fleet does NOT
  #       suppress the parametric fleet's penalty.
  {
    .sel_idx <- data_list$fleet_control$Selectivity_index
    .sel_typ <- data_list$flt_sel_type
    .lead <- rep(1L, length(.sel_typ))
    for(.i in seq_along(.sel_typ)) {
      if(.i > 1 && any(.sel_idx[seq_len(.i - 1)] == .sel_idx[.i] &
                       .sel_typ[seq_len(.i - 1)] == .sel_typ[.i], na.rm = TRUE)) .lead[.i] <- 0L
    }
    data_list$flt_sel_lead <- as.integer(.lead)
  }

  data_list$flt_sel_dim <- data_list$fleet_control %>%
    dplyr::mutate(Selectivity_dimension = dplyr::case_when(
      .data$Selectivity_dimension == "Age" ~ 0,
      .data$Selectivity_dimension == "Length" ~ 1,
      .default = NA
    )) %>%
    dplyr::pull(.data$Selectivity_dimension) %>% as.integer()

  # - 5) Number of ages/lengths for non-parametric selectivity
  data_list$flt_n_sel_bins <- data_list$fleet_control %>%
    dplyr::mutate(N_sel_bins = ifelse(is.na(.data$N_sel_bins), -999, .data$N_sel_bins)) %>%
    dplyr::pull(.data$N_sel_bins) %>% as.integer()

  # - 6) Time-varying selectivity type.
  data_list$flt_varying_sel <- data_list$fleet_control %>%
    dplyr::pull(.data$Time_varying_sel) %>% as.integer()

  # - 7) First age selected
  data_list$bin_first_selected <- data_list$fleet_control %>%
    dplyr::mutate(Bin_first_selected = .data$Bin_first_selected - 1, # R to C++ indexing
                  Bin_first_selected = ifelse(is.na(.data$Bin_first_selected), 0, .data$Bin_first_selected)) %>%
    dplyr::pull(.data$Bin_first_selected) %>% as.integer()

  # - 8) Age of max selectivity (used for normalization). If NA, does not normalize
  data_list$sel_norm_bin1 <- data_list$fleet_control %>%
    dplyr::mutate(
      Sel_norm_bin1 = .data$Sel_norm_bin1 - data_list$minage[Species],
      Sel_norm_bin1 = ifelse(.data$Sel_norm_bin1 < 0, -99, .data$Sel_norm_bin1),         # Less than zero, normalize by max
      Sel_norm_bin1 = ifelse(is.na(.data$Sel_norm_bin1), -999, .data$Sel_norm_bin1)) %>% # NA, do not normalize (unless type = 2)
    dplyr::pull(.data$Sel_norm_bin1) %>% as.integer()

  # - 9) upper age of max selectivity (used for normalization). If NA, does not normalize
  data_list$sel_norm_bin2 <- data_list$fleet_control %>%
    dplyr::mutate(Sel_norm_bin2 = .data$Sel_norm_bin2 - data_list$minage[Species],
                  Sel_norm_bin2 = ifelse(is.na(.data$Sel_norm_bin2), -999, .data$Sel_norm_bin2)) %>%
    dplyr::pull(.data$Sel_norm_bin2) %>% as.integer()

  # - 9b) Per-fleet selectivity start year (0-based from styr). Selectivity
  #       penalties begin the year after this (excludes pre-survey years + the
  #       start-year boundary). NA -> styr (index 0). Used by LogisticPM (type 11).
  data_list$flt_sel_start_yr <- data_list$fleet_control %>%
    dplyr::mutate(Sel_start_year = ifelse(is.na(.data$Sel_start_year), data_list$styr, .data$Sel_start_year),
                  Sel_start_year = .data$Sel_start_year - data_list$styr) %>%
    dplyr::pull(.data$Sel_start_year) %>% as.integer()

  # - 9c) First age (0-based) for the non-parametric shape penalty. NA -> -999
  #       (cpp falls back to bin_first_selected). Lets the ascending/descending
  #       shape constraint span a narrower age-range than the first selected age
  #       (e.g. ATS: age-1 is selected but the shape penalty starts at mina_ats).
  data_list$flt_sel_pen_first_age <- data_list$fleet_control %>%
    dplyr::mutate(Sel_pen_first_age = .data$Sel_pen_first_age - data_list$minage[Species],
                  Sel_pen_first_age = ifelse(is.na(.data$Sel_pen_first_age), -999, .data$Sel_pen_first_age)) %>%
    dplyr::pull(.data$Sel_pen_first_age) %>% as.integer()

  # - 9d) Last (left) age (0-based) of the shape-penalty adjacent pairs. NA -> -999
  #       (cpp falls back to nages-2). With 9c, bounds the shape penalty to a sub-
  #       range (e.g. RTMB "rpm" fishery smoothness over ages 6-11, ATS 5-7).
  data_list$flt_sel_pen_last_age <- data_list$fleet_control %>%
    dplyr::mutate(Sel_pen_last_age = .data$Sel_pen_last_age - data_list$minage[Species],
                  Sel_pen_last_age = ifelse(is.na(.data$Sel_pen_last_age), -999, .data$Sel_pen_last_age)) %>%
    dplyr::pull(.data$Sel_pen_last_age) %>% as.integer()

  # - 9e) Non-parametric shape-penalty mode: "Directional" (0; sign of Sel_curve_pen1
  #       -> one-sided decreasing/increasing, ADMB/AMAK) or "Smooth" (1; two-sided
  #       d^2 over adjacent ages, RTMB "rpm"). NA/other -> 0.
  data_list$flt_sel_shape_mode <- data_list$fleet_control %>%
    dplyr::mutate(Sel_shape_mode = dplyr::case_when(
      .data$Sel_shape_mode %in% c(1, "Smooth", "smooth") ~ 1L,
      .default = 0L)) %>%
    dplyr::pull(.data$Sel_shape_mode) %>% as.integer()

  # - 9f) NonParametricRPM (type 9) age cap (0-based): the realized selectivity is
  #       held flat at/after this age (RTMB cap_old_age). NA -> -999 (no cap).
  data_list$flt_sel_cap_age <- data_list$fleet_control %>%
    dplyr::mutate(Sel_cap_age = .data$Sel_cap_age - data_list$minage[Species],
                  Sel_cap_age = ifelse(is.na(.data$Sel_cap_age), -999, .data$Sel_cap_age)) %>%
    dplyr::pull(.data$Sel_cap_age) %>% as.integer()

  # - 10) Index indicating whether to do dirichlet multinomial or a multinomial
  data_list$comp_ll_type <- data_list$fleet_control %>%
    dplyr::pull(.data$Comp_loglike) %>% as.integer()
  data_list$caal_ll_type <- data_list$fleet_control %>%
    dplyr::pull(.data$CAAL_loglike) %>% as.integer()
  data_list$diet_ll_type <- data_list$Diet_loglike %>%
    as.integer()

  # - 11) Index units (1 = weight, 2 = numbers)
  data_list$flt_units <- data_list$fleet_control %>%
    dplyr::pull(.data$Weight1_Numbers2) %>% as.integer()

  # - 12) Dim1 of weight (what weight-at-age data set)
  data_list$flt_wt_index <- data_list$fleet_control %>%
    dplyr::pull(.data$Weight_index) %>% as.integer() - 1

  # - 13) Dim3 of age transition matrix (what ALK to use)
  data_list$flt_age_transition_index <- data_list$fleet_control %>%
    dplyr::pull(.data$Age_transition_index) %>% as.integer() - 1

  # - 14) Parametric form of q
  data_list$est_index_q <- data_list$fleet_control %>%
    dplyr::pull(.data$Catchability) %>% as.integer()

  # - 15) Time varying q type
  data_list$index_varying_q <- data_list$fleet_control %>%
    dplyr::pull(.data$Time_varying_q) %>% as.integer()

  # - 16) Whether to estimate standard deviation of index time series
  data_list$est_sigma_index <- data_list$fleet_control %>%
    dplyr::pull(.data$Estimate_index_sd) %>% as.integer()

  # - 17) Whether to estimate standard deviation of fishery time series
  data_list$est_sigma_fsh <- data_list$fleet_control %>%
    dplyr::pull(.data$Estimate_catch_sd) %>% as.integer()

  data_list$index_log_q_prior <- log(data_list$fleet_control$Q_prior)

  # Species names
  data_list$spnames <- NULL

  # * F reference points ----
  data_list$Ftarget_percent <- data_list$Ftarget
  data_list$Flimit_percent <- data_list$Flimit


  # 2 -  Index data ----
  # - Seperate index metadata from observation
  data_list$index_ctl <- data_list$index_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$index_n <- as.matrix(data_list$index_data[,c("Month")])

  data_list$index_obs <- data_list$index_data %>%
    dplyr::select(Observation, Log_sd) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()


  # 3 -  Catch data ----
  # - Seperate catch metadata from observation
  data_list$catch_ctl <- data_list$catch_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$catch_n <- as.matrix(data_list$catch_data[,c("Month")])

  data_list$catch_obs <- data_list$catch_data %>%
    dplyr::select(Catch, Log_sd) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()


  # 4 -  Comp data ----
  # - Seperate comp metadata from observation
  data_list$comp_ctl <- data_list$comp_data %>%
    dplyr::select(Fleet_code, Species, Sex, Age0_Length1, Year) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$comp_n <- data_list$comp_data %>%
    dplyr::select(Month, Sample_size) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$comp_obs <- data_list$comp_data %>%
    dplyr::select(contains("Comp_")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  if(nrow(data_list$comp_obs) > 0){ #FIXME: probably cleaner way to deal with this
    data_list$comp_obs <- t(apply(data_list$comp_obs, 1, function(x) as.numeric(x) / sum(as.numeric(x), na.rm = TRUE))) # Normalize
    data_list$comp_obs[is.infinite(data_list$comp_obs)] <- 0
    data_list$comp_obs[is.na(data_list$comp_obs)] <- 0
  }
  data_list <- check_composition_data(data_list)

  # 5 - CAAL data ----
  data_list$caal_ctl <- data_list$caal_data %>%
    dplyr::mutate(Length_bin = factor(.data$Length)) %>%
    dplyr::select(Fleet_code, Species, Sex, Year, Length_bin) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$caal_n <- data_list$caal_data %>%
    dplyr::select(Sample_size) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$caal_obs <- data_list$caal_data %>%
    dplyr::select(contains("CAAL_")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  if(nrow(data_list$caal_obs) > 0){#FIXME: probably cleaner way to deal with this
    data_list$caal_obs <- t(apply(data_list$caal_obs, 1, function(x) as.numeric(x) / sum(as.numeric(x), na.rm = TRUE))) # Normalize
    data_list$caal_obs[is.infinite(data_list$caal_obs)] <- 0
    data_list$caal_obs[is.na(data_list$caal_obs)] <- 0
  }
  data_list <- check_caal_data(data_list)

  # * Get lengths from CAAL ----
  data_list$lengths <- matrix(rep(1:max(data_list$nlengths), data_list$nspp), nrow = data_list$nspp, byrow = TRUE)

  if(nrow(data_list$caal_data) > 0){
    caal_lengths <- data_list$caal_data %>%
      dplyr::distinct(.data$Species, .data$Length) %>%
      dplyr::arrange(.data$Species, .data$Length) %>%
      dplyr::group_by(.data$Species) %>%
      dplyr::mutate(Bin = paste0("Bin", 1:dplyr::n())) %>%
      tidyr::pivot_wider(names_from = "Bin", values_from = "Length")

    data_list$lengths[caal_lengths$Species, 1:ncol(caal_lengths[,-1])] <- as.matrix(caal_lengths[,-1])

    # - Check to make sure n-bins in data match nlengths
    n_bins_data <- data_list$caal_data %>%
      dplyr::distinct(.data$Species, .data$Length) %>%
      dplyr::count(.data$Species)

    if(any(data_list$nlengths[n_bins_data$Species] != n_bins_data$n)){
      stop("Number of length bins in CAAL data does not match nlengths in control file.
           If some lengths are missing from CAAL data, please add them to the data
           as rows of all 1s with Sample_size = 0")
    }
  }


  # 6 -  Diet data ----
  # - Seperate diet metadata from observation
  if(!is.null(data_list$diet_data)){ # Add a check in case there's no diet data

    # Create a temporary diet data frame to work with
    diet_dat <- data_list$diet_data

    # Add n_stomach_obs and the stomach_id vector to the main data_list
    if(length(diet_dat$stomach_id) == 0){
      data_list$n_stomach_obs = 0
    } else{
      data_list$n_stomach_obs <- max(diet_dat$stomach_id) + 1
    }
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
      Wt_index = as.numeric(as.character(.data$Wt_index)),
      Species = as.numeric(as.character(.data$Species)),
      Sex = as.numeric(as.character(.data$Sex)),
      Year = as.numeric(as.character(.data$Year)),
      Year = ifelse(.data$Year == 0,
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
    dplyr::mutate(Species = as.numeric(as.character(.data$Species)),
           Sex = as.numeric(as.character(.data$Sex)),
           Year = as.numeric(as.character(.data$Year)) - data_list$styr + 1)

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
    dplyr::mutate(Species = as.numeric(as.character(.data$Species)),
                  Sex = as.numeric(as.character(.data$Sex)),
                  Year = as.numeric(as.character(.data$Year)),
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


  # 15 - Remove species column from maturity and sex_ratio ----
  data_list$sex_ratio <- data_list$sex_ratio %>%
    dplyr::select(contains("Age")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$maturity <- data_list$maturity %>%
    dplyr::select(contains("Age")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()


  # 16 - Make data.frames into matrices ----
  df_to_mat <- which(sapply(data_list, function(x) class(x)[1]) == "data.frame")
  data_list[df_to_mat] <- lapply(data_list[df_to_mat], as.matrix)

  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "catch_data",    "index_data", "comp_data", "caal_data", "env_data", "spnames",
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
#' data_list <- list(
#'   comp_obs = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2),
#'   comp_data = data.frame(Species = c(1, 2), Sex = c(1, 1), Age0_Length1 = c(0, 0)),
#'   nsex = c(1, 1),
#'   nages = c(3, 3),
#'   nlengths = c(3, 3)
#' )
#' cleaned_data_list <- check_composition_data(data_list)
#'
#' @export
check_composition_data <- function(data_list) {

  # If no data, convert to empty matrix
  if(is.null(dim(data_list$comp_obs))){
    data_list$comp_obs <- matrix(NA, ncol = 10, nrow = 0)
  } else if(nrow(data_list$comp_obs) > 0){

    # # Check for zero sum rows in composition data
    # if (any(rowSums(data_list$comp_obs, na.rm = TRUE) == 0)) {
    #   stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")
    # }

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
      message(sprintf("Composition data have NAs in row(s): %s. Converting to 0s.",
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
  if(is.null(dim(data_list$caal_obs)) | nrow(data_list$caal_obs) == 0){
    data_list$caal_obs <- matrix(NA, ncol = 10, nrow = 0)
  } else{

    # # Check for zero sum rows in composition data
    # if (any(rowSums(data_list$caal_obs, na.rm = TRUE) == 0)) {
    #   stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")
    # }

    # Calculate adjustments for joint sex and CAAL data
    joint_adjust <- ifelse(data_list$nsex[data_list$caal_data$Species] == 2 &
                             data_list$caal_data$Sex == 3, 2, 1)

    col_adjust <- data_list$nages[data_list$caal_data$Species]
    col_adjust <- col_adjust * joint_adjust

    # Check for NAs in composition data and convert them to 0
    na_check <- apply(cbind(col_adjust, data_list$caal_obs), 1, function(x) sum(x[2:(x[1] + 1)]))

    if (any(is.na(na_check))) {
      na_rows <- which(is.na(na_check))
      message(sprintf("Composition data have NAs in row(s): %s. Converting to 0s.",
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

#' Convert intuitive text strings to integer switches for TMB
#'
#' @param data_list Rceattle data list
#'
#' @importFrom rlang .data
convert_switches <- function(data_list) {

  # Helper: convert a single string value using a map, pass integers through unchanged
  .conv_single <- function(x, map) {
    if (is.character(x) && x %in% names(map)) unname(map[[x]]) else x
  }
  .conv <- Vectorize(.conv_single, vectorize.args = "x", USE.NAMES = FALSE)

  # Fleet controls ----
  # If vector is a string that exists in our map, replace it with the integer.
  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Catchability = .conv(.data$Catchability, q_map),
      Time_varying_q = .conv(.data$Time_varying_q, tv_q_map),
      Comp_loglike = .conv(.data$Comp_loglike, comp_loglike_map),
      CAAL_loglike = .conv(.data$CAAL_loglike, comp_loglike_map)
    ) %>%
    # CRITICAL: Force columns back to integers so TMB doesn't crash expecting ints but getting chars
    dplyr::mutate(
      Fleet_type = as.integer(.data$Fleet_type),
      Selectivity = as.integer(.data$Selectivity),
      Time_varying_sel = as.integer(.data$Time_varying_sel),
      Catchability = as.integer(.data$Catchability),
      Time_varying_q = as.integer(.data$Time_varying_q),
      Comp_loglike = as.integer(.data$Comp_loglike),
      CAAL_loglike = as.integer(.data$CAAL_loglike)
    )

  # Pop dy controls ----
  data_list$initMode <- as.integer(.conv(data_list$initMode, initMode_map))
  data_list$HCR <- as.integer(.conv(data_list$HCR, hcr_map))

  return(data_list)
}
