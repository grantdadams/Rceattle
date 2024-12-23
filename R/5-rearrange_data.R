#' Rearrange
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB CEATTLE
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @export
rearrange_dat <- function(data_list){
  '%!in%' <- function(x,y)!('%in%'(x,y))

  # Step 0 - Max functions on data
  data_list$max_bin <- max(data_list$nlengths) # Integer of maximum number of length/age bins.
  data_list$max_nsex <- max(data_list$nsex)
  data_list$max_nages <- max(data_list$nages)

  # Step 1 - remove non-integer objects from control
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
                  Weight_index,         # 12) Dim1 of wt (what weight-at-age data set)
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

  # - Input missing nselages and age first selected (use age range)
  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(Nselages = ifelse(is.na(Nselages), -999, Nselages),
                  Age_max_selected = ifelse(is.na(Age_max_selected), -999, Age_max_selected),
                  Age_first_selected = ifelse(is.na(Age_first_selected), data_list$minage[Species], Age_first_selected)
    )

  # Step 2 -  Seperate survey biomass info from observation
  data_list$index_ctl <- data_list$index_data[,c("Fleet_code", "Species", "Year")]
  data_list$index_n <- as.matrix(data_list$index_data[,c("Month")])
  data_list$index_obs <- data_list$index_data[,c("Observation", "Log_sd")]

  # Step 3 -  Seperate catch biomass info from observation
  data_list$catch_ctl <- data_list$catch_data[,c("Fleet_code", "Species", "Year")]
  data_list$catch_n <- as.matrix(data_list$catch_data[,c("Month")])
  data_list$catch_obs <- data_list$catch_data[,c("Catch", "Log_sd")]

  # Step 4 -  Seperate survey comp info from observation
  data_list$comp_ctl <- data_list$comp_data[,c("Fleet_code", "Species", "Sex", "Age0_Length1", "Year")]
  data_list$comp_n <- data_list$comp_data[,c("Month", "Sample_size")]
  data_list$comp_obs <- data_list$comp_data[,grep("Comp_", colnames(data_list$comp_data))]

  if(sum(rowSums(data_list$comp_obs, na.rm = TRUE) %in% 0) > 0){stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")}

  # - NA to 0 if in obs
  joint_adjust <- ifelse(data_list$nsex[data_list$comp_data$Species] == 2 & data_list$comp_data$Sex == 0, 2, 1)
  col_adjust <- ifelse(data_list$comp_data$Age0_Length1 == 0,
                       data_list$nages[data_list$comp_data$Species],
                       data_list$nlengths[data_list$comp_data$Species])
  col_adjust <- col_adjust * joint_adjust
  na_check <- apply(cbind(col_adjust, data_list$comp_obs), 1, function(x) sum(x[2:(x[1] + 1)]))
  if(is.na(sum(na_check))){
    warning(paste0("Composition data have NAs in row ", paste(which(is.na(na_check)), collapse = ", "), ". Converting to 0s"))
  }
  data_list$comp_obs[is.na(data_list$comp_obs)] <- 0

  # Step 5 -  Seperate uobs info from observation
  data_list$stom_prop_ctl <- data_list$stom_prop_data[,c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age", "Year")]
  data_list$stom_prop_obs <- data_list$stom_prop_data[,c("Sample_size", "Stomach_proportion_by_weight")]

  # Step 6 -  Seperate survey empirical selectivity info from observation
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

  data_list$emp_sel_ctl <- as.matrix(data_list$emp_sel[,c("Fleet_code", "Species", "Sex", "Year")])
  data_list$emp_sel_obs <- matrix(as.numeric(unlist(data_list$emp_sel[,grep("Comp_", colnames(data_list$emp_sel))])), nrow = nrow(data_list$emp_sel_ctl))


  # Step 7 - Rearrange age-transition matrix ----
  age_trans_matrix <- data_list$age_trans_matrix
  unique_age_transition <- unique(as.character(age_trans_matrix$Age_transition_index))
  if(sum(data_list$pop_age_transition_index %!in% unique_age_transition) > 0){
    stop("Check population age_transition index, not in age_transition file")
  }
  age_transition <- array(0, dim = c(length(unique_age_transition), 2, max(data_list$nages, na.rm = T), max(data_list$nlengths, na.rm = T)))


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

    # Assign
    if(sex == 0){ sex = c(1, 2)}
    for(j in 1:length(sex)){
      age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]] <- as.numeric(as.character(age_trans_matrix[i, (1:data_list$nlengths[sp]) + 5]))

      # Normalize
      age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]] <- age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]] / sum(age_transition[age_transition_ind, sex[j], age, 1:data_list$nlengths[sp]], na.rm = TRUE)
    }
  }
  data_list$age_trans_matrix <- age_transition


  # Step 8 - Rearrange age_error matrices ----
  arm <- array(0, dim = c(data_list$nspp, max(data_list$nages, na.rm = T), max(data_list$nages, na.rm = T)))

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


  # Step 10 - Normalize comp data ----
  data_list$comp_obs <- t(apply(data_list$comp_obs, 1, function(x) as.numeric(x) / sum(as.numeric(x), na.rm = TRUE)))


  # Step 11 - Set up weight-at-age ----
  # - Convert to array
  data_list$wt <- data_list$wt %>%
    mutate(
      Wt_index = as.numeric(as.character(Wt_index)),
      Species = as.numeric(as.character(Species)),
      Sex = as.numeric(as.character(Sex)),
      Year = as.numeric(as.character(Year)) - data_list$styr + 1)

  unique_wt <- unique(as.numeric(data_list$wt$Wt_index))
  if(sum(data_list$pop_wt_index %!in% unique_wt) > 0){
    stop("Check population weight index, not in weight file")
  }

  if(sum(data_list$ssb_wt_index %!in% unique_wt) > 0){
    stop("Check SSB weight index, not in weight file")
  }

  wt <- array(0, dim = c(length(unique_wt), 2, max(data_list$nages, na.rm = T), length(data_list$styr:data_list$endyr)))

  for (i in 1:nrow(data_list$wt)) {

    wt_ind <- data_list$wt$Wt_index[i]
    sp <- data_list$wt$Species[i]
    sex <- data_list$wt$Sex[i]
    yr <- data_list$wt$Year[i]

    if(yr <= data_list$endyr - data_list$styr + 1){

      # If year == 0, set weight to all years
      if(yr == (-data_list$styr + 1)){
        yr = 1:length(data_list$styr:data_list$endyr)
      }

      if(sex == 0){ sex = c(1, 2)}
      for(j in 1:length(sex)){
        if(sum(grepl("[[:space:]]", as.character(data_list$wt[i, (1:data_list$nages[sp]) + 5])))){
          stop(paste("Space found in wt data: row", i))
        }
        wt[wt_ind, sex[j], 1:data_list$nages[sp], yr] <- as.numeric(as.character(data_list$wt[i, (1:data_list$nages[sp]) + 5]))
      }
    }
  }
  data_list$wt <- wt


  # Step 12 - Set up NByageFixed ----
  # - Convert to array
  data_list$NByageFixed <- data_list$NByageFixed %>%
    mutate(Species = as.numeric(as.character(Species)),
           Sex = as.numeric(as.character(Sex)),
           Year = as.numeric(as.character(Year)) - data_list$styr + 1)

  NByageFixed <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T), length(data_list$styr:data_list$projyr)))

  if(nrow(data_list$NByageFixed) > 0){
    for (i in 1:nrow(data_list$NByageFixed)) {

      sp <- data_list$NByageFixed$Species[i]
      sex <- data_list$NByageFixed$Sex[i]
      yr <- data_list$NByageFixed$Year[i]

      if(sex == 0){ sex = c(1, 2)}

      for(j in 1:length(sex)){
        if(yr > 0){
          NByageFixed[sp, sex[j], 1:max(data_list$nages, na.rm = T), yr] <- as.numeric(data_list$NByageFixed[i,c(1:max(data_list$nages, na.rm = T))+4])
        }
      }
    }
  }

  data_list$NByageFixed <- NByageFixed


  # Step 13 - Set up environmental indices
  # - Fill in missing years with column mean
  data_list$env_index <- merge(data_list$env_data, data.frame(Year = data_list$styr:data_list$projyr), all = TRUE)
  data_list$env_index <-  data_list$env_index %>%
    dplyr::select(-Year) %>%
    mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%
    as.matrix()

  # - Create matrix for srr curve
  if(sum(sapply(data_list$srr_env_indices, function(x) x>ncol(data_list$env_index))) > 0){stop("srr_env_indices greater than the number of indices included")}
  data_list$env_index_srr <-  data_list$env_index[, data_list$srr_env_indices] %>%
    as.matrix()


  # Step 14 - Set up pyrs ----
  # - Convert to array
  data_list$Pyrs <- data_list$Pyrs %>%
    mutate(Species = as.numeric(as.character(Species)),
           Sex = as.numeric(as.character(Sex)),
           Year = as.numeric(as.character(Year)) - data_list$styr + 1)

  Pyrs <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages), length(data_list$styr:data_list$endyr))) # FIXME: Change for forecast

  if(nrow(data_list$Pyrs)>0){
    for (i in 1:nrow(data_list$Pyrs)) {
      sp <- data_list$Pyrs$Species[i]
      sex <- data_list$Pyrs$Sex[i]
      yr <- data_list$Pyrs$Year[i]

      if(sex == 0){
        sex = c(1,2)
      }

      if(yr <= (data_list$endyr - data_list$styr + 1)){
        for(j in 1:length(sex)){
          Pyrs[sp, sex[j], 1:data_list$nages[sp], yr] <- as.numeric(as.character(data_list$Pyrs[i, (1:data_list$nages[sp]) + 3]))
        }
      }
    }
  }
  data_list$Pyrs <- Pyrs


  # Step 15 - Remove species column from alw, pmature, sex_ratio
  data_list$sex_ratio <- data_list$sex_ratio[,-1]
  data_list$pmature <- data_list$pmature[,-1]
  data_list$aLW <- data_list$aLW[,-1]


  # Step 16 - Make data.frames into matrices
  df_to_mat <- which(sapply(data_list, function(x) class(x)[1]) == "data.frame")
  data_list[df_to_mat] <- lapply(data_list[df_to_mat], as.matrix)

  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "catch_data",    "index_data", "comp_data", "env_data", "spnames",
                       "aLW", "NByageFixed", "estDynamics", "minNByage")
  data_list[items_to_remove] <- NULL

  return(data_list)
}
