



set_osa_obs = function(data_list){

  # Set up data frame to store vectorized data
  obs.colnames <- c("Year", "Fleet_code", "Comp", "Type", "Observation")
  obs <- as.data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 5 components: ecov, fleet catch (log), index catch (log), paa catch, paa index
  # do all Ecov first
  # 1. Ecov
  if(any(data_list$fleet_control$Catchability == 6, na.rm = TRUE)){
    #
    # x <- as.data_list.frame(data_list$Ecov_obs)
    # x[data_list$Ecov_use_obs==0] <- NA # only include index data_list to fit in obsvec
    # colnames(x) <- paste0("Ecov_", 1:data_list$n_Ecov)
    # x$year <- seq(from=data_list$year1_Ecov-data_list$year1_model+1, length.out=data_list$n_years_Ecov) # don't assume Ecov and model years are the same
    # tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
    # tmp <- tmp[complete.cases(tmp),]
    # tmp$age <- NA
    # tmp$type <- "Ecov"
    # obs <- rbind(obs, tmp[, obs.colnames])
  }

  # Find out if there are any non-parametric selectivity parameters set to 0 and not estimated
  # an OSA residual for ages with selectivity assumed to be zero will be useless and causes estimation problems with conditional construction of age comp likelihoods.
  # Selectivy types:
  # - 2 = non-parametric selecitivty sensu Ianelli et al 2018
  # - 5 = non-parametric selectivity sensu Taylor et al 2014 (Hake)
  data_list$fleet_control$Age_first_selected
  ages_omit <- rep(0, nrow(data_list$fleet_control))
  for(flt in 1:nrow(data_list$fleet_control)){

    flt = data_list$fleet_control$Fleet_code[i]
    spp <- data_list$fleet_control$Species[i]

    if(!is.na(data_list$fleet_control$Age_first_selected[i])){
      ages_omit[flt] <- data_list$fleet_control$Age_first_selected[i] - data_list$minage[spp] # How many initial ages to remove
    }
  }


  # 1) Index data ----
  index_fleets <- data_list$fleet_control$Fleet_code[data_list$fleet_control$Fleet_type == 2] # Index fleets
  x <- data_list$index_data %>%
    dplyr::filter(Year > 0 & Fleet_code %in% index_fleets) %>% # Data is used in likelihood
    mutate(Comp = NA,
           Type = "logindex",
           Observation = log(Observation)) %>%
    select(Year, Fleet_code, Comp, Type, Observation)
  obs <- rbind(obs, x)



  # 2) Catch data ----
  catch_fleets <- data_list$fleet_control$Fleet_code[data_list$fleet_control$Fleet_type == 1]
  x <- data_list$catch_data %>%
    dplyr::filter(Year > 0 & Fleet_code %in% catch_fleets) %>% # Data is used in likelihood
    mutate(Comp = NA,
           Type = "logcatch",
           Observation = log(Catch)) %>%
    select(Year, Fleet_code, Comp, Type, Observation)
  obs <- rbind(obs, x)

  # 3) Composition data ----
  # - Calculate adjustments for sex and age/length
  joint_adjust <- ifelse(data_list$nsex[data_list$comp_data$Species] == 2 &
                           data_list$comp_data$Sex == 3, 2, 1)

  ncol_comp <- ifelse(data_list$comp_data$Age0_Length1 == 0,
                      data_list$nages[data_list$comp_data$Species],
                      data_list$nlengths[data_list$comp_data$Species])

  ages_omit <- ages_omit[data_list$comp_data$Fleet_code]

  comp_select <- apply(cbind(ncol_comp, joint_adjust), 1, function(x) (1:(x[1] * x[2]))) # Columns to select by row
  comp_on <- apply(cbind(ages_omit, ncol_comp, joint_adjust), 1, function(x) rep(rep(c(0, 1), times = c(x[1],x[2]-x[1])), x[3])) # Columns to set to 0 by row

  # - Loop through comp data
  for(i in 1:nrow(data_list$comp_data)){
    data_list$comp_ctl <- data_list$comp_data[,c("Fleet_code", "Species", "Sex", "Age0_Length1", "Year")]
    data_list$comp_n <- data_list$comp_data[,c("Month", "Sample_size")]
    comp_tmp <- data_list$comp_data[i, grep("Comp_", colnames(data_list$comp_data))][comp_select[[i]]]
    comp_tmp[!comp_on[[i]]] <- 0 # Set ages with selectivity not estimated to zero
    comp_tmp[is.na(comp_tmp)] <- 0
    comp_tmp <- comp_tmp/sum(comp_tmp) # Rescale

    # - If not estimated set to zero
    if(data_list$comp_data$Year[i] < 0 | data_list$comp_data$Sample_size == 0){
      comp_tmp[] <- 0
    }

    #FIXME: may want situations wheere only on comp has data to set everything to 0

    # Add to data
    x <- data.frame(
      Year = data_list$comp_data$Year[i],
      Fleet_code = data_list$comp_data$Fleet_code[i],
      Comp = comp_select[[i]],
      Type = "comp",
      Observation = comp_tmp)

    obs <- rbind(obs, x)
  }


}

# calculate obsvec indices in keep arrays
obs$ind <- 1:dim(obs)[1]


data_list$keep_E <- matrix(NA, nrow=data_list$n_years_Ecov, ncol=data_list$n_Ecov)
for(i in 1:data_list$n_Ecov){
  ind <- which(data_list$Ecov_use_obs[,i]==1)
  data_list$keep_E[ind,i] <- subset(obs, type=='Ecov' & fleet == paste0("Ecov_",i))$ind
}
# subtract 1 bc TMB indexes from 0
data_list$keep_E <- data_list$keep_E - 1

data_list$keep_C <- matrix(NA, nrow=data_list$n_years_catch, ncol=data_list$n_fleets)
for(y in 1:data_list$n_years_model) for(i in 1:data_list$n_fleets){
  if(data_list$use_agg_catch[y,i]==1){
    data_list$keep_C[y,i] <- subset(obs, year == y & type=='logcatch' & fleet == paste0("fleet_",i))$ind
  }
}
# subtract 1 bc TMB indexes from 0
data_list$keep_C <- data_list$keep_C - 1

data_list$keep_I <- matrix(NA, nrow=data_list$n_years_indices, ncol=data_list$n_indices)
for(y in 1:data_list$n_years_model) for(i in 1:data_list$n_indices){
  if(data_list$use_indices[y,i]==1){
    data_list$keep_I[y,i] <- subset(obs, year == y & type=='logindex' & fleet == paste0("index_",i))$ind
  }
}
# subtract 1 bc TMB indexes from 0
data_list$keep_I <- data_list$keep_I - 1

data_list$condition_no_osa = NULL #to condition on age comps for likelihoods we don't have osa residuals set up for.
data_list$subset_discrete_osa = NULL #age comp obs for likelihoods we need to specify as discrete obs.

data_list$keep_Cpaa <- array(NA, dim=c(data_list$n_fleets, data_list$n_years_model, 2))
for(i in 1:data_list$n_fleets) {
  for(y in 1:data_list$n_years_model) if(data_list$use_catch_paa[y,i]==1){
    tmp = subset(obs, year == y & type=='catchpaa' & fleet==paste0("fleet_",i))
    if(length(tmp$ind)) #should always be TRUE because use_paa changed above
    {
      data_list$keep_Cpaa[i,y,1:2] <- c(tmp$ind[1], length(tmp$ind))
      #if(data_list$age_comp_model_fleets[i] %in% 1:2) data_list$subset_discrete_osa = c(data_list$subset_discrete_osa, tmp$ind)
      #subset for oneStepPredict can't include these
      if(data_list$age_comp_model_fleets[i] %in% 8:10) data_list$condition_no_osa = c(data_list$condition_no_osa, tmp$ind)
    }
  }
}
# subtract 1 bc TMB indexes from 0
data_list$keep_Cpaa[,,1] <- data_list$keep_Cpaa[,,1] - 1

data_list$keep_Ipaa <- array(NA, dim=c(data_list$n_indices, data_list$n_years_model, 2))
for(i in 1:data_list$n_indices) {
  for(y in 1:data_list$n_years_model) if(data_list$use_index_paa[y,i]==1){
    tmp = subset(obs, year == y & type=='indexpaa' & fleet==paste0("index_",i))
    if(length(tmp$ind)) #should always be TRUE because use_paa changed above
    {
      data_list$keep_Ipaa[i,y,1:2] <- c(tmp$ind[1], length(tmp$ind))
      #if(data_list$age_comp_model_indices[i] %in% 1:2) data_list$subset_discrete_osa = c(data_list$subset_discrete_osa, tmp$ind)
      #subset for oneStepPredict can't include these
      if(data_list$age_comp_model_indices[i] %in% 8:10) data_list$condition_no_osa = c(data_list$condition_no_osa, tmp$ind)
    }
  }
}
# subtract 1 bc TMB indexes from 0
data_list$keep_Ipaa[,,1] <- data_list$keep_Ipaa[,,1] - 1

obs$cohort = as.numeric(obs$year) - as.numeric(obs$age)   #make cohort column. could be useful for analyzing age comp OSA residuals
data_list$obs <- obs
data_list$obsvec <- obs$val
data_list$agesvec <- obs$age #potentially needed for AR1 sigma correlation of logistic-normal paa obs.
data_list$do_osa = 0 #this will be changed when TMB::oneStepPredict is called by fit_wham
data_list$do_post_samp = rep(0,5) #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated

data_list$data_list = data_list
return(data_list)
}
