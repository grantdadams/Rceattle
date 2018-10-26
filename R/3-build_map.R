# Function to construct the TMB map argument for CEATTLE
# Grant Adams June 2018

build_map <- function(data_list, params, debug = FALSE, random_rec = FALSE) {

  map_list <- params

  # STEP 1 -- Convert map_list to seperate parameters
  for(i in 1:length(map_list)){
    map_list[[i]] <- replace(map_list[[i]], values = c(1:length(map_list[[i]])))
  }


  # STEP 2 -- NA out parameters not to be estimated
  # Initial population deviations
  for(i in 1:nrow(map_list$init_dev)){
    if((data_list$nages[i]-1) < ncol(map_list$init_dev)){
      map_list$init_dev[i,  (data_list$nages[i]) : ncol(map_list$init_dev) ] <- NA
    }
  }

  # Survey selectivity coefficients
  for( i in 1: nrow(map_list$srv_sel_coff)){
    if(data_list$logist_sel_phase[i] > 0){
      map_list$srv_sel_coff[i,] <- replace(map_list$srv_sel_coff[i,], values = rep(NA, length(map_list$srv_sel_coff[i,])))
    }
    if(data_list$logist_sel_phase[i] < 0){
      map_list$srv_sel_slp[i] <- NA
      map_list$srv_sel_inf[i] <- NA
    }
  }

  # Catchability of surveys
  map_list$log_srv_q <- replace(map_list$log_srv_q, values = rep(NA, length(map_list$log_srv_q)))

  # # Recruitment deviation sigmas - turn off if not estimating
  # if(random_rec == FALSE){
  #   map_list$ln_rec_sigma <- replace(map_list$ln_rec_sigma, values = rep(NA, length(map_list$ln_rec_sigma)))
  # }


  # STEP 3 - set up debug - I.E. turn off all parameters besides dummy
  map_list$dummy <- NA
  if(debug == TRUE){
    for(i in 1:length(map_list)){
      map_list[[i]] <- replace(map_list[[i]], values = rep(NA, length(map_list[[i]])))
    }

    map_list$dummy = 1

  }


  # STEP 4 -- Conver to factor
  for(i in 1:length(map_list)){
    map_list[[i]] <- factor(map_list[[i]])
  }


  return(map_list)
}
