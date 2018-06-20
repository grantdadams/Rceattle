# Function to construct the TMB map argument for CEATTLE
# Grant Adams June 2018

build_map <- function(data_list, params) {

  map_list <- params

  # STEP 1 -- Convert map_list to seperate parameters
  for(i in 1:length(map_list)){
    map_list[[i]] <- replace(map_list[[i]], values = c(1:length(map_list[[i]])))
  }


  # STEP 2 -- NA out parameters not to be estimated
  # Initial population deviations
  for(i in 1:nrow(map_list$init_dev)){
    if((data_list$nages[i] + 1) < ncol(map_list$init_dev)){
      map_list$init_dev[i,  (data_list$nages[i] + 1) : ncol(map_list$init_dev) ] <- NA
    }
  }

  # Survey selectivity coefficients
  for( i in 1: nrow(map_list$srv_sel_coff)){
    if(data_list$logist_sel_phase[i] > 0){
      map_list$srv_sel_coff[i,3:ncol(map_list$srv_sel_coff)] <- NA
    }

  }

  # STEP 3 -- Conver to factor
  for(i in 1:length(map_list)){
    map_list[[i]] <- factor(map_list[[i]])
  }

  return(map_list)
}
