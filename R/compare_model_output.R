# Function to compare output of models to check if there are within a certain threshold of relative error

compare_output <- function( tmb = rep, admb = tmp, cut_off = 0.001){
  param_check <- list()
  tmb_names <- c("fsh_age_obs","fsh_age_hat", "F", "pmature","tc_biom_hat", "biomass", "srv_bio_hat", "NByage", "R", "S", "srv_sel", "biomassSSB", "srv_age_obs", "srv_age_hat", "biomassByage", "fsh_sel", "M1", "tc_hat", "biomassSSB")
  admb_names <- c("fsh_age_obs","fsh_age_hat", "F", "pmature","obs_catch_hat", "Biomass", "srv_bio_hat", "NByage", "R", "S", "srv_sel", "BiomassSSB", "srv_age_obs", "srv_age_hat", "biomassByage", "fsh_sel", "M1", "tc_hat", "BiomassSSB")

  # Survey selectivity
  for( j in 1:length(tmb_names)){
    param <- tmb_names[j]
    admb_param <- admb_names[j]
    tmb_params <- tmb[[param]]
    param_check[[param]] <- tmb[[param]]
    if(class(tmb_params) == "matrix"){
      for(i in 1:nrow(tmb_params)){
        param_check[[param]][i,] <- replace(param_check[[param]][i,], values = rep(NA, length(param_check[[param]][i,])))
        admb_params <- admb[[paste(admb_param, i, sep = "_")]]
        if(class(admb_params) == "matrix" ){
          admb_params <- admb_params[1,]
        }
        diff <- (tmb_params[i,1:length(admb_params)] - admb_params) / admb_params < cut_off
        param_check[[param]][i, 1:length(admb_params)] <- replace( param_check[[param]][i,1:length(admb_params)], values = diff)
      }
    }
    if( class(tmb_params) == "array"){
      param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
      for(i in 1:dim(tmb_params)[3]){

        admb_params <- admb[[paste(param, i, sep = "_")]]
        diff <- (tmb_params[1:nrow(admb_params), 1:ncol(admb_params), i] - admb_params) / admb_params < cut_off
        param_check[[param]][1:nrow(admb_params), 1:ncol(admb_params), i] <- replace( param_check[[param]][1:nrow(admb_params), 1:ncol(admb_params), i], values = diff)
      }
    }
  }

  # EIT selectivity
  param <- "eit_hat"

  param_check[[param]] <- tmb[[param]]
  tmb_params <- tmb[[param]]
  param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
  admb_params <- admb[[param]]
  tmb_params <- tmb_params[which(tmb_params > 0)]
  diff <- (tmb_params[1:length(admb_params)] - admb_params) / admb_params < cut_off
  param_check[[param]][1:length(admb_params)] <- replace( param_check[1:length(admb_params)], values = diff)

  # Do a summary of how many items are wrong
  check_summary <- data.frame( parameter = names(param_check), n_wrong = rep(NA, length(names(param_check))))
  for(i in 1:length(names(param_check))){
    check_summary[i,2] = length(which(param_check[[i]] == 0))
  }

  return(list(param_check, check_summary))
}


res <- compare_output(rep, tmp)
res[[2]]
