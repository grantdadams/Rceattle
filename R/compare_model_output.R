# Function to compare output of models to check if there are within a certain threshold of relative error

compare_output <- function( rep = rep, tmp = tmp, cut_off = 0.01){
  param_check <- list()
  tmb_names <- c("fsh_age_obs","fsh_age_hat", "F", "pmature","tc_biom_hat", "biomass", "srv_bio_hat", "NByage", "R", "S", "srv_sel", "biomassSSB", "srv_age_obs", "srv_age_hat", "biomassByage", "fsh_sel", "M1", "tc_hat", "biomassSSB")
  admb_names <- c("fsh_age_obs","fsh_age_hat", "F", "pmature","obs_catch_hat", "Biomass", "srv_bio_hat", "NByage", "R", "S", "srv_sel", "BiomassSSB", "srv_age_obs", "srv_age_hat", "biomassByage", "fsh_sel", "M1", "tc_hat", "BiomassSSB")

  # Survey selectivity
  for( j in 1:length(tmb_names)){
    param <- tmb_names[j]
    admb_param <- admb_names[j]
    tmb_params <- rep[[param]]
    param_check[[param]] <- rep[[param]]
    if(class(tmb_params) == "matrix"){
      for(i in 1:nrow(tmb_params)){
        param_check[[param]][i,] <- replace(param_check[[param]][i,], values = rep(NA, length(param_check[[param]][i,])))
        admb_params <- tmp[[paste(admb_param, i, sep = "_")]]
        if(class(admb_params) == "matrix" ){
          admb_params <- admb_params[1,]
        }
        diff <- abs((tmb_params[i,1:length(admb_params)] - admb_params)) / admb_params < cut_off
        param_check[[param]][i, 1:length(admb_params)] <- replace( param_check[[param]][i,1:length(admb_params)], values = diff)
      }
    }
    if( class(tmb_params) == "array"){
      param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
      for(i in 1:dim(tmb_params)[3]){

        admb_params <- tmp[[paste(param, i, sep = "_")]]
        diff <- abs((tmb_params[1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params)), i] - admb_params[1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params))]) / admb_params[1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params))]) < cut_off
        param_check[[param]][1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params)), i] <- replace( param_check[[param]][1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params)), i], values = diff)
      }
    }
  }

  # EIT age hat
  param <- "eit_age_comp"
  admb_param <- "obs_eit_age"
  tmb_params <- rep[[param]]
  param_check[[param]] <- rep[[param]]
  param_check[[param]][] <- NA
  admb_params <- tmp[[admb_param]]
  diff <- abs((tmb_params - admb_params) / admb_params) < cut_off
  param_check[[param]] <- replace( param_check[[param]], values = diff)

  # EIT selectivity
  param <- "eit_hat"

  param_check[[param]] <- rep[[param]]
  tmb_params <- rep[[param]]
  param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
  admb_params <- tmp[[param]]
  tmb_params <- tmb_params[which(tmb_params > 0)]
  diff <- abs((tmb_params[1:length(admb_params)] - admb_params) / admb_params) < cut_off
  param_check[[param]][1:length(admb_params)] <- replace( param_check[1:length(admb_params)], values = diff)



  # Do a summary of how many items are wrong
  check_summary <- data.frame( parameter = names(param_check), n_wrong = rep(NA, length(names(param_check))))
  for(i in 1:length(names(param_check))){
    check_summary[i,2] = length(which(param_check[[i]] == 0))
  }

  # Check obejctive functions
  tmb_like <- rep$jnll_comp
  admb_like <- tmb_like
  admb_like[] <- 0
  admb_like[1,] <- c(tmp$srv_bio_like_1, tmp$srv_bio_like_2, tmp$srv_bio_like_3)
  admb_like[2,] <- c(tmp$srv_age_like_1, tmp$srv_age_like_2, tmp$srv_age_like_3)
  admb_like[3,] <- c(tmp$eit_srv_like, 0, 0)
  admb_like[4,] <- c(tmp$eit_age_like, 0, 0)
  admb_like[5,] <- c(tmp$fsh_cat_like_1, tmp$fsh_cat_like_2, tmp$fsh_cat_like_3)
  admb_like[6,] <- c(tmp$fsh_age_like_1, tmp$fsh_age_like_2, tmp$fsh_age_like_3)
  admb_like[7,] <- c(tmp$fsh_sel_like_1, tmp$fsh_sel_like_2, tmp$fsh_sel_like_3)
  admb_like[9,] <- c(tmp$srv_sel_like_1, tmp$srv_sel_like_2, tmp$srv_sel_like_3)

  like_check <- abs((tmb_like - admb_like) / admb_like) < cut_off
  like_check <- as.data.frame(like_check)
  admb_like <- round(admb_like,3)
  tmb_like <- round(tmb_like,3)
  admb_like <- as.data.frame(admb_like)
  tmb_like <- as.data.frame(tmb_like)


  like_names <- c("srv_bio_like", "srv_age_like", "eit_srv_like" , "eit_age_like", "fsh_cat_like", "fsh_age_like", "fsh_sel_like", NA , "srv_sel_like", NA, NA, NA, NA)
  like_check$names <- like_names
  tmb_like$names <- like_names
  admb_like$names <- like_names

  return(list(param_check, check_summary, like_check, tmb_like, admb_like))
}


# Note: TC-Hat and Fsh_age_hat will be off, but like is good
res <- compare_output(rep, tmp)

res[[2]]
res[[3]]
res[[4]]
res[[5]]


opt$objective
rep$jnll
tmp$obj_fun
