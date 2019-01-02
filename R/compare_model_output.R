#' Compare TMB vs ADMB
#' @description Function to compare output of TMB and ADMB based CEATTLE models to check if there are within a certain threshold of relative error. 0 means it is not within the tolerance, 1 means it is.
#'
#' @param rep TMB CEATTLE estimated quantities from exported from \code{\link{Rceattle}}
#' @param tmp list of ADMB CEATTLE estimated quantities from saved .Rdata file
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param rel_error relative error tolerance
compare_output <- function( rep = rep, tmp = tmp, data_list = data_list, rel_error = 0.01){

  param_check <- list()
  tmb_names <- c("fsh_age_obs","fsh_age_hat", "F", "pmature","tc_biom_hat", "biomass", "srv_bio_hat", "NByage", "R", "S", "srv_sel", "biomassSSB", "srv_age_obs", "srv_age_hat", "biomassByage", "fsh_sel", "M1", "tc_hat", "biomassSSB", "ration2Age", "AvgN", "M2", "avail_food", "suit_other", "of_stomKir", "stomKir", "stom_div_bio2")
  admb_names <- c("fsh_age_obs","fsh_age_hat", "F", "pmature","obs_catch_hat", "Biomass", "srv_bio_hat", "NByage", "R", "S", "srv_sel", "BiomassSSB", "srv_age_obs", "srv_age_hat", "biomassByage", "fsh_sel", "M1", "tc_hat", "BiomassSSB", "ration2", "AvgN", "M2", "avail_food", "suit_other", "of_stomKir", "stomKir", "stom_div_bio2")

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
        diff <- abs((tmb_params[i,1:length(admb_params)] - admb_params)) / admb_params < rel_error
        param_check[[param]][i, 1:length(admb_params)] <- replace( param_check[[param]][i,1:length(admb_params)], values = diff)
      }
    }
    if( class(tmb_params) == "array"){

      # 3D Array
      if(length(dim(tmb_params)) == 3){
        param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
        for(i in 1:dim(tmb_params)[3]){
          admb_params <- tmp[[paste(admb_param, i, sep = "_")]]

          diff <- abs((tmb_params[1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params)), i] - admb_params[1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params))]) / admb_params[1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params))]) < rel_error

          param_check[[param]][1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params)), i] <- replace( param_check[[param]][1:min(nrow(tmb_params[,,i]) , nrow(admb_params)), 1:min(ncol(tmb_params[,,i]) , ncol(admb_params)), i], values = diff)
        }
      }
#
#       # 5D Array
#       if(length(dim(tmb_params)) == 5){
#         param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
#
#         for(sp in 1:dim(tmb_params)[1]){
#           for(sp_age in 1:data_list$nages[sp]){
#             for(yr in 1:data_list$nyrs){
#               admb_params <- tmp[[paste(admb_param, yr, sp, sp_age, sep = "_")]]
#               tmb_tmp <- tmb_params[sp, , sp_age, ,yr]
#               diff <- abs(tmb_tmp -  admb_params) / (admb_params + 1e-16) < rel_error
#               param_check[[param]][sp, , sp_age, ,yr] <- replace( param_check[[param]][sp, , sp_age, ,yr], values = diff)
#             }
#           }
#         }
#       }
    }
  }

  # EIT age hat
  param <- "eit_age_comp"
  admb_param <- "obs_eit_age"
  tmb_params <- rep[[param]]
  param_check[[param]] <- rep[[param]]
  param_check[[param]][] <- NA
  admb_params <- tmp[[admb_param]]
  diff <- abs((tmb_params - admb_params) / admb_params) < rel_error
  param_check[[param]] <- replace( param_check[[param]], values = diff)

  # EIT selectivity
  param <- "eit_hat"
  param_check[[param]] <- rep[[param]]
  tmb_params <- rep[[param]]
  param_check[[param]] <- replace(param_check[[param]], values = rep(NA, length(param_check[[param]])))
  admb_params <- tmp[[param]]
  tmb_params <- tmb_params[which(tmb_params > 0)]
  diff <- abs((tmb_params[1:length(admb_params)] - admb_params) / admb_params) < rel_error
  param_check[[param]][1:length(admb_params)] <- replace( param_check[[param]][1:length(admb_params)], values = diff)

  # EIT age comp hat
  param <- "eit_age_comp_hat" #FIXME - get eit age hat from


  # Do a summary of how many items are wrong
  check_summary <- data.frame(matrix(NA, ncol = tmp$nspp + 1,nrow = length(param_check)))
  colnames(check_summary) <- c("Parameter", paste0("sp_", 1:tmp$nspp))
  check_summary$Parameter <- names(param_check)
  for(i in 1:length(names(param_check))){
    for(j in 1:tmp$nspp){

      # Vectors
      if(names(param_check)[i] == "eit_hat"){
        check_summary[i,2] = length(which(param_check[[i]] == 0))
      }

      # Matrices
      if(length(dim(param_check[[i]])) == 2){
        check_summary[i,j+1] = length(which(param_check[[i]][j,] == 0))
      }

      # 3D Arrays
      if(length(dim(param_check[[i]])) == 3){
        check_summary[i,j+1] = length(which(param_check[[i]][,,j] == 0))
      }

      # 5D Arrays
      if(length(dim(param_check[[i]])) == 5){
        check_summary[i,j+1] = length(which(param_check[[i]][j,,,,] == 0))
      }
    }
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

  like_check <- abs((tmb_like - admb_like) / admb_like) < rel_error
  like_check <- as.data.frame(like_check)
  admb_like <- round(admb_like,3)
  tmb_like <- round(tmb_like,3)
  admb_like <- as.data.frame(admb_like)
  tmb_like <- as.data.frame(tmb_like)


  like_names <- c("srv_bio_like", "srv_age_like", "eit_srv_like" , "eit_age_like", "fsh_cat_like", "fsh_age_like", "fsh_sel_like", NA , "srv_sel_like", NA, NA, NA, NA)
  like_check$names <- like_names
  tmb_like$names <- like_names
  admb_like$names <- like_names

  return(list(summary = check_summary, parameters = param_check, likelihood = like_check, tmb_like = tmb_like, admb_like = admb_like))
}


