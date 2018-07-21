

age_comp_like_comparison <- function(data_list, species, rep, tmp, ADMB_TMB){

  MNConst = 0.0001

  # TMB DATA
  if(ADMB_TMB == 1){
    fsh_age_obs <- rep$eit_age_comp
    fsh_age_hat <- rep$eit_age_hat

  }

  # ADMB DATA
  if(ADMB_TMB == 2){
    fsh_age_obs <- tmp[["obs_eit_age"]]
    fsh_age_hat <- tmp[["eit_age_hat"]]
  }

  offset_fsh <- 0

  for(y in 1:data_list$n_eit){
    for(j in 1: data_list$nages[1]){
      offset_fsh = offset_fsh - (data_list$eit_age_n[y] * (fsh_age_obs[y, j] + MNConst) * log(fsh_age_obs[y, j] + MNConst));
    }
  }
  offset_fsh
  nll_comp = 0

  for (y in 1:data_list$n_eit){
    for(j in 1:data_list$nages[1]){
      print(paste(y,j,(data_list$eit_age_n[y] * (fsh_age_obs[y, j] + MNConst) * log(fsh_age_hat[y, j] + MNConst )), sep = " "))
      nll_comp = nll_comp - (data_list$eit_age_n[y] * (fsh_age_obs[y, j] + MNConst) * log(fsh_age_hat[y, j] + MNConst ));
    }
  }
  nll_comp = nll_comp - offset_fsh

  round(nll_comp)

  target <- tmp[[paste("eit_age_like", species, sep = "_")]]
  my_nll <- rep$jnll_comp[4, 1]
  my_nll
  target
}
