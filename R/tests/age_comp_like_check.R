

age_comp_like_comparison <- function(data_list, species, rep, tmp, ADMB_TMB){

  MNConst = 0.0001
  tau = 200

  # TMB DATA
  if(ADMB_TMB == 1){
    fsh_age_obs <- rep$fsh_age_obs[,,species]
    fsh_age_hat <- rep$fsh_age_hat[,,species]

  }

  # ADMB DATA
  if(ADMB_TMB == 2){
    fsh_age_obs <- tmp[[paste("fsh_age_obs", species, sep = "_")]]
    fsh_age_hat <- tmp[[paste("fsh_age_hat", species, sep = "_")]]
  }

  offset_fsh <- 0

  for(y in 1:data_list$nyrs_fsh_comp[species]){
    for(j in 1: data_list$fsh_age_bins[species]){
      offset_fsh = offset_fsh - (tau * (fsh_age_obs[y, j] + MNConst) * log(fsh_age_obs[y, j] + MNConst));
    }
  }
  offset_fsh
  nll_comp = 0

    for (y in 1:data_list$nyrs_fsh_comp[species]){
      for(j in 1: data_list$fsh_age_bins[species]){
      fsh_yr_ind = data_list$yrs_fsh_comp[species, y] - data_list$styr + 1;
      nll_comp = nll_comp - tau * (fsh_age_obs[y, j] + MNConst) * log(fsh_age_hat[y, j] + MNConst );
      }
      nll_comp = nll_comp - offset_fsh
    }

  round(nll_comp)

  target <- tmp[[paste("fsh_age_like", species, sep = "_")]]
  my_nll <- rep$jnll_comp[6, species]
  my_nll
  target
}
