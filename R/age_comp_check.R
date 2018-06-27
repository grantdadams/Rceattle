# Function to return relative error of survey age comps using data from TMB or ADMB
rep <- Rceattle()
species <- 1

age_comp_comparison <- function( species, rep, tmp, ADMB_TMB = 1){

  # TMB DATA
  if(ADMB_TMB == 1){
    NByage <- rep$NByage[,,species]
    F = rep$F[,,species]
    M1 = rep$M1[species,]
    srv_sel = rep$srv_sel[species,]
    ALK <- data_list$age_trans_matrix[,,species]
  }

  # ADMB DATA
  if(ADMB_TMB == 2){
    NByage <- tmp[[paste("NByage", species, sep = "_")]]
    F = tmp[[paste("F", species, sep = "_")]]
    M1 = tmp[[paste("M1", species, sep = "_")]]
    srv_sel = tmp[[paste("srv_sel", species, sep = "_")]]
    ALK <- data_list$age_trans_matrix[,,species]
  }

  srv_age_hat <- sweep(NByage, 2, srv_sel, "*") # sweep(NByage * exp(-0.5 * sweep(F, 2, M1, "+")), 2, srv_sel, "*")
  if(species == 1){
    srv_age_com <- (srv_age_hat/ rowSums(srv_age_hat, na.rm = T))
    srv_age_com <- srv_age_com[1:36, 1:12]
  }
  if(species > 1){
    srv_len_hat <- matrix(NA, nrow = 36, ncol = 25)
    for(i in 1:36){
      srv_len_hat[i,] <- srv_age_hat[i,1:c(12,21)[species-1]] %*% ALK[1:c(12,21)[species-1], 1:25]
    }
    srv_age_com <- (srv_len_hat/ rowSums(srv_len_hat, na.rm = T))
  }

  # srv_age_com == tmp[[paste("srv_age_hat", species, sep = "_")]]
  rel_error <- (srv_age_com - tmp[[paste("srv_age_hat", species, sep = "_")]] )/ tmp[[paste("srv_age_hat", species, sep = "_")]]
  return(rel_error)
}

results_tmb <- age_comp_comparison( species = 1, rep, tmp, ADMB_TMB = 1)
results_admb <- age_comp_comparison( species = 1, rep, tmp, ADMB_TMB = 2)
