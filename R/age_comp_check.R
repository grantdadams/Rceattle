rep <- Rceattle()

species <- 2
NByage <- rep$NByage[,,species]
F = rep$F[,,species]
M1 = rep$M1[species,]
srv_sel = rep$srv_sel[species,]
ALK <- data_list$age_trans_matrix[,,species]


srv_age_hat <- sweep(NByage * exp(-0.5 * sweep(F, 2, M1, "+")), 2, srv_sel, "*")
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

srv_age_com == rep$srv_age_hat[1:36,,2]
