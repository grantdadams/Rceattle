sp = 3
offset_srv <- rep(0, 3)
srv_age_like <- rep(0, 3)
for(sp in 1:3){
obs1 <- tmp[[paste0("srv_age_obs_", sp)]]
hat1 <- tmp[[paste0("srv_age_hat_", sp)]]


for (y in 1: data_list$nyrs_srv_age[sp]){
  offset_srv[sp] = offset_srv[sp]  - 100 * sum((obs1[y,] + 1e-4) * log(obs1[y,] + 1e-4 )) ;
}


for (y in 1: data_list$nyrs_srv_age[sp]){
  srv_age_like[sp] = srv_age_like[sp] - 100 * sum((obs1[y,]+1e-4) * log(hat1[y,]+1e-4));
}
srv_age_like - offset_srv
}

tmp$srv_age_like_3
