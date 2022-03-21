spr_fun <- function(Rceattle){
  
  for(sp in 1:Rceattle$data_list$nspp){
    nages <- Rceattle$data_list$nages[sp]
    ages <- (1:nages)-1+Rceattle$data_list$minage[sp]
    wt <- Rceattle$data_list$wt[which(Rceattle$data_list$wt$Wt_index == Rceattle$data_list$ssb_wt_index[sp]),5+1:nages]
    wt <- c(0.011, 0.101, 0.267, 0.499, 0.638, 0.857, 1.043, 1.198, 1.409, 1.435)
    mat <- Rceattle$data_list$pmature[which(Rceattle$data_list$pmature$Species == sp), 1+1:nages]
    sexratio = Rceattle$data_list$sex_ratio[which(Rceattle$data_list$sex_ratio$Species == sp), 1+1:nages]
    M1 <- Rceattle$quantities$M1[sp, 1, 1:nages]
    spawn_mo <- Rceattle$data_list$spawn_month[sp]
    
    N = 1
    for(age in 2:(nages-1)){
      N[age] = N[age-1] * exp(-M1[age-1])
    }
    N[nages] = N[nages-1]*exp(-M1[nages-1])/(1-exp(-M1[nages]))
    SP0 = sum( N * wt * mat * sexratio * exp(-M1 * spawn_mo/12)) * exp(Rceattle$estimated_params$ln_mean_rec[sp])
  }
  return(sprF/spr0)
}