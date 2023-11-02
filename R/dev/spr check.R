
M <- ss_run$data_list$M1_base[1,3:12]
wt <- as.numeric(ss_run$data_list$wt[which(ss_run$data_list$wt$Wt_index==6)[42],6:15])
sex_ratio <- as.numeric(ss_run$data_list$sex_ratio[1,2:11])
pmature <- as.numeric(ss_run$data_list$pmature[1,2:11])
spawn_month <- ss_run$data_list$spawn_month[1]

natage = rep(1, 10)
for(age in 2:9){
  natage[age] <- natage[age-1] * as.numeric(exp(-M[age-1]))
}

natage[10] <-  natage[9] * as.numeric(exp(-M[9])) / as.numeric(1-exp(-M[10]))


ssb_at_age <- natage * wt * pmature * sex_ratio * exp(-M*spawn_month)
sum(ssb_at_age)

ss_run$quantities$SPR0
