
tc_hat <- tmp$tc_hat_1[1,]
fsh_age_hat <- tmp$fsh_age_hat_1

ch_atm_prod <- sweep(fsh_age_hat,MARGIN=1,tc_hat,"*")

age_trans_matrix <- as.matrix(read.csv("data/age_trans_matrix_1.csv"))

# fsh_age_hat(sp,yr)  = catch_hat(sp,yr)*age_trans_matrix(sp) / tc_hat(sp,yr);

 age_trans_matrix %*% ch_atm_prod


catch_hat <- runif(12, 100, 1000)

age_trans_matrix <- as.matrix(read.csv("data/age_trans_matrix_2.csv", header = F))

catch_hat * age_trans_matrix
