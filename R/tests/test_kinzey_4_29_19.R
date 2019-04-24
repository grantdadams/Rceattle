library(Rceattle)
data("BS2017MS")
data("BS2017SS")


ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            silent = TRUE)


ms_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = ss_run$estimated_params, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = 1, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 2, # Single species mode
                            niter = 2,
                            silent = TRUE)
ms_run$quantities$jnll_comp
sum(is.na(ms_run$quantities$Q_ksum_u))

  ms_run$quantities$M2[,,1]


comp <- c(0,0,0)
for (rsp in 1:3) {
  for (r_age in 1:ms_run$data_list$nages[rsp]) {
    for (ksp in 1:3) {
      for (k_age in 1:ms_run$data_list$nages[ksp]) {
        # if (UobsWtAge(rsp, ksp, r_age, k_age) > 0) { // (rsp, ksp, a, r_ln, yr)
          comp[rsp] = comp[rsp] - 20 * (ms_run$quantities$mn_UobsWtAge[rsp, ksp, r_age, k_age] + 0.001) * (log(ms_run$quantities$mn_UobsWtAge_hat[rsp, ksp, r_age, k_age] + 0.001))
        # }
      }
    }
  }
}

