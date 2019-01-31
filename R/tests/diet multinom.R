nspp <- BS2017MS$nspp
nages <- BS2017MS$nages
UobsWtAge  <- BS2017MS$UobsWtAge
mn_UobsWtAge_hat <-  ms_run_diet$quantities$mn_UobsWtAge_hat

jnll_comp <- c()
for (rsp in 1: nspp) {
  jnll_comp[rsp] = 0;
  for (r_age =in 1:nages[rsp]) {
    for (ksp in 1: nspp) {
      for (k_age in 1:nages[ksp]) {
        if (UobsWtAge[rsp, ksp, r_age, k_age] > 0) {
          jnll_comp[rsp] = jnll_comp[rsp] - 20 * (UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst) * (log(UobsWtAge_hat(rsp, ksp, r_age, k_age, yr) + MNConst))
        }
      }
    }
  }
  jnll_comp[rsp] = jnll_comp[rsp] - ms_run_diet$quantities$offset_uobsagewt[rsp]
}
