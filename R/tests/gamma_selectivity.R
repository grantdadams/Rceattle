# Predator–prey age selectivity. Eq 18 Kinzey and Punt

gamma_selectivity <- function(alpha = 2, beta = 1){
  g_hat <- (alpha - 1) * beta # Value of logrithm of ratio of predator/prey lengths which predator selectivity is 1
  print(exp(g_hat))

  curve((log(x)/g_hat)^(alpha - 1) * exp(-(log(x) - g_hat) / beta), from = 1, to = exp(g_hat) , xlab = "Ratio of predator length to prey length")
}

gamma_selectivity(3,1)


ksp_type <- c()
rk_sp <- c()

ind = 1
for(p in 1:3){
  for(k in 1:3){
    ksp_type[ind] <- (p - 1) * (4) + k
    rk_sp[ind] <- (p - 1) * 3 + k
    ind <- ind + 1

  }
}
ksp_type
