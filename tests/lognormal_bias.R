

# Simulate data
log_sigma_re <- 2
log_sigma <- 3
n_group <- 20
meanR <- 5

# Simulate random effects
re <- c()
x <- c()
group <- c()
ind <- 1
for(i in 1:n_group){
  re[i] <- rnorm(1, 0, exp(log_sigma_re))
  n_obs <- round(runif(1, min = 5, max = 100))
  # Simulate observations given random effects
  x <- c(x, rlnorm(n_obs, meanR + re[i], exp(log_sigma)))
  group <- c(group, rep(i-1, n_obs))
}

# Setup TMB
library(TMB)
library(TMBhelper)
version <- "lognormal_bias"
cpp_directory <- "tests"
cpp_file <- paste0(cpp_directory, "/", version)
TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)))
params <- list(log_re = rep(0, n_group), log_sigma = 2, log_sigma_re = 2, log_meanR = 50)

# Mod 1 - not in like
datalist <- list(x = x, group = group, n_group = n_group, mod = 1)


# Estimate mean, deviations, and variance as fixed effects
obj = TMB::MakeADFun(datalist, parameters = params,  DLL = version, silent = TRUE, random = "log_re")
opt = TMBhelper::Optimize( obj )
rep <- sdreport(obj)

