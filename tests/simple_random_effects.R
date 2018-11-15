

# Simulate data
log_sigma_re <- 2
log_sigma <- 3
n_group <- 20
meanR <- 100

# Simulate random effects
re <- c()
x <- c()
group <- c()
ind <- 1
for(i in 1:n_group){
  re[i] <- rnorm(1, 0, exp(log_sigma_re))
  n_obs <- round(runif(1, min = 5, max = 100))
  # Simulate observations given random effects
  x <- c(x, rnorm(n_obs, meanR + re[i], exp(log_sigma)))
  group <- c(group, rep(i, n_obs))
}

# Plot function
plot_fun <- function(obs, est, CI, x, obs_var, est_var, var_CI, file_name = NULL, mean_var = NULL){

  tiff( file = paste0("tests/",file_name,".tiff") , width=169 / 25.4, height = 100 / 25.4, family = "Helvetica", units = "in", res = 100)
  par( mar=c(3, 3 , 1 , 0.25) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
  layout(matrix(1:2, ncol = 2), widths = c(3,1))

  # Plot recruitment
  plot(NA, NA, xlab = "Year", ylab = "Recruitment", xlim = c(0, n_group), ylim = c(75, 125))
  points(y = est, x = x, pch = 16, cex = 2, col = "#5ab4ac")
  arrows(x, est - CI, x, est + CI, length=0.05, angle=90, code=3, col = "#5ab4ac", lwd = 3)
  points(y = obs, x = x + 0.1, pch = 16, cex = 2, col = "#d8b365")
  legend("topleft", c("True value", "Estimate"), pch = c(16, 16), col = c("#d8b365", "#5ab4ac"), bty = "n")

  if(!is.null(mean_var)){
    abline(h = mean_var, lwd = 3, col = "#d8b365")
  }

  # Plot variance
  min_y <- min( c(0, est_var - var_CI, est_var + var_CI, obs_var, est_var), na.rm = T)
  max_y <- max( c(0, est_var - var_CI, est_var + var_CI, obs_var, est_var), na.rm = T) + .75
  plot(NA, NA, xlab = NA, ylab = "Variance in recruitment", ylim = c(min_y, max_y), xlim = c(0.35, 0.75), xaxt = "n")
  points(y = est_var, x = 0.5, pch = 16, cex = 2, col = "#5ab4ac")
  arrows(0.5,  est_var - var_CI, 0.5,  est_var + var_CI, length=0.05, angle=90, code=3, col = "#5ab4ac", lwd = 3)
  points(y = obs_var, x = 0.5 + 0.1, pch = 16, cex = 2, col = "#d8b365")
  dev.off()
}


# Setup TMB
library(TMB)
library(TMBhelper)
version <- "simple_random_effects"
cpp_directory <- "tests"
cpp_file <- paste0(cpp_directory, "/", version)
TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)))
params <- list(re = rep(0, n_group), log_sigma = 2, log_sigma_re = 2, meanR = 50)
datalist <- list(x = x, group = group, n_group = n_group)


plot_fun(obs = meanR + re, est= rep(NA, 20), CI = NA, x = c(1:n_group), obs_var = exp(log_sigma_re), est_var = NA, var_CI = NA, file_name = "Annual coefs", mean_var = meanR)


# Strategy 1 - Estimate annual bits
mod1 <- lm(x ~ as.factor(group) - 1)
sum1 <- summary(mod1)
plot_fun(obs = meanR + re, est = mod1$coefficients, CI = sum1$coefficients[,2] * 1.92, x = c(1:n_group), obs_var = exp(log_sigma_re), est_var = sd(mod1$coefficients), var_CI = NA, file_name = "Strategy 1 - Annual coefficients")


# Strategy 2 - Estimate mean, deviations, and variance as fixed effects
obj = TMB::MakeADFun(datalist, parameters = params,  DLL = version, silent = TRUE)
opt = TMBhelper::Optimize( obj )
rep <- sdreport(obj)
estR <- obj$report()$rec
CI <- rep$sd[1:length(estR)] * 1.92
var_CI <- rep$sd[length(estR)+1] * 1.92
plot_fun(obs = meanR + re, est = estR, CI = CI, x = c(1:n_group), obs_var = exp(log_sigma_re), est_var = obj$report()$sigma_re, var_CI = NA, file_name = "Strategy 3 - Estimate all as fixed effects")


# Likelihood profile
log_sigma_re_vec <- log(seq(0.0001, 20, by = 0.01))
jnll <- c()
for(i in 1:length(log_sigma_re_vec)){
  params$log_sigma_re <- log_sigma_re_vec[i]
  map <- list(log_sigma_re = as.factor(NA))
  obj = TMB::MakeADFun(datalist, parameters = params,  DLL = version, map = map, silent = TRUE)
  opt = TMBhelper::Optimize( obj )
  jnll[i] <- opt$objective
}
jnll <- jnll - min(jnll)

tiff( file = paste0("tests/like_profile_pen_likelihood",".tiff") , width=169 / 25.4, height = 100 / 25.4, family = "Helvetica", units = "in", res = 100)
par( mar=c(3, 3 , 1 , 0.25) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
plot(x = exp(log_sigma_re_vec), y = jnll, xlab = "Sigma rec", ylab = "NLL", type = "l", lwd = 3)
abline(v = exp(log_sigma_re), lwd = 3, col = "#d8b365" , lty = 2)
abline(v = exp(log_sigma_re_vec)[which(jnll ==  min(jnll))], lwd = 3, col = "#5ab4ac", lty = 2)
legend("topright", c("True value", "Estimate"), lty = c(1, 2), col = c("#d8b365", "#5ab4ac"), lwd = rep(2, 2), bty = "n")
dev.off()


# Strategy 4 - Use penalized likelihood
tiff( file = paste0("tests/like_profile_pen_likelihood2",".tiff") , width=169 / 25.4, height = 100 / 25.4, family = "Helvetica", units = "in", res = 100)
par( mar=c(3, 3 , 1 , 0.25) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
plot(x = exp(log_sigma_re_vec), y = jnll, xlab = "Variance in recruitment", ylab = "Negative log-likelihood", type = "l", lwd = 3)
abline(v = exp(log_sigma_re), lwd = 3, col = "#d8b365" , lty = 2)
abline(v = exp(log_sigma_re_vec)[which(jnll ==  min(jnll))], lwd = 3, col = "#5ab4ac", lty = 2)
local_minima <- exp(log_sigma_re_vec)[which(jnll ==  min(jnll[ which(exp(log_sigma_re_vec) > 2)]) )]
abline(v = local_minima, lwd = 3, col = "#5ab4ac", lty = 3)
legend("bottomright", c("True value", "Estimate", "Local minima"), lty = c(1, 2, 3), col = c("#d8b365", "#5ab4ac", "#5ab4ac"), lwd = rep(2, 3), bty = "n")
dev.off()


params$log_sigma_re <- log(local_minima)
map <- list(log_sigma_re = as.factor(NA))
obj = TMB::MakeADFun(datalist, parameters = params,  DLL = version, map = map, silent = TRUE)
opt = TMBhelper::Optimize( obj )
rep <- sdreport(obj)
estR <- obj$report()$rec
CI <- rep$sd[1:length(estR)] * 1.92
plot_fun(obs = meanR + re, est = estR, CI = NA, x = c(1:n_group), obs_var = exp(log_sigma_re), est_var = local_minima, var_CI = NA, file_name = "Strategy 4 - Fix variance")


# Strategy 5 - Random effects
random <- "re"
obj = TMB::MakeADFun(datalist, parameters = params,  DLL = version, random = random, silent = TRUE)
opt = TMBhelper::Optimize( obj )
rep <- sdreport(obj)
estR <- obj$report()$rec
CI <- rep$sd[1:length(estR)] * 1.92

sigma_re <- rep$value[which(names(rep$value) == "sigma_re")]
sigma_CI <- rep$sd[which(names(rep$value) == "sigma_re")] * 1.92
plot_fun(obs = meanR + re, est = estR, CI = CI, x = c(1:n_group), obs_var = exp(log_sigma_re), est_var =  sigma_re, var_CI = sigma_CI, file_name = "Strategy 5 - Random effects")



