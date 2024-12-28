


library(TMB)
compile("linreg.cpp")
dyn.load(dynlib("linreg"))

f <- function(parms) {
  Y <- data$Y
  x <- data$x
  a <- parms$a
  b <- parms$b
  logSigma <- parms$logSigma
  ADREPORT(exp(2*logSigma));
  nll = -sum(dnorm(Y, a+b*x, exp(logSigma), TRUE))
  nll
}

parameters <- list(a=0, b=0, logSigma=0)
library(RTMB)

niter <- 1000
samp_mult <- 200
rtmb_time_vec <- time_vec <- matrix(0, samp_mult, niter)


for(j in 1:samp_mult){

  set.seed(1)
  nsamp = 25 * j

  for(i in 1:niter){

    # Sim
    data <- list(Y = rnorm(nsamp) * 5 + 2, x=1:nsamp)

    # TMB
    st <- Sys.time()
    obj <- TMB::MakeADFun(data, parameters, DLL="linreg", silent = TRUE)
    obj$hessian <- TRUE
    opt <- do.call("optim", obj)
    sd <- sdreport(obj)
    time_vec[j,i] <- Sys.time() - st

    # RTMB
    #print(i)
    st <- Sys.time()
    obj <- RTMB::MakeADFun(f, parameters, silent = TRUE, fast = TRUE)

    obj$hessian <- TRUE
    opt <- do.call("optim", obj)

    rep <- sdreport(obj)
    rtmb_time_vec[j,i] <- Sys.time() - st
  }
}

library(ggplot2)
library(tidyr)

time_vec <- as.data.frame(time_vec)
rtmb_time_vec <- as.data.frame(rtmb_time_vec)
colnames(time_vec) <- colnames(rtmb_time_vec) <- paste0("Iter", 1:niter)
rownames(time_vec) <- rownames(rtmb_time_vec) <- paste0("Samp", 1:samp_mult)

time_vec$N_samples <- rtmb_time_vec$N_samples <- 1:samp_mult*25

time_long <- time_vec %>%
  pivot_longer(cols = starts_with("Iter"), names_to = "Iter", values_to = "time")

rtmb_time_long <- rtmb_time_vec %>%
  pivot_longer(cols = starts_with("Iter"), names_to = "Iter", values_to = "time")

time_long$Framework <- "TMB"
rtmb_time_long$Framework <- "RTMB"

rbind(time_long, rtmb_time_long) %>%
  ggplot(aes(x = N_samples, y = time, group = Framework, color = Framework)) +
  # geom_line() +
  stat_smooth(aes(x = N_samples, y = time), method = "lm",
              formula = y ~ poly(x, 21), se = FALSE) +
  ylab("Average time (seconds)") +
  xlab("Sample size") +
  annotate("text", x = 1500, y = 0.03, label = "Y~mx+b over 1,000 simulations") +
  theme_classic()
