# Set up models ----
library(TMB)
library(RTMB)
library(microbenchmark)
library(ggplot2)
library(tidyr)
library(dplyr)

# TMB ----
# https://github.com/kaskr/adcomp/blob/master/tmb_examples/linreg.cpp
linreg <- "
// Simple linear regression.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));
  Type nll = -sum(dnorm(Y, a+b*x, exp(logSigma), true));
  return nll;
}
"

write(linreg, file = "linreg.cpp")
compile("linreg.cpp", framework = "TMBad")
dyn.load(dynlib("linreg"))

# RTMB ----
# https://github.com/kaskr/RTMB/blob/master/tmb_examples/linreg.R
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


# Parameters ----
parameters <- list(a=0, b=0, logSigma=0)


# Set up simulation dimensions ----
times <- 100        # number of microbenchmark iterations
samp_mult <- 200    # range of sample size 1:200 * 25

mat <- rep(0, samp_mult)
names(mat) <- paste0("Samp", 1:samp_mult)

rtmb_make_vec <- tmb_make_vec <- mat         # make time
rtmb_fit_vec <- tmb_fit_vec <- mat           # fit time
rtmb_sdreport_vec <- tmb_sdreport_vec <- mat # sdreport time
tmb_fn_vec <- rtmb_fn_vec <- mat             # fn() time
tmb_gr_vec <- rtmb_gr_vec <- mat             # gr() time


# Run simulation/experiment ----
for(j in 1:samp_mult){

  # Sample size
  set.seed(1)
  nsamp = 25 * j

  # Simulate data
  data <- list(Y = rnorm(nsamp) * 5 + 2, x=1:nsamp)

  # TMB
  tmb_make_vec[j] <-  mean(microbenchmark(
    obj <- TMB::MakeADFun(data, parameters, DLL="linreg", silent = TRUE)
    , times = times)$time) * 1e-6

  tmb_fit_vec[j] <- mean(microbenchmark(
    opt <- do.call("optim", obj)
    , times = times)$time) * 1e-6

  tmb_sdreport_vec[j] <- mean(microbenchmark(
    sd <- TMB::sdreport(obj)
    , times = times)$time) * 1e-6

  tmb_fn_vec[j] <- mean(microbenchmark(
    fn <- obj$fn()
    , times = times)$time) * 1e-6

  tmb_gr_vec[j] <- mean(microbenchmark(
    gr <- obj$gr()
    , times = times)$time) * 1e-6

  # RTMB
  rtmb_make_vec[j] <- mean(microbenchmark(
    obj <- RTMB::MakeADFun(f, parameters, silent = TRUE, fast = TRUE)
    , times = times)$time) * 1e-6

  rtmb_fit_vec[j] <- mean(microbenchmark(
    opt <- do.call("optim", obj)
    , times = times)$time) * 1e-6

  rtmb_sdreport_vec[j] <- mean(microbenchmark(
    sd <- RTMB::sdreport(obj)
    , times = times)$time) * 1e-6

  rtmb_fn_vec[j] <- mean(microbenchmark(
    fn <- obj$fn()
    , times = times)$time) * 1e-6

  rtmb_gr_vec[j] <- mean(microbenchmark(
    gr <- obj$gr()
    , times = times)$time) * 1e-6

}

# Plot results ----

# - Reformat
rtmb_results <-
  cbind(rtmb_make_vec,
        rtmb_fit_vec,
        rtmb_sdreport_vec,
        rtmb_gr_vec,
        rtmb_fn_vec)
rtmb_results <- as.data.frame(cbind(rtmb_results, rowSums(rtmb_results))) # Total
colnames(rtmb_results) <- c("MakeADFun", "Fit", "sdreport", "gr", "fn", "Total (Make, Fit, SDreport)")
rtmb_results$N_samples <- 1:samp_mult*25

rtmb_results <- rtmb_results %>%
  pivot_longer(!N_samples, names_to = "Action", values_to = "Seconds") %>%
 mutate(Framework = "RTMB")

tmb_results <-
  cbind(tmb_make_vec,
        tmb_fit_vec,
        tmb_sdreport_vec,
        tmb_gr_vec,
        tmb_fn_vec)
tmb_results <- as.data.frame(cbind(tmb_results, rowSums(tmb_results))) # Total
colnames(tmb_results) <- c("MakeADFun", "Fit", "sdreport", "gr", "fn", "Total (Make, Fit, SDreport)")
tmb_results$N_samples <- 1:samp_mult*25

tmb_results <- tmb_results %>%
  pivot_longer(!N_samples, names_to = "Action", values_to = "Seconds") %>%
  mutate(Framework = "TMB: TMBad")


results <- rbind(tmb_results, rtmb_results) %>%
  mutate(group = factor(paste(Action, Framework)),
         Action = factor(Action, levels = c("Total (Make, Fit, SDreport)", "MakeADFun", "Fit", "sdreport", "gr", "fn"))
  )


# - Plot
results %>%
  ggplot(aes(x = N_samples, y = time, color = Framework, group = group)) +
  # geom_line() +
  stat_smooth(aes(x = N_samples, y = Seconds), method = "lm",
              formula = y ~ poly(x, 21), se = FALSE) +
  ylab("Average time (seconds)") +
  facet_wrap(~Action, scales = "free_y") +
  xlab("Sample size") +
  theme_classic()
