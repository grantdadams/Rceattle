

# Set up models
library(TMB)
library(RTMB)
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

# Set up dimensions
niter <- 1000
samp_mult <- 200

mat <- matrix(0, samp_mult, niter)
colnames(mat) <- paste0("Iter", 1:niter)
rownames(mat) <- paste0("Samp", 1:samp_mult)

rtmb_make_vec <- tmb_make_vec <- mat
rtmb_fit_vec <- tmb_fit_vec <- mat
rtmb_sdreport_vec <- tmb_sdreport_vec <- mat
tmb_fn_vec <- rtmb_fn_vec <- tmb_gr_vec <- rtmb_gr_vec <- mat

# Run simulation/experiment
for(j in 1:samp_mult){

  # Sample size
  set.seed(1)
  nsamp = 25 * j

  for(i in 1:niter){

    # Simulate data
    data <- list(Y = rnorm(nsamp) * 5 + 2, x=1:nsamp)

    # TMB
    st <- Sys.time()
    obj <- TMB::MakeADFun(data, parameters, DLL="linreg", silent = TRUE)
    obj$hessian <- TRUE
    tmb_make_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    opt <- do.call("optim", obj)
    tmb_fit_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    sd <- sdreport(obj)
    tmb_sdreport_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    fn <- obj$fn()
    tmb_fn_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    gr <- obj$gr()
    tmb_gr_vec[j,i] <- Sys.time() - st

    # RTMB
    st <- Sys.time()
    obj <- RTMB::MakeADFun(f, parameters, silent = TRUE, fast = TRUE)
    obj$hessian <- TRUE
    rtmb_make_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    opt <- do.call("optim", obj)
    rtmb_fit_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    sd <- sdreport(obj)
    rtmb_sdreport_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    fn <- obj$fn()
    rtmb_fn_vec[j,i] <- Sys.time() - st

    st <- Sys.time()
    gr <- obj$gr()
    rtmb_gr_vec[j,i] <- Sys.time() - st
  }
}

# Plot results
library(ggplot2)
library(tidyr)

# - Reformat
results_list <- list(rtmb_make_vec, tmb_make_vec,
                     rtmb_fit_vec, tmb_fit_vec,
                     rtmb_sdreport_vec, tmb_sdreport_vec,
                     rtmb_gr_vec, tmb_gr_vec,
                     rtmb_fn_vec, tmb_fn_vec)
results_list <- lapply(results_list, as.data.frame)
results_list <- lapply(results_list, function(x) { x$N_samples <- 1:samp_mult*25;return(x)})
results_list[c(1,3,5,7,9)] <- lapply(results_list[c(1,3,5,7,9)], function(x) { x$Framework <- "RTMB";return(x)})
results_list[c(2,4,6,8,10)] <- lapply(results_list[c(2,4,6,8,10)], function(x) { x$Framework <- "TMB";return(x)})

results_list[c(1,2)] <- lapply(results_list[c(1,2)], function(x) { x$Action <- "MakeADFun";return(x)})
results_list[c(3,4)] <- lapply(results_list[c(3,4)], function(x) { x$Action <- "Fit";return(x)})
results_list[c(5,6)] <- lapply(results_list[c(5,6)], function(x) { x$Action <- "sdreport";return(x)})
results_list[c(7,8)] <- lapply(results_list[c(7,8)], function(x) { x$Action <- "gr";return(x)})
results_list[c(9,10)] <- lapply(results_list[c(9,10)], function(x) { x$Action <- "fn";return(x)})

results <- do.call("rbind", results_list)

results <- results %>%
  pivot_longer(cols = starts_with("Iter"), names_to = "Iter", values_to = "time")

results_total <- results %>%
  dplyr::filter(Action %in% c("MakeADFun", "Fit", "sdreport")) %>%
  group_by(Framework, N_samples, Iter) %>%
  summarise(time = sum(time)) %>%
  mutate(Action = "Total (Make, Fit, SDreport)") %>%
  select(N_samples, Framework, Action, Iter, time)


results <- rbind(results, results_total) %>%
  mutate(group = factor(paste(Action, Framework)),
         Action = factor(Action, levels = c("Total (Make, Fit, SDreport)", "MakeADFun", "Fit", "sdreport", "gr", "fn"))
  )


# - Plot
results %>%
  ggplot(aes(x = N_samples, y = time, color = Framework, group = group)) +
  # geom_line() +
  stat_smooth(aes(x = N_samples, y = time), method = "lm",
              formula = y ~ poly(x, 21), se = FALSE) +
  ylab("Average time (seconds)") +
  facet_wrap(~Action, scales = "free_y") +
  xlab("Sample size") +
  theme_classic()
