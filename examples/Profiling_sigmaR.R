# Code to profile over sigma-R when rec devs are treated as random effects vs fixed effects
library(Rceattle)


# ESTIMATION ----
# * EBS Combined ----
data("BS2017SS")
ebs_run <- Rceattle::fit_mod(data_list = BS2017SS,
                             inits = NULL, # Initial parameters = 0
                             file = NULL, # Don't save
                             estimateMode = 0, # Estimate
                             random_rec = FALSE, # No random recruitment
                             msmMode = 0, # Single species mode
                             fit_control = fit_control(
                               phase = TRUE,
                               verbose = 1))

# -- Treat recruitment as random effects
ebs_run_re <- Rceattle::fit_mod(data_list = BS2017SS,
                                inits = ebs_run$estimated_params, # Initial parameters from previous
                                file = NULL, # Don't save
                                estimateMode = 0, # Estimate
                                random_rec = TRUE, # Random recruitment
                                msmMode = 0, # Single species mode
                                fit_control = fit_control(
                                  phase = FALSE,
                                  verbose = 1))

# * EBS with Ricker ----
alpha = exp(c(4.121, 2.119, 1.553))
ebs_ricker_run <- Rceattle::fit_mod(
  data_list = BS2017SS,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  M1Fun = build_M1(M1_model = 0,
                   M1_use_prior = FALSE,
                   M2_use_prior = FALSE),
  recFun = build_srr(srr_fun = 0,
                     srr_pred_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior = alpha,
                     srr_prior_sd = 0.2),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  initMode = "NonEquilibrium", # Start at fished equilibrium (biases alpha and beta otherwise)
  fit_control = fit_control(phase = TRUE, verbose = 1))

# -- Treat recruitment as random effects
ebs_ricker_run_re <- Rceattle::fit_mod(
  data_list = BS2017SS,
  inits = ebs_ricker_run$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  M1Fun = build_M1(M1_model = 0,
                   M1_use_prior = FALSE,
                   M2_use_prior = FALSE),
  recFun = build_srr(srr_fun = 0,
                     srr_pred_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior = alpha,
                     srr_prior_sd = 0.2),
  random_rec = TRUE, # Random recruitment
  msmMode = 0, # Single species mode
  initMode = "NonEquilibrium",
  fit_control = fit_control(phase = FALSE, verbose = 1))


# * GOA Pollock ----
data("GOApollock")
GOApollock$styr = 1977 # The SAFE model starts at 1970, so change styr to 1970 to run the full time series model (data is in there). I start them all at 1977 because thats the years with overlap.
pollock_model <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = FALSE, # No random recruitment
  msmMode = 0,
  fit_control = fit_control(
    verbose = 1, # Silence optimization output
    phase = TRUE)) # Use default phasing

pollock_model_re <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = pollock_model$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = TRUE, # Random recruitment
  msmMode = 0,
  fit_control = fit_control(
    verbose = 1, # Silence optimization output
    phase = FALSE))


# * GOA Arrowtooth flounder ----
data("GOAatf")
GOAatf$styr = 1977 # The SAFE model starts at 1961, so change styr to 1961 to run the full time series model (data is in there). I start them all at 1977 because thats the years with overlap.
atf_model <- Rceattle::fit_mod(
  data_list = GOAatf,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = FALSE, # No random recruitment
  msmMode = 0,
  fit_control = fit_control(
    verbose = 1, # Silence optimization output
    phase = TRUE)) # Use default phasing

atf_model_re <- Rceattle::fit_mod(
  data_list = GOAatf,
  inits = atf_model$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = TRUE, # Random recruitment
  msmMode = 0,
  fit_control = fit_control(
    verbose = 1, # Silence optimization output
    phase = FALSE))


# * GOA Cod ----
data("GOAcod")
GOAcod$pmature[1,2:13] <- 2 # Spawn wt from SS model includes sex-ratio and maturity already, so setting Pmature (age-at-maturity) to 2 to have CEATTLE calculations be the same
cod_model <- Rceattle::fit_mod(
  data_list = GOAcod,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = FALSE, # No random recruitment
  msmMode = 0,
  fit_control = fit_control(
    verbose = 1, # Silence optimization output
    phase = TRUE)) # Use default phasing

cod_model_re <- Rceattle::fit_mod(
  data_list = GOAcod,
  inits = cod_model$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec =TRUE, # Random recruitment
  msmMode = 0,
  fit_control = fit_control(
    verbose = 1, # Silence optimization output
    phase = FALSE))


# PROFILE ----
# All profiles below use Rceattle::profile_param(), which generalises the
# original profile_rsigma() helper to any parameter slot and supports
# arbitrary N-D cross-profiles. The legacy helper is kept below the new
# calls for reference / reproducibility of older runs.
rsigma_vec <- seq(from = 0.05, to = 2, by = 0.05)

prof_sigmaR <- function(model, species){
  # Natural-scale alias: pass sigmaR values directly (no manual log()).
  profile_param(
    Rceattle = model,
    param    = "sigmaR",
    slots    = list(species),
    values   = list(rsigma_vec)
  )
}


# * Run profile ----
# - EBS
ebs_list1 <- prof_sigmaR(ebs_run, species = 1)
ebs_list2 <- prof_sigmaR(ebs_run, species = 2)
ebs_list3 <- prof_sigmaR(ebs_run, species = 3)

ebs_re_list1 <- prof_sigmaR(ebs_run_re, species = 1)
ebs_re_list2 <- prof_sigmaR(ebs_run_re, species = 2)
ebs_re_list3 <- prof_sigmaR(ebs_run_re, species = 3)

# - EBS w/ Ricker
ebsr_list1 <- prof_sigmaR(ebs_ricker_run, species = 1)
ebsr_list2 <- prof_sigmaR(ebs_ricker_run, species = 2)
ebsr_list3 <- prof_sigmaR(ebs_ricker_run, species = 3)

ebsr_re_list1 <- prof_sigmaR(ebs_ricker_run_re, species = 1)
ebsr_re_list2 <- prof_sigmaR(ebs_ricker_run_re, species = 2)
ebsr_re_list3 <- prof_sigmaR(ebs_ricker_run_re, species = 3)


# - GOA
goa_list1 <- prof_sigmaR(pollock_model, species = 1)
goa_list2 <- prof_sigmaR(atf_model,     species = 1)
goa_list3 <- prof_sigmaR(cod_model,     species = 1)

goa_re_list1 <- prof_sigmaR(pollock_model_re, species = 1)
goa_re_list2 <- prof_sigmaR(atf_model_re,     species = 1)
goa_re_list3 <- prof_sigmaR(cod_model_re,     species = 1)

goa_list <- list(pollock_model_re, cod_model_re, atf_model_re)

# * Combine ----
ebs_jnll <- list(ebs_list1, ebs_list2, ebs_list3)
ebs_re_jnll <- list(ebs_re_list1, ebs_re_list2, ebs_re_list3)

ebsr_jnll <- list(ebsr_list1, ebsr_list2, ebsr_list3)
ebsr_re_jnll <- list(ebsr_re_list1, ebsr_re_list2, ebsr_re_list3)

goa_jnll <- list(goa_list1, goa_list3, goa_list2)
goa_re_jnll <- list(goa_re_list1, goa_re_list3, goa_re_list2)


# PLOT ----
par(mfrow = c(3,3))

# -- EBS
for(i in 1:3){
  y = ebs_jnll[[i]]$nll
  y = y - min(y, na.rm = TRUE)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = paste("EBS", ebs_run$data_list$spnames[i]), col = "red", ylim = c(0,10))


  y = ebs_re_jnll[[i]]$nll
  y = y - min(y, na.rm = TRUE)
  lines(y = y, x = rsigma_vec, col = 1)


  abline(v = exp(ebs_run_re$estimated_params$R_log_sd[i]), lty = 2)
}

legend("topright", c("Penalized likelihood", "Random effects", "Minima"), col = c(2,1,1), lty = c(1,1,2), bty = "n")


# w/ Ricker
for(i in 1:3){
  y = ebsr_jnll[[i]]$nll
  y = y - min(y, na.rm = TRUE)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = paste("EBS-Ricker", ebs_ricker_run$data_list$spnames[i]), col = "red", ylim = c(0,10))


  y = ebsr_re_jnll[[i]]$nll
  y = y - min(y, na.rm = TRUE)
  lines(y = y, x = rsigma_vec, col = 1)

  abline(v = exp(ebs_ricker_run_re$estimated_params$R_log_sd[i]), lty = 2)
}


# -- GOA
for(i in 1:3){
  y = goa_jnll[[i]]$nll
  y = y - min(y, na.rm = TRUE)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = paste("GOA", goa_list[[i]]$data_list$spnames[1]), col = "red", ylim = c(0,10))


  y = goa_re_jnll[[i]]$nll
  y = y - min(y, na.rm = TRUE)
  lines(y = y, x = rsigma_vec, col = 1)

  abline(v = exp(goa_list[[i]]$estimated_params$R_log_sd[1]), lty = 2)
}


# PROFILE OTHER PARAMETERS ----
# profile_param() generalises to any parameter slot. The "natural-scale"
# aliases (sigmaR, M1, R0, alpha, beta) take values in natural units and
# log() them internally. For rec_pars aliases the column is inferred from
# the alias name, so `slots` only needs the species index.

# 1-D profile of SRR alpha for species 1 of the EBS-Ricker run
alpha_prof <- profile_param(
  Rceattle = ebs_ricker_run,
  param    = "alpha",
  slots    = list(1),
  values   = list(seq(20, 100, length.out = 20))
)

# 2-D cross-profile: M1 across sex for species 1 of the GOA cod run
# (males = sex 1, females = sex 2; profiled at age 1)
M_sex_prof <- profile_param(
  Rceattle = cod_model,
  param    = "M1",
  slots    = list(c(1, 1, 1), c(1, 2, 1)),
  values   = list(seq(0.20, 0.60, by = 0.05),
                  seq(0.20, 0.60, by = 0.05))
)

# 2-D cross-profile of sigmaR across species 1 and 2 in the EBS run
sigmaR_cross <- profile_param(
  Rceattle = ebs_run,
  param    = "sigmaR",
  slots    = list(1, 2),
  values   = list(seq(0.2, 1.5, by = 0.1),
                  seq(0.2, 1.5, by = 0.1))
)

# Equivalent raw-form call (operates directly on log-scale R_log_sd)
sigmaR_cross_raw <- profile_param(
  Rceattle  = ebs_run,
  param     = "R_log_sd",
  slots     = list(1, 2),
  values    = list(log(seq(0.2, 1.5, by = 0.1)),
                   log(seq(0.2, 1.5, by = 0.1))),
  transform = "identity"
)
