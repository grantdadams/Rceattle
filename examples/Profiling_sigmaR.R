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
rsigma_vec <- seq(from = 0.05, to = 2, by = 0.05)


profile_rsigma <- function(model = NULL, rsigma_vec = NULL, species = NULL){
  ### Set up parallel processing
  library(foreach)
  library(doParallel)

  cores = detectCores() - 6
  registerDoParallel(cores)

  # Loop through Rsigma
  profile_list <- foreach(i = 1:length(rsigma_vec)) %dopar% {
    library(Rceattle)
    library(dplyr)

    # Update sigmaR
    inits <- model$estimated_params
    inits$R_ln_sd[species] <- log(rsigma_vec[i])

    # Build map
    data_list <- model$data_list
    # data_list$estDynamics <- rep(1, data_list$nspp)
    # data_list$estDynamics[species] <- 0
    map <- Rceattle::build_map(data_list, params = inits, debug = FALSE, random_rec = FALSE)

    # Estimate
    mod_prof <- fit_mod(
      data_list = data_list,
      inits = inits,
      map =  map,
      bounds = NULL,
      file = NULL,
      estimateMode = 1,
      HCR = build_hcr(HCR = data_list$HCR,
                      DynamicHCR = data_list$DynamicHCR,
                      Ftarget = data_list$Ftarget,
                      Flimit = data_list$Flimit,
                      Ptarget = data_list$Ptarget,
                      Plimit = data_list$Plimit,
                      Alpha = data_list$Alpha,
                      Pstar = data_list$Pstar,
                      Sigma = data_list$Sigma,
                      Fmult = data_list$Fmult,
                      HCRorder = data_list$HCRorder
      ),
      recFun = build_srr(srr_fun = data_list$srr_fun,
                         srr_pred_fun  = data_list$srr_pred_fun ,
                         proj_mean_rec  = data_list$proj_mean_rec ,
                         srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr), # Update end year if less than srr_mse_switchyr
                         srr_hat_styr = data_list$srr_hat_styr,
                         srr_hat_endyr = data_list$srr_hat_endyr,
                         srr_est_mode  = data_list$srr_est_mode ,
                         srr_prior  = data_list$srr_prior,
                         srr_prior_sd   = data_list$srr_prior_sd,
                         Bmsy_lim = data_list$Bmsy_lim,
                         srr_indices = data_list$srr_indices),
      M1Fun =     build_M1(M1_model = data_list$M1_model,
                           M1_re = data_list$M1_re,
                           updateM1 = FALSE,  # Dont update M1 from data
                           M1_use_prior = data_list$M1_use_prior,
                           M2_use_prior = data_list$M2_use_prior,
                           M_prior = data_list$M_prior,
                           M_prior_sd = data_list$M_prior_sd,
                           M1_indices = data_list$M1_indices),
      growthFun = build_growth(fun = data_list$growth_fun,
                               linkages = data_list$growth_linkages),
      random_rec = data_list$random_rec,
      niter = data_list$niter,
      msmMode = data_list$msmMode,
      avgnMode = data_list$avgnMode,
      suitMode = data_list$suitMode,
      suit_styr = data_list$suit_styr,
      suit_endyr = min(data_list$suit_endyr, data_list$endyr),   # Update to end year if less than suit_endyr
      initMode = data_list$initMode,
      fit_control = fit_control(
        phase = FALSE,
        loopnum = data_list$loopnum,
        getsd = TRUE,
        verbose = 0))
    mod_prof
  }

  closeAllConnections()
  gc()

  return(profile_list)
}


# * Run profile ----
# - EBS
ebs_list1 <- profile_rsigma(model = ebs_run, rsigma_vec, species = 1)
ebs_list2 <- profile_rsigma(model = ebs_run, rsigma_vec, species = 2)
ebs_list3 <- profile_rsigma(model = ebs_run, rsigma_vec, species = 3)

ebs_re_list1 <- profile_rsigma(model = ebs_run_re, rsigma_vec, species = 1)
ebs_re_list2 <- profile_rsigma(model = ebs_run_re, rsigma_vec, species = 2)
ebs_re_list3 <- profile_rsigma(model = ebs_run_re, rsigma_vec, species = 3)

# - EBS w/ Ricker
ebsr_list1 <- profile_rsigma(model = ebs_ricker_run, rsigma_vec, species = 1)
ebsr_list2 <- profile_rsigma(model = ebs_ricker_run, rsigma_vec, species = 2)
ebsr_list3 <- profile_rsigma(model = ebs_ricker_run, rsigma_vec, species = 3)

ebsr_re_list1 <- profile_rsigma(model = ebs_ricker_run_re, rsigma_vec, species = 1)
ebsr_re_list2 <- profile_rsigma(model = ebs_ricker_run_re, rsigma_vec, species = 2)
ebsr_re_list3 <- profile_rsigma(model = ebs_ricker_run_re, rsigma_vec, species = 3)


# - GOA
goa_list1 <- profile_rsigma(model = pollock_model, rsigma_vec, species = 1)
goa_list2 <- profile_rsigma(model = atf_model, rsigma_vec, species = 1)
goa_list3 <- profile_rsigma(model = cod_model, rsigma_vec, species = 1)

goa_re_list1 <- profile_rsigma(model = pollock_model_re, rsigma_vec, species = 1)
goa_re_list2 <- profile_rsigma(model = atf_model_re, rsigma_vec, species = 1)
goa_re_list3 <- profile_rsigma(model = cod_model_re, rsigma_vec, species = 1)

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
  y = sapply(ebs_jnll[[i]], function(x) x$opt$objective)
  y = y-min(y)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = paste("EBS", ebs_run$data_list$spnames[i]), col = "red", ylim = c(0,10))


  y = sapply(ebs_re_jnll[[i]], function(x) x$opt$objective)
  y = y-min(y)
  lines(y = y, x = rsigma_vec, col = 1)


  abline(v = exp(ebs_run_re$estimated_params$ln_rec_sigma[i]), lty = 2)
}

legend("topright", c("Penalized likelihood", "Random effects", "Minima"), col = c(2,1,1), lty = c(1,1,2), bty = "n")


# w/ Ricker
for(i in 1:3){
  y = sapply(ebsr_jnll[[i]], function(x) x$opt$objective)
  y = y-min(y)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = paste("EBS-Ricker", ebs_ricker_run$data_list$spnames[i]), col = "red", ylim = c(0,10))


  y = sapply(ebsr_re_jnll[[i]], function(x) x$opt$objective)
  y = y-min(y)
  lines(y = y, x = rsigma_vec, col = 1)

  abline(v = exp(ebs_ricker_run_re$estimated_params$ln_rec_sigma[i]), lty = 2)
}


# -- GOA
for(i in 1:3){
  y = sapply(goa_jnll[[i]], function(x) x$opt$objective)
  y = y-min(y)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = paste("GOA", goa_list[[i]]$data_list$spnames[1]), col = "red", ylim = c(0,10))


  y = sapply(goa_re_jnll[[i]], function(x) x$opt$objective)
  y = y-min(y)
  lines(y = y, x = rsigma_vec, col = 1)

  abline(v = exp(goa_list[[i]]$estimated_params$ln_rec_sigma[1]), lty = 2)
}
