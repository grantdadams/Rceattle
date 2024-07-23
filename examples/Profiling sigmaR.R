# Code to profile over sigma-R when rec devs are treated as random effects vs fixed effects
library(Rceattle)


# ESTIMATION ----
# * EBS Combined ----
data("BS2017SS")
BS2017SS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ebs_run <- Rceattle::fit_mod(data_list = BS2017SS,
                             inits = NULL, # Initial parameters = 0
                             file = NULL, # Don't save
                             estimateMode = 0, # Estimate
                             random_rec = FALSE, # No random recruitment
                             msmMode = 0, # Single species mode
                             phase = "default",
                             verbose = 1)

# -- Treat recruitment as random effects
ebs_run_re <- Rceattle::fit_mod(data_list = BS2017SS,
                                inits = ebs_run$estimated_params, # Initial parameters = 0
                                file = NULL, # Don't save
                                estimateMode = 0, # Estimate
                                random_rec = TRUE, # Random recruitment
                                msmMode = 0, # Single species mode
                                phase = NULL,
                                verbose = 1)

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
                     srr_prior_mean = alpha,
                     srr_prior_sd = 0.2),
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # Single species mode
  phase = "default",
  verbose = 1,
  initMode = 2)

ebs_ricker_run_re <- Rceattle::fit_mod(
  data_list = BS2017SS,
  inits = ebs_ricker_run$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 1, # Estimate hindcast only
  M1Fun = build_M1(M1_model = 0,
                   M1_use_prior = FALSE,
                   M2_use_prior = FALSE),
  recFun = build_srr(srr_fun = 0,
                     srr_pred_fun = 3,
                     proj_mean_rec = FALSE,
                     srr_est_mode = 1,
                     srr_prior_mean = alpha,
                     srr_prior_sd = 0.2),
  random_rec = TRUE, # No random recruitment
  msmMode = 0, # Single species mode
  phase = NULL,
  verbose = 1,
  initMode = 2)


# * GOA Combined ----
# data("GOA2018SS")
# GOA2018SS$fleet_control$proj_F_prop <- rep(1, nrow(GOA2018SS$fleet_control))
# goa_run <- Rceattle::fit_mod(data_list = GOA2018SS,
#                              inits = NULL, # Initial parameters = 0
#                              file = NULL, # Don't save
#                              estimateMode = 0, # Estimate
#                              random_rec = FALSE, # No random recruitment
#                              msmMode = 0, # Single species mode
#                              phase = "default",
#                              verbose = 1)
#
# # -- Treat recruitment as random effects
# goa_run_re <- Rceattle::fit_mod(data_list = GOA2018SS,
#                                 inits = goa_run$estimated_params, # Initial parameters = 0
#                                 file = NULL, # Don't save
#                                 estimateMode = 0, # Estimate
#                                 random_rec = TRUE, # Random recruitment
#                                 msmMode = 0, # Single species mode
#                                 phase = NULL,
#                                 getsd = FALSE,
#                                 verbose = 1)

# * GOA Pollock ----
data("GOApollock")
GOApollock$styr = 1977 # The SAFE model starts at 1970, so change styr to 1970 to run the full time series model (data is in there). I start them all at 1977 because thats the years with overlap.
pollock_model <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Use default phasing

pollock_model_re <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = pollock_model$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = TRUE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = NULL) # Use default phasing


# * GOA Arrowtooth flounder ----
data("GOAatf")
GOAatf$styr = 1977 # The SAFE model starts at 1961, so change styr to 1961 to run the full time series model (data is in there). I start them all at 1977 because thats the years with overlap.
atf_model <- Rceattle::fit_mod(
  data_list = GOAatf,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Use default phasing

atf_model_re <- Rceattle::fit_mod(
  data_list = GOAatf,
  inits = atf_model$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = TRUE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = NULL) # Use default phasing


# * GOA Cod ----
data("GOAcod")
GOAcod$pmature[1,2:13] <- 2 # Spawn wt from SS model includes sex-ratio and maturity already, so setting Pmature (age-at-maturity) to 2 to have CEATTLE calculations be the same
cod_model <- Rceattle::fit_mod(
  data_list = GOAcod,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = "default") # Use default phasing

cod_model_re <- Rceattle::fit_mod(
  data_list = GOAcod,
  inits = cod_model$estimated_params, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = TRUE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = NULL) # Use default phasing


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
    inits$ln_rec_sigma[species] <- log(rsigma_vec[i])

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
      HCR = build_hcr(HCR = model$data_list$HCR, # Tier3 HCR
                      DynamicHCR = model$data_list$DynamicHCR,
                      FsprTarget = model$data_list$FsprTarget,
                      FsprLimit = model$data_list$FsprLimit,
                      Ptarget = model$data_list$Ptarget,
                      Plimit = model$data_list$Plimit,
                      Alpha = model$data_list$Alpha,
                      Pstar = model$data_list$Pstar,
                      Sigma = model$data_list$Sigma,
                      Fmult = model$data_list$Fmult,
                      HCRorder = model$data_list$HCRorder
      ),
      recFun = build_srr(srr_fun = model$data_list$srr_fun,
                         srr_pred_fun  = model$data_list$srr_pred_fun ,
                         proj_mean_rec  = model$data_list$proj_mean_rec ,
                         srr_meanyr = model$data_list$srr_meanyr,
                         R_hat_yr = model$data_list$R_hat_yr,
                         srr_est_mode  = model$data_list$srr_est_mode ,
                         srr_prior_mean  = model$data_list$srr_prior_mean,
                         srr_prior_sd   = model$data_list$srr_prior_sd,
                         Bmsy_lim = model$data_list$Bmsy_lim,
                         srr_env_indices = model$data_list$srr_env_indices),
      M1Fun = build_M1(M1_model= model$data_list$M1_model,
                       updateM1 = FALSE,
                       M1_use_prior = model$data_list$M1_use_prior,
                       M2_use_prior = model$data_list$M2_use_prior,
                       M1_prior_mean = model$data_list$M1_prior_mean,
                       M1_prior_sd = model$data_list$M1_prior_sd),
      random_rec = model$data_list$random_rec,
      niter = model$data_list$niter,
      msmMode = model$data_list$msmMode,
      avgnMode = model$data_list$avgnMode,
      suitMode = model$data_list$suitMode,
      suit_meanyr = model$data_list$suit_meanyr,
      initMode = model$data_list$initMode,
      phase = NULL,
      loopnum = 1,
      getsd = FALSE,
      verbose = 0)

    # Get JNLL by species
    # jnll_data <- data.frame(Species = mod_prof$data_list$fleet_control$Species, JNLL = colSums(mod_prof$quantities$jnll_comp[c(1:9, 13),]))
    # jnll_penalties <- data.frame(Species = 1:3,  JNLL = colSums(mod_prof$quantities$jnll_comp[10:12,1:3] ))
    # mod_prof$jnll <- rbind(jnll_data, jnll_penalties) %>%
    #   group_by(Species) %>%
    #   summarise(jnll = sum(JNLL)) %>%
    #   arrange(Species) %>%
    #   mutate(Species_name = mod_prof$data_list$spnames,
    #          sigmaR = rsigma_vec[i])
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
