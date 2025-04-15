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
                             phase = TRUE,
                             verbose = 1)

# -- Treat recruitment as random effects
ebs_run_re <- Rceattle::fit_mod(data_list = BS2017SS,
                                inits = ebs_run$estimated_params, # Initial parameters from previous
                                file = NULL, # Don't save
                                estimateMode = 0, # Estimate
                                random_rec = TRUE, # Random recruitment
                                msmMode = 0, # Single species mode
                                phase = FALSE,
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
  phase = TRUE,
  verbose = 1,
  initMode = 2) # Start at fished equilibrium (biases alpha and beta otherwise)

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
                     srr_prior_mean = alpha,
                     srr_prior_sd = 0.2),
  random_rec = TRUE, # Random recruitment
  msmMode = 0, # Single species mode
  phase = FALSE,
  verbose = 1,
  initMode = 2)


# * EBS Yellowfin sole ----
mydata_yfs <- Rceattle::read_data( file = "C:/Users/grant.adams/GitHub/yfs_ss3/Rceattle runs/Data/yfs_single_species_2022.xlsx")
mydata_yfs$estDynamics = 0
mydata_yfs$srv_biom$Log_sd <- mydata_yfs$srv_biom$Log_sd/mydata_yfs$srv_biom$Observation

mydata_yfs$fsh_biom$Catch <- mydata_yfs$fsh_biom$Catch*1000


yfs_model <- Rceattle::fit_mod(data_list = mydata_yfs,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 0, # Estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 0, # Single species mode
                               verbose = 1,
                               phase = TRUE,
                               initMode = 2)

yfs_model_re <- Rceattle::fit_mod(data_list = mydata_yfs,
                                  inits = yfs_model$estimated_params, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = TRUE, # Random recruitment
                                  msmMode = 0, # Single species mode
                                  verbose = 1,
                                  phase = FALSE,
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
#                              phase = TRUE,
#                              verbose = 1)
#
# # -- Treat recruitment as random effects
# goa_run_re <- Rceattle::fit_mod(data_list = GOA2018SS,
#                                 inits = goa_run$estimated_params, # Initial parameters from previous
#                                 file = NULL, # Don't save
#                                 estimateMode = 0, # Estimate
#                                 random_rec = TRUE, # Random recruitment
#                                 msmMode = 0, # Single species mode
#                                 phase = FALSE,
#                                 getsd = FALSE,
#                                 verbose = 1)

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
  verbose = 1, # Silence optimization output
  phase = TRUE) # Use default phasing

pollock_model_re <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = pollock_model$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = TRUE, # Random recruitment
  msmMode = 0,
  verbose = 1, # Silence optimization output
  phase = FALSE)


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
  verbose = 1, # Silence optimization output
  phase = TRUE) # Use default phasing

atf_model_re <- Rceattle::fit_mod(
  data_list = GOAatf,
  inits = atf_model$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec = TRUE, # Random recruitment
  msmMode = 0,
  verbose = 1, # Silence optimization output
  phase = FALSE)


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
  verbose = 1, # Silence optimization output
  phase = TRUE) # Use default phasing

cod_model_re <- Rceattle::fit_mod(
  data_list = GOAcod,
  inits = cod_model$estimated_params, # Initial parameters from previous
  file = NULL, # Don't save
  estimateMode = 0,
  random_rec =TRUE, # Random recruitment
  msmMode = 0,
  verbose = 1, # Silence optimization output
  phase = FALSE)


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
      phase = FALSE,
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


# EXPERIMENT ----
fix_sigmaR <- function(model = NULL, fix_sigmaR = TRUE, bias.correct = FALSE){

  # Build map (used to fix sigmaR)
  inits <- model$estimated_params
  data_list <- model$data_list
  map <- Rceattle::build_map(data_list, params = inits, debug = FALSE, random_rec = !fix_sigmaR)

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
    bias.correct = bias.correct,
    phase = FALSE,
    loopnum = 3,
    getsd = TRUE,
    verbose = 0)

  return(mod_prof)
}



# * EBS ----
# -- No stock recruit
ebs_run_est <- fix_sigmaR(ebs_run, fix_sigmaR = FALSE)
ebs_run_fixed <- ebs_run
ebs_run_fixed$estimated_params$ln_rec_sigma <- ebs_run_est$estimated_params$ln_rec_sigma
ebs_run_fixed <- fix_sigmaR(ebs_run_fixed, fix_sigmaR = TRUE)

ebs_run_re_fixed <- fix_sigmaR(ebs_run_re, fix_sigmaR = TRUE, bias.correct = FALSE)
ebs_run_re_est <- fix_sigmaR(ebs_run_re, fix_sigmaR = FALSE, bias.correct = FALSE)
ebs_run_re_fixedBC <- fix_sigmaR(ebs_run_re, fix_sigmaR = TRUE, bias.correct = TRUE)
ebs_run_re_estBC <- fix_sigmaR(ebs_run_re, fix_sigmaR = FALSE, bias.correct = TRUE)

# -- Ricker
ebs_ricker_run_est <- fix_sigmaR(ebs_ricker_run, fix_sigmaR = FALSE)
ebs_ricker_run_fixed <- ebs_ricker_run
ebs_ricker_run_fixed$estimated_params$ln_rec_sigma <- ebs_ricker_run_est$estimated_params$ln_rec_sigma
ebs_ricker_run_fixed <- fix_sigmaR(ebs_ricker_run_fixed, fix_sigmaR = TRUE)

ebs_ricker_run_re_fixed <- fix_sigmaR(ebs_ricker_run_re, fix_sigmaR = TRUE, bias.correct = FALSE)
ebs_ricker_run_re_est <- fix_sigmaR(ebs_ricker_run_re, fix_sigmaR = FALSE, bias.correct = FALSE)
ebs_ricker_run_re_fixedBC <- fix_sigmaR(ebs_ricker_run_re, fix_sigmaR = TRUE, bias.correct = TRUE)
ebs_ricker_run_re_estBC <- fix_sigmaR(ebs_ricker_run_re, fix_sigmaR = FALSE, bias.correct = TRUE)

# * YFS ----
yfs_model_est <- fix_sigmaR(yfs_model, fix_sigmaR = FALSE)
yfs_model_fixed <- yfs_model
yfs_model_fixed$estimated_params$ln_rec_sigma <- yfs_model_est$estimated_params$ln_rec_sigma
yfs_model_fixed <- fix_sigmaR(yfs_model_fixed, fix_sigmaR = TRUE)

yfs_model_re_fixed <- fix_sigmaR(yfs_model_re, fix_sigmaR = TRUE, bias.correct = FALSE)
yfs_model_re_est <- fix_sigmaR(yfs_model_re, fix_sigmaR = FALSE, bias.correct = FALSE)
yfs_model_re_fixedBC <- fix_sigmaR(yfs_model_re, fix_sigmaR = TRUE, bias.correct = TRUE)
yfs_model_re_estBC <- fix_sigmaR(yfs_model_re, fix_sigmaR = FALSE, bias.correct = TRUE)


# * GOA ----
# -- Pollock
pollock_model_est <- fix_sigmaR(pollock_model, fix_sigmaR = FALSE)
pollock_model_fixed <- pollock_model
pollock_model_fixed$estimated_params$ln_rec_sigma <- pollock_model_est$estimated_params$ln_rec_sigma
pollock_model_fixed <- fix_sigmaR(pollock_model_fixed, fix_sigmaR = TRUE)

pollock_model_re_fixed <- fix_sigmaR(pollock_model_re, fix_sigmaR = TRUE, bias.correct = FALSE)
pollock_model_re_est <- fix_sigmaR(pollock_model_re, fix_sigmaR = FALSE, bias.correct = FALSE)
pollock_model_re_fixedBC <- fix_sigmaR(pollock_model_re, fix_sigmaR = TRUE, bias.correct = TRUE)
pollock_model_re_estBC <- fix_sigmaR(pollock_model_re, fix_sigmaR = FALSE, bias.correct = TRUE)

# -- Cod
cod_model_est <- fix_sigmaR(cod_model, fix_sigmaR = FALSE)
cod_model_fixed <- cod_model
cod_model_fixed$estimated_params$ln_rec_sigma <- cod_model_est$estimated_params$ln_rec_sigma
cod_model_fixed <- fix_sigmaR(cod_model_fixed, fix_sigmaR = TRUE)

cod_model_re_fixed <- fix_sigmaR(cod_model_re, fix_sigmaR = TRUE, bias.correct = FALSE)
cod_model_re_est <- fix_sigmaR(cod_model_re, fix_sigmaR = FALSE, bias.correct = FALSE)
cod_model_re_fixedBC <- fix_sigmaR(cod_model_re, fix_sigmaR = TRUE, bias.correct = TRUE)
cod_model_re_estBC <- fix_sigmaR(cod_model_re, fix_sigmaR = FALSE, bias.correct = TRUE)

# -- ATF
# atf_model_est <- fix_sigmaR(atf_model, fix_sigmaR = FALSE)
# atf_model$estimated_params$ln_rec_sigma <- atf_model_est$estimated_params$ln_rec_sigma
atf_model_fixed <- fix_sigmaR(atf_model, fix_sigmaR = TRUE)

atf_model_re_fixed <- fix_sigmaR(atf_model_re, fix_sigmaR = TRUE, bias.correct = FALSE)
atf_model_re_est <- fix_sigmaR(atf_model_re, fix_sigmaR = FALSE, bias.correct = FALSE)
atf_model_re_fixedBC <- fix_sigmaR(atf_model_re, fix_sigmaR = TRUE, bias.correct = TRUE)
atf_model_re_estBC <- fix_sigmaR(atf_model_re, fix_sigmaR = FALSE, bias.correct = TRUE)


bias_correct_list <- list(
  ebs_run = ebs_run,
  ebs_run_est = ebs_run_est,
  ebs_run_fixed = ebs_run_fixed,
  ebs_run_re_fixed = ebs_run_re_fixed,
  ebs_run_re_est = ebs_run_re_est,
  ebs_run_re_fixedBC = ebs_run_re_fixedBC,
  ebs_run_re_estBC = ebs_run_re_estBC,

  ebs_ricker_run = ebs_ricker_run,
  ebs_ricker_run_est = ebs_ricker_run_est,
  ebs_ricker_run_fixed = ebs_ricker_run_fixed,
  ebs_ricker_run_re_fixed = ebs_ricker_run_re_fixed,
  ebs_ricker_run_re_est = ebs_ricker_run_re_est,
  ebs_ricker_run_re_fixedBC = ebs_ricker_run_re_fixedBC,
  ebs_ricker_run_re_estBC = ebs_ricker_run_re_estBC,

  yfs_model = yfs_model,
  yfs_model_est = yfs_model_est,
  yfs_model_fixed = yfs_model_fixed,
  yfs_model_re_fixed = yfs_model_re_fixed,
  yfs_model_re_est = yfs_model_re_est,
  yfs_model_re_fixedBC = yfs_model_re_fixedBC,
  yfs_model_re_estBC = yfs_model_re_estBC,

  pollock_model = pollock_model,
  pollock_model_est = pollock_model_est,
  pollock_model_fixed = pollock_model_fixed,
  pollock_model_re_fixed = pollock_model_re_fixed,
  pollock_model_re_est = pollock_model_re_est,
  pollock_model_re_fixedBC = pollock_model_re_fixedBC,
  pollock_model_re_estBC = pollock_model_re_estBC,

  cod_model = cod_model,
  cod_model_est = cod_model_est,
  cod_model_fixed = cod_model_fixed,
  cod_model_re_fixed = cod_model_re_fixed,
  cod_model_re_est = cod_model_re_est,
  cod_model_re_fixedBC = cod_model_re_fixedBC,
  cod_model_re_estBC = cod_model_re_estBC,

  #atf_model = atf_model,
  #atf_model_est = atf_model_est,
  atf_model_fixed = atf_model_fixed,
  atf_model_re_fixed = atf_model_re_fixed,
  atf_model_re_est = atf_model_re_est,
  atf_model_re_fixedBC = atf_model_re_fixedBC,
  atf_model_re_estBC = atf_model_re_estBC
)


# * Get quantities ----
data_tmp <- list()
ind = 1
for(i in 1:length(bias_correct_list)){
  for(sp in 1:bias_correct_list[[i]]$data_list$nspp){
    for(j in 1:3){ # SSB, R, or B

      # Get years
      yrs <- bias_correct_list[[i]]$data_list$styr:bias_correct_list[[i]]$data_list$endyr

      # Get output
      quantity <- bias_correct_list[[i]]$quantities[[c("R", "biomassSSB", "biomass")[j]]]

      # Get SD
      sd_temp <- which(names(bias_correct_list[[i]]$sdrep$value) == c("R", "biomassSSB", "biomass")[j])
      sd_temp <- bias_correct_list[[i]]$sdrep$sd[sd_temp]
      se_temp <- quantity
      se_temp <- replace(se_temp, values = sd_temp)

      # - Get unbiased
      if(!is.null(bias_correct_list[[i]]$sdrep$unbiased)){
        name_loc <- which(names(bias_correct_list[[i]]$sdrep$unbiased$value) == c("R", "biomassSSB", "biomass")[j])
        bs_temp <- bias_correct_list[[i]]$sdrep$unbiased$value[name_loc]
        quantity <- quantity
        quantity <- replace(quantity, values = bs_temp)

        sd_temp <- bias_correct_list[[i]]$sdrep$unbiased$sd[name_loc]
        se_temp <- quantity
        se_temp <- replace(se_temp, values = sd_temp)
      }

      # Cut to years and species
      quantity <- quantity[sp,1:length(yrs)]
      se_temp <- se_temp[sp,1:length(yrs)]

      # SE of sigmaR
      sigmaR_se <- which(names(bias_correct_list[[i]]$sdrep$par.fixed) == "ln_rec_sigma")
      sigmaR_se <- sqrt(diag(bias_correct_list[[i]]$sdrep$cov.fixed)[sigmaR_se])

      # Assign to df
      data_tmp[[ind]] <- data.frame(
        Model = names(bias_correct_list)[i],
        Ricker = bias_correct_list[[i]]$data_list$srr_pred_fun > 0,
        LaplaceApprox = bias_correct_list[[i]]$data_list$random_rec,
        EstSigmaR = !is.na(bias_correct_list[[i]]$map$mapList$ln_rec_sigma[sp]),
        BiasCorrect = !is.null(bias_correct_list[[i]]$opt$SD$unbiased),
        Species = bias_correct_list[[i]]$data_list$spnames[sp],
        sigmaR = exp(bias_correct_list[[i]]$estimated_params$ln_rec_sigma[sp]),
        ln_sigmaR_SE = sigmaR_se[sp],
        Year = yrs,
        Quantity = c("R", "SSB", "Biomass")[j],
        Value = quantity,
        SE = se_temp)
      ind = ind+1
    }
  }
}

exp_results <- do.call("rbind", data_tmp)
write.csv(exp_results, file = "Results/sigmaR_experiments.csv")

# - System
exp_results$System <- ifelse(grepl("ebs_", exp_results$Model), "EBS", "GOA")

# Figures ----
library(ggplot2)
library(cowplot)
species_df <- exp_results %>%
  distinct(System, Ricker, Species)

for(i in 1:nrow(species_df)){
  exp_results_spp <- exp_results %>%
    filter(System == species_df$System[i] &
             Ricker == species_df$Ricker[i] &
             Species == species_df$Species[i])

  quantities <- c("R", "SSB", "Biomass")
  p1 <- list()
  p2 <- list()
  for(k in 1:3){
    exp_results_tmp <- exp_results_spp %>%
      filter(Quantity == quantities[k])


    exp_results_tmp$Group = NA
    for(j in 1:nrow( exp_results_tmp)){
      exp_results_tmp$Group[j] <- paste0(
        ifelse(as.logical(exp_results_tmp$LaplaceApprox[j]), "MML;", "PML;"),
        ifelse(as.logical(exp_results_tmp$EstSigmaR[j]), " Est;", " Fix;"),
        ifelse(as.logical(exp_results_tmp$BiasCorrect[j]), " Bias.c", " No Bias.c"),
        "; SigmaR = ", round(exp_results_tmp$sigmaR[j], 3), "(",round(exp_results_tmp$ln_sigmaR_SE[j], 3), ")" )
    }


    p1[[k]] <- exp_results_tmp %>%
      ggplot(aes(x = Year, y = Value, colour = Group)) +
      geom_line(linewidth = 1.1) +
      ylab(quantities[k]) +
      theme_classic()

    p2[[k]] <- exp_results_tmp %>%
      ggplot(aes(x = Year, y = SE, colour = Group)) +
      geom_line(linewidth = 1.1) +
      ylab(paste(quantities[k], "SE")) +
      theme_classic() +
      theme(legend.position="none")
  }



  g1 <- plot_grid(p1[[1]] + theme(legend.position="none"), p2[[1]],
                  p1[[2]] + theme(legend.position="none"), p2[[2]],
                  p1[[3]] + theme(legend.position="none"), p2[[3]],
                  nrow = 3)

  grobs <- ggplotGrob(p1[[3]])$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

  g2 <- plot_grid(g1, legend, nrow = 2, rel_heights = c(0.75,0.25))
  ggsave(g2,
         filename = paste("Results/SigmaR experiment-", species_df$System[i], ifelse(species_df$Ricker[i], "Ricker", "MeanR"), species_df$Species[i], ".png"),
         width = 7, height = 10)
}
