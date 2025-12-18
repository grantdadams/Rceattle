test_that("Dynamics match CEATTLE classic", {

  # Load old model
  load("inst/extdata/CEATTLE_classic_ss.Rdata")


  # Load data and set up inits
  library(Rceattle)
  data(BS2017SS) # ?BS2017SS for more information on the data

  BS2017SS$srr_prior_mean <- 9
  BS2017SS$initMode  <- 1
  inits <- build_params(BS2017SS)

  # - Update population dynamics from previous parameters
  inits$init_dev <- CEATTLE_classic_SS$estimated_params$init_dev
  inits$rec_dev[,1:39] <- CEATTLE_classic_SS$estimated_params$rec_dev
  inits$rec_pars[,1] <- CEATTLE_classic_SS$estimated_params$ln_mn_rec
  inits$ln_F[1:3, 1:39] <- CEATTLE_classic_SS$estimated_params$F_dev[,1:39] + CEATTLE_classic_SS$estimated_params$ln_mean_F
  inits$sel_coff[1:3,1,] <- CEATTLE_classic_SS$estimated_params$fsh_sel_coff


  ss_run_old_params <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = inits, # Initial parameters = 0
                                         file = NULL, # Don't save
                                         estimateMode = 3, # Estimate
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0, # Single species mode
                                         phase = TRUE,
                                         verbose = 1)

  # - Estimate instead
  ss_run_new <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = NULL, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 0, # Single species mode
                                  phase = TRUE,
                                  verbose = 1)

  # Previous time series
  ss_run_old <- ss_run_new
  ss_run_old$quantities$R[,1:39] <- CEATTLE_classic_SS$quantities$R[,1:39]
  ss_run_old$quantities$biomass[,1:39] <- CEATTLE_classic_SS$quantities$biomass
  ss_run_old$quantities$ssb[,1:39] <- CEATTLE_classic_SS$quantities$biomassSSB


  # Plot
  plot_recruitment(list(ss_run_old_params, ss_run_old), model_names = 1:2)
  plot_biomass(list(ss_run_old_params, ss_run_old), model_names = 1:2)
  plot_ssb(list(ss_run_old_params, ss_run_old), model_names = 1:2)
  plot_index(list(ss_run_old_params, ss_run_old), model_names = 1:2)

  plot_recruitment(list(ss_run_new, ss_run_old_params, ss_run_old), model_names = 1:3)
  plot_biomass(list(ss_run_new, ss_run_old_params, ss_run_old), model_names = 1:3)
  plot_ssb(list(ss_run_new, ss_run_old_params, ss_run_old), model_names = 1:3)
})
#> Test passed 😀
#>
#>
#>
#>
test_that("Dynamics match multi-species CEATTLE classic", {

  # Load old model
  library(Rceattle)
  load("inst/extdata/CEATTLE_classic_MS.Rdata")

  # Load data and set up inits
  data("BS2017SS")
  data("BS2017MS")

  BS2017MS$srr_prior_mean <- 9
  BS2017MS$initMode  <- 2
  # BS2017MS$M1_base[,-c(1:2)] <- CEATTLE_classic_MS$quantities$M1 # + 0.0001 is in the old one

  inits <- build_params(BS2017MS)

  # - Update population dynamics from old parameters
  inits$init_dev[,1:20] <- CEATTLE_classic_MS$estimated_params$init_dev
  inits$rec_dev[,1:39] <- CEATTLE_classic_MS$estimated_params$rec_dev
  inits$rec_pars[,1] <- CEATTLE_classic_MS$estimated_params$ln_mn_rec
  inits$ln_F[1:3, 1:39] <- CEATTLE_classic_MS$estimated_params$F_dev[,1:39] + CEATTLE_classic_MS$estimated_params$ln_mean_F

  # -- Sel
  inits$sel_coff[1:3,1,] <- CEATTLE_classic_MS$estimated_params$fsh_sel_coff
  inits$sel_inf[1,4:6,1] <- CEATTLE_classic_MS$estimated_params$srv_sel_inf[1,]
  inits$ln_sel_slp[1,4:6,1] <- log(CEATTLE_classic_MS$estimated_params$srv_sel_slp[1,])
  inits$index_ln_q[4:6] <- CEATTLE_classic_MS$estimated_params$log_srv_q
  inits$index_ln_q[7] <- CEATTLE_classic_MS$estimated_params$log_eit_q # Need to scale by max sel because selectivity is rescaled to max = 1 in the CPP

  # - Update diet data
  BS2017MS_new <- BS2017MS
  BS2017MS_new$diet_data <- as.data.frame(BS2017MS_new$diet_data)
  for(i in 1:nrow(BS2017MS_new$diet_data)){
    BS2017MS_new$diet_data$Stomach_proportion_by_weight[i] <- CEATTLE_classic_MS$data_list$UobsAge[
      BS2017MS_new$diet_data$Pred[i], BS2017MS_new$diet_data$Prey[i],
      BS2017MS_new$diet_data$Pred_age[i], BS2017MS_new$diet_data$Prey_age[i]
    ]
  }

  # - Update foraging data
  BS2017MS_new$ration_data <- as.data.frame(BS2017MS_new$ration_data)
  BS2017MS_new$ration_data[,4:ncol(BS2017MS_new$ration_data)] <- rbind(CEATTLE_classic_MS$data_list$Pyrs[-40,,1],
                                                         CEATTLE_classic_MS$data_list$Pyrs[-40,,2],
                                                         CEATTLE_classic_MS$data_list$Pyrs[-40,,3])



  ms_run_old_params <- Rceattle::fit_mod(
    data_list = BS2017MS_new,
    inits = inits, # Initial parameters from old model
    M1Fun = build_M1(M1_model = 0,
                     updateM1 = TRUE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE),
    file = NULL, # Don't save
    estimateMode = 3, # Fix at parameter values
    niter = 10, # iterations around population and predation dynamics
    random_rec = FALSE, # No random recruitment
    msmMode = 1, # MSVPA based
    suitMode = 0, # empirical suitability
    verbose = 1)


  #  Fit models
  ss_run_new <- Rceattle::fit_mod(
    data_list = BS2017SS,
    inits = NULL, # Initial parameters = 0
    file = NULL, # Don't save
    estimateMode = 0, # Estimate
    random_rec = FALSE, # No random recruitment
    msmMode = 0, # Single species mode
    phase = TRUE,
    verbose = 1)


  # - Multi-species
  ms_run_new <- Rceattle::fit_mod(
    data_list = BS2017MS_new,
    inits = ss_run_new$estimated_params, # Initial parameters from single species ests
    M1Fun = build_M1(M1_model = 0,
                     updateM1 = TRUE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE),
    file = NULL, # Don't save
    estimateMode = 0, # Estimate
    niter = 10, # 3 iterations around population and predation dynamics
    random_rec = FALSE, # No random recruitment
    msmMode = 1, # MSVPA based
    suitMode = 0, # empirical suitability
    verbose = 1)



  # Previous time series
  ms_run_old <- ms_run_new
  ms_run_old$quantities$B_eaten_as_prey[,1,,1:39] <- CEATTLE_classic_MS$quantities$B_eaten
  ms_run_old$quantities$R[,1:39] <- CEATTLE_classic_MS$quantities$R[,1:39]
  ms_run_old$quantities$biomass[,1:39] <- CEATTLE_classic_MS$quantities$biomass
  ms_run_old$quantities$ssb[,1:39] <- CEATTLE_classic_MS$quantities$biomassSSB
  ms_run_old$quantities$index_hat <- c(CEATTLE_classic_MS$quantities$srv_bio_hat[1,1:36],
                                       CEATTLE_classic_MS$quantities$srv_bio_hat[2,1:36],
                                       CEATTLE_classic_MS$quantities$srv_bio_hat[3,1:36],
                                       CEATTLE_classic_MS$quantities$eit_hat[1:20])
  fsh_bio <- ms_run_old$data_list$catch_data
  fsh_bio$Pred <- ms_run_old$quantities$catch_hat
  fsh_bio$CatchOld <- NA
  fsh_bio$CatchOld[which(fsh_bio$Year < 2018)] <- c(CEATTLE_classic_MS$quantities$tc_biom_hat[1,],
                                                    CEATTLE_classic_MS$quantities$tc_biom_hat[2,],
                                                    CEATTLE_classic_MS$quantities$tc_biom_hat[3,])
  ms_run_old$quantities$catch_hat <- fsh_bio$CatchOld

  # Sel
  round(ms_run_old_params$quantities$sel[,1,,1], 5)
  round(CEATTLE_classic_MS$quantities$fsh_sel, 5)
  round(CEATTLE_classic_MS$quantities$srv_sel, 5)

  ms_run_old_params$quantities$sel[1:3,1,,1] - CEATTLE_classic_MS$quantities$fsh_sel
  ms_run_old_params$quantities$sel[4:6,1,,1] - CEATTLE_classic_MS$quantities$srv_sel

  # Plot
  plot_recruitment(list(ms_run_old_params, ms_run_old))
  plot_biomass(list(ms_run_old_params, ms_run_old))
  plot_ssb(list(ms_run_old_params, ms_run_old))
  plot_index(list(ms_run_old_params, ms_run_old), model_names = c(1,2))
  plot_catch(list(ms_run_old_params, ms_run_old))

  # ms_run_old$quantities$M1[,1,]
  # CEATTLE_classic_MS$quantities$M1
  # CEATTLE_classic_MS$quantities$ration2Age[1,1:10,1:39]
  # ms_run_old_params$quantities$ration[1,1,1:10,1:39]
  #
  #
  # CEATTLE_classic_MS$quantities$ConsumAge[1,1:10,1:39] - ms_run_old_params$quantities$ConsumAge[1,1,1:10,1:39]
  #
  # CEATTLE_classic_MS$quantities$Pyrs[1:39,1:12,1] - t(ms_run_old_params$quantities$Pyrs[1,1,1:12,1:39])
  #
  # head(t(CEATTLE_classic_MS$quantities$ConsumAge[1,1:12,1:39]))
  # head(t(ms_run_old_params$quantities$ConsumAge[1,1,1:10,1:39]))
  #
  #
  # ms_run_old_params$data_list$pop_wt_index
  #
  # # Weight
  # CEATTLE_classic_MS$data_list$wt[,1:12,1] - ms_run_old_params$data_list$wt %>%
  #   filter(Wt_index == 1) %>%
  #   dplyr::select(paste0("Age", 1:12)) %>%
  #   as.matrix()
  #
  #
  # CEATTLE_classic_MS$quantities$fT[,1:39] == ms_run_old_params$quantities$fT[,1:39]
  #
  # CEATTLE_classic_MS$data_list$CA == ms_run_old_params$data_list$CA
  #
  # CEATTLE_classic_MS$data_list$CB == ms_run_old_params$data_list$CB

  plot_recruitment(list(ms_run_new, ms_run_old_params, ms_run_old), model_names = 1:3)
  plot_biomass(list(ms_run_new, ms_run_old_params, ms_run_old), model_names = 1:3)
  plot_ssb(list(ms_run_new, ms_run_old_params, ms_run_old), model_names = 1:3)
  plot_b_eaten(list(ms_run_new, ms_run_old_params, ms_run_old), model_names = 1:3)

})
#> Test passed 😀
#>

