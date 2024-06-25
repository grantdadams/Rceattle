test_that("Dynamics match CEATTLE classic", {

  # Load old model
  load("~/GitHub/Rceattle/inst/extdata/CEATTLE_classic_SS.Rdata")
  load("C:/Users/Grant Adams/Documents/CEATTLE Runs/CEATTLE_classic_SS.Rdata")


  # Load data and set up inits
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$srr_prior_mean <- 9
  BS2017SS$initMode  <- 1
  inits <- build_params(BS2017SS)

  # - Update population dynamics from previous parameters
  inits$init_dev <- ss_run$estimated_params$init_dev
  inits$rec_dev[,1:39] <- ss_run$estimated_params$rec_dev
  inits$ln_mean_rec <- ss_run$estimated_params$ln_mn_rec
  inits$F_dev[1:3, 1:39] <- ss_run$estimated_params$F_dev[,1:39]
  inits$ln_mean_F[1:3] <- ss_run$estimated_params$ln_mean_F
  inits$sel_coff[1:3,1,] <- ss_run$estimated_params$fsh_sel_coff


  ss_run_old_params <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = inits, # Initial parameters = 0
                                         file = NULL, # Don't save
                                         estimateMode = 3, # Estimate
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0, # Single species mode
                                         phase = "default",
                                         verbose = 1)

  # - Estimate instead
  ss_run_new <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = NULL, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 0, # Single species mode
                                  phase = "default",
                                  verbose = 1)

  # Previous time series
  ss_run_old <- ss_run_new
  ss_run_old$quantities$R[,1:39] <- ss_run$quantities$R[,1:39]
  ss_run_old$quantities$biomass[,1:39] <- ss_run$quantities$biomass
  ss_run_old$quantities$biomassSSB[,1:39] <- ss_run$quantities$biomassSSB


  # Plot
  plot_recruitment(list(ss_run_new, ss_run_old_params, ss_run_old))
  plot_biomass(list(ss_run_new, ss_run_old_params, ss_run_old), model_names = 1:3)
})
#> Test passed 😀
#>
#>
#>
#>
test_that("Dynamics match multi-species CEATTLE classic", {

  # Load old model
  # load("~/GitHub/Rceattle/inst/extdata/CEATTLE_classic_MS.Rdata")
  load("C:/Users/Grant Adams/Documents/CEATTLE Runs/CEATTLE_classic_MS.Rdata")

  # Load data and set up inits
  library(Rceattle)
  data("BS2017SS")
  data("BS2017MS")
  BS2017MS$srr_prior_mean <- 9
  BS2017MS$initMode  <- 1
  BS2017MS$M1_base[,-c(1:2)] <- ms_run$quantities$M1 # + 0.0001 is in the old one
  inits <- build_params(BS2017MS)

  # - Update population dynamics from old parameters
  inits$init_dev[,1:20] <- ms_run$estimated_params$init_dev
  inits$rec_dev[,1:39] <- ms_run$estimated_params$rec_dev
  inits$ln_mean_rec <- ms_run$estimated_params$ln_mn_rec
  inits$F_dev[1:3, 1:39] <- ms_run$estimated_params$F_dev[,1:39]
  inits$ln_mean_F[1:3] <- ms_run$estimated_params$ln_mean_F
  inits$sel_coff[1:3,1,] <- ms_run$estimated_params$fsh_sel_coff

  # - Update stomach data
  check <- as.data.frame(BS2017MS$UobsWtAge)
  check$OldStomData <- NA
  check$OldStom <- NA
  check$NewStom <- NA

  for(i in 1:nrow(check)){
    check$OldStomData[i] <- ms_run$data_list$UobsAge[check$Pred[i], check$Prey[i],check$Pred_age[i], check$Prey_age[i]]
    check$OldStom[i] <- ms_run$quantities$stomKir[check$Pred[i], check$Prey[i],check$Pred_age[i], check$Prey_age[i], 1]
    # check$NewStom[i] <- ms_run_old_params$quantities$diet_prop_weight[check$Pred[i], check$Prey[i],check$Pred_age[i], check$Prey_age[i], 1]
  }

  BS2017MS$UobsWtAge$Stomach_proportion_by_weight <- check$OldStomData

  # - Update Pyrs data
  check <- as.data.frame(BS2017MS$Pyrs)

  for(i in 1:nrow(check)){
    check[i,3:23] <- ms_run$data_list$Pyrs[check$Year[i]-1978, , check$Species[i]]
  }

  BS2017MS$Pyrs <- check

  # - Update dynamics
  ms_run_old_params <- Rceattle::fit_mod(
    data_list = BS2017MS,
    inits = inits, # Initial parameters from old model
    file = NULL, # Don't save
    estimateMode = 4, # Estimate
    niter = 10, # 3 iterations around population and predation dynamics
    random_rec = FALSE, # No random recruitment
    msmMode = 1, # MSVPA based
    suitMode = 0, # empirical suitability
    verbose = 1)
#
#   ms_run$quantities$stomKir[3,1,,1:10,1]
#   ms_run_old_params$quantities$diet_prop_weight[3,1,,1:10,1]
#
#   # Previous time series
#   ms_run_old <- ms_run_old_params
#   ms_run_old$quantities$R[,1:39] <- ms_run$quantities$R[,1:39]
#   ms_run_old$quantities$biomass[,1:39] <- ms_run$quantities$biomass
#   ms_run_old$quantities$biomassSSB[,1:39] <- ms_run$quantities$biomassSSB
#
#
#   # Plot
#   plot_recruitment(list(ms_run_old_params, ms_run_old))
#   plot_biomass(list(ms_run_old_params, ms_run_old))
#
#
#   # - Check diet components
#   sum(ms_run$quantities$stomKir[1:3,1:3,,,1:39] != ms_run_old_params$quantities$diet_prop_weight[1:3,1:3,,,1:39]) # Stomacc is good
#   sum(ms_run$quantities$of_stomKir[1:3,,1:39] != ms_run_old_params$quantities$other_food_diet_prop_weight[1:3,1,,1:39]) # Other food is good
#   sum(abs(ms_run$quantities$AvgN[1:3,,1:39] - ms_run_old_params$quantities$AvgN[1:3,1,,1:39]) > 1e-6) # AvgN is NOT good
#   sum(abs(ms_run$quantities$NByage[1:3,,1:39] - ms_run_old_params$quantities$NByage[1:3,1,,1:39]) > 1e-8) # AvgN is NOT good
#   sum(abs(ms_run$quantities$S[1:3,,1:39] - ms_run_old_params$quantities$S[1:3,1,,1:39]) > 1e-8) # S is good
#   sum(abs(ms_run$quantities$Zed[1:3,,1:39] - ms_run_old_params$quantities$Zed[1:3,1,,1:39]) > 1e-8) # Zed is good
#
#   sum(abs(ms_run$quantities$stom_div_bio2[1:3,1:3,,,1:39] - ms_run_old_params$quantities$stom_div_bio2[1:3,1:3,,,1:39]) > 1e-8) # stom_div_bio2 is good
#
#   sum(abs(ms_run$quantities$suma_suit[1:3,,1:39] - ms_run_old_params$quantities$suma_suit[1:3,1,,1:39]) > 1e-8) # suma_suit is good
#
#   sum(abs(ms_run$quantities$suit_main[1:3,1:3,,,1:39] - ms_run_old_params$quantities$suit_main[1:3,1:3,,,1:39]) > 1e-8) # suit_main is good
#
#   sum(abs(ms_run$quantities$suit_other[1:3,] - ms_run_old_params$quantities$suit_other[1:3,1,,1]) > 1e-8) # suma_suit is good
#
#   sum(abs(ms_run$quantities$M2[1:3,,1:39] - ms_run_old_params$quantities$M2[1:3,1,,1:39]) > 0.1) # M2 is NOT good
#   sum(abs(ms_run$quantities$ration2Age[1:3,,1:39] - ms_run_old_params$quantities$ration[1:3,1,,1:39]) > 1e-1) # Ration is NOT good
#   sum(abs(ms_run$quantities$avail_food[1:3,,1:39] - ms_run_old_params$quantities$avail_food[1:3,1,,1:39]) > 1e-7) # suma_suit is good
#
#   sum(abs(ms_run$quantities$fT[1:3,] - ms_run_old_params$quantities$fT[1:3,1:39]) > 1e-8) # suma_suit is good
#   sum(abs(ms_run$quantities$ConsumAge[1:3,,] - ms_run_old_params$quantities$ConsumAge[1:3,1,,1:39]) > 1) # suma_suit is good
#
#   sum(abs(ms_run$data_list$Pyrs[1:3,,] - ms_run_old_params$data_list$Pyrs[1:3,1,,1:39]) > 1) # suma_suit is good



  #  Fit models
  ss_run_new <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = NULL, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 0, # Single species mode
                                  phase = "default",
                                  verbose = 1)


  # - Multi-species
  BS2017MS$est_M1 <- c(0,0,0) # Do not estimate residual M
  BS2017MS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
  ms_run_new <- Rceattle::fit_mod(data_list = BS2017MS,
                                  inits = ss_run_new$estimated_params, # Initial parameters from single species ests
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  niter = 10, # 3 iterations around population and predation dynamics
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 1, # MSVPA based
                                  suitMode = 0, # empirical suitability
                                  verbose = 1)



  # Previous time series
  ms_run_old <- ms_run_new
  ms_run_old$quantities$R[,1:39] <- ms_run$quantities$R[,1:39]
  ms_run_old$quantities$biomass[,1:39] <- ms_run$quantities$biomass
  ms_run_old$quantities$biomassSSB[,1:39] <- ms_run$quantities$biomassSSB


  # Plot
  plot_recruitment(list(ms_run_old_params, ms_run_old))
  plot_biomass(list(ms_run_old_params, ms_run_old))
  plot_biomass(list(ms_run_new, ms_run_old_params, ms_run_old))
  plot_ssb(list(ms_run_new, ms_run_old_params, ms_run_old))
  plot_ssb(list(ms_run_new, ms_run_old_params, ms_run_old))
  plot_f(list(ms_run_new, ms_run_old_params, ms_run_old))
})
#> Test passed 😀
#>

