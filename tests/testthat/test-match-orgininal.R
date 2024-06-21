test_that("Dynamics match CEATTLE classic", {

  # Load old model
  load("~/GitHub/Rceattle/inst/extdata/CEATTLE_classic_SS.Rdata")


  # Load data and set up inits
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$srr_prior_mean <- 9
  BS2017SS$initMode  <- 1
  inits <- build_params(BS2017SS)

  # - Update population dynamics from previous parameters
  inits$init_dev <- ss_run$estimated_params$init_dev
  inits$rec_dev[,1:39] <- ss_run$estimated_params$rec_dev
  inits$rec_pars[,1] <- ss_run$estimated_params$ln_mn_rec
  inits$F_dev[1:3, 1:39] <- ss_run$estimated_params$F_dev[,1:39]
  inits$ln_mean_F[1:3] <- ss_run$estimated_params$ln_mean_F
  inits$sel_coff[1:3,1,] <- ss_run$estimated_params$fsh_sel_coff


  ss_run_old_params <- Rceattle::fit_mod(data_list = mydata,
                                         inits = inits, # Initial parameters = 0
                                         file = NULL, # Don't save
                                         estimateMode = 3, # Estimate
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0, # Single species mode
                                         phase = "default",
                                         verbose = 1)

  # - Estimate instead
  ss_run_new <- Rceattle::fit_mod(data_list = mydata,
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
  plot_biomass(list(ss_run_new, ss_run_old_params, ss_run_old))
})
#> Test passed 😀
#>
#>
#>
#>
test_that("Dynamics match multi-species CEATTLE classic", {

  # Load old model
  load("~/GitHub/Rceattle/inst/extdata/CEATTLE_classic_MS.Rdata")

  # Load data and set up inits
  data("BS2017SS")
  data("BS2017MS")
  BS2017MS$srr_prior_mean <- 9
  BS2017MS$initMode  <- 1
  BS2017MS$M1_base[,-c(1:2)] <- ms_run$quantities$M1 # + 0.0001 is in the old one
  inits <- build_params(BS2017MS)

  # - Update population dynamics from old parameters
  inits$init_dev[,1:20] <- ms_run$estimated_params$init_dev
  inits$rec_dev[,1:39] <- ms_run$estimated_params$rec_dev
  inits$rec_pars[,1] <- ms_run$estimated_params$ln_mn_rec
  inits$F_dev[1:3, 1:39] <- ms_run$estimated_params$F_dev[,1:39]
  inits$ln_mean_F[1:3] <- ms_run$estimated_params$ln_mean_F
  inits$sel_coff[1:3,1,] <- ms_run$estimated_params$fsh_sel_coff


  ms_run_old_params <- Rceattle::fit_mod(
    data_list = BS2017MS,
    inits = inits, # Initial parameters from old model
    M1Fun = build_M1(M1_model = 0,
                     updateM1 = TRUE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE),
    file = NULL, # Don't save
    estimateMode = 3, # Estimate
    niter = 10, # 3 iterations around population and predation dynamics
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
    phase = "default",
    verbose = 1)


  # - Multi-species
  ms_run_new <- Rceattle::fit_mod(
    data_list = BS2017MS,
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
  ms_run_old$quantities$R[,1:39] <- ms_run$quantities$R[,1:39]
  ms_run_old$quantities$biomass[,1:39] <- ms_run$quantities$biomass
  ms_run_old$quantities$biomassSSB[,1:39] <- ms_run$quantities$biomassSSB


  # Plot
  plot_recruitment(list(ms_run_old_params, ms_run_old))
  plot_biomass(list(ms_run_new, ms_run_old_params, ms_run_old))
  plot_ssb(list(ms_run_new, ms_run_old_params, ms_run_old))


})
#> Test passed 😀
#>

