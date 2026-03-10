

testthat::test_that("Diet proportion multinomial likelihood (jnll_comp) matches R math", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip()


  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  inits <- suppressMessages( build_params(simData) )
  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  mod <- Rceattle::fit_mod(data_list = simData,
                           inits = inits, # Initial parameters from inits
                           file = NULL, # Don't save
                           estimateMode = 3, # Don't estimate
                           random_rec = FALSE, # No random recruitment
                           phase = FALSE,
                           msmMode = 1,
                           suitMode = 4,
                           niter = 5,
                           initMode = 2,
                           verbose = 0)

  # Extract the NLL calculated by C++
  tmb_diet_nll <- mod$quantities$jnll_comp[19, 1:nspp]


  # RECONSTRUCT DIET LIKELIHOOD IN R
  r_diet_nll <- rep(0, nspp)

  # Rceattle stores stomach_id in 0-indexed format for TMB, adjust to 1-indexed for R
  re_data <-Rceattle::rearrange_dat(mod$data_list)
  stomach_ids <-re_data$stomach_id + 1
  unique_stomachs <- unique(stomach_ids)

  for (i in unique_stomachs) {

    # Get all rows associated with this specific stomach/group
    idx <- which(stomach_ids == i)

    # Extract predator species (diet_ctl column 1 is 1-indexed for species)
    rsp <- re_data$diet_ctl[idx[1], 1]

    # Only calculate if suitability is turned on for this species
    if (mod$data_list$suitMode[rsp] <= 0) next

    # Extract Sample Size (N_s)
    N_s <- re_data$diet_obs[idx[1], 1]

    # Extract observed and predicted proportions
    obs_prop <- re_data$diet_obs[idx, 2]
    pred_prop <- mod$quantities$diet_hat[idx, 2]

    # --- 1. Add "Other prey" category mathematically ---
    sum_obs_p <- sum(obs_prop)
    obs_prop <- c(obs_prop, 1.0 - sum_obs_p)

    sum_est_p <- sum(pred_prop)
    pred_prop <- c(pred_prop, 1.0 - sum_est_p)

    # --- 2. Apply C++ Offsets & Normalize ---
    obs_prop <- obs_prop + 0.00001
    pred_prop <- pred_prop + 0.00001

    obs_prop <- obs_prop / sum(obs_prop)
    pred_prop <- pred_prop / sum(pred_prop)

    # --- 3. Calculate Likelihood ---
    obs_counts <- obs_prop * N_s
    stomach_nll <- calc_multinom_nll(obs_counts, pred_prop)

    # Apply diet_comp_weight for the specific predator
    weight <- mod$data_list$Diet_comp_weights[rsp]

    # Accumulate species-specific NLL
    r_diet_nll[rsp] <- r_diet_nll[rsp] + (weight * stomach_nll)
  }

  # Verify that the R calculation perfectly matches the C++ tape
  testthat::expect_equal(as.numeric(tmb_diet_nll), as.numeric(r_diet_nll), tolerance = 1e-5)
})

