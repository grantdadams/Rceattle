testthat::test_that("Composition likelihoods match (Multinomial and Dirichlet-Multinomial)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip()

  # ========================================================================
  # SIMULATE DATA
  # ========================================================================
  set.seed(123)
  sim <- make_msm_test_data(nspp = 2, years = 1:10) # 1 species, 10 years to keep it fast
  simData <- sim$data_list

  # Test the likelihood for the first survey fleet
  flt_idx <- 1

  # Find the data rows for this fleet
  comp_ind <- which(simData$comp_data$Fleet_code == flt_idx)


  # ========================================================================
  # TEST 1: MULTINOMIAL (comp_ll_type = 0)
  # ========================================================================
  simData_mn <- simData
  simData_mn$fleet_control$Comp_loglike <- 0 # 0 = Multinomial
  simData_mn$fleet_control$Comp_weights <- 0.75

  inits_mn <- build_params(simData_mn)
  # Set a specific likelihood weight for testing (e.g., 0.75)
  inits_mn$comp_weights[flt_idx] <- 0.75

  mod_mn <- Rceattle::fit_mod(data_list = simData_mn,
                              inits = inits_mn,
                              estimateMode = 3, # Just evaluate NLL, do not optimize
                              phase = FALSE,
                              verbose = 0)

  # Extract the NLL calculated by C++ (Slot 2 in C++ is Row 3 in R)
  tmb_nll_mn <- mod_mn$quantities$jnll_comp[3, flt_idx]

  # Reconstruct the NLL in R
  r_nll_mn <- 0
  for(i in comp_ind){
    row_data <- mod_mn$data_list$comp_data[i,]
    n_ages <- mod_mn$data_list$nages[1]

    # Raw proportions from data
    obs_prop <- as.numeric(row_data[grep("Comp_", colnames(row_data))])[1:n_ages]
    obs_prop <- obs_prop/sum(obs_prop) # Normalize
    samp_size <- row_data$Sample_size

    # Expected proportions from TMB
    hat_prop <- mod_mn$quantities$comp_hat[i, 1:n_ages]

    # Apply C++ math offsets exactly as written in ceattle_v01_11.cpp
    obs_prop_offset <- obs_prop + 0.00001
    hat_prop_offset <- hat_prop + 0.00001

    obs_prop_offset <- obs_prop_offset/sum(obs_prop_offset)
    hat_prop_offset <- hat_prop_offset/sum(hat_prop_offset)

    # Convert proportion to numbers
    obs_num <- obs_prop_offset * samp_size

    r_nll_mn <- r_nll_mn + (inits_mn$comp_weights[flt_idx] * calc_multinom_nll(obs_num, hat_prop_offset))
  }

  testthat::expect_equal(tmb_nll_mn, as.numeric(r_nll_mn), tolerance = 1e-5)


  # ========================================================================
  # TEST 2: DIRICHLET-MULTINOMIAL (comp_ll_type = 1)
  # ========================================================================
  simData_dm <- simData
  simData_dm$fleet_control$Comp_loglike <- 1 # 1 = Dirichlet-Multinomial

  inits_dm <- build_params(simData_dm)

  # In Dirichlet-Multinomial, 'comp_weights' acts as log(theta).
  inits_dm$comp_weights[flt_idx] <- log(5.0)

  mod_dm <- Rceattle::fit_mod(data_list = simData_dm,
                              inits = inits_dm,
                              estimateMode = 3,
                              phase = FALSE,
                              verbose = 0)

  # Extract NLL
  tmb_nll_dm <- mod_dm$quantities$jnll_comp[3, flt_idx]

  # Reconstruct the NLL in R
  r_nll_dm <- 0
  theta <- exp(inits_dm$comp_weights[flt_idx]) # DM_pars_comp

  for(i in comp_ind){
    row_data <- mod_dm$data_list$comp_data[i,]
    n_ages <- mod_dm$data_list$nages[1]

    obs_prop <- as.numeric(row_data[grep("Comp_", colnames(row_data))])[1:n_ages]
    obs_prop <- obs_prop/sum(obs_prop)
    samp_size <- row_data$Sample_size
    hat_prop <- mod_dm$quantities$comp_hat[i, 1:n_ages]

    # Apply C++ offsets
    obs_prop_offset <- obs_prop + 0.00001
    hat_prop_offset <- hat_prop + 0.00001

    obs_prop_offset <- obs_prop_offset/sum(obs_prop_offset)
    hat_prop_offset <- hat_prop_offset/sum(hat_prop_offset)

    # Calculate Dirichlet Alphas
    obs_num <- obs_prop_offset * samp_size
    alpha <- samp_size * hat_prop_offset * theta

    r_nll_dm <- r_nll_dm + calc_dirmultinom_nll(obs_num, alpha)
  }

  testthat::expect_equal(tmb_nll_dm, r_nll_dm, tolerance = 1e-5)
})
