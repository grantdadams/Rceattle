# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Composition likelihoods match (Multinomial and Dirichlet-Multinomial)", {
  # DM likelihood drifts ~1e-5 under covr's -O0 instrumentation vs the -O2 build.
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")


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
  simData_mn$fleet_control$Comp_distribution <- 0 # 0 = Multinomial
  simData_mn$fleet_control$Comp_weights <- 0.75

  # * Inits ----
  ss_run <- Rceattle::fit_mod(data_list = simData_mn,
                              estimateMode = 3,
                              fit_control = fit_control(
                                phase = FALSE,
                                verbose = 0))
  inits_mn <- ss_run$estimated_params

  # Set a specific likelihood weight for testing (e.g., 0.75)
  inits_mn$comp_weights[flt_idx] <- 0.75

  mod_mn <- Rceattle::fit_mod(data_list = simData_mn,
                              inits = inits_mn,
                              estimateMode = 3, # Just evaluate NLL, do not optimize
                              fit_control = fit_control(
                                phase = FALSE,
                                verbose = 0))

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

    # C++ adds the comp_offset (default 1e-5) to obs and pred, then the OSA
    # conditional-binomial multinomial (dmultinom_osa) normalizes the predicted
    # proportions to sum to 1 (matching WHAM). Reproduce that: offset both, then
    # normalize the predicted proportions before the log term.
    obs_prop_offset <- obs_prop + 0.00001
    hat_prop_offset <- hat_prop + 0.00001
    hat_prop_offset <- hat_prop_offset / sum(hat_prop_offset)

    obs_num <- obs_prop_offset * samp_size

    ll <- lgamma(sum(obs_num) + 1) - sum(lgamma(obs_num + 1)) +
      sum(obs_num * log(hat_prop_offset))
    r_nll_mn <- r_nll_mn + (inits_mn$comp_weights[flt_idx] * (-ll))
  }

  testthat::expect_equal(tmb_nll_mn, as.numeric(r_nll_mn), tolerance = 1e-5)


  # ========================================================================
  # TEST 2: DIRICHLET-MULTINOMIAL (comp_ll_type = 1)
  # ========================================================================
  simData_dm <- simData
  simData_dm$fleet_control$Comp_distribution <- 1 # 1 = Dirichlet-Multinomial
  # In Dirichlet-Multinomial, 'comp_weights' acts as log(theta). Set in
  # fleet_control so the value survives the second fit_mod() call.
  simData_dm$fleet_control$Comp_weights[flt_idx] <- log(5.0)

  # * Inits ----
  ss_run <- Rceattle::fit_mod(data_list = simData_dm,
                              estimateMode = 3,
                              fit_control = fit_control(
                                phase = FALSE,
                                verbose = 0))
  inits_dm <- ss_run$estimated_params
  inits_dm$comp_weights[flt_idx] <- log(5.0)

  mod_dm <- Rceattle::fit_mod(data_list = simData_dm,
                              inits = inits_dm,
                              estimateMode = 3,
                              fit_control = fit_control(
                                phase = FALSE,
                                verbose = 0))

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

    # C++ adds the 0.00001 offset but does NOT renormalize before alphas.
    obs_prop_offset <- obs_prop + 0.00001
    hat_prop_offset <- hat_prop + 0.00001

    obs_num <- obs_prop_offset * samp_size
    alpha <- samp_size * hat_prop_offset * theta

    r_nll_dm <- r_nll_dm + calc_dirmultinom_nll(obs_num, alpha)
  }

  testthat::expect_equal(tmb_nll_dm, r_nll_dm, tolerance = 1e-5)
})


testthat::test_that("Invalid comp_loglike", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Comp_distribution <- 9 # Not in scope

  # Run
  testthat::expect_error(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

})

testthat::test_that("Invalid caal_loglike", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$CAAL_distribution <- 9 # Not in scope

  # Run
  testthat::expect_error(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

})


testthat::test_that("Integer to text multinomial comp_loglike", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Comp_distribution <- 0

  # Run
  mod <- Rceattle::fit_mod(data_list = GOA2018SS,
                           estimateMode = 3, # Don't estimate
                           fit_control = fit_control(
                             verbose = 0))

  testthat::expect_equal(mod$data_list$fleet_control$Comp_distribution,
                         rep("Multinomial", nrow(mod$data_list$fleet_control))
  )

})

testthat::test_that("Integer to text DM comp_loglike", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Comp_distribution <- 1

  # Run
  mod <- Rceattle::fit_mod(data_list = GOA2018SS,
                           estimateMode = 3, # Don't estimate
                           fit_control = fit_control(
                             verbose = 0))

  testthat::expect_equal(mod$data_list$fleet_control$Comp_distribution,
                         rep("DirichletMultinomial", nrow(mod$data_list$fleet_control))
  )

  nflt <- nrow(GOA2018SS$fleet_control)
  fleets <- 1:nflt
  no_comp_data <- !(fleets %in% mod$data_list$comp_data$Fleet_code)
  fleets[mod$data_list$fleet_control$Fleet_type == "Off"] <- NA
  fleets[no_comp_data] <- NA

  testthat::expect_equal(as.numeric(mod$map$mapList$comp_weights), fleets)
})



