testthat::test_that("Composition likelihoods match (Multinomial and Dirichlet-Multinomial)", {
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
  simData_mn$fleet_control$Comp_loglike <- 0 # 0 = Multinomial
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

    # C++ adds the 0.00001 offset but does NOT renormalize, and TMB's dmultinom
    # does not renormalize p — so reproduce that here with unnormalized vectors.
    obs_prop_offset <- obs_prop + 0.00001
    hat_prop_offset <- hat_prop + 0.00001

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
  simData_dm$fleet_control$Comp_loglike <- 1 # 1 = Dirichlet-Multinomial
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


testthat::test_that("SS3Robust comp likelihood (case 2) matches closed-form kernel", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # SS3 Method-5 robust multinomial kernel form:
  #   NLL = weight * N * sum_a (obs_s_a * log(obs_s_a / hat_s_a))
  # where obs_s = (obs + ac)/(1 + n*ac) and hat_s = (hat + ac)/(1 + n*ac).
  # In the cpp, hat is pre-smoothed in the pred_comp section, so the loop
  # uses comp_hat directly as hat_s. (Mirrors SS3 SS_objfunc.tpl:430:
  #   age_like = -offset_a + Comp_logL_multinomial(N, obs, exp))

  set.seed(123)
  sim <- make_msm_test_data(nspp = 2, years = 1:10)
  simData <- sim$data_list

  flt_idx <- 1
  comp_ind <- which(simData$comp_data$Fleet_code == flt_idx)

  # Enable SS3Robust + a non-trivial addtocomp
  ac_val <- 1e-4
  simData$fleet_control$Comp_loglike      <- 2  # SS3Robust
  simData$fleet_control$Comp_weights      <- 0.6
  simData$fleet_control$Comp_addtocomp    <- ac_val

  # Build + run
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              estimateMode = 3,
                              fit_control = fit_control(phase = FALSE, verbose = 0))
  inits <- ss_run$estimated_params
  inits$comp_weights[flt_idx] <- 0.6

  mod <- Rceattle::fit_mod(data_list = simData,
                           inits = inits,
                           estimateMode = 3,
                           fit_control = fit_control(phase = FALSE, verbose = 0))

  tmb_nll <- mod$quantities$jnll_comp[3, flt_idx]

  # Reconstruct in R using the same kernel as the cpp case 2
  r_nll <- 0
  for (i in comp_ind) {
    row_data <- mod$data_list$comp_data[i, ]
    n_ages   <- mod$data_list$nages[1]
    n_bins   <- if (row_data$Age0_Length1 == 1) mod$data_list$nlengths[1] else n_ages

    obs_prop <- as.numeric(row_data[grep("Comp_", colnames(row_data))])[1:n_bins]
    obs_prop <- obs_prop / sum(obs_prop)
    samp     <- row_data$Sample_size
    hat_s    <- mod$quantities$comp_hat[i, 1:n_bins]   # already smoothed in cpp

    denom <- 1 + n_bins * ac_val
    obs_s <- (obs_prop + ac_val) / denom

    r_nll <- r_nll + 0.6 * samp * sum(obs_s * log(obs_s / pmax(hat_s, 1e-300)))
  }

  testthat::expect_equal(tmb_nll, as.numeric(r_nll), tolerance = 1e-5)
})


# NOTE: SS3Robust for CAAL uses the IDENTICAL kernel form as for comp
# (cpp case 2 in both the comp_ll_type and caal_ll_type switches:
#   jnll += weight * N * obs_s * log(obs_s / hat_s)
# with obs_s = (obs + ac)/(1 + n*ac), hat_s = comp/caal_hat (pre-smoothed)).
# The comp test above exercises this kernel directly. A separate synthetic
# CAAL test would require populating CAAL obs, setting parametric growth
# (so the cpp populates growth_matrix and pred_CAAL is non-zero), and
# matching the cpp's single-pop-bin selector logic -- all fragile in
# synthetic data. The kernel form is validated empirically on the GOA Pcod
# 2024 model: R-side reconstruction = 374.58 vs SS3 official 374.10.


testthat::test_that("Invalid comp_loglike", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Comp_loglike <- 9 # Not in scope

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
  GOA2018SS$fleet_control$CAAL_loglike <- 9 # Not in scope

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
  GOA2018SS$fleet_control$Comp_loglike <- 0

  # Run
  mod <- Rceattle::fit_mod(data_list = GOA2018SS,
                           estimateMode = 3, # Don't estimate
                           fit_control = fit_control(
                             verbose = 0))

  testthat::expect_equal(mod$data_list$fleet_control$Comp_loglike,
                         rep("Multinomial", nrow(mod$data_list$fleet_control))
  )

})

testthat::test_that("Integer to text DM comp_loglike", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Comp_loglike <- 1

  # Run
  mod <- Rceattle::fit_mod(data_list = GOA2018SS,
                           estimateMode = 3, # Don't estimate
                           fit_control = fit_control(
                             verbose = 0))

  testthat::expect_equal(mod$data_list$fleet_control$Comp_loglike,
                         rep("DirichletMultinomial", nrow(mod$data_list$fleet_control))
  )

  nflt <- nrow(GOA2018SS$fleet_control)
  fleets <- 1:nflt
  no_comp_data <- !(fleets %in% mod$data_list$comp_data$Fleet_code)
  fleets[mod$data_list$fleet_control$Fleet_type == "Off"] <- NA
  fleets[no_comp_data] <- NA

  testthat::expect_equal(as.numeric(mod$map$mapList$comp_weights), fleets)
})



