
testthat::test_that("Estimated composition data matches expected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))
  nspp = 2
  ages = 1:15
  years = 1:30
  sim <- make_msm_test_data(
    nspp = nspp,
    ages = ages,
    years = years,
    sigma_R = 0, # No rec devs
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-2.5 * (ages - 4)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(0.2, nrow = nspp, ncol = length(years), byrow = TRUE),

    # Multispecies bits
    niter = 1,
    log_phi =  matrix(c(-Inf,-Inf,-Inf,-Inf), nspp, nspp, byrow = TRUE) # No species interactions
  )

  simData <- sim$data_list

  # Build parameter object
  inits <- suppressMessages(build_params(simData))
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(0.2)
  inits$ln_F[4,] <- log(0.2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$ln_M1[,1,] <- log(sim$model_quantities$M)

  # Run
  fit <- suppressMessages(
    Rceattle::fit_mod(data_list = simData,
                      inits = inits,
                      estimateMode = 3,
                      phase = FALSE,
                      loopnum = 1,
                      verbose = 0)
  )

  # Expected
  # * Comp data
  # - Index
  index_comp <- sim$model_quantities$NAA
  for(sp in 1:nspp){
    for(y in years){
      index_comp[sp, , y] <- index_comp[sp, , y] * sim$model_quantities$srv_sel[sp,]
    }
  }
  index_comp <- do.call(rbind, lapply(seq_len(dim(index_comp)[1]), function(i) t(index_comp[i, , ])))
  colnames(index_comp) <- paste0("Comp_", ages)

  # - Fishery
  fishery_comp <- sim$model_quantities$CAA
  fishery_comp <- do.call(rbind, lapply(seq_len(dim(fishery_comp)[1]), function(i) t(fishery_comp[i, , ])))
  colnames(fishery_comp) <- paste0("Comp_", ages)
  expected <- rbind(index_comp, fishery_comp)
  expected <- expected/ rowSums(expected)

  # Test
  testthat::expect_equal(c(fit$quantities$comp_hat), c(expected))
}
)
