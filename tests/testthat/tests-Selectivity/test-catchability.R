testthat::test_that("Fixed catchability", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Catchability <- 0 # Fixed

  # Set params
  inits <- suppressMessages( build_params(GOA2018SS) )
  inits$index_ln_q[] <- 0

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      verbose = 1)
  )

  # Map
  nflt <- nrow(GOA2018SS$fleet_control)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_ln_q), as.numeric(rep(NA, nflt)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_ln_q), as.numeric(rep(NA, nflt)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_beta), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_beta))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev_ln_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev_ln_sd))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_ln_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_ln_sd))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_rho), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_rho))))

  # Check q
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$index_q), as.numeric(rep(1, length(ss_run$quantities$index_q))))
})

testthat::test_that("Estimated catchability", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nflt <- nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Catchability <- 1 # Estimated
  GOA2018SS$fleet_control$Q_index <- 1:nflt

  # Set params
  inits <- suppressMessages( build_params(GOA2018SS) )
  inits$index_ln_q[] <- 0

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      verbose = 1)
  )

  # Map
  fleets <- 1:nflt
  fleets[GOA2018SS$fleet_control$Fleet_type != 2] <- NA
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_ln_q), fleets)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_beta), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_beta))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev_ln_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev_ln_sd))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_ln_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_ln_sd))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_rho), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_rho))))

  # Check q
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$index_q), as.numeric(rep(1, length(ss_run$quantities$index_q))))
})

testthat::test_that("Invalid catchability", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nflt <- nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Catchability <- 9 # Not in scope
  GOA2018SS$fleet_control$Q_index <- 1:nflt

  # Set params
  inits <- suppressMessages( build_params(GOA2018SS) )
  inits$index_ln_q[] <- 0

  # Run
  testthat::expect_error(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      verbose = 0)
  )

})

