# Tests for various forms of logistic selectivity
# - Test correct value
# - Test correct mapping

testthat::test_that("Sex-specific age-based logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Selectivity <- 1 # Age-based logistic
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA

  # Specify logistic selectivity
  inf = 13; alpha = 0.2
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  inits <- ss_run$estimated_params
  inits$ln_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf     # Females
  inits$sel_inf[1,9:11,2] <- inf + 1 # Males

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits,
                      map = ss_run$map,
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      verbose = 1)
  )

  # Map
  # - Slope
  testthat::expect_equal(c(ss_run$map$mapList$ln_sel_slp[1,,]), c(1:6, NA, 7:8,10, 12, 14:18, rep(NA, 8), 9, 11, 13, rep(NA, 5)))
  # - Desc slope
  testthat::expect_equal(c(ss_run$map$mapList$ln_sel_slp[2,,]), as.numeric(rep(NA, 16*2)))

  # - Asympt
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf[1,,]), c(1:6, NA, 7:8,10, 12, 14:18, rep(NA, 8), 9, 11, 13, rep(NA, 5)))
  # - Desc Asympt
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf[2,,]), as.numeric(rep(NA, 16*2)))

  # - Devs
  testthat::expect_equal(c(ss_run$map$mapList$ln_sel_slp_dev), as.numeric(rep(NA, 2688)))
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf_dev), as.numeric(rep(NA, 2688)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_ln_sd), as.numeric(rep(NA, 16))) # Dev sigma


  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10], tolerance = 0.0001))

  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21], tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21], tolerance = 0.0001))
})


testthat::test_that("Sex-specific age-based time-varying logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Selectivity <- 1 # Age-based logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 10
  inf_dev <- rnorm(nyrs)
  ln_slp_dev <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(ln_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-alpha*exp(x[1]) * (ages - inf - x[2]))))
  sel2 <- apply(cbind(ln_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-(alpha+1)*exp(x[1]) * (ages - inf - x[2]))))


  # Set params to logistic
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  inits <- ss_run$estimated_params
  inits$ln_sel_slp[1,,1] <- log(alpha) # Females
  inits$ln_sel_slp[1,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$ln_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$ln_sel_slp_dev[1,,,])[2]){
      inits$ln_sel_slp_dev[1,i,j,] <- ln_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits,
                      map = ss_run$map,
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      verbose = 0)
  )

  # Map
  # - Slope
  testthat::expect_equal(c(!is.na(ss_run$map$mapList$ln_sel_slp[1,,])), !is.na(c(1:6, NA, 7:8,10, 12, 14:18, rep(NA, 8), 9, 11, 13, rep(NA, 5))))
  # - Desc slope
  testthat::expect_equal(c(ss_run$map$mapList$ln_sel_slp[2,,]), as.numeric(rep(NA, 16*2)))

  # - Asympt
  testthat::expect_equal(c(!is.na(ss_run$map$mapList$sel_inf[1,,])), !is.na(c(1:6, NA, 7:8,10, 12, 14:18, rep(NA, 8), 9, 11, 13, rep(NA, 5))))
  # - Desc Asympt
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf[2,,]), as.numeric(rep(NA, 16*2)))

  # - Devs
  # -- Ascending
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$ln_sel_slp_dev[1,,,]))), 18 * nyrs) # slope
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[1,,,]))), 18 * nyrs) # asymptote

  testthat::expect_equal(length(unique(c(ss_run$map$mapList$ln_sel_slp_dev[1,,,]))), 18 * nyrs + 1) # slope
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$sel_inf_dev[1,,,]))), 18 * nyrs + 1) # asymptote

  testthat::expect_equal(c(ss_run$map$mapList$ln_sel_slp_dev[2,,,]), as.numeric(rep(NA, 2688/2))) # Descending slope
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf_dev[2,,,]), as.numeric(rep(NA, 2688/2))) # Descending asymptote

  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_ln_sd), as.numeric(rep(NA, 16))) # Dev sigma turned off


  # Check selectivity
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[1,1,1:10,1:42]), as.numeric(sel[1:10,]), tolerance = 0.0001)

  # - ATF
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,1,1:21,1:42]), as.numeric(sel[1:21,]), tolerance = 0.0001) # Females
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,2,1:21,1:42]), as.numeric(sel2[1:21,]), tolerance = 0.0001) # Males
})


testthat::test_that("Sex-specific age-based time-varying logistic selectivity not normalized (devs as random effects)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  GOA2018SS$fleet_control$Selectivity <- 1 # Age-based logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = NULL,
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      verbose = 0)
  )

  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_ln_sd), c(1:6, NA, 8:16)) # Dev sigma turned on except for not estimated fleet

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  2 * 18 * nyrs)
})


testthat::test_that("Time-varying logistic selectivity likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Selectivity <- 1 # Age-based logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1 # Note that inf does 4*sd prior, should adapt to scale-invariant sd
  #FIXME: move to scale invariant setup
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 10
  inf_dev <- rnorm(nyrs)
  ln_slp_dev <- rnorm(nyrs)

  # Set params to logistic
  alpha = 0.5
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  inits <- ss_run$estimated_params
  inits$ln_sel_slp[1,,1] <- log(alpha) # Females
  inits$ln_sel_slp[1,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$ln_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$ln_sel_slp_dev[1,,,])[2]){
      inits$ln_sel_slp_dev[1,i,j,] <- ln_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits,
                      map = ss_run$map,
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      verbose = 0)
  )

  # Nll
  rcnll <- sum(ss_run$quantities$jnll_comp[6,])
  single_dv_nll <- sum(-dnorm(inf_dev, 0, 1, log = TRUE) - dnorm(ln_slp_dev, 0, 4, log = TRUE))

  testthat::expect_equal(rcnll, single_dv_nll * 18) # 18 selectivities
})


