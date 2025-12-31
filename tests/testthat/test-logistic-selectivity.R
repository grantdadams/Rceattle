test_that("logistic selectivity divided by max sel", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  inits <- build_params(GOA2018SS)

  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  inits$ln_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf     # Females
  inits$sel_inf[1,9:11,2] <- inf + 1 # Males
  GOA2018SS$fleet_control$Selectivity <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- -999

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/max(sel[1:10]), tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/max(c(sel, sel2)), tolerance = 0.0001))
  apply(ss_run$quantities$sel[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21]/max(c(sel, sel2)), tolerance = 0.0001))
})


test_that("logistic selectivity not normalized", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  inits <- build_params(GOA2018SS)

  # Specify logistic selectivity
  inf = 13; alpha = 0.2
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  inits$ln_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf     # Females
  inits$sel_inf[1,9:11,2] <- inf + 1 # Males
  GOA2018SS$fleet_control$Selectivity <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10], tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21], tolerance = 0.0001))
  apply(ss_run$quantities$sel[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21], tolerance = 0.0001))
})


test_that("logistic selectivity divided by sel-at-age", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  inits <- build_params(GOA2018SS)

  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  inits$ln_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf
  GOA2018SS$fleet_control$Selectivity <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- 7

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/sel[7], tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/sel[7], tolerance = 0.0001))
  apply(ss_run$quantities$sel[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/sel[7], tolerance = 0.0001))
})



test_that("logistic selectivity divided by sel-at-age-RANGE", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  inits <- build_params(GOA2018SS)

  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  inits$ln_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf
  GOA2018SS$fleet_control$Selectivity <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- 7
  GOA2018SS$fleet_control$Sel_norm_bin2 <- 9

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/mean(sel[7:9]), tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/mean(sel[7:9]), tolerance = 0.0001))
  apply(ss_run$quantities$sel[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/mean(sel[7:9]), tolerance = 0.0001))
})



test_that("time-varying logistic selectivity divided by max sel", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  inits <- build_params(GOA2018SS)

  # Specify logistic selectivity
  nyrs <- GOA2018SS$styr:GOA2018SS$endyr
  inf = 10
  inf_dev <- rnorm(nyrs)
  ln_slp_dev <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(ln_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-alpha*exp(x[1]) * (ages - inf - x[2]))))

  # Set params to logistic
  inits$ln_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$ln_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$ln_sel_slp_dev[1,,,])[2]){
      inits$ln_sel_slp_dev[1,i,j,] <- ln_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
    }
  }
  GOA2018SS$fleet_control$Selectivity <- 1
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Sel_sd_prior <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1 <- -999

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Check selectivity
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$sel[1,1,1:10,1:42]), as.numeric(apply(sel[1:10,], 2, function(x) x/max(x))), tolerance = 0.0001)


  # - ATF
  testthat::expect_equal(as.numeric(ss_run$quantities$sel[9,1,1:21,1:42]), as.numeric(apply(sel[1:21,], 2, function(x) x/max(x))), tolerance = 0.0001)
  testthat::expect_equal(as.numeric(ss_run$quantities$sel[9,2,1:21,1:42]), as.numeric(apply(sel[1:21,], 2, function(x) x/max(x))), tolerance = 0.0001)
})

