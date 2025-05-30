test_that("M1 input and output the same", {
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
  GOA2018SS$fleet_control$Age_max_selected <- NA

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
  apply(ss_run$quantities$sel[1,1,1:10,], 2, function(x) expect_equal(as.numeric(x), sel[1:10]/max(sel[1:10]), tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel[9,1,1:21,], 2, function(x) expect_equal(as.numeric(x), sel[1:21]/max(sel[1:21]), tolerance = 0.0001))
  apply(ss_run$quantities$sel[9,2,1:21,], 2, function(x) expect_equal(as.numeric(x), sel[1:21]/max(sel[1:21]), tolerance = 0.0001))
})
