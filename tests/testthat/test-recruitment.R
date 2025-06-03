test_that("mean recruitment and devs", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  inits <- build_params(GOA2018SS)

  # Specify logistic selectivity
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Set params
  inits$rec_pars[,1] <- R0
  inits$rec_dev[,1:nyrs] <- rep(Rdev, each = 3)
  inits$R_ln_sd <- rep(0, 3)

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Check R
  for(sp in 1:3){
    expect_equal(as.numeric(ss_run$quantities$R[sp,1:nyrs]), exp(R0[sp] + Rdev), tolerance = 0.0001)
  }

  # JNLL
  expect_equal(as.numeric(ss_run$quantities$jnll_comp[12,1:3]), rep(sum(-dnorm(Rdev,
                                                                               mean = 1/2, # lognormal bias correction to center around 0
                                                                               sd = 1, log = TRUE)), 3), tolerance = 0.0001)
})

test_that("ssb under mean recruitment", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify recruitment
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Set params
  inits <- build_params(GOA2018SS)
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- rep(0, 3)
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 2,
                              verbose = 1)
  # Check ssb
  for(yr in 1:nyrs){
    expect_equal(as.numeric(ss_run$quantities$ssb[,yr]),  exp(R0)*1/(1-exp(-0.2))*0.5, tolerance = 0.0001)
  }
})


test_that("ssb and ricker recruitment", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2

  # Specify logistic selectivity
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Set params
  GOA2018SS$srr_fun <- 4
  GOA2018SS$initMode <- 1
  inits <- build_params(GOA2018SS)
  alpha = 0.3
  beta = 1e-6
  inits$rec_pars[,2] <- alpha
  inits$rec_pars[,3] <- beta
  inits$ln_F[] <- -999 # No fishing
  curve(alpha*x*exp(-beta*x), from = 0, to = 1e7)

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode

                              recFun = build_srr(srr_fun = 4,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1
                              ),
                              initMode = 2,
                              verbose = 1)

  # Check ssb
  for(yr in 1:nyrs){
    expect_equal(as.numeric(ss_run$quantities$ssb[,yr]),  exp(R0)*1/(1-exp(-0.2))*0.5, tolerance = 0.0001)
  }

})

