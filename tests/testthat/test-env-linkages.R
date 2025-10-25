
test_that("test recruitment linkeage and mean rec", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs))

  # Set params
  GOA2018SS$srr_fun <- 1
  GOA2018SS$initMode <- 1
  inits <- build_params(GOA2018SS)
  alpha = 0.4
  beta = 1e-6
  inits$rec_pars[,1] <- R0
  inits$beta_rec_pars[,] <- 1:3
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode

                              recFun = build_srr(srr_fun = 1,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1,
                                                 srr_indices = 1
                              ),
                              initMode = 2,
                              verbose = 1)


  # Check ssb
  for(sp in 1:3){
    expect_equal(as.numeric(ss_run$quantities$R[sp,1:nyrs]),  exp(R0[sp] + GOA2018SS$env_data$EnvIndex * sp), tolerance = 0.0001)
  }

})




test_that("test multiple recruitment linkeages and mean rec", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs), EnvIndex2 = seq(0,1, length.out = nyrs), EnvIndex3 = seq(0,1, length.out = nyrs))

  # Set params
  GOA2018SS$srr_fun <- 1
  GOA2018SS$initMode <- 1
  inits <- build_params(GOA2018SS)
  alpha = 0.4
  beta = 1e-6
  inits$rec_pars[,1] <- R0
  inits$beta_rec_pars[,] <- 1:9
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode

                              recFun = build_srr(srr_fun = 1,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1,
                                                 srr_indices = c(1,2,3)
                              ),
                              initMode = 1,
                              verbose = 1)


  # Check ssb
  for(sp in 1:3){
    expect_equal(as.numeric(ss_run$quantities$R[sp,1:nyrs]),  as.numeric(exp(R0[sp] + as.matrix(GOA2018SS$env_data[,-1]) %*% inits$beta_rec_pars[sp,])), tolerance = 0.0001)
  }
})



test_that("test multiple M linkeages", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs), EnvIndex2 = seq(0,1, length.out = nyrs), EnvIndex3 = seq(0,1, length.out = nyrs))

  # Set params
  GOA2018SS$srr_fun <- 0
  GOA2018SS$M1_model <- 4
  GOA2018SS$initMode <- 1
  inits <- build_params(GOA2018SS)
  alpha = 0.4
  beta = 1e-6
  inits$ln_M1[] <- log(0.2)
  inits$M1_beta[] <- 1
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              M1Fun = build_M1(M1_model = 4,
                                               M1_indices = c(1,2,3)
                              ),
                              initMode = 1,
                              verbose = 1)



  expect_equal(as.numeric(ss_run$quantities$M_at_age),  as.numeric(ss_run$quantities$M1_at_age), tolerance = 0.0001)

  # Check M
  for(sp in 1:3){
    expect_equal(as.numeric(ss_run$quantities$M_at_age[sp,1,1,1:nyrs]),  as.numeric(0.2*exp(as.matrix(GOA2018SS$env_data[,-1]) %*% inits$M1_beta[sp,1,])), tolerance = 0.0001)
    expect_equal(as.numeric(ss_run$quantities$M_at_age[sp,1,5,1:nyrs]),  as.numeric(0.2*exp(as.matrix(GOA2018SS$env_data[,-1]) %*% inits$M1_beta[sp,1,])), tolerance = 0.0001)
  }
})



test_that("test single M, multiple M/sex linkeages, M both-sex linkage", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs), EnvIndex2 = seq(0,1, length.out = nyrs), EnvIndex3 = seq(0,1, length.out = nyrs))

  # Set params
  GOA2018SS$srr_fun <- 0
  GOA2018SS$M1_model <- c(1,5,4)
  GOA2018SS$initMode <- 1
  inits <- build_params(GOA2018SS)
  alpha = 0.4
  beta = 1e-6
  inits$ln_M1[] <- log(0.2)
  inits$M1_beta[,1,] <- 1
  inits$M1_beta[,2,] <- 0.5
  inits$M1_beta[1,,] <- 0
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              M1Fun = build_M1(M1_model = c(1,5,4),
                                               M1_indices = c(1,2,3)
                              ),
                              initMode = 1,
                              verbose = 1)



  expect_equal(as.numeric(ss_run$quantities$M_at_age),  as.numeric(ss_run$quantities$M1_at_age), tolerance = 0.0001)

  # Check M
  # - Pollock fixed M
  expect_equal(as.numeric(ss_run$quantities$M_at_age[1,1,1:10,]),  rep(0.2, length(ss_run$quantities$M_at_age[1,1,1:10,])), tolerance = 0.0001)

  # - ATF sex-time env varying
  expect_equal(as.numeric(ss_run$quantities$M_at_age[2,1,5,1:nyrs]),  as.numeric(0.2*exp(as.matrix(GOA2018SS$env_data[,-1]) %*% inits$M1_beta[2,1,])), tolerance = 0.0001)
  expect_equal(as.numeric(ss_run$quantities$M_at_age[2,2,2,1:nyrs]),  as.numeric(0.2*exp(as.matrix(GOA2018SS$env_data[,-1]) %*% inits$M1_beta[2,2,])), tolerance = 0.0001)



  # - Cod env varying
  expect_equal(as.numeric(ss_run$quantities$M_at_age[3,1,5,1:nyrs]),  as.numeric(0.2*exp(as.matrix(GOA2018SS$env_data[,-1]) %*% inits$M1_beta[3,1,])), tolerance = 0.0001)
  expect_equal(as.numeric(ss_run$quantities$M_at_age[3,2,2,1:nyrs]),  rep(0, length(ss_run$quantities$M_at_age[3,2,2,1:nyrs])), tolerance = 0.0001)
})


