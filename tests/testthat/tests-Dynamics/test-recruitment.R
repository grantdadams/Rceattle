testthat::test_that("mean recruitment and devs", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nyrs = 20
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)
  R0 = 11
  Rdev <- rnorm(nyrs)

  # Set params
  dat$srr_fun = 0 # Set to mean R plus devs

  # inits <- Rceattle::build_params(dat)
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  inits$rec_pars[1,1] <- R0
  inits$x_tj[1:nyrs,1] <- Rdev
  inits$R_ln_sd <- 0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters from input
                              file = NULL, # Don't save
                              map = map,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  # Check R
  testthat::expect_equal(as.numeric(ss_run$quantities$R[1,1:nyrs]), exp(R0 + Rdev), tolerance = 0.0001)

  # JNLL
  # testthat::expect_equal(as.numeric(ss_run$quantities$jnll_dsem),
  #                        sum(-dnorm(Rdev,
  #                                   mean = 0, # FIXME: Doesn't match, and no lognormal bias correction to center around 0
  #                                   sd = 1, log = TRUE)), tolerance = 0.0001)
})

testthat::test_that("ssb under mean recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
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
  # inits <- suppressMessages(build_params(GOA2018SS))
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 2,
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- rep(0, 3)
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters from input
                              file = NULL, # Don't save
                              map = map,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 2,
                              verbose = 0)
  # Check ssb
  for(yr in 1:nyrs){
    testthat::expect_equal(as.numeric(ss_run$quantities$ssb[,yr]),  exp(R0)*1/(1-exp(-0.2))*0.5, tolerance = 0.0001)
  }
})


testthat::test_that("ssb and beverton recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify dims
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)

  # Set params
  GOA2018SS$srr_fun <- 4
  GOA2018SS$initMode <- 1
  # inits <- suppressMessages(build_params(GOA2018SS))
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              recFun = build_srr(srr_fun = 2,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1
                              ),
                              initMode = 2,
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  alpha = 0.4
  beta = 1e-6
  inits$rec_pars[,2] <- log(alpha)
  inits$rec_pars[,3] <- log(beta)
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters from input
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode

                              recFun = build_srr(srr_fun = 2,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1
                              ),
                              initMode = 2,
                              verbose = 0)
  # Calculate SPR
  M <- 0.2
  wt <- 1
  sex_ratio <- 0.5
  pmature <- 1

  natage = rep(1, 10)
  for(age in 2:9){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[10] <-  natage[9] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  ssb_at_age <- natage * wt * pmature * sex_ratio
  SPR0 <- sum(ssb_at_age)
  R0 <- (alpha-1/SPR0)/beta

  # Check ssb
  for(yr in 1:nyrs){
    for(sp in 1:ss_run$data_list$nspp){
      testthat::expect_equal(as.numeric(ss_run$quantities$ssb[sp,yr]),  R0/(1-exp(-0.2))*0.5, tolerance = 0.0001)
    }
  }

})


testthat::test_that("ssb and ricker recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)
  GOA2018SS$srr_fun <- 4
  GOA2018SS$initMode <- 1

  # Specify dims
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)

  # Set params
  # inits <- suppressMessages(build_params(GOA2018SS))
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              estimateMode = 3, # Don't estimate
                              msmMode = 0, # Single species mode
                              recFun = build_srr(srr_fun = 4,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1
                              ),
                              initMode = 2,
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  alpha = 0.4
  beta = 1e-6
  inits$rec_pars[,2] <- log(alpha)
  inits$rec_pars[,3] <- log(beta)
  inits$ln_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters from input
                              file = NULL, # Don't save
                              map = map,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              recFun = build_srr(srr_fun = 4,
                                                 proj_mean_rec = FALSE,
                                                 srr_est_mode = 1
                              ),
                              initMode = 2,
                              verbose = 0)
  # Calculate SPR
  M <- 0.2
  wt <- 1
  sex_ratio <- 0.5
  pmature <- 1

  natage = rep(1, 10)
  for(age in 2:9){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[10] <-  natage[9] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  ssb_at_age <- natage * wt * pmature * sex_ratio
  SPR0 <- sum(ssb_at_age)
  R0 <- log(alpha * SPR0) / (beta * SPR0) * 1e6 # Times 1e6 because beta is divided by 1e6 in C++

  # Check SSB
  for(yr in 1:nyrs){
    for(sp in 1:ss_run$data_list$nspp){
      testthat::expect_equal(as.numeric(ss_run$quantities$ssb[sp,yr]),  R0/(1-exp(-0.2))*0.5, tolerance = 0.0001)
    }
  }

})

