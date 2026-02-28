testthat::test_that("Basic index and index likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 3.5
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  # Set params
  inits <- build_params(dat)
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- 0
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1

  # Set logistic params
  inits$ln_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Calculate SPR
  M <- dat$M1_base$Age1
  wt <- dat$weight$Age1
  sex_ratio <- dat$sex_ratio$Age1
  pmature <- dat$maturity$Age1

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0)

  # Check derived index
  testthat::expect_equal(ss_run$quantities$index_hat, 0.5 * SB0[ss_run$data_list$index_data$Species])

  # Check index sd
  testthat::expect_equal(as.numeric(ss_run$quantities$ln_index_sd), rep(0.1, nyrs))

  # Check observed
  testthat::expect_equal(dat$index_data$Observation, rep(100, nyrs))

  # Check observed in model
  rearranged_check <- rearrange_dat(dat)
  testthat::expect_equal(rearranged_check$index_obs[,1], rep(100, nyrs))

  # Check Loglikelood
  exp_calc <- -sum(dnorm(log(dat$index_data$Observation), log(0.5 * SB0) - (0.1^2)/2, 0.1, log = TRUE))
  exp_mod <- ss_run$quantities$jnll_comp[1,1]
  testthat::expect_equal(exp_mod, exp_calc)

})
