testthat::skip_on_cran()
testthat::test_that("Fit sanity model: key quantities match baseline", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 1)

  # Ensure compiled TMB exists
  compile_tmb_if_needed()

  # Run a short/safe fit
  fit <- Rceattle::fit_mod(data_list = dat,
                           inits = NULL,
                           estimateMode = 0,
                           phase = FALSE,
                           loopnum = 1,
                           getsd = FALSE,
                           use_gradient = FALSE,
                           control = list(eval.max = 200, iter.max = 200),
                           verbose = 0)

  testthat::expect_s3_class(fit, "Rceattle")
  testthat::expect_true(!is.null(fit$quantities))

  # Collect a small set of key quantities to compare
  got <- list(
    R0 = as.numeric(fit$quantities$R0),
    ssb1 = as.numeric(fit$quantities$ssb[1,]),
    biomass1 = as.numeric(fit$quantities$biomass[1,])
  )

  baseline_path <- file.path("tests", "testthat", "fixtures", "fit_baseline.rds")

  # Regenerate baseline explicitly when requested by env var
  if (identical(tolower(Sys.getenv("REGEN_TEST_BASELINES")), "true")) {
    dir.create(dirname(baseline_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(got, baseline_path)
    testthat::skip("Baseline regenerated; rerun test without REGEN_TEST_BASELINES to validate.")
  }

  testthat::expect_true(file.exists(baseline_path), info = paste0("Baseline file missing: ", baseline_path, ". To create it run with REGEN_TEST_BASELINES=true."))
  expected <- readRDS(baseline_path)

  # Compare shapes
  testthat::expect_equal(length(got$R0), length(expected$R0))
  testthat::expect_equal(length(got$ssb1), length(expected$ssb1))
  testthat::expect_equal(length(got$biomass1), length(expected$biomass1))

  # Compare numeric values with tolerant thresholds
  testthat::expect_equal(got$R0, expected$R0, tolerance = 1e-6)
  testthat::expect_equal(got$ssb1, expected$ssb1, tolerance = 1e-6)
  testthat::expect_equal(got$biomass1, expected$biomass1, tolerance = 1e-6)
})


test_that("Index, biomass, and catch = 0 match expected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  # Set params
  inits <- build_params(dat)
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- 0
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1


  # Specify logistic selectivity
  inf = 0; alpha = 0
  ages <- 1:5
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)


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

  # Check biomass
  testthat::expect_equal(as.numeric(ss_run$quantities$biomass[1,1:nyrs]), rep(SB0, nyrs), tolerance = 0.0001)

  # Check derived index
  testthat::expect_equal(ss_run$quantities$index_hat, 0.5 * SB0[ss_run$data_list$index_data$Species])

  # Check catch
  testthat::expect_equal(ss_run$quantities$catch_hat[1:nyrs], rep(0, nyrs))
})
