# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test no multi-species data", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  M1rho = 0.9
  M1sd = 0.1
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    log_phi = log_phi # set to -Inf so no predation
  )

  # sum(sim$model_quantities$B_eaten_as_prey)

  # Set up Rceattle data
  simData <- sim$base_data


  testthat::expect_no_error(fit_mod(data_list = simData,
                                              inits = NULL,
                                              estimateMode = 0,
                                              random_rec = FALSE,
                                              msmMode = 0,
                                              initMode = "NonEquilibrium",
                                              fit_control = fit_control(phase = FALSE, verbose = 0))
  )
}
)
