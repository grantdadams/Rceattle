testthat::test_that("Invalid initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  # Run
  testthat::expect_error( Rceattle::fit_mod(data_list = dat,
                                              initMode = 5,
                                              estimateMode = 3)
  )

  # Run
  testthat::expect_error( Rceattle::fit_mod(data_list = dat,
                                            initMode = "wrong",
                                            estimateMode = 3)
  )
})

testthat::test_that("Free parameter initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 0,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9), rep(1, nages-1)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "FreeParams",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9), rep(1, nages-1)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))
})

testthat::test_that("Equilibrium initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 1,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), as.numeric(rep(NA, length(mod$map$mapList$init_dev))))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "Equilibrium",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), as.numeric(rep(NA, length(mod$map$mapList$init_dev))))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))
})


testthat::test_that("NonEquilibrium initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 2,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "NonEquilibrium",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1
})


testthat::test_that("FishedNonEquilibrium initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 3,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9-(0.2+Finit) * 0:3), exp(9) * exp(-(0.2+Finit)*4)/(1-exp(-0.2 - Finit))))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "FishedNonEquilibrium",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9-(0.2+Finit) * 0:3), exp(9) * exp(-(0.2+Finit)*4)/(1-exp(-0.2 - Finit))))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1
})


testthat::test_that("NonEquilibriumScaled initMode (=4)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 4,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9), exp(9-(0.2) * 1:3)* exp(-Finit), exp(9) * exp(-(0.2)*4)/(1-exp(-0.2 - Finit))* exp(-Finit)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1


  # Run string (renamed; previously "FishedNonEquilibriumScaled")
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "NonEquilibriumScaled",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9), exp(9-(0.2) * 1:3)* exp(-Finit), exp(9) * exp(-(0.2)*4)/(1-exp(-0.2 - Finit))* exp(-Finit)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1
})


testthat::test_that("EquilibriumScaled initMode (=5): init_dev mapped out, Finit alive", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Mode 5 keeps Finit alive but maps init_dev OFF (mirrors SS3 SR_regime)
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 5,
                           estimateMode = 3)

  # init_dev mapped entirely to NA
  testthat::expect_true(all(is.na(mod$map$mapList$init_dev)))
  # log_Finit kept estimated (= integer 1)
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # N_init formula reduces to pure regime equilibrium:
  #   N_init[a] = R_init * exp(-Finit) * exp(-sum(M1[0..a-1]))
  # with plus-group geometric correction. Same closed form as mode 4 with
  # init_dev = 0 (because mode 5 maps init_dev out, cpp reads 0).
  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(
    as.numeric(natage1),
    c(exp(9),
      exp(9 - (0.2) * 1:3) * exp(-Finit),
      exp(9) * exp(-(0.2) * 4) / (1 - exp(-0.2 - Finit)) * exp(-Finit))
  )

  # init_dev NLL slot should be zero in mode 5 (cpp gate skips the penalty)
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10, 1]), 0)

  # Run string
  mod2 <- Rceattle::fit_mod(data_list = dat,
                            initMode = "EquilibriumScaled",
                            estimateMode = 3)
  testthat::expect_true(all(is.na(mod2$map$mapList$init_dev)))
  testthat::expect_equal(as.numeric(mod2$map$mapList$log_Finit), 1)
})
