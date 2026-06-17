# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test multiple M linkeages", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

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
  Rdev <- stats::rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = yrs, EnvIndex = seq(0,1, length.out = nyrs), EnvIndex2 = seq(0,1, length.out = nyrs), EnvIndex3 = seq(0,1, length.out = nyrs))

  # Set params
  m1_spec <- Rceattle::build_M1(
    M1_model = "sex_age_invariant",
    linkages = list(
      M1 = Rceattle::linkage_spec(
        formula = ~ 0 + EnvIndex + EnvIndex2 + EnvIndex3,
        by = ~ species
      )
    )
  )
  mod0 <- suppressMessages(
    fit_mod(data_list = GOA2018SS,
            inits = NULL,
            estimateMode = 3,
            random_rec = FALSE,
            msmMode = 0,
            M1Fun = m1_spec,
            initMode = 1,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$log_M1[] <- log(0.2)
  tbl <- mod0$data_list$linkage_table
  is_M1 <- tbl$param == "M1"
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  inits$beta_linkage[is_M1] <- 1
  inits$log_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              M1Fun = m1_spec,
                              initMode = 1,
                              fit_control = fit_control(
                              verbose = 0))



  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age),  as.numeric(ss_run$quantities$M1_at_age), tolerance = 0.0001)

  # Check M
  for(sp in 1:3){
    testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[sp,1,1,1:nyrs]),
                           as.numeric(0.2 * exp(as.matrix(GOA2018SS$env_data[,-1]) %*% rep(1, 3))),
                           tolerance = 0.0001)
    testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[sp,1,5,1:nyrs]),
                           as.numeric(0.2 * exp(as.matrix(GOA2018SS$env_data[,-1]) %*% rep(1, 3))),
                           tolerance = 0.0001)
  }
})


testthat::test_that("Test single M, multiple M/sex linkeages, M both-sex linkage (hindcast and projection)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)


  # Specify R0
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  projyrs <- GOA2018SS$styr:GOA2018SS$projyr
  nyrs <- length(yrs)
  nyrsproj <- length(projyrs)
  R0 = 10:12
  Rdev <- stats::rnorm(nyrs)

  # Env data
  GOA2018SS$env_data <- data.frame(Year = projyrs,
                                   EnvIndex = seq(0,1, length.out = nyrsproj),
                                   EnvIndex2 = seq(0,1, length.out = nyrsproj),
                                   EnvIndex3 = seq(0,1, length.out = nyrsproj))

  # Set params
  m1_spec <- Rceattle::build_M1(
    M1_model = c("sex_age_invariant", "sex_specific", "sex_age_invariant"),
    linkages = list(
      M1 = list(
        Rceattle::linkage_spec(
          formula = ~ 0 + EnvIndex + EnvIndex2 + EnvIndex3,
          by = ~ species + sex,
          species = 2
        ),
        Rceattle::linkage_spec(
          formula = ~ 0 + EnvIndex + EnvIndex2,
          by = ~ species,
          species = 3
        )
      )
    )
  )

  # Inits
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              M1Fun = m1_spec,
                              initMode = 1,
                              fit_control = fit_control(
                                verbose = 0))

  inits <- ss_run$estimated_params
  inits$log_M1[] <- log(0.2)
  tbl <- ss_run$data_list$linkage_table
  env_cols <- c("EnvIndex", "EnvIndex2", "EnvIndex3")
  is_M1 <- tbl$param == "M1"
  inits$beta_linkage <- as.numeric(inits$beta_linkage)
  for (j in seq_along(env_cols)) {
    idx <- which(is_M1 & tbl$species == 2 & tbl$sex == 1 & tbl$design_col == env_cols[j])
    inits$beta_linkage[idx] <- 1
    idx <- which(is_M1 & tbl$species == 2 & tbl$sex == 2 & tbl$design_col == env_cols[j])
    inits$beta_linkage[idx] <- 0.5
    idx <- which(is_M1 & tbl$species == 3 & is.na(tbl$sex) & tbl$design_col == env_cols[j])
    inits$beta_linkage[idx] <- 1
  }
  inits$log_F[] <- -999 # No fishing

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              M1Fun = m1_spec,
                              initMode = 1,
                              fit_control = fit_control(
                              verbose = 0))



  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age),  as.numeric(ss_run$quantities$M1_at_age), tolerance = 0.0001)

  # Check M
  # - Pollock fixed M
  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[1,1,1:10,]),  rep(0.2, length(ss_run$quantities$M_at_age[1,1,1:10,])), tolerance = 0.0001)

  # - ATF sex-time env varying
  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[2,1,5,1:nyrsproj]),
                         as.numeric(0.2 * exp(as.matrix(GOA2018SS$env_data[,-1]) %*% rep(1, 3))),
                         tolerance = 0.0001)
  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[2,2,2,1:nyrsproj]),
                         as.numeric(0.2 * exp(as.matrix(GOA2018SS$env_data[,-1]) %*% rep(0.5, 3))),
                         tolerance = 0.0001)

  # - Cod env varying
  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[3,1,5,1:nyrsproj]),
                         as.numeric(0.2 * exp(as.matrix(GOA2018SS$env_data[,2:3]) %*% rep(1, 2))),
                         tolerance = 0.0001)
  testthat::expect_equal(as.numeric(ss_run$quantities$M_at_age[3,2,2,1:nyrsproj]),
                         rep(0, length(ss_run$quantities$M_at_age[3,2,2,1:nyrsproj])),
                         tolerance = 0.0001)
})
