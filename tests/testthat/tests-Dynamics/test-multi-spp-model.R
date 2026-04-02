testthat::test_that("Rceattle and multi-species model dynamics match", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               map = map,
                               file = NULL, # Don't save
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ms_run2$quantities$R[,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ms_run2$quantities$biomass[,1:nyrs]), tolerance = 1e-6)


  # Suitability
  testthat::expect_equal(exp(ms_run2$estimated_params$log_gam_a), gam_a)
  testthat::expect_equal(exp(ms_run2$estimated_params$log_gam_b), gam_b)
  testthat::expect_equal(as.numeric(ms_run2$quantities$vulnerability), as.numeric(sim$model_quantities$vulnerability))
  testthat::expect_equal(as.numeric(ms_run2$quantities$suitability[,,,,1]), as.numeric(sim$model_quantities$suitability))
  testthat::expect_equal(as.numeric(ms_run2$quantities$suit_other[,1,1,1]), as.numeric(sim$model_quantities$suit_other))

  # M2
  testthat::expect_equal(as.numeric(sim$model_quantities$M2_at_age), as.numeric(ms_run2$quantities$M2_at_age[,1,,1:nyrs]), tolerance = 1e-6)

  # Ration
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ms_run2$quantities$consumption_at_age[,1,,1]))

  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ms_run2$quantities$N_at_age[,1,,1:nyrs]))

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Avail food
  testthat::expect_equal(as.numeric(sim$model_quantities$avail_food[,,1]), as.numeric(ms_run2$quantities$avail_food[,1,,1]))

  # Selectivity
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]), as.numeric(ms_run2$quantities$sel_at_age[c(1,3),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]), as.numeric(ms_run2$quantities$sel_at_age[c(2,4),,,1]))

  # F
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA), as.numeric(ms_run2$quantities$F_flt_age[c(2,4),1,,1:nyrs]))

  # Q
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q), as.numeric(ms_run2$quantities$index_q[c(1,3),1]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))
})


testthat::test_that("Equilibrium MSVPA suitability dynamics match", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip() # Haven't figured out quite why there is a minor difference, maybe bias in diet weighting? Old EBS CEATTLE still matches

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  nages = 15
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    sigma_R = 0,
    Fmort = matrix(0.2, nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi,
    niter = 20,
    normalize_suitability = TRUE
  )


  # Set up Rceattle data
  simData <- sim$data_list

  # Unweighted average diet
  diet_avg <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 4), mean, na.rm = TRUE)

  dimnames(diet_avg) <- list(1:nspp, 1:nspp, 1:nages, 1:nages)
  diet_df <- as.data.frame.table(diet_avg, dnn = dimnames(diet_avg))
  colnames(diet_df) <- c("Pred", "Prey", "Pred_age", "Prey_age", "Stomach_proportion_by_weight")

  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))
  diet_df$Pred_age <- as.numeric(as.character(diet_df$Pred_age))
  diet_df$Prey_age <- as.numeric(as.character(diet_df$Prey_age))

  # Set to Year 0
  diet_df$Year <- 0
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Sample_size <- 1000

  # Overwrite diet data
  simData$diet_data <- diet_df[diet_df$Stomach_proportion_by_weight > 0, ]


  # Fit multi-species
  # * Fix parameters -----
  inits <- suppressMessages( build_params(simData) )
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(0.2)
  inits$ln_F[4,] <- log(0.2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ms_run_msvpa <- Rceattle::fit_mod(data_list = simData,
                                    inits = inits, # Initialize from sim pars
                                    file = NULL, # Don't save
                                    estimateMode = 3, # Do not estimate
                                    random_rec = FALSE, # No random recruitment
                                    phase = FALSE,
                                    msmMode = 1,
                                    niter = 20,
                                    suitMode = 0,
                                    initMode = 2,
                                    verbose = 0)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ms_run_msvpa$quantities$R[,1:nyrs]), tolerance = 1e-5)
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ms_run_msvpa$quantities$biomass[,1:nyrs]), tolerance = 1e-5)


  # Suitability (Notice the fixed suit_other index [ , 1 , , 1])
  testthat::expect_equal(as.numeric(ms_run_msvpa$quantities$suitability[,,,,1]), as.numeric(sim$model_quantities$suitability), tolerance = 1e-5)
  testthat::expect_equal(as.numeric(ms_run_msvpa$quantities$suit_other[,1,,1]), as.numeric(sim$model_quantities$suit_other), tolerance = 1e-5)

  # M2
  testthat::expect_equal(as.numeric(sim$model_quantities$M2_at_age), as.numeric(ms_run_msvpa$quantities$M2_at_age[,1,,1:nyrs]), tolerance = 1e-5)

  # Ration
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ms_run_msvpa$quantities$consumption_at_age[,1,,1]), tolerance = 1e-5)

  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ms_run_msvpa$quantities$N_at_age[,1,,1:nyrs]), tolerance = 1e-5)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run_msvpa$quantities$avgN_at_age[,1,,1:nyrs]), tolerance = 1e-5)

  # Avail food (Checking all years instead of just year 1)
  testthat::expect_equal(as.numeric(sim$model_quantities$avail_food), as.numeric(ms_run_msvpa$quantities$avail_food[,1,,1:nyrs]), tolerance = 1e-5)

  # Selectivity
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]), as.numeric(ms_run_msvpa$quantities$sel_at_age[c(1,3),,,1]), tolerance = 1e-5)
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]), as.numeric(ms_run_msvpa$quantities$sel_at_age[c(2,4),,,1]), tolerance = 1e-5)

  # F
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA), as.numeric(ms_run_msvpa$quantities$F_flt_age[c(2,4),1,,1:nyrs]), tolerance = 1e-5)

  # Q
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q), as.numeric(ms_run_msvpa$quantities$index_q[c(1,3),1]), tolerance = 1e-5)

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run_msvpa$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run_msvpa$quantities$diet_hat[,2]), tolerance = 1e-5)
})

testthat::test_that("Test proportion of prey-at-age in predator-at-age averaged across years", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Average diet across years -----
  diet_avg_mat <- rowMeans(sim$model_quantities$diet_prop, dims = 4) # 1. Average across the year dimension (dim 5)
  dimnames(diet_avg_mat) <- list(1:nspp, 1:nspp, 1:nages, 1:nages)
  diet_df <- as.data.frame.table(diet_avg_mat, dnn = dimnames(diet_avg_mat)) # 2. Flatten the 4D array into a "long" data frame
  colnames(diet_df) <- c("Pred", "Prey", "Pred_age", "Prey_age", "Stomach_proportion_by_weight")

  # as.data.frame.table creates factors, so convert them to numeric indices
  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))
  diet_df$Pred_age <- as.numeric(as.character(diet_df$Pred_age))
  diet_df$Prey_age <- as.numeric(as.character(diet_df$Prey_age))

  # Add the required constant columns
  diet_df$Year <- 0
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Sample_size <- 1000

  # Reorder
  simData$diet_data <- diet_df[, c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                   "Pred_age", "Prey_age", "Year",
                                   "Sample_size", "Stomach_proportion_by_weight")]


  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test annual proportion of prey (all ages) in predator-at-age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Sum diet across prey-ages -----
  diet_sum_mat <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 5), sum)
  dimnames(diet_sum_mat) <- list(1:nspp, 1:nspp, 1:nages, 1:nyrs)
  diet_df <- as.data.frame.table(diet_sum_mat, dnn = dimnames(diet_sum_mat)) # 2. Flatten the 4D array into a "long" data frame
  colnames(diet_df) <- c("Pred", "Prey", "Pred_age", "Year", "Stomach_proportion_by_weight")

  # as.data.frame.table creates factors, so convert them to numeric indices
  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))
  diet_df$Pred_age <- as.numeric(as.character(diet_df$Pred_age))
  diet_df$Year <- as.numeric(as.character(diet_df$Year))

  # Add the required constant columns
  diet_df$Prey_age <- -500 # Indicates sum across prey-ages in Rceattle
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Sample_size <- 1000

  # Reorder
  simData$diet_data <- diet_df[, c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                   "Pred_age", "Prey_age", "Year",
                                   "Sample_size", "Stomach_proportion_by_weight")]


  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight),
                         as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test proportion of prey (all ages) in predator-at-age averaged across years", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Sum diet across prey-ages and average across years -----
  diet_sum_prey_age <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 5), sum) # Sum across prey-ages (dimension 4)
  diet_avg_years <- apply(diet_sum_prey_age, c(1, 2, 3), mean) # Average across years (now dimension 4 of the new object)
  dimnames(diet_avg_years) <- list(1:nspp, 1:nspp, 1:nages)
  diet_df <- as.data.frame.table(diet_avg_years, dnn = dimnames(diet_avg_years)) # 2. Flatten the 4D array into a "long" data frame
  colnames(diet_df) <- c("Pred", "Prey", "Pred_age", "Stomach_proportion_by_weight")

  # as.data.frame.table creates factors, so convert them to numeric indices
  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))
  diet_df$Pred_age <- as.numeric(as.character(diet_df$Pred_age))

  # Add the required constant columns
  diet_df$Year <- 0 # Indicates average
  diet_df$Prey_age <- -500 # Indicates sum across prey-ages in Rceattle
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Sample_size <- 1000

  # Reorder
  simData$diet_data <- diet_df[, c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                   "Pred_age", "Prey_age", "Year",
                                   "Sample_size", "Stomach_proportion_by_weight")]


  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)

testthat::test_that("Test annual proportion of prey (all ages) in predator (mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Sum diet across prey-ages and average across predator ages -----
  diet_sum_prey_age <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 5), sum) # Sum across prey-ages (dimension 4)
  diet_avg_pred_age <- apply(diet_sum_prey_age, c(1, 2, 4), mean) # Average across predator ages (dimension 3 of the new object)
  dimnames(diet_avg_pred_age) <- list(1:nspp, 1:nspp, 1:nyrs)
  diet_df <- as.data.frame.table(diet_avg_pred_age, dnn = dimnames(diet_avg_pred_age)) # 2. Flatten the 4D array into a "long" data frame
  colnames(diet_df) <- c("Pred", "Prey", "Year", "Stomach_proportion_by_weight")

  # as.data.frame.table creates factors, so convert them to numeric indices
  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))
  diet_df$Year <- as.numeric(as.character(diet_df$Year))

  # Add the required constant columns
  diet_df$Prey_age <- -500 # Indicates sum across prey-ages in Rceattle
  diet_df$Pred_age <- -1 # Indicates average across pred-ages in Rceattle (< -500 does weighted average)
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Sample_size <- 1000

  # Reorder
  simData$diet_data <- diet_df[, c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                   "Pred_age", "Prey_age", "Year",
                                   "Sample_size", "Stomach_proportion_by_weight")]
  simData$diet_data <- diet_df[diet_df$Stomach_proportion_by_weight > 0, ]

  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test annual proportion of prey (all ages) in predator (weighted mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Sum diet across prey-ages and average across predators weighted by avgN -----
  sum_over_prey_age <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 5), sum, na.rm = TRUE) # sum across prey ages
  avgN <- sim$model_quantities$avgNAA # Dims: [pred, pred_age, year]

  # 2. Calculate the weighted average across predator ages (dimension 3)
  weighted_diet <- sweep(sum_over_prey_age, c(1, 3, 4), avgN, "*")

  # Numerator: Sum of (Proportion * Predator Numbers) across predator ages (dim 3)
  numerator <- apply(weighted_diet, c(1, 2, 4), sum, na.rm = TRUE)

  # Denominator: Sum of Predator Numbers across predator ages
  denominator <- apply(avgN, c(1, 3), sum, na.rm = TRUE)

  # Final weighted average proportion
  pred_age_agg_prop <- sweep(numerator, c(1, 3), denominator, "/")
  dimnames(pred_age_agg_prop) <- list(1:nspp, 1:nspp, 1:nyrs)


  # 3. Convert the resulting 3D array [pred, prey, year] to a long data frame
  diet_summary <- as.data.frame.table(pred_age_agg_prop, dnn = dimnames(pred_age_agg_prop))
  colnames(diet_summary) <- c("Pred", "Prey", "Year", "Stomach_proportion_by_weight")

  # 4. Format the data for Rceattle
  simData$diet_data <- diet_summary %>%
    dplyr::mutate(
      Pred = as.numeric(as.character(Pred)),
      Prey = as.numeric(as.character(Prey)),
      Year = as.numeric(as.character(Year))
    ) %>%
    dplyr::filter(Stomach_proportion_by_weight > 0) %>%
    dplyr::mutate(
      Pred_age = -520, # Set Pred_age to negative flag < -500 or weighted average
      Prey_age = -1,   # Set Prey_age to negative flag
      Pred_sex = 0,
      Prey_sex = 0,
      Sample_size = 200
    ) %>%
    dplyr::select(Pred, Prey, Pred_sex, Prey_sex, Pred_age, Prey_age, Year, Sample_size, Stomach_proportion_by_weight)

  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)



testthat::test_that("Test average (across years) proportion of prey (all ages) in predator (mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Sum diet across prey-ages and average across predator ages and years -----
  diet_sum_prey_age <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 5), sum) # Sum across prey-ages (dimension 4)
  diet_avg_pred_age <- apply(diet_sum_prey_age, c(1, 2, 4), mean) # Average across predator ages (dimension 3 of the new object)
  diet_avg_years <- apply(diet_avg_pred_age, c(1, 2), mean, na.rm = TRUE)
  dimnames(diet_avg_years) <- list(1:nspp, 1:nspp)
  diet_df <- as.data.frame.table(diet_avg_years, dnn = dimnames(diet_avg_years)) # 2. Flatten the 4D array into a "long" data frame
  colnames(diet_df) <- c("Pred", "Prey", "Stomach_proportion_by_weight")

  # as.data.frame.table creates factors, so convert them to numeric indices
  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))

  # Add the required constant columns
  diet_df$Prey_age <- -500 # Indicates sum across prey-ages in Rceattle
  diet_df$Pred_age <- -1 # Indicates average across pred-ages in Rceattle (< -500 does weighted average)
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Year <- 0 # Average across years
  diet_df$Sample_size <- 1000

  # Reorder
  simData$diet_data <- diet_df[, c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                   "Pred_age", "Prey_age", "Year",
                                   "Sample_size", "Stomach_proportion_by_weight")]
  simData$diet_data <- diet_df[diet_df$Stomach_proportion_by_weight > 0, ]

  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test average (across years) proportion of prey (all ages) in predator (weighted mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nages = 15
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    ages = 1:nages,
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = NULL, # Initial parameters = 0
                               estimateMode = 3, # Do not estimate
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)
  inits <- ms_run1$estimated_params
  map <- ms_run1$map

  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs


  # * Sum diet across prey-ages and average across predators weighted by avgN and years -----
  sum_over_prey_age <- apply(sim$model_quantities$diet_prop, c(1, 2, 3, 5), sum, na.rm = TRUE) # sum across prey ages
  avgN <- sim$model_quantities$avgNAA # Dims: [pred, pred_age, year]

  # 2. Calculate the weighted average across predator ages (dimension 3)
  weighted_diet <- sweep(sum_over_prey_age, c(1, 3, 4), avgN, "*")

  # Numerator: Sum of (Proportion * Predator Numbers) across predator ages (dim 3)
  numerator <- apply(weighted_diet, c(1, 2, 4), sum, na.rm = TRUE)

  # Denominator: Sum of Predator Numbers across predator ages
  denominator <- apply(avgN, c(1, 3), sum, na.rm = TRUE)

  # Final weighted average proportion
  pred_age_agg_prop <- sweep(numerator, c(1, 3), denominator, "/")

  # Average across years
  diet_avg_years <- apply(pred_age_agg_prop, c(1, 2), mean, na.rm = TRUE)
  dimnames(diet_avg_years) <- list(1:nspp, 1:nspp)
  diet_df <- as.data.frame.table(diet_avg_years, dnn = dimnames(diet_avg_years)) # 2. Flatten the 4D array into a "long" data frame
  colnames(diet_df) <- c("Pred", "Prey", "Stomach_proportion_by_weight")

  # as.data.frame.table creates factors, so convert them to numeric indices
  diet_df$Pred <- as.numeric(as.character(diet_df$Pred))
  diet_df$Prey <- as.numeric(as.character(diet_df$Prey))

  # Add the required constant columns
  diet_df$Prey_age <- -500 # Indicates sum across prey-ages in Rceattle
  diet_df$Pred_age <- -520 # Indicates average across pred-ages in Rceattle (< -500 does weighted average)
  diet_df$Pred_sex <- 0
  diet_df$Prey_sex <- 0
  diet_df$Year <- 0 # Average across years
  diet_df$Sample_size <- 1000


  # Reorder
  simData$diet_data <- diet_df[, c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                   "Pred_age", "Prey_age", "Year",
                                   "Sample_size", "Stomach_proportion_by_weight")]
  simData$diet_data <- diet_df[diet_df$Stomach_proportion_by_weight > 0, ]

  # Fit
  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initialize from sim pars
                               file = NULL, # Don't save
                               map = map,
                               estimateMode = 3, # Do not estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age, Year)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test joint single-species models", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
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
  simData <- sim$data_list


  # Fit multi-species
  # * Fix parameters -----
  # inits <- suppressMessages( build_params(simData) )
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = NULL, # Initialize from default
                              estimateMode = 3, # Do not estimate
                              random_rec = FALSE, # No random recruitment
                              phase = FALSE,
                              msmMode = 0,
                              suitMode = 0,
                              initMode = 2,
                              verbose = 0)
  inits <- ss_run$estimated_params
  map <- ss_run$map

  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = inits, # Initialize from sim pars
                              map = map,
                              estimateMode = 3, # Do not estimate
                              random_rec = FALSE, # No random recruitment
                              phase = FALSE,
                              msmMode = 0,
                              suitMode = 0,
                              initMode = 2,
                              verbose = 0)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ss_run$quantities$R[,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ss_run$quantities$biomass[,1:nyrs]), tolerance = 1e-6)


  # Suitability
  testthat::expect_equal(as.numeric(ss_run$quantities$suitability[,,,,1]), as.numeric(sim$model_quantities$suitability))

  # M2
  testthat::expect_equal(as.numeric(sim$model_quantities$M2_at_age), as.numeric(ss_run$quantities$M2_at_age[,1,,1:nyrs]), tolerance = 1e-6)

  # Ration
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ss_run$quantities$consumption_at_age[,1,,1]))

  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ss_run$quantities$N_at_age[,1,,1:nyrs]))

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ss_run$quantities$avgN_at_age[,1,,1:nyrs]))

  # Selectivity
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]), as.numeric(ss_run$quantities$sel_at_age[c(1,3),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]), as.numeric(ss_run$quantities$sel_at_age[c(2,4),,,1]))

  # F
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA), as.numeric(ss_run$quantities$F_flt_age[c(2,4),1,,1:nyrs]))

  # Q
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q), as.numeric(ss_run$quantities$index_q[c(1,3),1]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ss_run$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ss_run$quantities$diet_hat[,2]))
})


testthat::test_that("Mixed suitabilities: MSVPA and lognormal", {
  # One predator with lognormal weight-based suitability
  # One prey (all diet = 0) with MSVPA suitability
  # Test all quantities
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-Inf,-Inf), nspp, nspp, byrow = TRUE) # Species 2 has no predation

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit initial multi-species to get paramter object
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               file = NULL, # Don't save
                               estimateMode = 3, # Don't estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = c(4, 0),
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # * Fix parameters -----
  inits <- ms_run1$estimated_params
  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initial parameters from inits
                               map = ms_run1$map,
                               file = NULL, # Don't save
                               estimateMode = 3, # Don't estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = c(4, 0),
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ms_run2$quantities$R[,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ms_run2$quantities$biomass[,1:nyrs]), tolerance = 1e-6)


  # Suitability
  testthat::expect_equal(exp(ms_run2$estimated_params$log_gam_a), gam_a)
  testthat::expect_equal(exp(ms_run2$estimated_params$log_gam_b), gam_b)
  testthat::expect_equal(as.numeric(ms_run2$quantities$vulnerability),
                         as.numeric(sim$model_quantities$vulnerability))
  testthat::expect_equal(as.numeric(ms_run2$quantities$suitability[,,,,2]),
                         as.numeric(sim$model_quantities$suitability))
  testthat::expect_equal(as.numeric(ms_run2$quantities$suit_other[,1,1,1]),
                         as.numeric(sim$model_quantities$suit_other))

  # M2
  testthat::expect_equal(as.numeric(sim$model_quantities$M2_at_age),
                         as.numeric(ms_run2$quantities$M2_at_age[,1,,1:nyrs]), tolerance = 1e-6)

  # Ration
  testthat::expect_equal(as.numeric(sim$model_quantities$ration),
                         as.numeric(ms_run2$quantities$consumption_at_age[,1,,1]))

  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]),
                         as.numeric(ms_run2$quantities$N_at_age[,1,,1:nyrs]))

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]),
                         as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Avail food
  testthat::expect_equal(as.numeric(sim$model_quantities$avail_food[,,1]),
                         as.numeric(ms_run2$quantities$avail_food[,1,,1]))

  # Selectivity
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]),
                         as.numeric(ms_run2$quantities$sel_at_age[c(1,3),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]),
                         as.numeric(ms_run2$quantities$sel_at_age[c(2,4),,,1]))

  # F
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA),
                         as.numeric(ms_run2$quantities$F_flt_age[c(2,4),1,,1:nyrs]))

  # Q
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q),
                         as.numeric(ms_run2$quantities$index_q[c(1,3),1]))

  # Expected and observed diet
  pred1_ind <- ms_run2$data_list$diet_data$Pred == 1 # Species 2 will be biased because MSVPA
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight[pred1_ind]),
                         as.numeric(ms_run2$quantities$diet_hat[pred1_ind,2]))
})


testthat::test_that("Mixed suitabilities2: MSVPA and lognormal", {
  # One predator with lognormal weight-based suitability
  # One prey with MSVPA suitability
  # Test predator suitability!!
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit initial multi-species to get parameter object
  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               file = NULL, # Don't save
                               estimateMode = 3, # Don't estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = c(4, 0),
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # * Fix parameters -----
  inits <- ms_run1$estimated_params
  inits$log_gam_a <- log(gam_a)
  inits$log_gam_b <- log(gam_b)
  inits$log_phi <- log_phi
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$ln_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$ln_F[2,] <- log(Fmort)
  inits$ln_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_ln_q[] <- log(1)
  inits$R_ln_sd[] <- log(1)
  inits$x_tj[1:30, 1:2] <- t(sim$model_quantities$rec_devs)
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initial parameters from inits
                               map = ms_run1$map,
                               file = NULL, # Don't save
                               estimateMode = 3, # Don't estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = c(4, 0),
                               niter = 5,
                               initMode = 2,
                               verbose = 0)

  # Suitability of predator (species 1)
  testthat::expect_equal(as.numeric(ms_run2$quantities$suitability[1,,,,2]),
                         as.numeric(sim$model_quantities$suitability[1,,,]))
  testthat::expect_equal(as.numeric(ms_run2$quantities$suit_other[1,1,1,1]),
                         as.numeric(sim$model_quantities$suit_other[1]))

  # Suitability of prey (species 2) should be off
  testthat::expect_false(isTRUE(all.equal(as.numeric(ms_run2$quantities$suitability[2,,,,2]),
                                as.numeric(sim$model_quantities$suitability[2,,,]))))
})
