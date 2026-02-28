testthat::test_that("Rceattle and multi-species model dynamics match", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # Recruitment
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ms_run1$quantities$R[,1:nyrs]))
  testthat::expect_equal(as.numeric(sim$model_quantities$Total_Biom), as.numeric(ms_run1$quantities$biomass[,1:nyrs]), tolerance = 1e-6)


  # Suitability
  testthat::expect_equal(exp(ms_run1$estimated_params$log_gam_a), gam_a)
  testthat::expect_equal(exp(ms_run1$estimated_params$log_gam_b), gam_b)
  testthat::expect_equal(as.numeric(ms_run1$quantities$vulnerability), as.numeric(sim$model_quantities$vulnerability))
  testthat::expect_equal(as.numeric(ms_run1$quantities$suitability[,,,,1]), as.numeric(sim$model_quantities$suitability))
  testthat::expect_equal(as.numeric(ms_run1$quantities$suit_other[,1,1,1]), as.numeric(sim$model_quantities$suit_other))

  # M2
  testthat::expect_equal(as.numeric(sim$model_quantities$M2_at_age), as.numeric(ms_run1$quantities$M2_at_age[,1,,1:nyrs]), tolerance = 1e-6)

  # Ration
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ms_run1$quantities$consumption_at_age[,1,,1]))

  # N
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,,]), as.numeric(ms_run1$quantities$N_at_age[,1,,1:nyrs]))

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run1$quantities$avgN_at_age[,1,,1:nyrs]))

  # Avail food
  testthat::expect_equal(as.numeric(sim$model_quantities$avail_food[,,1]), as.numeric(ms_run1$quantities$avail_food[,1,,1]))

  # Selectivity
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_sel[,]), as.numeric(ms_run1$quantities$sel_at_age[c(1,3),,,1]))
  testthat::expect_equal(as.numeric(sim$model_quantities$fish_sel[,]), as.numeric(ms_run1$quantities$sel_at_age[c(2,4),,,1]))

  # F
  testthat::expect_equal(as.numeric(sim$model_quantities$FAA), as.numeric(ms_run1$quantities$F_flt_age[c(2,4),1,,1:nyrs]))

  # Q
  testthat::expect_equal(as.numeric(sim$model_quantities$srv_q), as.numeric(ms_run1$quantities$index_q[c(1,3),1]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run1$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run1$quantities$diet_hat[,2]))
})


testthat::test_that("Test proportion of prey-at-age in predator-at-age averaged across years", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test annual proportion of prey (all ages) in predator-at-age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test proportion of prey (all ages) in predator-at-age averaged across years", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)

testthat::test_that("Test annual proportion of prey (all ages) in predator (mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test annual proportion of prey (all ages) in predator (weighted mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)



testthat::test_that("Test average (across years) proportion of prey (all ages) in predator (mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test average (across years) proportion of prey (all ages) in predator (weighted mean across all ages)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

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
  inits <- build_params(simData)
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
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
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
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 3, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  # AvgN
  testthat::expect_equal(as.numeric(sim$model_quantities$avgNAA[,,]), as.numeric(ms_run2$quantities$avgN_at_age[,1,,1:nyrs]))

  # Expected and observed diet
  testthat::expect_equal(as.numeric(ms_run2$data_list$diet_data$Stomach_proportion_by_weight), as.numeric(ms_run2$quantities$diet_hat[,2]))

  # Diet data (no modifications when rearranged)
  diet_data1 <- ms_run2$data_list$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)

  diet_data2 <- simData$diet_data %>%
    arrange(Pred, Prey, Pred_age, Prey_age)
  testthat::expect_equal(as.numeric(diet_data1$Stomach_proportion_by_weight), as.numeric(diet_data2$Stomach_proportion_by_weight))
}
)


testthat::test_that("Test joint single-species models", {


})


