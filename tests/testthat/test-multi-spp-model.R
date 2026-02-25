test_that("Rceattle and simple multi-species model match", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nspp = 2
  nyrs = 30
  years <- 1:nyrs
  ages <- 1:15
  WAA <- 2 / (1 + exp(-0.8 * (ages - 3)))
  WAA2 <- 1.4 / (1 + exp(-1 * (ages - 3)))
  MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
  MatAA2 <- 1 / (1 + exp(-0.8 * (ages - 5)))
  sigma_R <- 0.3
  sigma_Catch <- 0.001
  sigma_SrvIdx <- 0.05
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)
  other_food <- rep(1e5, nspp)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    nspp = 2,
    years = years,
    ages = ages,
    WAA = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE),
    MatAA = matrix(c(MatAA, MatAA2), nspp, length(ages), byrow = TRUE),
    mean_Rec = c(1e2, 1e3),
    sigma_R = 1,
    sigma_catch = sigma_Catch,
    sigma_srv = sigma_SrvIdx,
    diet_ISS = 1e5,
    fish_ISS = 1e5,
    srv_ISS = 1e5,
    M = c(0.2, 0.3),
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-2.5 * (ages - 4)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(c(Fmort, Fmort2), nspp, length(years), byrow = TRUE),
    srv_q = rep(1, nspp),

    # Multispecies bits
    niter = 5,
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi,
    other_food = other_food,
    ration = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE) * 50
  )


  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # * Fix suitability -----
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
  testthat::expect_equal(as.numeric(sim$model_quantities$ration), as.numeric(ms_run1$quantities$ration[,1,,1]))

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
})


test_that("Simulated simple multi-species model the same", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nspp = 2
  nyrs = 30
  years <- 1:nyrs
  ages <- 1:15
  WAA <- 2 / (1 + exp(-0.8 * (ages - 3)))
  WAA2 <- 1.4 / (1 + exp(-1 * (ages - 3)))
  MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
  MatAA2 <- 1 / (1 + exp(-0.8 * (ages - 5)))
  sigma_R <- 0.3
  sigma_Catch <- 0.001
  sigma_SrvIdx <- 0.05
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)
  other_food <- rep(1e5, nspp)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    nspp = 2,
    years = years,
    ages = ages,
    WAA = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE),
    MatAA = matrix(c(MatAA, MatAA2), nspp, length(ages), byrow = TRUE),
    mean_Rec = c(1e2, 1e3),
    sigma_R = 1,
    sigma_catch = sigma_Catch,
    sigma_srv = sigma_SrvIdx,
    diet_ISS = 1e5,
    fish_ISS = 1e5,
    srv_ISS = 1e5,
    M = c(0.2, 0.3),
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-2.5 * (ages - 4)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(c(Fmort, Fmort2), nspp, length(years), byrow = TRUE),
    srv_q = rep(1, nspp),

    # Multispecies bits
    niter = 5,
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi,
    other_food = other_food,
    ration = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE) * 50
  )

  # * Plot
  par(mfrow = c(4,1), mar = c(4,4,0.1,0))
  plot(y = sim$model_quantities$Total_Biom[1,], x = years, type = "l", ylab = "Species 1 B")
  plot(y = sim$model_quantities$Total_Biom[2,], x = years, type = "l", ylab = "Species 2 B")
  plot(y = colSums(sim$model_quantities$B_eaten_as_prey[1,,]), x = years, type = "l", ylab = "Species 1 B consumed")
  plot(y = colSums(sim$model_quantities$B_eaten_as_prey[2,,]), x = years, type = "l", ylab = "Species 2 B consumed")

  # 2) Set up Rceattle data
  simData <- sim$data_list

  # 4)  Fit multi-species
  # * Fix suitability -----
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


  plot(x = sim$model_quantities$NAA[1,1,], y = ms_run1$quantities$R[1,1:nyrs]); abline(1,1)
  plot(x = sim$model_quantities$NAA[2,1,], y = ms_run1$quantities$R[2,1:nyrs]); abline(1,1)

  plot(x = sim$Total_Biom[1,], y = ms_run1$quantities$biomass[1,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom[2,], y = ms_run1$quantities$biomass[2,1:nyrs]); abline(1,1)

  # Suitability is OK
  exp(ms_run1$estimated_params$log_gam_a)-gam_a
  exp(ms_run1$estimated_params$log_gam_b)-gam_b

  sum(ms_run1$quantities$suitability[,,,,1] - sim$suitability)
  ms_run1$quantities$vulnerability-sim$vulnerability

  sim$suit_other -  ms_run1$quantities$vulnerability_other

  # M2
  sum(sim$M2_at_age - ms_run1$quantities$M2_at_age[,1,,1:nyrs])

  # Ration
  ms_run1$quantities$ration[,1,,1] - sim$ration

  # AvgN
  sum(ms_run1$quantities$avgN_at_age[,1,,1:nyrs] -sim$avgNAA[,,])

  # Avail food
  ms_run1$quantities$avail_food[,1,,1]-sim$avail_food[,,1]

  # Selectivity
  ms_run1$quantities$sel_at_age[c(1,3),,,1]-sim$srv_sel
  ms_run1$quantities$sel_at_age[c(2,4),,,1]-sim$fish_sel

  # F
  sum(ms_run1$quantities$F_flt_age[c(2,4),1,,1:nyrs] - sim$FAA)

  # Q
  ms_run1$quantities$index_q
  sim$srv_q


  # * Full  diet -----
  # Get all combinations of indices
  idx <- which(!is.na(sim$diet_prop), arr.ind = TRUE)  # or sim$diet_prop != 0 for nonzero only

  # Build the data frame directly
  simData$diet_data <- data.frame(
    Year = idx[, 5],
    Pred = idx[, 1],
    Prey = idx[, 2],
    Pred_sex = 0,
    Prey_sex = 0,
    Pred_age = idx[, 3],
    Prey_age = idx[, 4],
    Sample_size = 100,
    Stomach_proportion_by_weight = sim$diet_prop[idx]
  )
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))

  ss_run$estimated_params$log_gam_a <- log(gam_a)
  ss_run$estimated_params$log_gam_b <- log(gam_b)
  ss_run$estimated_params$log_phi <- log_phi

  ms_run1 <- Rceattle::fit_mod(data_list = simData,
                               inits = inits, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 0, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               suitMode = 4,
                               niter = 5,
                               initMode = 2,
                               verbose = 1)

  plot(x = sim$Total_Biom[1,], y = ms_run1$quantities$biomass[1,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom[2,], y = ms_run1$quantities$biomass[2,1:nyrs]); abline(1,1)
  exp(ms_run1$estimated_params$log_gam_a)
  gam_a
  exp(ms_run1$estimated_params$log_gam_b)
  gam_b

  ms_run1$quantities$vulnerability
  sim$vulnerability


  # * Average diet across years -----
  simData$diet_data[] <- NA
  ind = 1
  for(i in 1:dim(sim$diet_prop)[1]){ # pred
    for(j in 1:dim(sim$diet_prop)[2]){ # prey
      for(k in 1:dim(sim$diet_prop)[3]){ # pred age
        for(l in 1:dim(sim$diet_prop)[4]){ # prey-age
          simData$diet_data[ind,] <- NA
          simData$diet_data$Year[ind] <- 0 # average years
          simData$diet_data$Pred[ind] <- i
          simData$diet_data$Prey[ind] <- j
          simData$diet_data$Pred_sex[ind] <- 0
          simData$diet_data$Prey_sex[ind] <- 0
          simData$diet_data$Pred_age[ind] <- k
          simData$diet_data$Prey_age[ind] <- l
          simData$diet_data$Sample_size[ind] <- 200
          simData$diet_data$Stomach_proportion_by_weight[ind] <- mean(sim$diet_prop[i,j,k,l,])
          ind = ind+1
        }
      }
    }
  }
  simData$diet_data$Sample_size <- 1000
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))


  ms_run2 <- Rceattle::fit_mod(data_list = simData,
                               inits = ss_run$estimated_params, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 0, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               niter = 5,
                               suitMode = 4,
                               initMode = 2,
                               verbose = 1)

  plot(x = sim$Total_Biom[1,], y = ms_run2$quantities$biomass[1,1:nyrs]); abline(0,1)
  plot(x = sim$Total_Biom[2,], y = ms_run2$quantities$biomass[2,1:nyrs]); abline(0,1)

  plot(x = ms_run2$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run2$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Average sum diet")
  abline(0,1)



  # * Annual prey-spp diet -----
  simData$diet_data[] <- NA
  ind = 1
  for(i in 1:dim(sim$diet_prop)[1]){ # pred
    for(j in 1:dim(sim$diet_prop)[2]){ # prey
      for(k in 1:dim(sim$diet_prop)[3]){ # pred age
        for(m in 1:dim(sim$diet_prop)[5]){ # year
          simData$diet_data[ind,] <- NA
          simData$diet_data$Year[ind] <- m
          simData$diet_data$Pred[ind] <- i
          simData$diet_data$Prey[ind] <- j
          simData$diet_data$Pred_sex[ind] <- 0
          simData$diet_data$Prey_sex[ind] <- 0
          simData$diet_data$Pred_age[ind] <- k
          simData$diet_data$Prey_age[ind] <- -500  # sum across prey-ages
          simData$diet_data$Sample_size[ind] <- 200
          simData$diet_data$Stomach_proportion_by_weight[ind] <- sum(sim$diet_prop[i,j,k,,m])
          ind = ind+1
        }
      }
    }
  }
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))


  ms_run3 <- Rceattle::fit_mod(data_list = simData,
                               inits = ss_run$estimated_params, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 0, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               niter = 5,
                               suitMode = 4,
                               initMode = 2,
                               verbose = 1)

  plot(x = sim$Total_Biom[1,], y = ms_run3$quantities$biomass[1,1:nyrs]); abline(0,1)
  plot(x = sim$Total_Biom[2,], y = ms_run3$quantities$biomass[2,1:nyrs]); abline(0,1)

  plot(x = ms_run3$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run3$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Average sum diet")
  abline(0,1)


  # * Average prey-spp diet across years -----
  simData$diet_data[] <- NA
  ind = 1
  for(i in 1:dim(sim$diet_prop)[1]){ # pred
    for(j in 1:dim(sim$diet_prop)[2]){ # prey
      for(k in 1:dim(sim$diet_prop)[3]){ # pred age
        simData$diet_data[ind,] <- NA
        simData$diet_data$Year[ind] <- 0 # average of years
        simData$diet_data$Pred[ind] <- i
        simData$diet_data$Prey[ind] <- j
        simData$diet_data$Pred_sex[ind] <- 0
        simData$diet_data$Prey_sex[ind] <- 0
        simData$diet_data$Pred_age[ind] <- k
        simData$diet_data$Prey_age[ind] <- -500 # sum across prey-ages
        simData$diet_data$Sample_size[ind] <- 200
        simData$diet_data$Stomach_proportion_by_weight[ind] <- mean(rowSums(sim$diet_prop[i,j,k,,]))
        ind = ind+1
      }
    }
  }
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))


  ms_run4 <- Rceattle::fit_mod(data_list = simData,
                               inits = ss_run$estimated_params, # Initial parameters = 0
                               file = NULL, # Don't save
                               estimateMode = 0, # Estimate
                               random_rec = FALSE, # No random recruitment
                               phase = FALSE,
                               msmMode = 1,
                               niter = 5,
                               suitMode = 4,
                               initMode = 2,
                               verbose = 1)

  plot(x = sim$Total_Biom[1,], y = ms_run4$quantities$biomass[1,1:nyrs]); abline(0,1)
  plot(x = sim$Total_Biom[2,], y = ms_run4$quantities$biomass[2,1:nyrs]); abline(0,1)

  plot(x = ms_run4$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run4$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Average sum diet")
  abline(0,1)

})



test_that("Test joint single-species models", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

  # Set up simulation
  nspp = 2
  nyrs = 30
  years <- 1:nyrs
  ages <- 1:15
  WAA <- 2 / (1 + exp(-0.8 * (ages - 3)))
  WAA2 <- 1.4 / (1 + exp(-1 * (ages - 3)))
  MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
  MatAA2 <- 1 / (1 + exp(-0.8 * (ages - 5)))
  sigma_R <- 0.3
  sigma_Catch <- 0.001
  sigma_SrvIdx <- 0.2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    nspp = 2,
    years = years,
    ages = ages,
    WAA = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE),
    MatAA = matrix(c(MatAA, MatAA2), nspp, length(ages), byrow = TRUE),
    mean_Rec = c(1e2, 1e3),
    sigma_R = 1,
    sigma_catch = sigma_Catch,
    sigma_srv = sigma_SrvIdx,
    diet_ISS = 1e5,
    fish_ISS = 1e5,
    srv_ISS = 1e5,
    M = c(0.2, 0.3),
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-1 * (ages - 5)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(c(Fmort, Fmort2), nspp, length(years), byrow = TRUE),
    srv_q = rep(1, nspp),

    # Multispecies bits
    gam_a = c(1, 0.1),
    gam_b = rep(0.3, nspp),
    log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE),
    other_food = rep(1e5, nspp),
    ration = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE) * 50,
  )


  # Plot ------------
  par(mfrow = c(4,1), mar = c(4,4,0.1,0))
  plot(y = sim$Total_Biom[1,], x = years, type = "l", ylab = "Species 1 B")
  plot(y = sim$Total_Biom[2,], x = years, type = "l", ylab = "Species 2 B")
  plot(y = colSums(sim$B_eaten_as_prey[1,,]), x = years, type = "l", ylab = "Species 1 B consumed")
  plot(y = colSums(sim$B_eaten_as_prey[2,,]), x = years, type = "l", ylab = "Species 2 B consumed")


  # Fit single-species -------------------------------------------------------------
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              initMode = 2,
                              verbose = 1)

  plot(x = sim$SSB[1,], y = ss_run$quantities$ssb[1,1:nyrs]); abline(1,1)
  plot(x = sim$SSB[2,], y = ss_run$quantities$ssb[2,1:nyrs]); abline(1,1)

})



test_that("Test expected diet composition", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nspp = 2
  nyrs = 30
  years <- 1:nyrs
  ages <- 1:15
  WAA <- 2 / (1 + exp(-0.8 * (ages - 3)))
  WAA2 <- 1.4 / (1 + exp(-1 * (ages - 3)))
  MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
  MatAA2 <- 1 / (1 + exp(-0.8 * (ages - 5)))
  sigma_R <- 0.3
  sigma_Catch <- 0.001
  sigma_SrvIdx <- 0.05
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a = c(1, 0.1)
  gam_b = rep(0.15, nspp)
  log_phi = matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE)
  other_food <- rep(1e5, nspp)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    nspp = 2,
    years = years,
    ages = ages,
    WAA = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE),
    MatAA = matrix(c(MatAA, MatAA2), nspp, length(ages), byrow = TRUE),
    mean_Rec = c(1e2, 1e3),
    sigma_R = 1,
    sigma_catch = sigma_Catch,
    sigma_srv = sigma_SrvIdx,
    diet_ISS = 1e5,
    fish_ISS = 1e5,
    srv_ISS = 1e5,
    M = c(0.2, 0.3),
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-2.5 * (ages - 4)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(c(Fmort, Fmort2), nspp, length(years), byrow = TRUE),
    srv_q = rep(1, nspp),

    # Multispecies bits
    niter = 5,
    gam_a = gam_a,
    gam_b = gam_b,
    log_phi = log_phi,
    other_food = other_food,
    ration = matrix(c(WAA, WAA2), nspp, length(ages), byrow = TRUE) * 50
  )


  # * Plot
  par(mfrow = c(4,1), mar = c(4,4,0.1,0))
  plot(y = sim$model_quantities$Total_Biom[1,], x = years, type = "l", ylab = "Species 1 B")
  plot(y = sim$model_quantities$Total_Biom[2,], x = years, type = "l", ylab = "Species 2 B")
  plot(y = colSums(sim$model_quantities$B_eaten_as_prey[1,,]), x = years, type = "l", ylab = "Species 1 B consumed")
  plot(y = colSums(sim$model_quantities$B_eaten_as_prey[2,,]), x = years, type = "l", ylab = "Species 2 B consumed")


  # 4)  Fit multi-species
  simData <- sim$data_list

  # * Fix all parameters -----
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
  inits$rec_dev[,1:30] <- sim$rec_devs
  inits$init_dev[,1:14] <- sim$init_devs

  # * Full  diet -----
  # Get all combinations of indices
  idx <- which(!is.na(sim$diet_prop), arr.ind = TRUE)  # or sim$diet_prop != 0 for nonzero only

  # Build the data frame directly
  simData$diet_data <- data.frame(
    Year = idx[, 5],
    Pred = idx[, 1],
    Prey = idx[, 2],
    Pred_sex = 0,
    Prey_sex = 0,
    Pred_age = idx[, 3],
    Prey_age = idx[, 4],
    Sample_size = 1000,
    Stomach_proportion_by_weight = sim$diet_prop[idx]
  )
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))

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

  # * Plot output
  plot(x = sim$NAA[1,1,], y = ms_run1$quantities$R[1,1:nyrs]); abline(1,1)
  plot(x = sim$NAA[2,1,], y = ms_run1$quantities$R[2,1:nyrs]); abline(1,1)

  plot(x = sim$Total_Biom[1,], y = ms_run1$quantities$biomass[1,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom[2,], y = ms_run1$quantities$biomass[2,1:nyrs]); abline(1,1)


  # ** Diet check ----
  dev.off()
  plot(x = ms_run1$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run1$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Full diet")



  # * Average diet across years -----
  simData$diet_data[] <- NA
  ind = 1
  for(i in 1:dim(sim$diet_prop)[1]){ # pred
    for(j in 1:dim(sim$diet_prop)[2]){ # prey
      for(k in 1:dim(sim$diet_prop)[3]){ # pred age
        for(l in 1:dim(sim$diet_prop)[4]){ # prey-age
          simData$diet_data[ind,] <- NA
          simData$diet_data$Year[ind] <- 0 # average years
          simData$diet_data$Pred[ind] <- i
          simData$diet_data$Prey[ind] <- j
          simData$diet_data$Pred_sex[ind] <- 0
          simData$diet_data$Prey_sex[ind] <- 0
          simData$diet_data$Pred_age[ind] <- k
          simData$diet_data$Prey_age[ind] <- l
          simData$diet_data$Sample_size[ind] <- 200
          simData$diet_data$Stomach_proportion_by_weight[ind] <- mean(sim$diet_prop[i,j,k,l,])
          ind = ind+1
        }
      }
    }
  }
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))


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


  # ** Diet check ----
  plot(x = ms_run2$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run2$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Average diet")


  # * Annual prey-spp diet -----
  # Sum across prey-ages for each year
  simData$diet_data[] <- NA
  ind = 1
  for(i in 1:dim(sim$diet_prop)[1]){ # pred
    for(j in 1:dim(sim$diet_prop)[2]){ # prey
      for(k in 1:dim(sim$diet_prop)[3]){ # pred age
        for(m in 1:dim(sim$diet_prop)[5]){ # year
          simData$diet_data[ind,] <- NA
          simData$diet_data$Year[ind] <- m
          simData$diet_data$Pred[ind] <- i
          simData$diet_data$Prey[ind] <- j
          simData$diet_data$Pred_sex[ind] <- 0
          simData$diet_data$Prey_sex[ind] <- 0
          simData$diet_data$Pred_age[ind] <- k
          simData$diet_data$Prey_age[ind] <- -500  # sum across prey-ages
          simData$diet_data$Sample_size[ind] <- 200
          simData$diet_data$Stomach_proportion_by_weight[ind] <- sum(sim$diet_prop[i,j,k,,m])
          ind = ind+1
        }
      }
    }
  }
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))


  ms_run3 <- Rceattle::fit_mod(data_list = simData,
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


  plot(x = ms_run3$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run3$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Annual sum diet")


  # * Average prey-spp diet across years -----
  simData$diet_data[] <- NA
  ind = 1
  for(i in 1:dim(sim$diet_prop)[1]){ # pred
    for(j in 1:dim(sim$diet_prop)[2]){ # prey
      for(k in 1:dim(sim$diet_prop)[3]){ # pred age
        simData$diet_data[ind,] <- NA
        simData$diet_data$Year[ind] <- 0 # average of years
        simData$diet_data$Pred[ind] <- i
        simData$diet_data$Prey[ind] <- j
        simData$diet_data$Pred_sex[ind] <- 0
        simData$diet_data$Prey_sex[ind] <- 0
        simData$diet_data$Pred_age[ind] <- k
        simData$diet_data$Prey_age[ind] <- -500 # sum across prey-ages
        simData$diet_data$Sample_size[ind] <- 200
        simData$diet_data$Stomach_proportion_by_weight[ind] <- mean(rowSums(sim$diet_prop[i,j,k,,]))
        ind = ind+1
      }
    }
  }
  simData$diet_data <- simData$diet_data %>%
    filter(!is.na(Pred))


  ms_run4 <- Rceattle::fit_mod(data_list = simData,
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

  plot(x = ms_run4$data_list$diet_data$Stomach_proportion_by_weight, y = ms_run4$quantities$diet_hat[,2], xlab = "True diet", ylab ="Est diet", main = "Average sum diet"); abline(0,1)

})

