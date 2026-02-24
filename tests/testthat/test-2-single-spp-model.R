# TODO: fleet_code to fleet_index
# Get rid of sex ratio in control
# call "wt" "weight"
# call "pmature" "maturity"

test_that("Simulated simple multi-species model the same", {

  # Set up simulation -------------------------------------------------------------
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
  sim <- sim_msm_model(
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

  # Set up Rceattle data -------------------------------------------------------------
  library(Rceattle)
  data("GOAcod")
  simData <- GOAcod

  # * Data controls ----
  simData$nspp <- nspp
  simData$spnames <- paste0("Species",1:2)
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs+10
  simData$nsex <- rep(1, nspp)
  simData$spawn_month <- rep(0, nspp)
  simData$nages <- rep(15, nspp)
  simData$minage <- rep(1, nspp)
  simData$nlengths <- rep(15, nspp)
  simData$estDynamics <- rep(0, nspp)
  simData$pop_wt_index <- 1:nspp
  simData$ssb_wt_index <- 1:nspp
  simData$other_food <- rep(1e5, nspp)
  simData$sigma_rec_prior <- rep(0.707, nspp)
  simData$pop_age_transition_index <- rep(1, nspp)

  # * Fleet control ----
  simData$fleet_control <- simData$fleet_control[c(1,3,1,3),] # BT and Trawl Fishery are both simple logistic
  simData$fleet_control$Fleet_name <- c("Survey1", "Fishery1","Survey2", "Fishery2")
  simData$fleet_control$Fleet_code <- 1:4
  simData$fleet_control$Selectivity_index <- 1:4
  simData$fleet_control$Catchability
  simData$fleet_control$Species <- c(1,1,2,2)
  simData$fleet_control$Q_index <- c(1,NA,2,NA)
  simData$fleet_control$Weight1_Numbers2 <- 1
  simData$fleet_control$Weight_index <- c(1,1,2,2)

  # * Index data ----
  simData$index_data <- rbind(
    data.frame(Fleet_name = "Survey1",
               Fleet_code = 1,
               Species = 1,
               Year = 1:nyrs,
               Month = 0,
               Selectivity_block = 1,
               Q_block = 1,
               Observation = sim$SrvIdx[1,],
               Log_sd = sigma_SrvIdx),

    data.frame(Fleet_name = "Survey2",
               Fleet_code = 3,
               Species = 2,
               Year = 1:nyrs,
               Month = 0,
               Selectivity_block = 1,
               Q_block = 1,
               Observation = sim$SrvIdx[2,],
               Log_sd = sigma_SrvIdx)
  )

  # * Catch data ----
  simData$catch_data <- rbind(
    data.frame(Fleet_name = "Fishery1",
               Fleet_code = 2,
               Species = 1,
               Year = 1:nyrs,
               Month = 0,
               Selectivity_block = 1,
               Catch = sim$ObsCatch[1,],
               Log_sd = sigma_Catch),

    data.frame(Fleet_name = "Fishery2",
               Fleet_code = 4,
               Species = 2,
               Year = 1:nyrs,
               Month = 0,
               Selectivity_block = 1,
               Catch = sim$ObsCatch[2,],
               Log_sd = sigma_Catch)
  )

  # * Comp data ----
  # - Index
  tmp <- sim$ObsSrvAges[1,,]
  colnames(tmp) <- paste0("Comp_",ages)
  index_comp <- cbind(data.frame(Fleet_name = "Survey1",
                                 Fleet_code = 1,
                                 Species = 1,
                                 Sex = 0,
                                 Age0_Length1 = 0,
                                 Year = 1:nyrs,
                                 Month = 0,
                                 Sample_size = rowSums(tmp)),
                      tmp
  )

  tmp <- sim$ObsSrvAges[2,,]
  colnames(tmp) <- paste0("Comp_",ages)
  index_comp2 <- cbind(data.frame(Fleet_name = "Survey2",
                                  Fleet_code = 3,
                                  Species = 2,
                                  Sex = 0,
                                  Age0_Length1 = 0,
                                  Year = 1:nyrs,
                                  Month = 0,
                                  Sample_size = rowSums(tmp)),
                       tmp
  )

  # - Fishery
  tmp <- sim$ObsFishAges[1,,]
  colnames(tmp) <- paste0("Comp_",1:15)
  fishery_comp <- cbind(data.frame(Fleet_name = "Fishery1",
                                   Fleet_code = 2,
                                   Species = 1,
                                   Sex = 0,
                                   Age0_Length1 = 0,
                                   Year = 1:nyrs,
                                   Month = 0,
                                   Sample_size = rowSums(tmp)),
                        tmp
  )

  tmp <- sim$ObsFishAges[2,,]
  colnames(tmp) <- paste0("Comp_",1:15)
  fishery_comp2 <- cbind(data.frame(Fleet_name = "Fishery2",
                                    Fleet_code = 4,
                                    Species = 2,
                                    Sex = 0,
                                    Age0_Length1 = 0,
                                    Year = 1:nyrs,
                                    Month = 0,
                                    Sample_size = rowSums(tmp)),
                         tmp
  )

  simData$comp_data <- rbind(index_comp, index_comp2, fishery_comp, fishery_comp2)

  # * Empirical selectivity ----
  simData$emp_sel[] <- NA

  # * Fixed numbers ----
  simData$NByageFixed[] <- NA


  # * Age transition matrix ----
  tmp <- as.data.frame(diag(1,15))
  colnames(tmp) <- paste0("Length_",1:15)
  atf1 <- cbind(data.frame(Age_transition_name = "Base1",
                           Age_transition_index = 1,
                           Species = 1,
                           Sex = 0,
                           Age = 1:15),
                tmp
  )

  tmp <- as.data.frame(diag(1,15))
  colnames(tmp) <- paste0("Length_",1:15)
  atf2 <- cbind(data.frame(Age_transition_name = "Base2",
                           Age_transition_index = 2,
                           Species = 2,
                           Sex = 0,
                           Age = 1:15),
                tmp
  )

  simData$age_trans_matrix <- rbind(atf1, atf2)


  # * Age error ----
  tmp <- as.data.frame(diag(1,15))
  colnames(tmp) <- paste0("Obs_age",1:15)
  age_error <- cbind(data.frame(Species = 1,
                                True_age = 1:15),
                     tmp
  )

  tmp <- as.data.frame(diag(1,15))
  colnames(tmp) <- paste0("Obs_age",1:15)
  age_error2 <- cbind(data.frame(Species = 2,
                                 True_age = 1:15),
                      tmp
  )

  simData$age_error <- rbind(age_error, age_error2)


  # * Weight-at-age ----
  WAA <- as.data.frame(matrix(WAA, ncol = 15))
  colnames(WAA) <- paste0("Age",1:15)
  weight1 <- cbind(data.frame(Wt_name = "Base1",
                              Wt_index = 1,
                              Species = 1,
                              Sex = 0,
                              Year = 0),
                   WAA
  )

  WAA2 <- as.data.frame(matrix(WAA2, ncol = 15))
  colnames(WAA2) <- paste0("Age",1:15)
  weight2 <- cbind(data.frame(Wt_name = "Base2",
                              Wt_index = 2,
                              Species = 2,
                              Sex = 0,
                              Year = 0),
                   WAA2
  )

  simData$weight <- rbind(weight1, weight2)


  # * Maturity ----
  MatAA <- as.data.frame(matrix(MatAA, ncol = 15))
  colnames(MatAA) <- paste0("Age",1:15)
  maturity1 <- cbind(data.frame(Species = 1),
                     MatAA
  )

  MatAA2 <- as.data.frame(matrix(MatAA2, ncol = 15))
  colnames(MatAA2) <- paste0("Age",1:15)
  maturity2 <- cbind(data.frame(Species = 2),
                     MatAA2
  )

  simData$maturity <- rbind(maturity1, maturity2)


  # * Sex ratio ----
  sexratio <- as.data.frame(matrix(0.5, ncol = 15))
  colnames(sexratio) <- paste0("Age",1:15)
  simData$sex_ratio <- cbind(data.frame(Species = 1),
                             sexratio
  )
  simData$sex_ratio <- rbind(simData$sex_ratio,
                             cbind(data.frame(Species = 2),
                                   sexratio
                             )
  )

  # * Mortality ----
  mort <- as.data.frame(matrix(0.2, ncol = 15))
  colnames(mort) <- paste0("Age",1:15)
  M1_base <- cbind(data.frame(Species = 1,
                              Sex = 0),
                   mort
  )

  mort <- as.data.frame(matrix(0.3, ncol = 15))
  colnames(mort) <- paste0("Age",1:15)
  M1_base2 <- cbind(data.frame(Species = 2,
                               Sex = 0),
                    mort
  )
  simData$M1_base <- rbind(M1_base, M1_base2)

  # * Environmental data ----
  simData$env_data <- data.frame(Year = 1:nyrs,
                                 EnvData = 1)


  # * Relative foraging rate (days) ----
  WAA <- as.data.frame(matrix(WAA * 50, ncol = 15))
  colnames(WAA) <- paste0("Age",1:15)
  ration_data1 <- cbind(data.frame(Species = 1,
                            Sex = 0,
                            Year = 0),
                 WAA
  )

  WAA2 <- as.data.frame(matrix(WAA2 * 50, ncol = 15))
  colnames(WAA2) <- paste0("Age",1:15)
  ration_data2 <- cbind(data.frame(Species = 2,
                            Sex = 0,
                            Year = 0),
                 WAA2
  )

  simData$ration_data <- rbind(ration_data1, ration_data2)


  # * Bioenergetics ----
  simData$Ceq <- rep(4,nspp)
  simData$Cindex <- rep(simData$Cindex,nspp)
  simData$Pvalue <- rep(simData$Pvalue,nspp)
  simData$fday <- rep(simData$fday,nspp)
  simData$CA <- rep(simData$CA,nspp)
  simData$CB <- rep(simData$CB,nspp)
  simData$Qc <- rep(simData$Qc,nspp)
  simData$Tco <- rep(simData$Tco,nspp)
  simData$Tcm <- rep(simData$Tcm,nspp)
  simData$Tcl <- rep(simData$Tcl,nspp)
  simData$CK1 <- rep(simData$CK1,nspp)
  simData$CK4 <- rep(simData$CK4,nspp)
  simData$Diet_comp_weights <- rep(1,nspp)


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
