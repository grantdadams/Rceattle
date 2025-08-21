# TODO: fleet_code to fleet_index
# Get rid of sex ratio in control
# call "wt" "weight"
# call "pmature" "maturity"

test_that("Simulated simple multi-species model the same" {
  # Simulate Data -----------------------------------------------------------
  # Adapted from Matt Cheng
  sim_msm_model <- function(

    nspp = 2,
    years = 1:40,
    ages = 1:10,
    WAA = matrix(1:2, nspp, length(ages)),
    MatAA = matrix(1, nspp, length(ages)),
    mean_Rec = c(1e6, 1e5),
    sigma_R = 1,
    sigma_catch = 0.05,
    sigma_srv = 0.2,
    diet_ISS = 200,
    fish_ISS = 200,
    srv_ISS = 200,
    M = c(0.3, 0.2),
    fish_sel = matrix(1, nspp, length(ages)),
    srv_sel = matrix(1, nspp, length(ages)),
    Fmort = matrix(0.5, nspp, length(years)),
    srv_q = rep(1, nspp),

    # Multispecies bits
    iter = 10,
    gam_a = rep(0.1, nspp),
    gam_b = rep(0.1, nspp),
    log_phi = matrix(0, nspp, nspp),
    other_food = rep(1e5, nspp),
    ration = matrix(1, nspp, length(ages))

  ) {

    n_yrs <- length(years)
    n_ages <- length(ages)

    # Initialize arrays
    NAA <- array(0, dim = c(nspp, n_yrs + 1, n_ages))  # Numbers at age
    avgNAA <- array(0, dim = c(nspp, n_yrs + 1, n_ages))  # Numbers at age
    ZAA <- array(0, dim = c(nspp, n_yrs, n_ages))      # Total mortality
    CAA <- array(0, dim = c(nspp, n_yrs, n_ages))      # Catch at age
    FAA <- array(0, dim = c(nspp, n_yrs, n_ages))      # Fishing mortality at age

    # Population metrics
    SSB <- matrix(0, nspp, n_yrs)
    Total_Biom <- matrix(0, nspp, n_yrs)
    Catch <- matrix(0, nspp, n_yrs)

    # Generate recruitment deviations
    rec_devs <- matrix(rnorm(n_yrs * nspp, 0, sigma_R), nspp, n_yrs)
    init_devs <- matrix(rnorm((n_ages-1) * nspp, 0, sigma_R), nspp, n_ages-1)

    # Vulnerability
    suitability = array(0, dim = c(nspp, nspp, n_ages, n_ages))
    vulnerability <- matrix(0, nspp, nspp)
    suit_other <- c()

    for(sp in 1:nspp) {
      vulnerability[sp,] = exp(log_phi[sp,])/(1+sum(exp(log_phi[sp,]))) # multinomial logistic transformation
      suit_other[sp] = 1 - sum(vulnerability[sp,])
    }

    # Calculate suitability ----
    for(sp in 1:nspp) { # Predator
      for(r_age in 1:n_ages){  # Pred age
        for(ksp in 1:nspp) {   # Prey loop
          for(k_age in 1:n_ages){ # Prey age

            # Weight-based lognormal suitability
            log_size_ratio = log(WAA[sp,r_age] / WAA[ksp,k_age]) # Log ratio of weights

            if(log_size_ratio > 0){
              suitability[sp, ksp, r_age, k_age] = vulnerability[sp, ksp] * dnorm(log_size_ratio, gam_a[sp], gam_b[sp]) / dnorm(gam_a[sp], gam_a[sp], gam_b[sp]) # Divide by mode to scale to 1
            }
          }
        }
      }
    }

    # Initialize population
    for(sp in 1:nspp){
      init_age_idx <- 1:(n_ages - 2)
      NAA[sp, 1, init_age_idx + 1] <- mean_Rec[sp] * exp(- (init_age_idx * M[sp]))
      NAA[sp, 1, n_ages] <- mean_Rec[sp] * exp(-(n_ages - 1) * M[sp]) / (1 - exp(-M[sp]))
      NAA[sp,1,2:n_ages] <- NAA[sp,1,2:n_ages] * exp(init_devs[sp,])
    }

    # Project population forward
    M2_at_age <- array(0, dim = c(nspp, n_ages, n_yrs))
    for(iter in 1:10){

      for(y in 1:n_yrs) {
        for(sp in 1:nspp){
          # New recruits
          NAA[sp, y, 1] <- mean_Rec[sp] * exp(rec_devs[sp,y])

          # Calculate mortality
          FAA[sp, y,] <- Fmort[sp, y] * fish_sel[sp,]
          ZAA[sp, y,] <- FAA[sp,y,] + M[sp] + M2_at_age[sp, , y]

          # Calculate catch
          CAA[sp, y,] <- FAA[sp, y,] / ZAA[sp, y,] * NAA[sp, y,] * (1 - exp(-ZAA[sp, y,]))

          # Project survivors
          if(y < n_yrs) {
            for(a in 1:(n_ages-1)) {
              NAA[sp, y+1, a+1] <- NAA[sp, y, a] * exp(-ZAA[sp, y, a])
            }
            # Plus group
            NAA[sp, y+1, n_ages] <- NAA[sp, y+1, n_ages] + NAA[sp, y, n_ages] * exp(-ZAA[sp, y, n_ages])
          }
          avgNAA[sp,y,] = NAA[sp, y, ] * (1-exp(-ZAA[sp, y, ]))/ZAA[sp,y,]

          # Calculate annual metrics
          Total_Biom[sp, y] <- sum(NAA[sp, y,] * WAA[sp, ])
          SSB[sp, y] <- sum(NAA[sp, y,] * WAA[sp, ] * MatAA[sp, ]) * 0.5
          Catch[sp, y] <- sum(CAA[sp, y,] * WAA[sp, ])
        }
      }

      # Available food ----
      avail_food <- array(0, dim = c(nspp, n_ages, n_yrs))
      for(rsp in 1:nspp) {    # Predator species loop
        for(r_age in 1:n_ages) { # Predator age loop
          for(y in 1:n_yrs) {
            for(ksp in 1:nspp) {
              for(k_age in 1:n_ages) { # Prey age loop
                avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + suitability[rsp, ksp, r_age, k_age] * avgNAA[ksp, y, k_age] * WAA[ksp, k_age]
              }
            }
            # Other food
            avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + other_food[rsp] * suit_other[sp]
          }
        }
      }

      # Predation mortality ----
      M2_at_age[] <- 0
      B_eaten_as_prey <- array(0, dim = c(nspp, n_ages, n_yrs))
      diet_prop_hat <- array(0, dim = c(nspp, nspp, n_ages, n_ages, n_yrs))
      for(ksp in 1:nspp) {
        for(k_age in 1:n_ages) { # Prey age loop
          for(rsp in 1:nspp) {    # Predator species loop
            for(r_age in 1:n_ages) { # Predator age loop
              for(y in 1:n_yrs) {

                if(avail_food[rsp, r_age, y] > 0){

                  # MSVPA
                  # - M2
                  M2_at_age[ksp, k_age, y] = M2_at_age[ksp, k_age, y] + avgNAA[rsp, y, r_age] * ration[rsp, r_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]

                  # Biomass consumed as prey
                  B_eaten_as_prey[ksp, k_age, y] = B_eaten_as_prey[ksp, k_age, y] + avgNAA[ksp, y, k_age] * WAA[ksp, k_age] * avgNAA[rsp, y, r_age] * ration[rsp, r_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]

                  # Diet
                  diet_prop_hat[rsp, ksp, r_age, k_age, y] = avgNAA[ksp, y, k_age] * WAA[ksp, k_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]
                }
              }
            }
          }
        }
      }

      # End iterations
    }

    # Observation model
    ObsCatch <- Catch * rlnorm(n_yrs * nspp, 0, sigma_catch)

    # Survey observations
    SrvIdx <- srv_q * Total_Biom * rlnorm(n_yrs * nspp, 0, sigma_srv)

    # Age composition data (simplified multinomial)
    ObsFishAges <- array(0, dim=c(nspp, n_yrs, n_ages))
    ObsSrvAges <- array(0, dim=c(nspp, n_yrs, n_ages))

    for(sp in 1:nspp){
      for(y in 1:n_yrs) {
        # Fishery ages
        ObsFishAges[sp, y,] <- rmultinom(1, fish_ISS, CAA[sp, y,])
        # Survey ages
        ObsSrvAges[sp, y,] <- rmultinom(1, srv_ISS, NAA[sp, y,] * srv_sel[sp,])
      }
    }


    # Diet composition data
    ObsDiet <- array(0, dim=c(nspp, nspp, n_ages, n_yrs))
    for(rsp in 1:nspp){
      for(r_age in 1:n_ages){
        for(y in 1:n_yrs) {
          if(sum(diet_prop_hat[rsp, , r_age, , y]) > 0){
            ObsDiet[rsp, , r_age, y] = rowSums(diet_prop_hat[rsp, , r_age, , y]) # rmultinom(1, diet_ISS,  )[,1]
          }
        }
      }
    }

    # Return list of true and observed values
    return(list(
      NAA = NAA,
      CAA = CAA,
      FAA = FAA,
      SSB = SSB,
      Total_Biom = Total_Biom,
      B_eaten_as_prey = B_eaten_as_prey,
      Catch = Catch,
      ObsCatch = ObsCatch,
      SrvIdx = SrvIdx,
      ObsFishAges = ObsFishAges,
      ObsSrvAges = ObsSrvAges,
      fish_sel = fish_sel,
      srv_sel = srv_sel,
      WAA = WAA,
      MatAA = MatAA,
      M = M,
      srv_q = srv_q,
      rec_devs = rec_devs,
      init_devs = init_devs,
      ObsDiet = ObsDiet,
      vulnerability = vulnerability,
      suitability = suitability
    ))
  }

  # Set up simulation -------------------------------------------------------------
  nspp = 2
  nyrs = 20
  years <- 1:nyrs
  ages <- 1:15
  WAA <- 2 / (1 + exp(-0.8 * (ages - 3)))
  WAA2 <- 1.4 / (1 + exp(-1 * (ages - 3)))
  MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
  MatAA2 <- 1 / (1 + exp(-0.8 * (ages - 5)))
  sigma_R <- 0.3
  sigma_Catch <- 0.001
  sigma_SrvIdx <- 0.3
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
                        1 / (1 + exp(-1 * (ages - 3)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2)))), nspp, length(ages), byrow = TRUE),
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
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs+10
  simData$nsex <- rep(1, nspp)
  simData$nages <- rep(15, nspp)
  simData$minage <- rep(1, nspp)
  simData$nlengths <- rep(15, nspp)
  simData$pop_wt_index <- 1:nspp
  simData$ssb_wt_index <- 1:nspp
  simData$pop_age_transition_index <- rep(1, nspp)

  # * Fleet control ----
  simData$fleet_control <- simData$fleet_control[c(1,3,1,3),] # BT and Trawl Fishery are both simple logistic
  simData$fleet_control$Fleet_name <- c("Survey1", "Fishery1","Survey2", "Fishery2")
  simData$fleet_control$Fleet_code <- 1:4
  simData$fleet_control$Selectivity_index <- 1:4
  simData$fleet_control$Weight_index <- c(1,1,2,2)

  # * Index data ----
  simData$index_data <- rbind(
    data.frame(Fleet_name = "Survey1",
               Fleet_code = 1,
               Species = 1,
               Year = 1:20,
               Month = 0,
               Selectivity_block = 1,
               Q_block = 1,
               Observation = sim$SrvIdx[1,],
               Log_sd = sigma_SrvIdx),

    data.frame(Fleet_name = "Survey2",
               Fleet_code = 3,
               Species = 2,
               Year = 1:20,
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
               Year = 1:20,
               Month = 0,
               Selectivity_block = 1,
               Catch = sim$ObsCatch[1,],
               Log_sd = sigma_Catch),

    data.frame(Fleet_name = "Fishery2",
               Fleet_code = 4,
               Species = 2,
               Year = 1:20,
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
                                 Year = 1:20,
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
                                  Year = 1:20,
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
                                   Year = 1:20,
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
                                    Year = 1:20,
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

  simData$age_trans_matrix <- cbind(atf1, atf2)


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
  maturity2 <- cbind(data.frame(Species = 1),
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
  simData$env_data <- data.frame(Year = 1:20,
                                 EnvData = 1)


  # * Relative foraging rate (days) ----
  WAA <- as.data.frame(matrix(WAA * 50, ncol = 15))
  colnames(WAA) <- paste0("Age",1:15)
  Pyrs1 <- cbind(data.frame(Species = 1,
                            Sex = 0,
                            Year = 0),
                 WAA
  )

  WAA2 <- as.data.frame(matrix(WAA2 * 50, ncol = 15))
  colnames(WAA2) <- paste0("Age",1:15)
  Pyrs2 <- cbind(data.frame(Species = 2,
                            Sex = 0,
                            Year = 0),
                 WAA2
  )

  simData$Pyrs <- rbind(Pyrs1, Pyrs2)


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

  # * Diet -----
  ind = 1
  for(i in 1:dim(sim$ObsDiet)[1]){
    for(j in 1:dim(sim$ObsDiet)[2]){
      for(k in 1:dim(sim$ObsDiet)[3]){
        simData$diet_data[ind,] <- NA
        simData$diet_data$Pred[ind] <- i
        simData$diet_data$Prey[ind] <- i
        simData$diet_data$Pred_sex[ind] <- 0
        simData$diet_data$Prey_sex[ind] <- 0
        simData$diet_data$Pred_age[ind] <- k
        simData$diet_data$Prey_age[ind] <- -999
        simData$diet_data$Sample_size[ind] <- 200
        simData$diet_data$Stomach_proportion_by_weight[ind] <- mean(sim$ObsDiet[i,j,k,])
        ind = ind+1
      }
    }
  }


  # Fit Rceattle -------------------------------------------------------------
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

  plot(x = sim$SSB, y = ss_run$quantities$ssb[1,1:20]); abline(1,1)
  plot(x = sim$Total_Biom, y = ss_run$quantities$biomass[1,1:20], ylab = "Rceattle biomass", xlab = "True biomass"); abline(1,1)
  plot(x = sim$NAA[1:20,1], y = ss_run$quantities$R[1,1:20]); abline(1,1)

  expect_equal(sim$SSB[1], 45.80138, tolerance = 0.0001)
})
