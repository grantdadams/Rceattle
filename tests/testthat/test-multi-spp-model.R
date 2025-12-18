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
    niter = 10,
    gam_a = rep(0.1, nspp),
    gam_b = rep(0.1, nspp),
    log_phi = matrix(0, nspp, nspp),
    other_food = rep(1e5, nspp),
    ration = matrix(1, nspp, length(ages))

  ) {

    n_yrs <- length(years)
    n_ages <- length(ages)

    # Initialize arrays
    NAA <- array(0, dim = c(nspp, n_ages, n_yrs))  # Numbers at age
    avgNAA <- array(0, dim = c(nspp, n_ages, n_yrs))  # Numbers at age
    ZAA <- array(0, dim = c(nspp, n_ages, n_yrs))      # Total mortality
    MAA <- array(0, dim = c(nspp, n_ages, n_yrs))      # Total mortality
    M2_at_age <- array(0, dim = c(nspp, n_ages, n_yrs))
    CAA <- array(0, dim = c(nspp, n_ages, n_yrs))      # Catch at age
    FAA <- array(0, dim = c(nspp, n_ages, n_yrs))      # Fishing mortality at age

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
      NAA[sp, init_age_idx + 1, 1] <- mean_Rec[sp] * exp(- (init_age_idx * M[sp]))
      NAA[sp, n_ages, 1] <- mean_Rec[sp] * exp(-(n_ages - 1) * M[sp]) / (1 - exp(-M[sp]))
      NAA[sp,2:n_ages, 1] <- NAA[sp,2:n_ages, 1] * exp(init_devs[sp,])
    }

    # Project population forward
    for(iter in 1:niter){

      for(y in 1:n_yrs) {
        for(sp in 1:nspp){
          # New recruits
          NAA[sp, 1, y] <- mean_Rec[sp] * exp(rec_devs[sp,y])

          # Calculate mortality
          FAA[sp, ,y] <- Fmort[sp, y] * fish_sel[sp,]
          ZAA[sp, ,y] <- FAA[sp,,y] + M[sp] + M2_at_age[sp, , y]
          MAA[sp, ,y] <- M[sp] + M2_at_age[sp, ,y]

          # Calculate catch
          CAA[sp, ,y] <- FAA[sp, ,y] / ZAA[sp, ,y] * NAA[sp, ,y] * (1 - exp(-ZAA[sp, ,y]))

          # Project survivors
          if(y < n_yrs) {
            for(a in 1:(n_ages-1)) {
              NAA[sp, a+1, y+1] <- NAA[sp, a, y] * exp(-ZAA[sp, a, y])
            }
            # Plus group
            NAA[sp, n_ages, y+1] <- NAA[sp, n_ages, y+1] + NAA[sp, n_ages, y] * exp(-ZAA[sp, n_ages, y])
          }
          avgNAA[sp,,y] = NAA[sp, ,y] * (1-exp(-ZAA[sp, ,y]))/ZAA[sp,,y]

          # Calculate annual metrics
          Total_Biom[sp, y] <- sum(NAA[sp, ,y] * WAA[sp, ])
          SSB[sp, y] <- sum(NAA[sp, ,y] * WAA[sp, ] * MatAA[sp, ]) * 0.5
          Catch[sp, y] <- sum(CAA[sp, ,y] * WAA[sp, ])
        }
      }

      # Available food ----
      avail_food <- array(0, dim = c(nspp, n_ages, n_yrs))
      for(rsp in 1:nspp) {    # Predator species loop
        for(r_age in 1:n_ages) { # Predator age loop
          for(y in 1:n_yrs) {
            for(ksp in 1:nspp) {
              for(k_age in 1:n_ages) { # Prey age loop
                avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + suitability[rsp, ksp, r_age, k_age] * avgNAA[ksp, k_age, y] * WAA[ksp, k_age]
              }
            }
            # Other food
            avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + other_food[rsp] * suit_other[rsp]
          }
        }
      }

      # Predation mortality ----
      M2_at_age[] <- 0
      B_eaten_as_prey <- array(0, dim = c(nspp, n_ages, n_yrs))
      diet_prop <- array(0, dim = c(nspp, nspp, n_ages, n_ages, n_yrs))
      for(ksp in 1:nspp) {
        for(k_age in 1:n_ages) { # Prey age loop
          for(rsp in 1:nspp) {    # Predator species loop
            for(r_age in 1:n_ages) { # Predator age loop
              for(y in 1:n_yrs) {

                if(avail_food[rsp, r_age, y] > 0){

                  # MSVPA
                  # - M2
                  M2_at_age[ksp, k_age, y] = M2_at_age[ksp, k_age, y] + avgNAA[rsp, r_age, y] * ration[rsp, r_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]

                  # Biomass consumed as prey
                  B_eaten_as_prey[ksp, k_age, y] = B_eaten_as_prey[ksp, k_age, y] + avgNAA[ksp, k_age, y] * WAA[ksp, k_age] * avgNAA[rsp, r_age, y] * ration[rsp, r_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]

                  # Diet
                  diet_prop[rsp, ksp, r_age, k_age, y] = avgNAA[ksp, k_age, y] * WAA[ksp, k_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]
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
    ObsFishAges <- array(0, dim=c(nspp, n_ages, n_yrs))
    ObsSrvAges <- array(0, dim=c(nspp, n_ages, n_yrs))

    for(sp in 1:nspp){
      for(y in 1:n_yrs) {
        # Fishery ages
        ObsFishAges[sp, ,y] <- rmultinom(1, fish_ISS, CAA[sp, ,y])
        # Survey ages
        ObsSrvAges[sp, ,y] <- rmultinom(1, srv_ISS, NAA[sp, ,y] * srv_sel[sp,])
      }
    }


    # Diet composition data
    ObsDiet <- array(0, dim=c(nspp, nspp, n_ages, n_yrs))
    for(rsp in 1:nspp){
      for(r_age in 1:n_ages){
        for(y in 1:n_yrs) {
          if(sum(diet_prop[rsp, , r_age, , y]) > 0){
            ObsDiet[rsp, , r_age, y] = rowSums(diet_prop[rsp, , r_age, , y]) # rmultinom(1, diet_ISS,  )[,1]
          }
        }
      }
    }

    # Return list of true and observed values ----
    return(list(
      NAA = NAA,
      avgNAA = avgNAA,
      CAA = CAA,
      FAA = FAA,
      MAA = MAA,
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
      M2_at_age = M2_at_age,
      vulnerability = vulnerability,
      suitability = suitability,
      suit_other = suit_other,
      avail_food = avail_food,
      diet_prop = diet_prop,
      ration = ration
    ))
  }



  # 1) Set up simulation -------------------------------------------------------------
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

  sim$M2_at_age[2,,]

  # * Plot ------------
  par(mfrow = c(4,1), mar = c(4,4,0.1,0))
  plot(y = sim$Total_Biom[1,], x = years, type = "l", ylab = "Species 1 B")
  plot(y = sim$Total_Biom[2,], x = years, type = "l", ylab = "Species 2 B")
  plot(y = colSums(sim$B_eaten_as_prey[1,,]), x = years, type = "l", ylab = "Species 1 B consumed")
  plot(y = colSums(sim$B_eaten_as_prey[2,,]), x = years, type = "l", ylab = "Species 2 B consumed")

  # 2) Set up Rceattle data -------------------------------------------------------------
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
  simData$fleet_control$Estimate_q
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
  tmp <- t(sim$ObsSrvAges[1,,])
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

  tmp <- t(sim$ObsSrvAges[2,,])
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
  tmp <- t(sim$ObsFishAges[1,,])
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

  tmp <- t(sim$ObsFishAges[2,,])
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

  # * Emp sel and N ----
  simData$emp_sel[] <- NA
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


  # * Ration ----
  ration_data1 <- WAA * 50
  colnames(ration_data1) <- paste0("Age",1:15)
  ration_data1 <- cbind(data.frame(Species = 1,
                            Sex = 0,
                            Year = 0),
                 ration_data1
  )

  ration_data2 <- WAA2 * 50
  colnames(ration_data2) <- paste0("Age",1:15)
  ration_data2 <- cbind(data.frame(Species = 2,
                            Sex = 0,
                            Year = 0),
                 ration_data2
  )

  simData$ration_data <- rbind(ration_data1, ration_data2)

  simData$Ceq <- rep(4,nspp)
  simData$Cindex <- rep(1, nspp)
  simData$Pvalue <- rep(1, nspp)
  simData$fday <- rep(1, nspp)
  simData$CA <- rep(1, nspp)
  simData$CB <- rep(-1, nspp)
  simData$Qc <- rep(1,nspp)
  simData$Tco <- rep(1, nspp)
  simData$Tcm <- rep(1, nspp)
  simData$Tcl <- rep(1, nspp)
  simData$CK1 <- rep(1, nspp)
  simData$CK4 <- rep(1, nspp)
  simData$Diet_comp_weights <- rep(1,nspp)


  # 3) Fit single-species -------------------------------------------------------------
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              initMode = 2,
                              verbose = 1)

  ss_run$quantities$ration[2,1,,]
  plot(x = sim$SSB[1,], y = ss_run$quantities$ssb[1,1:nyrs]); abline(1,1)
  plot(x = sim$SSB[2,], y = ss_run$quantities$ssb[2,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom[1,], y = ss_run$quantities$biomass[1,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom[2,], y = ss_run$quantities$biomass[2,1:nyrs]); abline(1,1)
  plot(x = sim$NAA[1,1,], y = ss_run$quantities$R[1,1:nyrs]); abline(1,1)
  plot(x = sim$NAA[2,1,], y = ss_run$quantities$R[2,1:nyrs]); abline(1,1)





  # 4)  Fit multi-species -------------------------------------------------------------
  # * Fix suitability -----
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

  inits <- ss_run$estimated_params
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


  plot(x = sim$NAA[1,1,], y = ms_run1$quantities$R[1,1:nyrs]); abline(1,1)
  plot(x = sim$NAA[2,1,], y = ms_run1$quantities$R[2,1:nyrs]); abline(1,1)

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
  ms_run1$quantities$sel[c(1,3),,,1]-sim$srv_sel
  ms_run1$quantities$sel[c(2,4),,,1]-sim$fish_sel

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
