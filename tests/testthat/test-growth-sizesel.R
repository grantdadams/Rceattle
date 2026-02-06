

test_that("Simulated simple model with growth curve and size-based logistic selectivity" {
  # Simulate Data -----------------------------------------------------------
  # Adapted from Matt Cheng
  sim_pop_model <- function(years,
                            ages,
                            lengths,
                            WAA,
                            growth_matrix,
                            MatAA,
                            mean_Rec,
                            sigma_R,
                            sigma_catch,
                            sigma_srv,
                            fish_ISS,
                            srv_ISS,
                            fish_CAAL_ISS,
                            srv_CAAL_ISS,
                            M,
                            fish_size_sel,
                            srv_size_sel,
                            Fmort,
                            srv_q) {

    n_yrs <- length(years)
    n_ages <- length(ages)
    n_lengths <- length(lengths)

    # Initialize arrays
    NAA <- array(0, dim = c(n_yrs + 1, n_ages))  # Numbers at age
    ZAA <- array(0, dim = c(n_yrs, n_ages))      # Total mortality
    CAA <- array(0, dim = c(n_yrs, n_ages))      # Catch at age
    FAA <- array(0, dim = c(n_yrs, n_ages))      # Fishing mortality at age

    # Population metrics
    SSB <- numeric(n_yrs)
    Total_Biom <- numeric(n_yrs)
    Catch <- numeric(n_yrs)

    # Generate recruitment deviations
    rec_devs <- rnorm(n_yrs, 0, sigma_R)
    init_devs <- rnorm(n_ages - 2, 0, sigma_R)

    # Do any growth transformations
    fish_sel <- growth_matrix %*% fish_size_sel
    srv_sel <- growth_matrix %*% srv_size_sel

    # Initialize population
    init_age_idx <- 1:(n_ages - 2)
    NAA[1, init_age_idx + 1] <- mean_Rec * exp(init_devs - (init_age_idx * M))
    NAA[1, n_ages] <- mean_Rec * exp(-(n_ages - 1) * M) / (1 - exp(-M))

    # Project population forward
    for(y in 1:n_yrs) {
      # New recruits
      NAA[y, 1] <- mean_Rec * exp(rec_devs[y])

      # Calculate mortality
      FAA[y,] <- Fmort[y] * fish_sel
      ZAA[y,] <- FAA[y,] + M

      # Calculate catch
      CAA[y,] <- FAA[y,] / ZAA[y,] * NAA[y,] * (1 - exp(-ZAA[y,]))

      # Project survivors
      if(y < n_yrs) {
        for(a in 1:(n_ages-1)) {
          NAA[y+1, a+1] <- NAA[y, a] * exp(-ZAA[y, a])
        }
        # Plus group
        NAA[y+1, n_ages] <- NAA[y+1, n_ages] + NAA[y, n_ages] * exp(-ZAA[y, n_ages])
      }

      # Calculate annual metrics
      Total_Biom[y] <- sum(NAA[y,] * WAA)
      SSB[y] <- sum(NAA[y,] * WAA * MatAA) * 0.5
      Catch[y] <- sum(CAA[y,] * WAA)
    }

    # Observation model
    ObsCatch <- Catch * rlnorm(n_yrs, 0, sigma_catch)

    # Survey observations
    SrvIdx <- srv_q * Total_Biom * rlnorm(n_yrs, 0, sigma_srv)

    # Age composition data (simplified multinomial)
    ObsFishAges <- array(0, dim=c(n_yrs, n_ages))
    ObsSrvAges <- array(0, dim=c(n_yrs, n_ages))

    for(y in 1:n_yrs) {
      # Fishery ages
      ObsFishAges[y,] <- rmultinom(1, fish_ISS, CAA[y,])
      # Survey ages
      ObsSrvAges[y,] <- rmultinom(1, srv_ISS, NAA[y,] * srv_sel)
    }

    # Age composition data (simplified multinomial)
    ObsFishCAAL <- array(0, dim=c(n_yrs, n_ages, n_lengths))
    ObsSrvCAAL <- array(0, dim=c(n_yrs, n_ages, n_lengths))

    for(y in 1:n_yrs) {
      for(ln in 1:n_lengths) {
        # Fishery CAAL
        ObsFishCAAL[y,,ln] <- rmultinom(1, fish_CAAL_ISS, CAA[y,] * growth_matrix[, ln])
        # Survey CAAL
        ObsSrvCAAL[y,,ln] <- rmultinom(1, srv_CAAL_ISS, NAA[y,] * srv_sel * growth_matrix[, ln])
      }
    }

    # Return list of true and observed values
    return(list(
      NAA = NAA,
      CAA = CAA,
      FAA = FAA,
      SSB = SSB,
      Total_Biom = Total_Biom,
      Catch = Catch,
      ObsCatch = ObsCatch,
      SrvIdx = SrvIdx,
      ObsFishAges = ObsFishAges,
      ObsSrvAges = ObsSrvAges,
      ObsFishCAAL = ObsFishCAAL,
      ObsSrvCAAL = ObsSrvCAAL,
      WAA = WAA,
      MatAA = MatAA,
      M = M,
      fish_sel = fish_sel,
      srv_sel = srv_sel,
      srv_q = srv_q,
      rec_devs = rec_devs,
      init_devs = init_devs
    ))
  }

  # Set up simulation -------------------------------------------------------------
  years <- 1:30
  nyrs <- length(years)
  ages <- 1:15
  nages <- length(ages)
  MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
  sigma_R <- 0.4
  sigma_Catch <- 0.001
  sigma_SrvIdx <- 0.3
  Fmort <- c(seq(0.02, 0.4, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))


  # - Growth
  nlengths <- 20
  lengths <- 5 * 1:nlengths

  # - Growth params: (sex, yr, param) -> params: K, L1, Linf, m
  gp <- array(rep(c(0.3, 4.5, 90, 1.0), each = 1), dim=c(1, 1, 4))
  gsd <- array(log(3), dim=c(1, 2))

  # - Run Growth Matrix
  gm <- get_growth_matrix_r(fracyr=0,
                            nsex_sp=1,
                            nages_sp=nages,
                            nlengths_sp = nlengths,
                            nyrs = 1,
                            lengths_sp= lengths,
                            minage_sp= 1,
                            maxage_sp = nages,
                            growth_params_sp = gp,
                            growth_ln_sd_sp = gsd,
                            growth_model_sp = 1)

  plot(lengths, gm$growth_matrix[1, 5, , 1], type="h", main="Length distribution for Age 5")

  # 4. Calculate Weight-at-Age
  lw_p <- array(c(0.00001, 3.0), dim=c(1, 1, 2)) # a=0.00001, b=3.0
  waa <- get_weight_at_age_r(nsex_sp = 1, nages_sp = nages, nlengths_sp = nlengths, nyrs = 1,
                             lengths_sp = lengths, growth_matrix = gm$growth_matrix, lw_params = lw_p)


  # First, simulate some data for the model
  set.seed(123)
  sim <- sim_pop_model(years = years,
                       ages = ages,
                       lengths = lengths,
                       WAA = waa[1,,1],
                       MatAA = MatAA,
                       growth_matrix = gm$growth_matrix[1,,,1],
                       mean_Rec = 50,
                       sigma_R = sigma_R,
                       sigma_catch = sigma_Catch,
                       sigma_srv = sigma_SrvIdx,
                       fish_ISS = 1e5,
                       srv_ISS = 1e5,
                       fish_CAAL_ISS = 1e5,
                       srv_CAAL_ISS = 1e5,
                       M = 0.3,
                       fish_size_sel = 1 / (1 + exp(-0.5 * (lengths - 60))),
                       srv_size_sel = 1 / (1 + exp(-1 * (lengths - 30))),
                       Fmort = Fmort,
                       srv_q = 1)


  # Set up Rceattle data -------------------------------------------------------------
  library(Rceattle)
  data("GOAcod")
  simData <- GOAcod

  # * Data controls ----
  simData$nspp <- 1
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs + 10
  simData$nsex <- 1
  simData$nages <- nages
  simData$minage <- 1
  simData$nlengths <- nlengths
  simData$pop_wt_index <- 1
  simData$ssb_wt_index <- 1
  simData$alpha_wt_len <- 0.00001
  simData$beta_wt_len <- 3
  simData$pop_age_transition_index <- 1

  # * Fleet control ----
  simData$fleet_control <- simData$fleet_control[c(1,3),] # BT and Trawl Fishery are both simple logistic
  simData$fleet_control$Fleet_name <- c("Survey", "Fishery")
  simData$fleet_control$Fleet_code <- 1:2
  simData$fleet_control$Selectivity <- 6
  simData$fleet_control$Selectivity_index <- 1:2
  simData$fleet_control$Weight_index <- 1
  simData$fleet_control$Bin_first_selected <- 1

  # * Index data ----
  simData$index_data <- data.frame(Fleet_name = "Survey",
                                   Fleet_code = 1,
                                   Species = 1,
                                   Year = years,
                                   Month = 0,
                                   Selectivity_block = 1,
                                   Q_block = 1,
                                   Observation = sim$SrvIdx,
                                   Log_sd = sigma_SrvIdx)

  # * Catch data ----
  simData$catch_data <- data.frame(Fleet_name = "Fishery",
                                   Fleet_code = 2,
                                   Species = 1,
                                   Year = years,
                                   Month = 0,
                                   Selectivity_block = 1,
                                   Catch = sim$ObsCatch,
                                   Log_sd = sigma_Catch)

  # * Comp data ----
  # - Index
  tmp <- sim$ObsSrvAges
  colnames(tmp) <- paste0("Comp_",1:nages)
  index_comp <- cbind(data.frame(Fleet_name = "Survey",
                                 Fleet_code = 1,
                                 Species = 1,
                                 Sex = 0,
                                 Age0_Length1 = 0,
                                 Year = -(years), # Negative turns data off
                                 Month = 0,
                                 Sample_size = rowSums(tmp)),
                      tmp
  )

  # - Fishery
  tmp <- sim$ObsFishAges
  colnames(tmp) <- paste0("Comp_",1:nages)
  fishery_comp <- cbind(data.frame(Fleet_name = "Fishery",
                                   Fleet_code = 2,
                                   Species = 1,
                                   Sex = 0,
                                   Age0_Length1 = 0,
                                   Year = -(years), # Negative turns data off
                                   Month = 0,
                                   Sample_size = rowSums(tmp)),
                        tmp
  )

  simData$comp_data <- rbind(index_comp, fishery_comp)


  # * CAAL data ----
  # - Index
  index_caal <- list()
  for(yr in 1:nyrs){
    tmp <- t(sim$ObsSrvCAAL[yr,,])
    colnames(tmp) <- paste0("CAAL_",1:nages)
    index_caal[[yr]] <- cbind(data.frame(Fleet_name = "Survey",
                                         Fleet_code = 1,
                                         Species = 1,
                                         Sex = 0,
                                         Year = yr,
                                         Length = lengths,
                                         Sample_size = rowSums(tmp)),
                              tmp
    )
  }
  index_caal <- do.call(rbind, index_caal)

  # - Fishery

  fish_caal <- list()
  for(yr in 1:nyrs){
    tmp <- t(sim$ObsFishCAAL[yr,,])
    colnames(tmp) <- paste0("CAAL_",1:nages)
    fish_caal[[yr]] <- cbind(data.frame(Fleet_name = "Fishery",
                                        Fleet_code = 2,
                                        Species = 1,
                                        Sex = 0,
                                        Year = yr,
                                        Length = lengths,
                                        Sample_size = rowSums(tmp)),
                             tmp
    )
  }
  fish_caal <- do.call(rbind, fish_caal)

  simData$caal_data <- rbind(index_caal, fish_caal)

  # * Empirical selectivity ----
  simData$emp_sel[] <- NA

  # * Fixed numbers ----
  simData$NByageFixed[] <- NA


  # * Age transition matrix ----
  tmp <- as.data.frame(diag(1,nlengths))[1:15,]
  colnames(tmp) <- paste0("Length_",1:nlengths)
  simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                               Age_transition_index = 1,
                                               Species = 1,
                                               Sex = 0,
                                               Age = 1:nages),
                                    tmp
  )


  # * Age error ----
  tmp <- as.data.frame(diag(1,nages))
  colnames(tmp) <- paste0("Obs_age",1:nages)
  simData$age_error <- cbind(data.frame(Species = 1,
                                        True_age = 1:nages),
                             tmp
  )

  # * Weight-at-age ----
  WAA <- as.data.frame(matrix(waa[1,,1], ncol = nages))
  colnames(WAA) <- paste0("Age",1:nages)
  simData$weight <- cbind(data.frame(Wt_name = "Base",
                                     Wt_index = 1,
                                     Species = 1,
                                     Sex = 0,
                                     Year = 0),
                          WAA
  )

  # * Maturity ----
  MatAA <- as.data.frame(matrix(MatAA, ncol = nages))
  colnames(MatAA) <- paste0("Age",1:nages)
  simData$maturity <- cbind(data.frame(Species = 1),
                            MatAA
  )

  # * Sex ratio ----
  sexratio <- as.data.frame(matrix(0.5, ncol = nages))
  colnames(sexratio) <- paste0("Age",1:nages)
  simData$sex_ratio <- cbind(data.frame(Species = 1),
                             sexratio
  )

  # * Mortality ----
  mort <- as.data.frame(matrix(0.3, ncol = nages))
  colnames(mort) <- paste0("Age",1:nages)
  simData$M1_base <- cbind(data.frame(Species = 1,
                                      Sex = 0),
                           mort
  )

  # * Environmental data ----
  simData$env_data <- data.frame(Year = years,
                                 EnvData = 1)


  # * Ration data ----
  pyrs <- as.data.frame(matrix(1, nrow = nyrs, ncol = nages))
  colnames(pyrs) <- paste0("Age",1:nages)
  simData$ration_data <- cbind(data.frame(Species = 1,
                                          Sex = 0,
                                          Year = years),
                               pyrs
  )


  # Fit Rceattle -------------------------------------------------------------
  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Estimate
                              growthFun = build_growth(growth_model = 0),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = FALSE,
                              verbose = 1)
  ss_run$quantities$jnll_comp

  # Set up inits
  inits <- ss_run$estimated_params
  map = ss_run$map
  map$mapList$ln_growth_pars[1,1,] <- c(1, NA, 2, NA)
  map$mapFactor$ln_growth_pars <- factor(map$mapList$ln_growth_pars)
  inits$ln_growth_pars[1,1,] <- log(c(0.3, 4.5, 90, 1.0))
  inits$growth_ln_sd[1,1,] <- log(3)
  ss_run_init1 <- Rceattle::fit_mod(data_list = simData,
                                    inits = inits, # Initial parameters = 0
                                    map = map,
                                    file = NULL, # Don't save
                                    estimateMode = 0, # Estimate
                                    growthFun = build_growth(growth_model = 1),
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 0, # Single species mode
                                    phase = TRUE,
                                    verbose = 1)


  inits$ln_F[2,] <- log(c(seq(0.02, 0.3, length.out = 10), seq(0.3, 0.05, length.out = 10)))
  inits$rec_pars[1,1] <- log(50)
  inits$rec_dev[1,1:nyrs] <- sim$rec_devs
  inits$init_dev[1,1:(nages-2)] <- sim$init_devs
  inits$ln_sel_slp[1,,1] <- log(c(1, 0.5))
  inits$sel_inf[1,,1] <- c(30, 60)
  ss_run_init2 <- Rceattle::fit_mod(data_list = simData,
                                   inits = inits, # Initial parameters = 0
                                   map = ss_run$map,
                                   file = NULL, # Don't save
                                   estimateMode = 3, # Estimate
                                   growthFun = build_growth(growth_model = 1),
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 0, # Single species mode
                                   phase = TRUE,
                                   verbose = 1)

  ss_run_init1$quantities$jnll_comp
  ss_run_init2$quantities$jnll_comp
  sum(ss_run_init1$quantities$jnll_comp)
  sum(ss_run_init2$quantities$jnll_comp)


  plot(x = sim$SSB, y = ss_run_init1$quantities$ssb[1,1:nyrs]); abline(1,1)
  plot(x = sim$Total_Biom, y = ss_run_init$quantities$biomass[1,1:nyrs], ylab = "Rceattle biomass", xlab = "True biomass"); abline(1,1)
  plot(x = sim$NAA[1:nyrs,1], y = ss_run_init$quantities$R[1,1:nyrs]); abline(1,1)

  ss_run_init$quantities$weight_hat[4,1,,1]
  sim$WAA

  # 1. Check size-selectivity ----
  fish_size_sel = 1 / (1 + exp(-0.5 * (lengths - 60)))
  srv_size_sel = 1 / (1 + exp(-1 * (lengths - 30)))

  testthat::expect_equal(as.numeric(srv_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[1,1,,1]))
  testthat::expect_equal(as.numeric(fish_size_sel), as.numeric(ss_run_init$quantities$sel_at_length[2,1,,1]))


  # 2. Check age-selectivity ----
  testthat::expect_equal(as.numeric(sim$srv_sel), as.numeric(ss_run_init$quantities$sel_at_age[1,1,,1]))
  testthat::expect_equal(as.numeric(sim$fish_sel), as.numeric(ss_run_init$quantities$sel_at_age[2,1,,1]))

})
