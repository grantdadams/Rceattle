test_that("lognormal index", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Specify rec
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Set logistic sel
  GOA2018SS$fleet_control$Selectivity <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Age_max_selected <- NA
  GOA2018SS$index_data$Month <- 0

  # Set params
  GOA2018SS$srr_fun <- 0
  GOA2018SS$initMode <- 1
  inits <- build_params(GOA2018SS)
  inits$rec_pars[,1] <- R0
  inits$R_ln_sd <- rep(0, 3)
  inits$ln_F[] <- -999 # No fishing
  inits$index_ln_q[] <- 0 # Set q to 1
  inits$index_ln_q[] <- 0 # Set sex to 1


  # Specify logistic selectivity
  inf = 0; alpha = 0
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)


  # Set logistic params
  inits$ln_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

  # Calculate SPR
  M <- 0.2
  wt <- 1
  sex_ratio <- 0.5
  pmature <- 1

  natage = rep(1, 10)
  for(age in 2:9){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[10] <-  natage[9] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0)

  # Check biomass
  for(sp in 1:3){
    testthat::expect_equal(as.numeric(ss_run$quantities$biomass[sp,1:nyrs]), rep(SB0[sp], nyrs), tolerance = 0.0001)
  }

  # Check derived index
  testthat::expect_equal(ss_run$quantities$index_hat*2, SB0[ss_run$data_list$index_data$Species])


})
