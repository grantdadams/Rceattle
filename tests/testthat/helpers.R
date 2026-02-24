# Helpers for tests
# Minimal test data factory and small utilities used by testthat files

#' Make single-species test data set
#'
#' @param nyrs
#' @param nages
#' @param seed
make_test_data <- function(nyrs = 8, nages = 5, seed = NULL) {
  if (!requireNamespace("Rceattle", quietly = TRUE)) {
    stop("Rceattle package required for test helpers")
  }

  dat <- list()

  if (!is.null(seed)) set.seed(seed)
  nspp = 1
  dat$nspp <- nspp
  dat$styr <- 1
  dat$endyr <- nyrs
  dat$projyr <- nyrs + 2
  dat$spnames <- "Sim"
  dat$nsex = 1
  dat$spawn_month = 0
  dat$nages = nages
  dat$minage = 1
  dat$nlengths = nages
  dat$pop_wt_index = 1
  dat$ssb_wt_index = 1
  dat$alpha_wt_len = 0.0001
  dat$beta_wt_len = 3
  dat$pop_age_transition_index = 1
  dat$sigma_rec_prior = 1
  dat$other_food = 1e6
  dat$estDynamics = 0


  years <- seq(dat$styr, dat$endyr)

  # Fleet control
  # - Two simple fleets (survey + fishery)
  dat$fleet_control <- data.frame(
    Fleet_name = c("Survey", "Fishery"),
    Fleet_code = 1:2,
    Fleet_type = 2:1,
    Species = 1,
    Month = 0,
    Selectivity_index = 1:2,
    Selectivity = 1,
    N_sel_bins = NA,
    Sel_curve_pen1 = NA,
    Sel_curve_pen2 = NA,
    Time_varying_sel = 0,
    Time_varying_sel_sd_prior = 1,
    Bin_first_selected = 1,
    Sel_norm_bin1 = NA,
    Sel_norm_bin2 = NA,
    Comp_loglike = 0,
    Comp_weights = 1,
    CAAL_loglike = 0,
    Weight1_Numbers2 = 1,
    Weight_index = 1,
    Age_transition_index = 1,
    Q_index = c(1, NA),
    Catchability = c(0, NA),
    Q_prior = c(1, NA),
    Q_sd_prior = c(0.2, NA),
    Time_varying_q = c(0, NA),
    Time_varying_q_sd_prior = c(1, NA),
    Estimate_index_sd = c(0, NA),
    Index_sd_prior = c(1, NA),
    Estimate_catch_sd = c(NA, 0),
    Catch_sd_prior = c(NA, 1),
    proj_F_prop = c(NA, 1),
    CAAL_weights = 1,
    Est_weights_mcallister = 1
  )

  # Deterministic simple observations for fast tests
  total_biom <- rep(100.0, length(years))
  dat$index_data <- data.frame(
    Fleet_name = "Survey",
    Fleet_code = 1,
    Species = 1,
    Year = years,
    Month = 0,
    Selectivity_block = 1,
    Q_block = 1,
    Observation = total_biom,
    Log_sd = 0.1
  )

  dat$catch_data <- data.frame(
    Fleet_name = "Fishery",
    Fleet_code = 2,
    Species = 1,
    Year = years,
    Month = 0,
    Selectivity_block = 1,
    Catch = total_biom * 0.1,
    Log_sd = 0.05
  )

  # Minimal composition
  dat$comp_data <- data.frame(matrix(NA, nrow = 0, ncol = 8 + nages))
  colnames(dat$comp_data ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Age0_Length1", "Year", "Month", "Sample_size", paste("Comp_", 1:nages))

  # Minimal CAAL
  dat$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 7 + nages))
  colnames(dat$caal_data ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", paste("CAAL_", 1:nages))

  #  Empirical selectivity
  dat$emp_sel <- data.frame(matrix(NA, nrow = 0, ncol = 5 + nages))
  colnames(dat$emp_sel ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", paste("Comp_", 1:nages))

  # Input N-at-age
  dat$NByageFixed <- data.frame(matrix(NA, nrow = 0, ncol = 4 + nages))
  colnames(dat$NByageFixed ) = c("Species_name ", "Species", "Sex", "Year", paste("Age", 1:nages))

  # Age-transition
  age_transition <- as.data.frame(diag(nages))
  colnames(age_transition) <- paste0("Length_", seq_len(nages))
  dat$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                           Age_transition_index = 1,
                                           Species = 1,
                                           Sex = 0,
                                           Age = 1:nages),
                                age_transition)

  # Age-error
  age_error <- as.data.frame(diag(nages))
  colnames(age_error) <- paste0("Obs_age", seq_len(nages))
  dat$age_error <- cbind(data.frame(
    Species = 1,
    True_age = 1:nages),
    age_error)

  # Weight-at-age
  WAA <- rep(1, nages)
  WAA_df <- as.data.frame(matrix(WAA, nrow = 1))
  colnames(WAA_df) <- paste0("Age", seq_len(nages))
  dat$weight <- cbind(data.frame(Wt_name = "Base", Wt_index = 1, Species = 1, Sex = 0, Year = 0), WAA_df)

  # Maturity
  MatAA_df <- as.data.frame(matrix(1, nrow = 1, ncol = nages))
  colnames(MatAA_df) <- paste0("Age", seq_len(nages))
  dat$maturity <- cbind(data.frame(Species = 1), MatAA_df)

  # Sex ratio
  sexratio <- as.data.frame(matrix(0.5, nrow = 1, ncol = nages))
  colnames(sexratio) <- paste0("Age", seq_len(nages))
  dat$sex_ratio <- cbind(data.frame(Species = 1), sexratio)

  # Mortality
  mort <- as.data.frame(matrix(0.2, nrow = 1, ncol = nages))
  colnames(mort) <- paste0("Age", seq_len(nages))
  dat$M1_base <- cbind(data.frame(Species = 1, Sex = 0), mort)

  # Bioenergetics
  dat$Ceq = 1
  dat$Cindex = 1
  dat$Pvalue = 1
  dat$fday = 1
  dat$CA = 1
  dat$CB = 1
  dat$Qc = 1
  dat$Tco = 1
  dat$Tcm = 1
  dat$Tcl = 1
  dat$CK1 = 1
  dat$CK4 = 1
  dat$Diet_comp_weights = 1

  # Environmental data
  dat$env_data <- data.frame(
    Year = years,
    Index1 = rnorm(nyrs)
  )

  # Diet information Pyrs (relative foraging rate) ----
  dat$ration_data <- dat$weight %>%
    select(Species, Sex, Year, contains("Age"))


  # Diet proportion ----
  dat$diet_data <- as.data.frame(matrix(NA, nrow = 0, ncol = 9))
  colnames(dat$diet_data) = c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age",
                              "Year", "Sample_size", "Stomach_proportion_by_weight")


  # Clean and return
  dat <- Rceattle::clean_data(dat)
  return(dat)
}

compile_tmb_if_needed <- function() {
  compile_script <- file.path("src", "TMB", "compile.R")
  if (file.exists(compile_script)) {
    tryCatch(source(compile_script), error = function(e) stop("TMB compile failed: ", e$message))
  }
}

with_loaded_dll <- function(lib, code) {
  if (!file.exists(lib)) stop("DLL not found: ", lib)
  dll <- dyn.load(lib)
  on.exit(if (!is.null(dll)) dyn.unload(dll[["path"]]), add = TRUE)
  force(code)
}


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
  MAA <- array(0, dim = c(nspp, n_yrs, n_ages))      # Total mortality
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
        ZAA[sp, y,] <- FAA[sp,y,] + M[sp] # + M2_at_age[sp, , y]
        MAA[sp, y,] <- M[sp] + M2_at_age[sp, , y]

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
    diet_prop <- array(0, dim = c(nspp, nspp, n_ages, n_ages, n_yrs))
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
                diet_prop[rsp, ksp, r_age, k_age, y] = avgNAA[ksp, y, k_age] * WAA[ksp, k_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]
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
        if(sum(diet_prop[rsp, , r_age, , y]) > 0){
          ObsDiet[rsp, , r_age, y] = rowSums(diet_prop[rsp, , r_age, , y])
        }
      }
    }
  }

  # Return list of true and observed values ----
  return(list(
    NAA = NAA,
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
    diet_prop = diet_prop
  ))
}

