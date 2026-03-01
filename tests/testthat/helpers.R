# Helpers for tests
# Minimal test data factory and small utilities used by
calc_nll_ar1_1d <- function(x, sd, rho) {
  n <- length(x)
  Sigma_M <- sqrt(sd^2 / (1 - rho^2))

  # 2. Build the Correlation Matrix for an AR1 process
  times <- 1:n
  distance_matrix <- abs(outer(times, times, "-"))
  Correlation_Matrix <- rho^distance_matrix

  # 3. Scale Correlation to Covariance using the Marginal Variance
  Covariance_Matrix <- (Sigma_M^2) * Correlation_Matrix

  # 4. Calculate Negative Log-Likelihood
  # (mean is a vector of 0s because these are mean-zero deviations)
  nll <- -mvtnorm::dmvnorm(x, mean = rep(0, n), sigma = Covariance_Matrix, log = TRUE)

  return(nll)
}

#' Make single-species test data set
#'
#' @param nyrs
#' @param nages
#' @param seed
make_test_data <- function(nyrs = 8, nages = 5, seed = NULL) {
  if (!requireNamespace("Rceattle", quietly = TRUE)) {
    stop("Rceattle package required for test helpers")
  }

  simData <- list()

  if (!is.null(seed)) set.seed(seed)
  nspp = 1
  simData$nspp <- nspp
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs + 2
  simData$spnames <- "Sim"
  simData$nsex = 1
  simData$spawn_month = 0
  simData$nages = nages
  simData$minage = 1
  simData$nlengths = nages
  simData$pop_wt_index = 1
  simData$ssb_wt_index = 1
  simData$alpha_wt_len = 0.0001
  simData$beta_wt_len = 3
  simData$pop_age_transition_index = 1
  simData$sigma_rec_prior = 1
  simData$other_food = 1e6
  simData$estDynamics = 0


  years <- seq(simData$styr, simData$endyr)

  # Fleet control
  # - Two simple fleets (survey + fishery)
  simData$fleet_control <- data.frame(
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
  simData$index_data <- data.frame(
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

  simData$catch_data <- data.frame(
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
  simData$comp_data <- data.frame(matrix(NA, nrow = 0, ncol = 8 + nages))
  colnames(simData$comp_data ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Age0_Length1", "Year", "Month", "Sample_size", paste("Comp_", 1:nages))

  # Minimal CAAL
  simData$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 7 + nages))
  colnames(simData$caal_data ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", paste("CAAL_", 1:nages))

  #  Empirical selectivity
  simData$emp_sel <- data.frame(matrix(NA, nrow = 0, ncol = 5 + nages))
  colnames(simData$emp_sel ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", paste("Comp_", 1:nages))

  # Input N-at-age
  simData$NByageFixed <- data.frame(matrix(NA, nrow = 0, ncol = 4 + nages))
  colnames(simData$NByageFixed ) = c("Species_name ", "Species", "Sex", "Year", paste("Age", 1:nages))

  # Age-transition
  age_transition <- as.data.frame(diag(nages))
  colnames(age_transition) <- paste0("Length_", seq_len(nages))
  simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                               Age_transition_index = 1,
                                               Species = 1,
                                               Sex = 0,
                                               Age = 1:nages),
                                    age_transition)

  # Age-error
  age_error <- as.data.frame(diag(nages))
  colnames(age_error) <- paste0("Obs_age", seq_len(nages))
  simData$age_error <- cbind(data.frame(
    Species = 1,
    True_age = 1:nages),
    age_error)

  # Weight-at-age
  WAA <- rep(1, nages)
  WAA_df <- as.data.frame(matrix(WAA, nrow = 1))
  colnames(WAA_df) <- paste0("Age", seq_len(nages))
  simData$weight <- cbind(data.frame(Wt_name = "Base", Wt_index = 1, Species = 1, Sex = 0, Year = 0), WAA_df)

  # Maturity
  MatAA_df <- as.data.frame(matrix(1, nrow = 1, ncol = nages))
  colnames(MatAA_df) <- paste0("Age", seq_len(nages))
  simData$maturity <- cbind(data.frame(Species = 1), MatAA_df)

  # Sex ratio
  sexratio <- as.data.frame(matrix(0.5, nrow = 1, ncol = nages))
  colnames(sexratio) <- paste0("Age", seq_len(nages))
  simData$sex_ratio <- cbind(data.frame(Species = 1), sexratio)

  # Mortality
  mort <- as.data.frame(matrix(0.2, nrow = 1, ncol = nages))
  colnames(mort) <- paste0("Age", seq_len(nages))
  simData$M1_base <- cbind(data.frame(Species = 1, Sex = 0), mort)

  # Bioenergetics
  simData$Ceq = 1
  simData$Cindex = 1
  simData$Pvalue = 1
  simData$fday = 1
  simData$CA = 1
  simData$CB = 1
  simData$Qc = 1
  simData$Tco = 1
  simData$Tcm = 1
  simData$Tcl = 1
  simData$CK1 = 1
  simData$CK4 = 1
  simData$Diet_comp_weights = 1

  # Environmental data
  simData$env_data <- data.frame(
    Year = years,
    Index1 = rnorm(nyrs)
  )

  # Diet information Pyrs (relative foraging rate) ----
  simData$ration_data <- simData$weight |>
    dplyr::select(Species, Sex, Year, contains("Age"))


  # Diet proportion ----
  simData$diet_data <- as.data.frame(matrix(NA, nrow = 0, ncol = 9))
  colnames(simData$diet_data) = c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age",
                                  "Year", "Sample_size", "Stomach_proportion_by_weight")


  # Clean and return
  simData <- Rceattle::clean_data(simData)
  return(simData)
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
make_msm_test_data <- function(
    nspp = 2,
    years = 1:30,
    ages = 1:15,
    WAA = matrix(c(2 / (1 + exp(-0.8 * (ages - 3))),
                   1.4 / (1 + exp(-1 * (ages - 3)))),
                 nrow = nspp, ncol = length(ages), byrow = TRUE),
    MatAA = matrix(c(1 / (1 + exp(-1 * (ages - 5))),
                     1 / (1 + exp(-0.8 * (ages - 5)))),
                   nspp, length(ages), byrow = TRUE),
    mean_Rec = c(1e2, 1e3),
    sigma_R = 1,
    sigma_catch = 0.001,
    sigma_srv = 0.05,
    diet_ISS = 1e5,
    fish_ISS = 1e5,
    srv_ISS = 1e5,
    M = c(0.2, 0.3),
    rhoM = 0,
    sigmaM = 0,
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-2.5 * (ages - 4)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2),
                     seq(0.02, 0.3, length.out = nyrs)),
                   nrow = nspp, ncol = length(years), byrow = TRUE),
    srv_q = rep(1, nspp),

    # Multispecies bits
    niter = 5,
    gam_a = c(1, 0.1),
    gam_b =  rep(0.15, nspp),
    log_phi =  matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE),
    other_food =  rep(1e5, nspp),
    ration = WAA * 50,
    normalize_suitability = FALSE # NEW: Flag to force suitability to sum to 1.0

) {

  nyrs <- length(years)
  nages <- length(ages)

  # Initialize arrays
  NAA <- array(0, dim = c(nspp, nages, nyrs))  # Numbers at age
  avgNAA <- array(0, dim = c(nspp, nages, nyrs))  # Numbers at age
  ZAA <- array(0, dim = c(nspp, nages, nyrs))      # Total mortality
  MAA <- array(0, dim = c(nspp, nages, nyrs))      # Total mortality
  M2_at_age <- array(0, dim = c(nspp, nages, nyrs))
  CAA <- array(0, dim = c(nspp, nages, nyrs))      # Catch at age
  FAA <- array(0, dim = c(nspp, nages, nyrs))      # Fishing mortality at age

  # Population metrics
  SSB <- matrix(0, nspp, nyrs)
  Total_Biom <- matrix(0, nspp, nyrs)
  Catch <- matrix(0, nspp, nyrs)

  # Generate M deviations
  M_vec <- matrix(arima.sim(n = nyrs * nspp, list(order=c(1,0,0), ar=rhoM)
                            , sd = sigmaM), nspp, nyrs)
  Mtv = M * exp(M_vec)

  # Generate recruitment deviations
  rec_devs <- matrix(rnorm(nyrs * nspp, 0, sigma_R), nspp, nyrs)
  init_devs <- matrix(rnorm((nages-1) * nspp, 0, sigma_R), nspp, nages-1)

  # Vulnerability
  suitability = array(0, dim = c(nspp, nspp, nages, nages))
  vulnerability <- matrix(0, nspp, nspp)
  suit_other <- c()

  for(sp in 1:nspp) {
    vulnerability[sp,] = exp(log_phi[sp,])/(1+sum(exp(log_phi[sp,]))) # multinomial logistic transformation
    suit_other[sp] = 1 - sum(vulnerability[sp,])
  }

  # Calculate suitability ----
  for(sp in 1:nspp) { # Predator
    for(r_age in 1:nages){  # Pred age
      for(ksp in 1:nspp) {   # Prey loop
        for(k_age in 1:nages){ # Prey age

          # Weight-based lognormal suitability
          log_size_ratio = log(WAA[sp,r_age] / WAA[ksp,k_age]) # Log ratio of weights

          if(log_size_ratio > 0){
            suitability[sp, ksp, r_age, k_age] = vulnerability[sp, ksp] * dnorm(log_size_ratio, gam_a[sp], gam_b[sp]) / dnorm(gam_a[sp], gam_a[sp], gam_b[sp]) # Divide by mode to scale to 1
          }
        }
      }
    }
  }

  # NEW: Normalize suitability to sum to 1.0 across all prey for MSVPA testing
  suit_other_mat <- matrix(0, nspp, nages)
  for(sp in 1:nspp) suit_other_mat[sp, ] <- suit_other[sp]

  if(normalize_suitability) {
    for(sp in 1:nspp) {
      for(r_age in 1:nages) {
        S_tot = sum(suitability[sp, , r_age, ]) + suit_other_mat[sp, r_age]
        if(S_tot > 0) {
          suitability[sp, , r_age, ] = suitability[sp, , r_age, ] / S_tot
          suit_other_mat[sp, r_age] = suit_other_mat[sp, r_age] / S_tot
        }
      }
    }
  }

  # Initialize population
  for(sp in 1:nspp){
    init_age_idx <- 1:(nages - 2)
    NAA[sp, init_age_idx + 1, 1] <- mean_Rec[sp] * exp(- (init_age_idx * Mtv[sp, 1]))
    NAA[sp, nages, 1] <- mean_Rec[sp] * exp(-(nages - 1) * Mtv[sp, 1]) / (1 - exp(-Mtv[sp, 1]))
    NAA[sp,2:nages, 1] <- NAA[sp,2:nages, 1] * exp(init_devs[sp,])
  }

  # Project population forward
  for(iter in 1:niter){

    for(y in 1:nyrs) {
      for(sp in 1:nspp){
        # New recruits
        NAA[sp, 1, y] <- mean_Rec[sp] * exp(rec_devs[sp,y])

        # Calculate mortality
        FAA[sp, ,y] <- Fmort[sp, y] * fish_sel[sp,]
        ZAA[sp, ,y] <- FAA[sp,,y] + Mtv[sp, y] + M2_at_age[sp, , y]
        MAA[sp, ,y] <- Mtv[sp, y] + M2_at_age[sp, ,y]

        # Calculate catch
        CAA[sp, ,y] <- FAA[sp, ,y] / ZAA[sp, ,y] * NAA[sp, ,y] * (1 - exp(-ZAA[sp, ,y]))

        # Project survivors
        if(y < nyrs) {
          for(a in 1:(nages-1)) {
            NAA[sp, a+1, y+1] <- NAA[sp, a, y] * exp(-ZAA[sp, a, y])
          }
          # Plus group
          NAA[sp, nages, y+1] <- NAA[sp, nages, y+1] + NAA[sp, nages, y] * exp(-ZAA[sp, nages, y])
        }
        avgNAA[sp,,y] = NAA[sp, ,y] * (1-exp(-ZAA[sp, ,y]))/ZAA[sp,,y]

        # Calculate annual metrics
        Total_Biom[sp, y] <- sum(NAA[sp, ,y] * WAA[sp, ])
        SSB[sp, y] <- sum(NAA[sp, ,y] * WAA[sp, ] * MatAA[sp, ]) * 0.5
        Catch[sp, y] <- sum(CAA[sp, ,y] * WAA[sp, ])
      }
    }

    # Available food ----
    avail_food <- array(0, dim = c(nspp, nages, nyrs))
    for(rsp in 1:nspp) {    # Predator species loop
      for(r_age in 1:nages) { # Predator age loop
        for(y in 1:nyrs) {
          for(ksp in 1:nspp) {
            for(k_age in 1:nages) { # Prey age loop
              avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + suitability[rsp, ksp, r_age, k_age] * avgNAA[ksp, k_age, y] * WAA[ksp, k_age]
            }
          }
          # Other food (NEW logic handling normalization switch)
          if(normalize_suitability) {
            avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + other_food[rsp] * suit_other_mat[rsp, r_age]
          } else {
            avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + other_food[rsp] * suit_other[rsp]
          }
        }
      }
    }

    # Predation mortality ----
    M2_at_age[] <- 0
    B_eaten_as_prey <- array(0, dim = c(nspp, nages, nyrs))
    diet_prop <- array(0, dim = c(nspp, nspp, nages, nages, nyrs))
    for(ksp in 1:nspp) {
      for(k_age in 1:nages) { # Prey age loop
        for(rsp in 1:nspp) {    # Predator species loop
          for(r_age in 1:nages) { # Predator age loop
            for(y in 1:nyrs) {

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
  ObsCatch <- Catch * rlnorm(nyrs * nspp, 0, sigma_catch)
  SrvIdx <- srv_q * Total_Biom * rlnorm(nyrs * nspp, 0, sigma_srv)

  # Age composition data (simplified multinomial)
  ObsFishAges <- array(0, dim=c(nspp, nages, nyrs))
  ObsSrvAges <- array(0, dim=c(nspp, nages, nyrs))

  for(sp in 1:nspp){
    for(y in 1:nyrs) {
      ObsFishAges[sp, ,y] <- rmultinom(1, fish_ISS, CAA[sp, ,y])
      ObsSrvAges[sp, ,y] <- rmultinom(1, srv_ISS, NAA[sp, ,y] * srv_sel[sp,])
    }
  }

  # Export the array instead of the vector if normalized
  if(normalize_suitability) {
    suit_other <- suit_other_mat
  }

  # Set up Rceattle data -------------------------------------------------------------
  simData <- list()

  # * Data controls
  simData$nspp <- nspp
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs+10
  simData$spnames <- paste0("Species",1:nspp)
  simData$nsex <- rep(1, nspp)
  simData$spawn_month <- rep(0, nspp)
  simData$nages <- rep(nages, nspp)
  simData$minage <- rep(1, nspp)
  simData$nlengths <- rep(nages, nspp)
  simData$estDynamics <- rep(0, nspp)
  simData$pop_wt_index <- 1:nspp
  simData$ssb_wt_index <- 1:nspp
  simData$alpha_wt_len = rep(0.0001, nspp)
  simData$beta_wt_len = rep(3, nspp)
  simData$pop_age_transition_index <- rep(1, nspp)
  simData$sigma_rec_prior = rep(1, nspp)
  simData$other_food <- rep(1e5, nspp)
  simData$estDynamics = rep(0, nspp)

  # * Fleet control
  simData$fleet_control <- data.frame(
    Fleet_name = paste0(c("Survey", "Fishery"), rep(paste(" Species", 1:nspp), each = 2)),
    Fleet_code = rep(1:2, nspp),
    Fleet_type = rep(2:1, nspp),
    Species = rep(1:nspp, each = 2),
    Month = 0,
    Selectivity_index = 1:(nspp * 2),
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
    Weight_index = rep(1:nspp, each = 2),
    Age_transition_index = 1,
    Q_index = 1:(nspp * 2),
    Catchability = rep(c(0, NA), nspp),
    Q_prior = rep(c(1, NA), nspp),
    Q_sd_prior = rep(c(0.2, NA), nspp),
    Time_varying_q = rep(c(0, NA), nspp),
    Time_varying_q_sd_prior = rep(c(1, NA), nspp),
    Estimate_index_sd = rep(c(0, NA), nspp),
    Index_sd_prior = rep(c(1, NA), nspp),
    Estimate_catch_sd = rep(c(NA, 0), nspp),
    Catch_sd_prior = rep(c(NA, 1), nspp),
    proj_F_prop = rep(c(NA, 1), nspp),
    CAAL_weights = 1,
    Est_weights_mcallister = 1
  )

  # * Index data
  simData$index_data <- data.frame(
    Fleet_name = rep(paste("Survey Species", 1:nspp), each = nyrs),
    Fleet_code = rep(seq(1, nspp * 2, by = 2), each = nyrs),
    Species = rep(1:nspp, each = nyrs),
    Year = rep(years, nspp),
    Month = 0,
    Selectivity_block = 1,
    Q_block = 1,
    Observation = as.numeric(t(SrvIdx)),
    Log_sd = sigma_srv
  )

  # * Catch data
  simData$catch_data <- data.frame(
    Fleet_name = rep(paste("Fishery Species", 1:nspp), each = nyrs),
    Fleet_code = rep(seq(2, nspp * 2, by = 2), each = nyrs),
    Species = rep(1:nspp, each = nyrs),
    Year = rep(years, nspp),
    Month = 0,
    Selectivity_block = 1,
    Catch = as.numeric(t(ObsCatch)),
    Log_sd = sigma_catch
  )

  # * Comp data
  # - Index
  tmp <- do.call(rbind, lapply(seq_len(dim(ObsSrvAges)[1]), function(i) t(ObsSrvAges[i, , ])))
  colnames(tmp) <- paste0("Comp_",ages)
  index_comp <- cbind(
    data.frame(
      Fleet_name = rep(paste("Survey Species", 1:nspp), each = nyrs),
      Fleet_code = rep(seq(1, nspp * 2, by = 2), each = nyrs),
      Species = rep(1:nspp, each = nyrs),
      Sex = 0,
      Age0_Length1 = 0,
      Year = rep(years, nspp),
      Month = 0,
      Sample_size = rowSums(tmp)),
    tmp
  )

  # - Fishery
  tmp <- do.call(rbind, lapply(seq_len(dim(ObsSrvAges)[1]), function(i) t(ObsSrvAges[i, , ])))
  colnames(tmp) <- paste0("Comp_",ages)
  fishery_comp <- cbind(
    data.frame(
      Fleet_name = rep(paste("Fishery Species", 1:nspp), each = nyrs),
      Fleet_code = rep(seq(2, nspp * 2, by = 2), each = nyrs),
      Species = rep(1:nspp, each = nyrs),
      Sex = 0,
      Age0_Length1 = 0,
      Year = rep(years, nspp),
      Month = 0,
      Sample_size = rowSums(tmp)),
    tmp
  )

  simData$comp_data <- rbind(index_comp, fishery_comp)

  # Minimal CAAL
  simData$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 7 + nages))
  colnames(simData$caal_data ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", paste("CAAL_", 1:nages))

  #  Empirical selectivity
  simData$emp_sel <- data.frame(matrix(NA, nrow = 0, ncol = 5 + nages))
  colnames(simData$emp_sel ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", paste("Comp_", 1:nages))

  # Input N-at-age
  simData$NByageFixed <- data.frame(matrix(NA, nrow = 0, ncol = 4 + nages))
  colnames(simData$NByageFixed ) = c("Species_name ", "Species", "Sex", "Year", paste("Age", 1:nages))

  # * Age transition matrix
  age_trans_list <- lapply(1:nspp, function(sp) {
    tmp <- as.data.frame(diag(nages))
    colnames(tmp) <- paste0("Length_",1:nages)
    cbind(
      data.frame(Age_transition_name = paste("Species", sp),
                 Age_transition_index = sp,
                 Species = sp,
                 Sex = 0,
                 Age = 1:nages),
      tmp
    )
  })
  simData$age_trans_matrix <- do.call(rbind, age_trans_list)

  # * Age error
  age_error_list <- lapply(1:nspp, function(sp) {
    tmp <- as.data.frame(diag(1,nages))
    colnames(tmp) <- paste0("Obs_age",1:nages)
    cbind(data.frame(Species = sp,
                     True_age = 1:nages),
          tmp
    )
  })
  simData$age_error <- do.call(rbind, age_error_list)

  # * Weight-at-age
  weight_list <- lapply(1:nspp, function(sp) {
    weight_tmp <- as.data.frame(matrix(WAA[sp,], ncol = nages))
    colnames(weight_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(Wt_name = paste("Species", sp),
                     Wt_index = sp,
                     Species = sp,
                     Sex = 0,
                     Year = 0),
          weight_tmp
    )
  })
  simData$weight <- do.call(rbind, weight_list)

  # * Maturity
  maturity_list <- lapply(1:nspp, function(sp) {
    mat_tmp <- as.data.frame(matrix(MatAA[sp,], ncol = nages))
    colnames(mat_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(Species = sp),
          mat_tmp
    )
  })
  simData$maturity <- do.call(rbind, maturity_list)

  # * Sex ratio
  sexratio <- as.data.frame(matrix(0.5, nrow = nspp, ncol = nages))
  colnames(sexratio) <- paste0("Age",1:nages)
  simData$sex_ratio <- cbind(data.frame(Species = 1:nspp),
                             sexratio
  )

  # * Mortality
  mort_list <- lapply(1:nspp, function(sp) {
    mort_tmp <- as.data.frame(matrix(M[sp], nrow = 1, ncol = nages))
    colnames(mort_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(Species = sp,
                     Sex = 0),
          mort_tmp
    )
  })
  simData$M1_base <- do.call(rbind, mort_list)

  # * Environmental data ----
  simData$env_data <- data.frame(Year = 1:nyrs, EnvData = 1)

  # * Ration ----
  ration_list <- lapply(1:nspp, function(sp) {
    ration_tmp <- as.data.frame(matrix(WAA[sp,] * 50, ncol = nages))
    colnames(ration_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(
      Species = sp,
      Sex = 0,
      Year = 0),
      ration_tmp
    )
  })
  simData$ration_data <- do.call(rbind, ration_list)

  # * Bioenergetics
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

  # * Diet proportion ----
  # Get all combinations of indices
  idx <- which(!is.na(diet_prop), arr.ind = TRUE)  # or diet_prop != 0 for nonzero only

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
    Stomach_proportion_by_weight = diet_prop[idx]
  )
  simData$diet_data <- simData$diet_data |>
    dplyr::filter(!is.na(Pred))

  # Return list of true and observed values ----
  return(list(
    model_quantities = list(
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
      Mtv = Mtv,
      M_vec = M_vec,
      srv_q = srv_q,
      rec_devs = rec_devs,
      init_devs = init_devs,
      M2_at_age = M2_at_age,
      vulnerability = vulnerability,
      suitability = suitability,
      suit_other = suit_other, # Automatically outputs as matrix if normalized
      avail_food = avail_food,
      diet_prop = diet_prop,
      ration = ration
    ),
    data_list = simData
  ))
}


calc_nll_ar1_2d <- function(x_matrix, sigma_innov, rho_a, rho_y) {
  n_age <- nrow(x_matrix)
  n_yr  <- ncol(x_matrix)

  # 1. Calculate the Marginal Standard Deviation
  Sigma_M <- sqrt(sigma_innov^2 / ((1 - rho_a^2) * (1 - rho_y^2)))

  # 2. Build the Age Correlation Matrix (Rows)
  dist_age <- abs(outer(1:n_age, 1:n_age, "-"))
  R_age <- rho_a^dist_age

  # 3. Build the Year Correlation Matrix (Columns)
  dist_yr <- abs(outer(1:n_yr, 1:n_yr, "-"))
  R_yr <- rho_y^dist_yr

  # 4. Combine them using the Kronecker Product
  # TMB's SEPARABLE(AR1(age), AR1(yr)) assumes rows=age, cols=yr.
  Correlation_2D <- kronecker(R_yr, R_age)

  # 5. Scale to Covariance
  Covariance_2D <- (Sigma_M^2) * Correlation_2D

  # 6. Flatten the 2D matrix of deviations into a 1D vector
  x_vec <- as.vector(x_matrix)

  # 7. Calculate Negative Log-Likelihood
  nll <- -mvtnorm::dmvnorm(x_vec, mean = rep(0, length(x_vec)), sigma = Covariance_2D, log = TRUE)

  return(nll)
}

