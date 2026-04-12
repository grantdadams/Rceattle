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
make_test_data <- function(nyrs = 8, nprojyrs = 10, nages = 5, seed = NULL) {
  if (!requireNamespace("Rceattle", quietly = TRUE)) {
    stop("Rceattle package required for test helpers")
  }

  simData <- list()

  if (!is.null(seed)) set.seed(seed)
  nspp = 1
  simData$nspp <- nspp
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs + nprojyrs
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
    Selectivity = "Logistic",
    Selectivity_dimension = "Age",
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
  simData$Diet_loglike = 1
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


calc_multinom_nll <- function(obs_num, hat_prop) {
  p <- hat_prop / sum(hat_prop)
  # TMB uses the continuous lgamma instead of factorial: x! = gamma(x+1)
  ll <- lgamma(sum(obs_num) + 1) - sum(lgamma(obs_num + 1)) + sum(obs_num * log(p))
  return(-ll)
}

calc_dirmultinom_nll <- function(obs_num, alpha) {
  N <- sum(obs_num)
  sum_alpha <- sum(alpha)
  ll <- lgamma(N + 1) - sum(lgamma(obs_num + 1)) +
    lgamma(sum_alpha) - lgamma(N + sum_alpha) +
    sum(lgamma(obs_num + alpha) - lgamma(alpha))
  return(-ll)
}

