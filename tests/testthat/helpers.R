# Helpers for tests
# Minimal test data factory and small utilities used by testthat files

make_test_data <- function(nyrs = 8, nages = 5, nspp = 1, seed = NULL) {
  if (!requireNamespace("Rceattle", quietly = TRUE)) {
    stop("Rceattle package required for test helpers")
  }

  data("GOAcod", package = "Rceattle")
  dat <- GOAcod

  if (!is.null(seed)) set.seed(seed)

  dat$nspp <- nspp
  dat$nsex <- 1
  dat$styr <- 1
  dat$endyr <- nyrs
  dat$projyr <- nyrs + 2
  dat$nages <- nages
  dat$minage <- 1
  dat$nlengths <- nages
  dat$pop_wt_index <- 1
  dat$ssb_wt_index <- 1
  dat$pop_age_transition_index <- 1

  # Keep two simple fleets (survey + fishery)
  if (is.data.frame(dat$fleet_control) && nrow(dat$fleet_control) >= 3) {
    dat$fleet_control <- dat$fleet_control[c(1, 3), ]
  }
  dat$fleet_control$Fleet_name <- c("Survey", "Fishery")
  dat$fleet_control$Fleet_code <- 1:2
  dat$fleet_control$Selectivity_index <- 1:2
  dat$fleet_control$Weight_index <- 1

  years <- seq(dat$styr, dat$endyr)

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

  # Minimal composition / auxiliary tables
  dat$comp_data <- data.frame()
  if (!is.null(dat$emp_sel)) dat$emp_sel[] <- NA
  if (!is.null(dat$NByageFixed)) dat$NByageFixed[] <- NA

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
