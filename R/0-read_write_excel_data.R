

#' Write data file
#'
#' @param data_list Rceattle data_list object
#' @param file Filname to be used. Must end with '.xlsx'
#'
#' @seealso [build_data()] to assemble or edit a data list in R, [read_data()]
#'   to read one back. Note a [model_config()] slot is not written to the
#'   workbook and does not survive the round-trip (a warning is issued); persist
#'   it separately with [save_config()] / [load_config()].
#'
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' out_file <- file.path(tempdir(), "BS2017SS.xlsx")
#' write_data(data_list = BS2017SS, file = out_file)
#' file.remove(out_file)
write_data <- function(data_list, file = "Rceattle_data.xlsx") {

  # Upgrade any deprecated column / element names to canonical first, so the
  # written workbook uses the canonical names (the control + bioenergetics
  # sheets below are assembled from the canonical schema names).
  data_list$fleet_control <-
    .rce_upgrade_fleet_control_aliases(data_list$fleet_control)
  data_list <- .rce_upgrade_data_list_aliases(data_list)

  # A model_config slot is code-side model structure, not a workbook data sheet,
  # so it is not written and will not survive the xlsx round-trip. Warn rather
  # than drop it silently; re-attach it in code after read_data(), or persist it
  # separately with save_config() / load_config().
  if (!is.null(data_list$model_config)) {
    warning("data_list$model_config is not written to the xlsx workbook and ",
            "will be lost on read_data(); re-attach it in code with ",
            "model_config(), or persist it with save_config()/load_config().", call. = FALSE)
  }

  # Setup a workbook
  data_names <- names(data_list)
  names_used <- c()
  xcel_list <- list()

  # Metadata
  meta_filename <- system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")
  meta_data <- suppressMessages((suppressWarnings(as.data.frame(readxl::read_xlsx(meta_filename)))))
  xcel_list$meta_data <- meta_data


  # Control-sheet objects (names + order) come from the canonical schema
  # (R/0-column_schema.R) instead of a hand-maintained vector. read_data parses
  # the control sheet by object NAME, so the row order is not load-bearing. The
  # four model-wide dimensions carry a single meaningful value in column 1 and
  # NA in the rest; every other object is a per-species vector (length nspp) or a
  # single value recycled across species (rbind's default) -- matching the prior
  # hand-written assembly exactly.
  row_labels <- .rce_schema_names("control")
  .model_dims <- c("nspp", "styr", "endyr", "projyr")
  control <- do.call(rbind, lapply(row_labels, function(nm) {
    v <- data_list[[nm]]
    if (nm %in% .model_dims) c(v[1], rep(NA, data_list$nspp - 1)) else v
  }))

  # Move object names into a column
  control <- data.frame(Object = row_labels, control, row.names = NULL)
  colnames(control) <- c("Object", data_list$spnames)

  names_used <- c(names_used, as.character(control$Object))
  xcel_list$control <- control


  # Fleet control ----
  xcel_list$fleet_control <- as.data.frame(data_list$fleet_control)
  names_used <- c(names_used, "fleet_control")


  # Composition, fleet control, fixed selectivity, n-at-age ---
  matrix_data <- c("index_data", "catch_data", "comp_data",  "caal_data", "emp_sel", "NByageFixed", "age_trans_matrix")
  for (i in 1:length(matrix_data)) {
    xcel_list[[matrix_data[i]]] <- data_list[[matrix_data[i]]]
  }
  names_used <- c(names_used, matrix_data)


  # Ageing error ----
  xcel_list$age_error <- as.data.frame(data_list$age_error)
  names_used <- c(names_used, "age_error")


  # Weight-at-age ----
  xcel_list$weight <- data_list$weight
  names_used <- c(names_used, "weight")


  # Maturity-at-age ----
  colnames(data_list$maturity) <- c("Species", paste0("Age", 1:max(data_list$nages)))
  xcel_list$maturity <- as.data.frame(data_list$maturity)
  names_used <- c(names_used, "maturity")


  # Sex ratio-at-age (proportion F) ----
  colnames(data_list$sex_ratio) <- c("Species", paste0("Age", 1:max(data_list$nages)))
  xcel_list$sex_ratio <- as.data.frame(data_list$sex_ratio)
  names_used <- c(names_used, "sex_ratio")


  # Residual natural mortality ----
  xcel_list$M1_base <- as.data.frame(data_list$M1_base)
  names_used <- c(names_used, "M1_base")


  # # Length-weight parameters ----
  # xcel_list$aLW <- as.data.frame(data_list$aLW)
  # names_used <- c(names_used, "aLW")


  # Bioenergetics specifications  ----
  if(is.null(data_list$Diet_comp_weights)){
    data_list$Diet_comp_weights <- rep(1, data_list$nspp)
  }
  if(is.null(data_list$Diet_distribution)){
    data_list$Diet_distribution <- rep(0, data_list$nspp)
  }
  bioenergetics_control <- matrix(NA, ncol = data_list$nspp, nrow = 14)
  bioenergetics_control[1, ] <- data_list$Ceq
  bioenergetics_control[2, ] <- data_list$Cindex
  bioenergetics_control[3, ] <- data_list$Pvalue
  bioenergetics_control[4, ] <- data_list$fday
  bioenergetics_control[5, ] <- data_list$CA
  bioenergetics_control[6, ] <- data_list$CB
  bioenergetics_control[7, ] <- data_list$Qc
  bioenergetics_control[8, ] <- data_list$Tco
  bioenergetics_control[9, ] <- data_list$Tcm
  bioenergetics_control[10, ] <- data_list$Tcl
  bioenergetics_control[11, ] <- data_list$CK1
  bioenergetics_control[12, ] <- data_list$CK4
  bioenergetics_control[13, ] <- data_list$Diet_distribution
  bioenergetics_control[14, ] <- data_list$Diet_comp_weights

  bioenergetics_control <- as.data.frame(bioenergetics_control)

  # Object names + order for the bioenergetics sheet come from the schema; the
  # matrix rows above are assembled by hard-coded index in this exact order, so
  # assert the schema order matches (guards against a silent mislabel if the
  # schema's bioenergetics rows are ever reordered).
  bio_labels <- .rce_schema_names("bioenergetics_control")
  stopifnot(identical(bio_labels,
    c("Ceq", "Cindex", "Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm",
      "Tcl", "CK1", "CK4", "Diet_distribution", "Diet_comp_weights")))
  bioenergetics_control <- cbind(bio_labels, bioenergetics_control)
  colnames(bioenergetics_control) <- c("Object", data_list$spnames)


  xcel_list$bioenergetics_control <- as.data.frame(bioenergetics_control)
  names_used <- c(names_used, as.character(bioenergetics_control$Object))


  # Environmental data ----
  xcel_list$env_data <- data_list$env_data
  names_used <- c(names_used, c("env_data"))


  # Diet information Pyrs (relative foraging rate) ----
  xcel_list$ration_data <- as.data.frame(data_list$ration_data)
  names_used <- c(names_used, "ration_data")


  # Diet proportion ----
  # - Proportion of prey-at-age in the stomach of a predator-at-age
  xcel_list$diet_data <- as.data.frame(data_list$diet_data)
  names_used <- c(names_used, "diet_data")


  # Survey-index covariance ----
  # - One sheet per MVN/MVNORM fleet, named "index_cov_<Fleet_name>", holding that
  #   fleet's variance-covariance matrix. Sheet names are capped at 31 characters
  #   by the xlsx format, so long fleet names are truncated.
  if(!is.null(data_list$index_cov) && length(data_list$index_cov) > 0){
    for(nm in names(data_list$index_cov)){
      sheet <- substr(paste0("index_cov_", nm), 1, 31)
      xcel_list[[sheet]] <- as.data.frame(as.matrix(data_list$index_cov[[nm]]))
      names_used <- c(names_used, sheet)
    }
  }

  # data_names[!data_names %in% names_used]

  # Write the data
  writexl::write_xlsx(xcel_list, file)
}


#' Write a minimal starter CEATTLE data workbook
#'
#' Emits a small, structurally complete single-species workbook -- one survey
#' and one fishery, flat placeholder data -- that a user can open and edit as a
#' starting point. fleet_control is built on the canonical column names with
#' schema defaults filled by \code{\link{switch_check}}, so the template is
#' always in sync with the current schema. The template round-trips through
#' \code{\link{read_data}} and \code{\link{data_check}} and builds under
#' \code{fit_mod(estimateMode = 3)}; replace the placeholder observations with
#' real data before fitting.
#'
#' @param file Output path for the \code{.xlsx} workbook.
#' @param nages Number of ages to template (default 10).
#' @param nyrs Number of hindcast years (default 30).
#' @param nprojyrs Number of projection years beyond the hindcast (default 12).
#' @param minage Minimum age (recruitment age; default 1).
#'
#' @return Invisibly, the minimal \code{data_list} that was written.
#' @export
write_template <- function(file = "Rceattle_data_template.xlsx",
                           nages = 10, nyrs = 30, nprojyrs = 12, minage = 1) {
  ages  <- seq_len(nages)
  years <- seq_len(nyrs)

  d <- list(
    nspp = 1, styr = 1, endyr = nyrs, projyr = nyrs + nprojyrs,
    spnames = "Species_1", nsex = 1, spawn_month = 0, nages = nages,
    minage = minage, nlengths = nages, pop_wt_index = 1, ssb_wt_index = 1,
    alpha_wt_len = 1e-4, beta_wt_len = 3, pop_age_transition_index = 1,
    sigma_rec = 1, other_food = 1e6, estDynamics = 0)

  # fleet_control: one survey + one fishery on canonical column names. The full
  # column set is written explicitly (some columns -- e.g. Time_varying_q -- are
  # read by convert_switches()/rearrange_data() and are not schema-defaulted).
  d$fleet_control <- data.frame(
    Fleet_name = c("Survey", "Fishery"), Fleet_code = 1:2,
    Fleet_type = c("Survey", "Fishery"), Species = 1L, Month = 0L,
    Selectivity_index = 1:2, Selectivity = "Logistic",
    Selectivity_dimension = "Age", N_sel_bins = NA,
    Sel_curve_pen1 = NA, Sel_curve_pen2 = NA, Time_varying_sel = 0L,
    Time_varying_sel_sd = 1, Bin_first_selected = 1L,
    Sel_norm_bin = NA, Sel_norm_bin_upper = NA,
    Comp_distribution = "Multinomial", Comp_weights = 1,
    CAAL_distribution = 0L, CAAL_weights = 1, Observation_units = 1L,
    Weight_index = 1L, Age_transition_index = 1L,
    Catchability_index = c(1L, NA), Catchability = c("Fixed", NA),
    Catchability_init = c(1, NA), Catchability_prior_sd = c(0.2, NA),
    Time_varying_q = c(0L, NA), Time_varying_q_sd = c(1, NA),
    Estimate_index_sd = c(0L, NA), Index_sd = c(1, NA),
    Estimate_catch_sd = c(NA, 0L), Catch_sd = c(NA, 1),
    Proj_F_proportion = c(NA, 1), stringsAsFactors = FALSE)

  # Placeholder observations (flat) -- replace with real data before fitting.
  d$index_data <- data.frame(Fleet_name = "Survey", Fleet_code = 1L, Species = 1L,
    Year = years, Month = 0L, Selectivity_block = 1L, Observation = 100,
    Log_sd = 0.2)
  d$catch_data <- data.frame(Fleet_name = "Fishery", Fleet_code = 2L, Species = 1L,
    Year = years, Month = 0L, Selectivity_block = 1L, Catch = 10, Log_sd = 0.05)

  # Minimal flat age composition (one row per fleet) so a fleet with estimated
  # selectivity has comp data to fit.
  comp_cols <- c("Fleet_name", "Fleet_code", "Species", "Sex", "Age0_Length1",
                 "Year", "Month", "Sample_size", paste0("Comp_", ages))
  flat <- as.data.frame(matrix(1 / nages, nrow = 2, ncol = nages))
  d$comp_data <- stats::setNames(
    cbind(data.frame(Fleet_name = c("Survey", "Fishery"), Fleet_code = 1:2,
                     Species = 1L, Sex = 0L, Age0_Length1 = 0L, Year = years[1],
                     Month = 0L, Sample_size = 1L, stringsAsFactors = FALSE), flat),
    comp_cols)

  # Empty optional matrices (typed, so read/write round-trip cleanly).
  d$caal_data   <- stats::setNames(data.frame(matrix(NA_real_, 0, 7 + nages)),
    c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length",
      "Sample_size", paste0("CAAL_", ages)))
  d$emp_sel     <- stats::setNames(data.frame(matrix(NA_real_, 0, 5 + nages)),
    c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", paste0("Comp_", ages)))
  d$NByageFixed <- stats::setNames(data.frame(matrix(NA_real_, 0, 4 + nages)),
    c("Species_name", "Species", "Sex", "Year", paste0("Age", ages)))

  # Weight / maturity / sex-ratio / mortality (flat placeholders).
  waa <- stats::setNames(as.data.frame(matrix(1, 1, nages)), paste0("Age", ages))
  d$weight   <- cbind(data.frame(Wt_name = "Base", Wt_index = 1L, Species = 1L,
                                 Sex = 0L, Year = 0L), waa)
  d$maturity <- cbind(data.frame(Species = 1L),
                      stats::setNames(as.data.frame(matrix(1, 1, nages)), paste0("Age", ages)))
  d$sex_ratio <- cbind(data.frame(Species = 1L),
                       stats::setNames(as.data.frame(matrix(0.5, 1, nages)), paste0("Age", ages)))
  d$M1_base <- cbind(data.frame(Species = 1L, Sex = 0L),
                     stats::setNames(as.data.frame(matrix(0.2, 1, nages)), paste0("Age", ages)))

  # Age-transition (identity) + ageing-error (identity).
  atm <- stats::setNames(as.data.frame(diag(nages)), paste0("Length_", ages))
  d$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
    Age_transition_index = 1L, Species = 1L, Sex = 0L,
    Age = minage:(minage + nages - 1L)), atm)
  aerr <- stats::setNames(as.data.frame(diag(nages)), paste0("Obs_age", ages))
  d$age_error <- cbind(data.frame(Species = 1L,
    True_age = minage:(minage + nages - 1L)), aerr)

  # Bioenergetics scalars (single-species defaults) + a minimal environmental
  # index (its column count sizes the M1 environmental-linkage parameter) +
  # empty ration / diet frames.
  d$Ceq <- 1L; d$Cindex <- 0L; d$Pvalue <- 1; d$fday <- 365
  d$CA <- 0; d$CB <- 0; d$Qc <- 0; d$Tco <- 0; d$Tcm <- 0; d$Tcl <- 0
  d$CK1 <- 0; d$CK4 <- 0; d$Diet_distribution <- 0L; d$Diet_comp_weights <- 1
  d$env_data  <- data.frame(Year = years, BottomTemp = 0)
  d$ration_data <- d$weight[, c("Species", "Sex", "Year", paste0("Age", ages))]
  d$diet_data <- stats::setNames(data.frame(matrix(NA_real_, 0, 9)),
    c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age",
      "Year", "Sample_size", "Stomach_proportion_by_weight"))

  d <- clean_data(d)
  d <- suppressMessages(switch_check(d))
  write_data(d, file)
  invisible(d)
}







#' Read a CEATTLE excel data file
#'
#' @param file Filname to be used. Must end with '.xlsx'
#'
#' @seealso [build_data()] to assemble or edit a data list in R (and to read a
#'   workbook then override blocks in one call, `build_data(file = ...)`),
#'   [data_requirements()] to see which inputs a configuration needs.
#'
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' out_file <- file.path(tempdir(), "BS2017SS.xlsx")
#' write_data(data_list = BS2017SS, file = out_file)
#' data_list <- read_data(file = out_file)
#' # Or read + override blocks in one call:
#' data_list <- build_data(file = out_file, projyr = 2060)
#' file.remove(out_file)
read_data <- function(file = "Rceattle_data.xlsx") {

  # Setup a list object
  data_list <- list()

  # Sheet inventory, hoisted so every read below can check presence. The two
  # required sheets error clearly if absent; every other sheet is optional and
  # simply skipped when missing (defaulted downstream by clean_data() /
  # switch_check()), so a minimal single-species workbook reads cleanly instead
  # of failing with a cryptic readxl error.
  sheetnames <- readxl::excel_sheets(file)
  for (req in c("control", "fleet_control")) {
    if (!req %in% sheetnames)
      stop("Required sheet '", req, "' not found in workbook '", file,
           "'. Sheets present: ", paste(sheetnames, collapse = ", "), ".",
           call. = FALSE)
  }

  # Control (model dimensions) ----
  # 1. Read and Transpose
  sheet1 <- readxl::read_xlsx(file, sheet = "control")
  control <- as.data.frame(t(sheet1[, -1]))
  colnames(control) <- sheet1$Object

  # 2. Automatically convert numeric columns (cleans up as.numeric calls)
  control[] <- lapply(control, function(x) type.convert(as.character(x), as.is = TRUE))

  # 3. Extract Scalars (just take the first value)
  scalar_vars <- c("nspp", "styr", "endyr", "projyr")
  data_list[scalar_vars] <- lapply(control[1, scalar_vars], as.numeric)

  # 4. Extract Species Names
  data_list$spnames <- rownames(control)[1:data_list$nspp]

  # 5. Extract Vectors
  # vec_vars <- c("nsex", "spawn_month", "nages", "minage", "nlengths",
  #               "pop_wt_index", "ssb_wt_index", "alpha_wt_len", "beta_wt_len",
  #               "pop_age_transition_index", "sigma_rec", "other_food", "estDynamics")
  # data_list[vec_vars] <- lapply(control[1:data_list$nspp, vec_vars], as.numeric)

  # Coerce a per-species control/bioenergetics row to numeric, but error on a
  # non-empty cell that is not a number (a typo like "O.5" or a stray label)
  # rather than silently turning it into NA -- which used to propagate a wrong
  # default deep into the fit. Genuinely empty / NA cells stay NA.
  checked_numeric <- function(x, object, sheet) {
    chr <- trimws(as.character(x))
    num <- suppressWarnings(as.numeric(chr))
    # Flag only a genuinely present, non-empty, non-numeric cell. A blank Excel
    # cell reads back as NA (never ""), and the literal strings "NA"/"NaN" are
    # legitimate ways to say "no value" -- all of these stay NA, not an error.
    bad <- !is.na(chr) & nzchar(chr) & is.na(num) &
      !(tolower(chr) %in% c("na", "nan"))
    if (any(bad)) {
      stop("Non-numeric value(s) ",
           paste(sprintf("'%s'", chr[bad]), collapse = ", "),
           " for '", object, "' in the '", sheet, "' sheet; expected numbers.",
           call. = FALSE)
    }
    num
  }

  for (i in 5:nrow(sheet1)) {
    data_list[[sheet1$Object[i]]] <-
      checked_numeric(sheet1[i, ((1:data_list$nspp) + 1)], sheet1$Object[i], "control")
  }


  # Composition, fleet control, fixed selectivity, n-at-age. fleet_control is
  # required (guarded above); comp_data / emp_sel / NByageFixed are optional --
  # read only if present.
  matrix_data <- c("fleet_control" , "comp_data", "emp_sel", "NByageFixed")
  for (i in 1:length(matrix_data)) {
    if (!matrix_data[i] %in% sheetnames) next
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = matrix_data[i]))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list[[matrix_data[i]]] <- sheet
  }

  # -- Upgrade any deprecated fleet_control column names to canonical, from the
  # single schema-driven migration (the `aliases` field). The same call
  # switch_check() uses, so the xlsx-read and direct-data-list paths upgrade
  # identically from one source (no hand-maintained cascade to drift). Legacy
  # names accepted: Estimate_survey_sd, Survey_sd_prior/Index_sd_prior, Nselages,
  # Sel_sd_prior, Time_varying_{q,sel}_sd_prior, Q_prior, Catch_sd_prior,
  # Estimate_q, Age_first_selected, Age_max_selected(_upper).
  data_list$fleet_control <-
    .rce_upgrade_fleet_control_aliases(data_list$fleet_control)

  # Default a missing `Month` column.
  if(is.null(data_list$fleet_control$Month)){
    data_list$fleet_control <- data_list$fleet_control |>
      dplyr::mutate(Month = 0)
    message("Adding 'Month' column to 'fleet_control' with default value of 0")
  }


  if(length(data_list$fleet_control$Nselages) > 0){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::rename(N_sel_bins = Nselages)
    print("Renaming 'Nselages' to 'N_sel_bins'")
  }


  if(length(data_list$fleet_control$Age_first_selected) > 0){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::rename(Bin_first_selected = Age_first_selected)
    print("Renaming 'Age_first_selected' to 'Bin_first_selected'")
  }

  if(length(data_list$fleet_control$Age_max_selected) > 0){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::rename(Sel_norm_bin1 = Age_max_selected)
    print("Renaming 'Age_max_selected' to 'Sel_norm_bin1'")
  }

  if(length(data_list$fleet_control$Age_max_selected_upper) > 0){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::rename(Sel_norm_bin2 = Age_max_selected_upper)
    print("Renaming 'Age_max_selected_upper' to 'Sel_norm_bin2'")
  }

  if(!is.null(data_list$fleet_control$Month)){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::mutate(Month = 0)
    print("Adding 'Month' column to 'fleet_control' with default value of 0")
  }

  # Index and catch data ----
  # (sheetnames was hoisted to the top of read_data)

  # -- Catch
  if("catch_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "catch_data"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$catch_data <- sheet
  }

  if("fsh_biom" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "fsh_biom"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$catch_data <- sheet
    message("Renaming 'fsh_biom' to 'catch_data'")
  }

  # -- Index
  if("index_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "index_data"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$index_data <- sheet
  }

  if("srv_biom" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "srv_biom"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$index_data <- sheet
    message("Renaming 'srv_biom' to 'index_data'")
  }

  # -- CAAL
  if("caal_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "caal_data"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$caal_data <- sheet
  } else{
    data_list$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 11))
    colnames(data_list$caal_data) <- c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", "CAAL_1", "CAAL_2", "CAAL_3", "CAAL_4")
  }

  # -- CAAL
  if("caal_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "caal_data"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$caal_data <- sheet
  } else{
    data_list$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 11))
    colnames(data_list$caal_data) <- c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", "CAAL_1", "CAAL_2", "CAAL_3", "CAAL_4")
  }


  # Age transition matrix (age to length) ----
  if("age_trans_matrix" %in% sheetnames){
    data_list$age_trans_matrix <-
      as.data.frame(readxl::read_xlsx(file, sheet = "age_trans_matrix"))
  }


  # Ageing error ----
  if("age_error" %in% sheetnames){
    age_error <- as.data.frame(readxl::read_xlsx(file, sheet = "age_error"))
    age_error <- age_error[rowSums(is.na(age_error)) != ncol(age_error), ] # Remove rows with all NA's
    data_list$age_error <- age_error
  }


  # Weight-at-age ----
  if("weight" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "weight"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$weight <- sheet
  }

  if("wt" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "wt"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$weight <- sheet
    message("Renaming 'wt' to 'weight'")
  }


  # Maturity-at-age ----
  if("pmature" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "pmature"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$maturity <- sheet
    message("Renaming 'pmature' to 'maturity'")
  }

  if("maturity" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "maturity"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$maturity <- sheet
  }


  # Sex ratio-at-age (proportion F) ----
  if("sex_ratio" %in% sheetnames){
    data_list$sex_ratio <-
      suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "sex_ratio")))
  }


  # Residual natural mortality ----
  if("M1_base" %in% sheetnames){
    data_list$M1_base <-
      suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "M1_base")))
  }


  # # Length-weight parameters ----
  # aLW <- as.data.frame(readxl::read_xlsx(file, sheet = "aLW"))
  # data_list$aLW <- aLW


  # Bioenergetics specifications (optional -- multispecies only) ----
  if("bioenergetics_control" %in% sheetnames){
    bioenergetics_control <- as.data.frame(readxl::read_xlsx(file, sheet = "bioenergetics_control"))

    for (i in 1:nrow(bioenergetics_control)) {
      data_list[[bioenergetics_control$Object[i]]] <-
        checked_numeric(bioenergetics_control[i, ((1:data_list$nspp) + 1)],
                        bioenergetics_control$Object[i], "bioenergetics_control")
    }
  }


  # Environmental data ----
  if("env_data" %in% sheetnames){
    data_list$env_data <- as.data.frame(readxl::read_xlsx(file, sheet = "env_data"))
  }


  # Consumption information  ----
  # - Pyrs (relative foraging rate)
  # - or input consumption-at-age (kg/individual)
  if("ration_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "ration_data"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$ration_data <- sheet
  }

  if("Pyrs" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "Pyrs"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ] # Keep rows with data
    data_list$ration_data <- sheet
    message("Renaming 'Pyrs' to 'ration_data'")
  }


  # Diet proportion ----
  # - Proportion of prey-at-age in the stomach of a predator-at-age
  if("diet_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "diet_data"))
    data_list$diet_data <- sheet
  }

  if("UobsWtAge" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "UobsWtAge"))
    data_list$diet_data <- sheet
    message("Renaming 'UobsWtAge' to 'diet_data'")
  }

  if("stom_prop_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "stom_prop_data"))
    data_list$diet_data <- sheet
    message("Renaming 'stom_prop_data' to 'diet_data'")
  }

  # Survey-index covariance ----
  # - One "index_cov_<Fleet_name>" sheet per MVN/MVNORM fleet (see write_data).
  cov_sheets <- grep("^index_cov_", sheetnames, value = TRUE)
  if(length(cov_sheets) > 0){
    data_list$index_cov <- stats::setNames(
      lapply(cov_sheets, function(sh)
        as.matrix(as.data.frame(readxl::read_xlsx(file, sheet = sh)))),
      sub("^index_cov_", "", cov_sheets))
  }

  # Trim trailing all-NA age padding from NByageFixed. Older writers (e.g. the
  # Hake multispecies workbooks) pad the age columns to a fixed width wider than
  # max(nages); those extra columns are entirely empty and would otherwise trip
  # the column-count check in data_check().
  if(!is.null(data_list$NByageFixed) && !is.null(data_list$nages) &&
     ncol(data_list$NByageFixed) > 4 + max(data_list$nages)){
    keep <- 4 + max(data_list$nages)
    extra <- data_list$NByageFixed[, (keep + 1):ncol(data_list$NByageFixed), drop = FALSE]
    if(all(is.na(extra))){
      data_list$NByageFixed <- data_list$NByageFixed[, seq_len(keep), drop = FALSE]
      message("Trimming ", ncol(extra), " empty trailing age column(s) from 'NByageFixed'")
    }
  }

  # Upgrade any deprecated control / bioenergetics element names read from an
  # older workbook (e.g. sigma_rec_prior -> sigma_rec, Diet_loglike ->
  # Diet_distribution) now that every sheet has been read. fleet_control column
  # aliases were already upgraded above.
  data_list <- .rce_upgrade_data_list_aliases(data_list)

  return(data_list)
}
