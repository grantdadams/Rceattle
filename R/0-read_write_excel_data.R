

#' Write data file
#'
#' @param data_list Rceattle data_list object
#' @param file Filname to be used. Must end with '.xlsx'
#'
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' write_data(data_list = BS2017SS, file = 'BS2017SS.xlsx')
write_data <- function(data_list, file = "Rceattle_data.xlsx") {

  # Setup a workbook
  data_names <- names(data_list)
  names_used <- c()
  xcel_list <- list()

  # Metadata
  meta_filename <- system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")
  meta_data <- suppressMessages((suppressWarnings(as.data.frame(readxl::read_xlsx(meta_filename)))))
  xcel_list$meta_data <- meta_data


  # Define the row names in order
  row_labels <- c(
    "nspp", "styr", "endyr", "projyr", "nsex", "spawn_month", "nages", "minage",
    "nlengths", "pop_wt_index", "ssb_wt_index", "alpha_wt_len", "beta_wt_len",
    "pop_age_transition_index", "sigma_rec_prior", "other_food", "estDynamics"
  )

  # Combine the data elements into a matrix
  control <- rbind(
    nspp                     = c(data_list$nspp, rep(NA, data_list$nspp - 1)),
    styr                     = c(data_list$styr, rep(NA, data_list$nspp - 1)),
    endyr                    = c(data_list$endyr, rep(NA, data_list$nspp - 1)),
    projyr                   = c(data_list$projyr, rep(NA, data_list$nspp - 1)),
    nsex                     = data_list$nsex,
    spawn_month              = data_list$spawn_month,
    nages                    = data_list$nages,
    minage                   = data_list$minage,
    nlengths                 = data_list$nlengths,
    pop_wt_index             = data_list$pop_wt_index,
    ssb_wt_index             = data_list$ssb_wt_index,
    alpha_wt_len             = data_list$alpha_wt_len,
    beta_wt_len              = data_list$beta_wt_len,
    pop_age_transition_index = data_list$pop_age_transition_index,
    sigma_rec_prior          = data_list$sigma_rec_prior,
    other_food               = data_list$other_food,
    estDynamics              = data_list$estDynamics
  )

  # Convert to data frame and move row names to a column
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
  bioenergetics_control <- matrix(NA, ncol = data_list$nspp, nrow = 13)
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
  bioenergetics_control[13, ] <- data_list$Diet_comp_weights

  bioenergetics_control <- as.data.frame(bioenergetics_control)

  bioenergetics_control <- cbind(c("Ceq", "Cindex","Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4", "Diet_comp_weights"), bioenergetics_control)
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


  # data_names[!data_names %in% names_used]

  # Write the data
  writexl::write_xlsx(xcel_list, file)
}







#' Read a CEATTLE excel data file
#'
#' @param file Filname to be used. Must end with '.xlsx'
#'
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' write_data(data_list = BS2017SS, file = 'BS2017SS.xlsx')
#' data_list <- read_data(file = 'BS2017SS.xlsx')
read_data <- function(file = "Rceattle_data.xlsx") {

  # Setup a list object
  data_list <- list()

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
  vec_vars <- c("nsex", "spawn_month", "nages", "minage", "nlengths",
                "pop_wt_index", "ssb_wt_index", "alpha_wt_len", "beta_wt_len",
                "pop_age_transition_index", "sigma_rec_prior", "other_food", "estDynamics")
  data_list[vec_vars] <- lapply(control[1:data_list$nspp, vec_vars], as.numeric)


  # Composition, fleet control, fixed selectivity, n-at-age
  matrix_data <- c("fleet_control" , "comp_data", "emp_sel", "NByageFixed")
  for (i in 1:length(matrix_data)) {
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = matrix_data[i]))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list[[matrix_data[i]]] <- sheet
  }

  # Transpose data list if necessary
  data_list <- Rceattle::transpose_fleet_control(data_list)

  # -- Update names if necessary
  if(length(data_list$fleet_control$Survey_sd_prior) > 0){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::rename(Estimate_index_sd = Estimate_survey_sd,
                    Index_sd_prior = Survey_sd_prior)
    print("Renaming 'Estimate_survey_sd' to 'Estimate_index_sd' and 'Survey_sd_prior' to 'Index_sd_prior'")
  }


  if(length(data_list$fleet_control$Nselages) > 0){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::rename(N_sel_bins = Nselages,
                    Sel_norm_bin1 = Age_max_selected,
                    Sel_norm_bin2 = Age_max_selected_upper)
    print("Renaming 'Nselages' to 'N_sel_bins', 'Age_max_selected' to 'Sel_norm_bin1', and 'Age_max_selected_upper' to 'Sel_norm_bin2'")
  }

  if(!exists(data_list$fleet_control$Month)){
    data_list$fleet_control <- data_list$fleet_control %>%
      dplyr::mutate(Month = 0)
    print("Adding 'Month' column to 'fleet_control' with default value of 0")
  }

  # Index and catch data ----
  sheetnames <- readxl::excel_sheets(file)

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
    print("Renaming 'fsh_biom' to 'catch_data'")
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
    print("Renaming 'srv_biom' to 'index_data'")
  }

  # -- CAAL
  if("caal_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "caal_data"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$caal_data <- sheet
  } else{
    data_list$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 10))
    colnames(data_list$caal_data) <- c("Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", "CAAL_1", "CAAL_2", "CAAL_3", "CAAL_4")
  }


  # Age transition matrix (age to length) ----
  age_trans_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "age_trans_matrix"))
  data_list$age_trans_matrix <- age_trans_matrix


  # Ageing error ----
  age_error <- as.data.frame(readxl::read_xlsx(file, sheet = "age_error"))
  age_error <- age_error[rowSums(is.na(age_error)) != ncol(age_error), ] # Remove rows with all NA's
  data_list$age_error <- age_error


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
    print("Renaming 'wt' to 'weight'")
  }


  # Maturity-at-age ----
  if("pmature" %in% sheetnames){ # Old name
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "pmature"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$maturity <- sheet
    print("Renaming 'pmature' to 'maturity'")
  }

  if("maturity" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "maturity"))
    sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]
    data_list$maturity <- sheet
  }


  # Sex ratio-at-age (proportion F) ----
  data_list$sex_ratio <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "sex_ratio")))


  # Residual natural mortality ----
  M1_base <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "M1_base")))
  data_list$M1_base <- M1_base


  # # Length-weight parameters ----
  # aLW <- as.data.frame(readxl::read_xlsx(file, sheet = "aLW"))
  # data_list$aLW <- aLW


  # Bioenergetics specifications  ----
  bioenergetics_control <- as.data.frame(readxl::read_xlsx(file, sheet = "bioenergetics_control"))

  for (i in 1:nrow(bioenergetics_control)) {
    data_list[[bioenergetics_control$Object[i]]] <- suppressWarnings(as.numeric(as.character(bioenergetics_control[i, ((1:data_list$nspp) + 1)])))
  }


  # Environmental data ----
  env_data <- as.data.frame(readxl::read_xlsx(file, sheet = "env_data"))
  data_list$env_data <- env_data


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
    print("Renaming 'Pyrs' to 'ration_data'")
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
    print("Renaming 'UobsWtAge' to 'diet_data'")
  }

  if("stom_prop_data" %in% sheetnames){
    sheet <- as.data.frame(readxl::read_xlsx(file, sheet = "stom_prop_data"))
    data_list$diet_data <- sheet
    print("Renaming 'stom_prop_data' to 'diet_data'")
  }


  # write the data
  return(data_list)
}
