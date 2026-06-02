#' Function to clean data prior to Rceattle runs
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
clean_data <- function(data_list){

  # --- 0. Default-fill optional data.frame fields ----
  # These fields can be NULL when the user is not using the corresponding feature
  # (e.g., comp_data in a model without composition likelihoods, ration_data /
  # diet_data in single-species mode, NByageFixed when estDynamics == 0).
  # Filling with empty data.frames that carry the metadata columns the
  # downstream code expects lets rearrange_data / build_params use uniform
  # `nrow > 0` checks without separate NULL guards. Whether the missing data is
  # actually a problem is enforced by data_check() based on the model
  # configuration (msmMode, growth_model, estDynamics, Selectivity, etc.).
  default_dfs <- list(
    comp_data   = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Age0_Length1 = integer(),
                             Year = integer(), Month = integer(),
                             Sample_size = numeric()),
    caal_data   = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Year = integer(),
                             Length = numeric(), Sample_size = numeric()),
    emp_sel     = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Year = integer()),
    NByageFixed = data.frame(Species = integer(), Sex = integer(),
                             Year = integer()),
    ration_data = data.frame(Species = integer(), Sex = integer(),
                             Year = integer()),
    diet_data   = data.frame(Pred = integer(), Prey = integer(),
                             Pred_sex = integer(), Prey_sex = integer(),
                             Pred_age = integer(), Prey_age = integer(),
                             Year = integer(),
                             Sample_size = numeric(),
                             Stomach_proportion_by_weight = numeric())
  )
  for(df_name in names(default_dfs)){
    if(is.null(data_list[[df_name]])){
      data_list[[df_name]] <- default_dfs[[df_name]]
    }
  }

  # env_data: default to a Year-only data.frame so downstream code that does
  # `ncol(env_data) - 1` (number of indices) gets 0 when no environmental data
  # are supplied, and `merge(env_data, ...)` in rearrange_data() works.
  # data_check() still errors when a model feature (env-q catchability,
  # temperature-dependent consumption, env linkages, srr/M1 indices) needs
  # actual indices.
  if(is.null(data_list$env_data)){
    yrs <- if(!is.null(data_list$styr) && !is.null(data_list$projyr))
             data_list$styr:data_list$projyr else integer(0)
    data_list$env_data <- data.frame(Year = yrs)
  }

  # --- 1. Filter Data by Year ----
  # Data in likelihood (use absolute Year)
  abs_year_data <- c("index_data", "catch_data", "comp_data", "caal_data")
  for(df_name in abs_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- data_list[[df_name]] |>
        dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
    }
  }

  # Fixed data (allow Year == 0)
  fixed_year_data <- c("diet_data", "weight", "emp_sel", "NByageFixed", "ration_data")
  for(df_name in fixed_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- as.data.frame(data_list[[df_name]]) |>
        dplyr::filter((Year >= data_list$styr & Year <= data_list$projyr) | Year == 0)
    }
  }

  # --- 2. Add temp multi-species SB0 ----
  #FIXME:may be redundant now?
  if(is.null(data_list$MSSB0)){
    data_list$MSSB0 <- rep(999, data_list$nspp)
    data_list$MSB0 <- rep(999, data_list$nspp)
  }


  # --- 3. Extend catch data to proj year for projections ----
  if(data_list$projyr > data_list$endyr){
    for(flt in unique(data_list$catch_data$Fleet_code)){
      catch_data_sub <- data_list$catch_data |> dplyr::filter(Fleet_code == flt)

      yrs_proj <- (data_list$endyr + 1):data_list$projyr
      yrs_proj <- yrs_proj[which(!yrs_proj %in% catch_data_sub$Year)]
      nyrs_proj <- length(yrs_proj)

      if(nyrs_proj > 0) {
        proj_catch_data <- data.frame(
          Fleet_name = rep(catch_data_sub$Fleet_name[1], nyrs_proj),
          Fleet_code = rep(flt, nyrs_proj),
          Species = rep(catch_data_sub$Species[1], nyrs_proj),
          Year = yrs_proj,
          Month = rep(dplyr::last(catch_data_sub$Month), nyrs_proj),
          Selectivity_block = rep(dplyr::last(catch_data_sub$Selectivity_block), nyrs_proj),
          Catch = rep(NA, nyrs_proj),
          Log_sd = rep(dplyr::last(catch_data_sub$Log_sd), nyrs_proj)
        )
        data_list$catch_data <- rbind(data_list$catch_data, proj_catch_data)
      }
    }
  }


  # --- 4. Column names ----
  expected_cols <- c("Species", paste0("Age", 1:max(data_list$nages)))

  if(any(!colnames(data_list$sex_ratio) %in% expected_cols)){
    colnames(data_list$sex_ratio) <- expected_cols
    message("Renaming column names in 'sex_ratio' data to 'Age1', 'Age2', ....")
  }

  if(any(!colnames(data_list$maturity) %in% expected_cols)){
    colnames(data_list$maturity) <- expected_cols
    message("Renaming column names in 'maturity' data to 'Age1', 'Age2', ....")
  }


  # --- 5. Arrange diet data ----
  if(!is.null(data_list$diet_data)){
    data_list$diet_data <- data_list$diet_data |>
      dplyr::arrange(Pred, Pred_sex, Pred_age, Prey, Prey_sex, Prey_age, Year) |>
      dplyr::mutate(stratum_id = paste(Pred, Pred_sex, Pred_age, Year, sep = "_"),
                    stomach_id = as.numeric(as.factor(stratum_id)) - 1) |>
      dplyr::arrange(stomach_id)
  }

  return(data_list)
}
