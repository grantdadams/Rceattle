#' Build parameter list from cpp file
#'
#' @description Function to read a TMB cpp file and construct parameter list object for Rceattle
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#'
#' @return a list of map arguments for each parameter
#' @export
build_params <- function(data_list) {

  # - Dimensions
  param_list <- list()

  max_age <- max(data_list$nages, na.rm = T)
  max_sex <- 2 # max(data_list$nsex, na.rm = T)
  sex_labels <- c("Sex combined or females", "males")
  if(max_sex == 1){
    sex_labels <- "Sex combined"
  }
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 1. Population dynamics parameters ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  param_list$dummy = 0  # Variable to test derived quantities given input parameters; n = [1]

  # * 1.0. Population scalar ----
  param_list$ln_pop_scalar = matrix(0, nrow = data_list$nspp, ncol = max_age,
                                    dimnames = list(data_list$spnames, paste0("Age", 1:max_age)))

  # * 1.1. Recruitment parameters ----
  # - Stock recruit parameters
  param_list$rec_pars = matrix(9, nrow = data_list$nspp, ncol = 3,
                               dimnames = list(data_list$spnames, c("R0", "Alpha", "Beta")))  # col 1 = mean rec, col 2 = alpha from srr curve, col 3 = beta from srr curve
  param_list$rec_pars[,3] <- log(3) # Starting low here for beta
  param_list$rec_pars[,2] <- 3

  if(!is.null(data_list$srr_prior)){
    param_list$rec_pars[,2] <- log(data_list$srr_prior)
  } else{
    print("Warnings: alpha was not initialized to `srr_prior` from `build_srr`")
  }

  # - Rec devs
  param_list$rec_dev = matrix(0, nrow = data_list$nspp, ncol = nyrs_proj,
                              dimnames = list(data_list$spnames, yrs_proj))  # Annual recruitment deviation; n = [nspp, nyrs_hind]

  param_list$R_ln_sd = log(as.numeric(data_list$sigma_rec_prior))  # Standard deviation of recruitment deviations; n = [1, nspp]
  names(param_list$R_ln_sd) <- data_list$spnames

  # - Env regression parameters for recruitment
  param_list$beta_rec_pars <- matrix(0, nrow = data_list$nspp, ncol = length(data_list$srr_env_indices))


  # * 1.3. Initial age-structure parameters ----
  param_list$init_dev = matrix(0, nrow = data_list$nspp, ncol = max_age,
                               dimnames = list(data_list$spnames, paste0("Age", 1:max_age)))

  if(!is.null(data_list$initMode)){
    for(sp in 1:data_list$nspp){

      # Fill in ages above max age with -999
      if(data_list$nages[sp] != max_age){
        param_list$init_dev[sp,(data_list$nages[sp]+1):max_age] = -999
      }

      # Estimate as devs (fill in ages above max age w/ -999)
      if(data_list$initMode > 0){
        param_list$init_dev[sp,data_list$nages[sp]] = -999
      }
    }
  }


  # * 1.4. Residual natural mortality ----
  m1 <- array(1, dim = c(data_list$nspp, max_sex, max_age),
              dimnames = list(data_list$spnames, sex_labels, paste0("Age", 1:max_age))) # Set up array

  # Initialize from inputs
  for (i in 1:nrow(data_list$M1_base)) {
    sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
    sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))


    # Handle sex == 0 case for 2-sex species
    sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

    # Fill in M1 array from fixed values for each sex
    for(j in 1:length(sex_values)){
      m1[sp, sex_values[j], 1:max_age] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = T)) + 2])
    }
  }
  param_list$ln_M1 <- log(m1)
  # param_list$sex_ratio_ln_sd = log(data_list$sex_ratio_sigma)


  # * 1.5. fishing mortality parameters ----

  # Future fishing mortality limit
  param_list$ln_Flimit = rep(0, data_list$nspp)
  names(param_list$ln_Flimit) <- data_list$spnames

  # - Future fishing mortality target
  param_list$ln_Ftarget = rep(0, data_list$nspp)
  names(param_list$ln_Ftarget) <- data_list$spnames

  # - Initial F when population is not at equilibrium
  param_list$ln_Finit = rep(-10, data_list$nspp)
  names(param_list$ln_Finit) <- data_list$spnames

  # - Proportion of future fishing mortality for projections for each fleet
  param_list$proj_F_prop = data_list$fleet_control$proj_F_prop
  names(param_list$proj_F_prop) <- data_list$fleet_control$Fleet_name

  # - Annual fishing mortality deviations
  param_list$ln_F = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind,
                           dimnames = list(data_list$fleet_control$Fleet_name, yrs_hind))

  # -- Make ln_F very low if the fleet is turned off or not a fishery
  for (i in 1:nrow(data_list$fleet_control)) {
    if (data_list$fleet_control$Fleet_type[i] %in% c(0,2)) {
      param_list$ln_F[i,] <- -999
    }
  }

  # -- Set Fdev for years with 0 catch to very low number
  zero_catch <- data_list$catch_data %>%
    dplyr::filter(Year <= data_list$endyr &
                    Catch == 0) %>%
    dplyr::mutate(Year = Year - data_list$styr + 1) %>%
    select(Fleet_code, Year) %>%
    as.matrix()
  param_list$ln_F[zero_catch] <- -999


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 2. Observation model parameters ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  # * 2.1. Catchability parameters ----
  # - Catchability on log scale
  param_list$index_ln_q = log(data_list$fleet_control$Q_prior)
  names(param_list$index_ln_q) <- data_list$fleet_control$Fleet_name

  # - Regression coefficients for environment-q linkage
  param_list$index_q_beta = matrix(0, nrow = nrow(data_list$fleet_control), ncol = ncol(data_list$env_data) - 1,
                                   dimnames = list(data_list$fleet_control$Fleet_name, colnames(data_list$env_data)[-1]))

  # - Rho for environment-q linkage (sensu GOA Pollock)
  param_list$index_q_rho = rep(0, nrow(data_list$fleet_control))
  names(param_list$index_q_rho) <- data_list$fleet_control$Fleet_name

  # param_list$index_q_pow = rep(0, nrow(data_list$fleet_control))

  # - Annual index catchability deviations
  param_list$index_q_dev = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind,
                                  dimnames = list(data_list$fleet_control$Fleet_name, yrs_hind))

  # - Log standard deviation prior on Q (maybe should be data...)
  param_list$index_q_ln_sd <- log(data_list$fleet_control$Q_sd_prior)
  names(param_list$index_q_ln_sd) <- data_list$fleet_control$Fleet_name

  # - Log standard deviation for survey selectivity random walk - used for logistic
  param_list$index_q_dev_ln_sd <- log(data_list$fleet_control$Time_varying_q_sd_prior)
  names(param_list$index_q_dev_ln_sd) <- data_list$fleet_control$Fleet_name


  # * 2.2. Selectivity parameters ----
  n_selectivities <- nrow(data_list$fleet_control)
  max_sel_ages <- suppressWarnings(max(1, as.numeric(c(data_list$fleet_control$Nselages)), na.rm = T))

  # - Non-parametric selectivity coefficients
  param_list$sel_coff =  array(0, dim = c(n_selectivities, max_sex, max_sel_ages),
                               dimnames = list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_sel_ages)))

  # - Non-parametric selectivity penalties (sensu Ianelli)
  param_list$sel_curve_pen = matrix( c(data_list$fleet_control$Time_varying_sel, data_list$fleet_control$Sel_sd_prior), nrow = n_selectivities, ncol = 2)

  # - Non-parametric selectivity coef annual deviates
  param_list$sel_coff_dev = array(0, dim = c(n_selectivities, max_sex, max(1, as.numeric(c(data_list$fleet_control$Nselages) ), na.rm = T), nyrs_hind),
                                  dimnames = list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_sel_ages), yrs_hind))

  # - Selectivity slope parameters for logistic
  param_list$ln_sel_slp = array(0.5, dim = c(2, n_selectivities, max_sex),
                                dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels))

  # - Selectivity asymptotic parameters for logistic
  param_list$sel_inf = array(0, dim = c(2, n_selectivities, max_sex),
                             dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels))
  param_list$sel_inf[1,,] <- 0
  param_list$sel_inf[2,,] <- 10

  # - Annual selectivity slope deviation for logistic
  param_list$ln_sel_slp_dev = array(0, dim = c(2, n_selectivities, max_sex, nyrs_hind),
                                    dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels, yrs_hind))

  # - Annual selectivity asymptotic deviations for logistic
  param_list$sel_inf_dev = array(0, dim = c(2, n_selectivities, max_sex, nyrs_hind),
                                 dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels, yrs_hind))

  # - Log standard deviation for selectivity random walk - used for logistic
  param_list$sel_dev_ln_sd <- log(data_list$fleet_control$Sel_sd_prior)
  names(param_list$sel_dev_ln_sd) <- data_list$fleet_control$Fleet_name


  # * 2.3. Variance of survey and fishery time series ----
  # - Log standard deviation of survey index time-series
  param_list$index_ln_sd = log(data_list$fleet_control$Index_sd_prior)
  names(param_list$index_ln_sd) <- data_list$fleet_control$Fleet_name

  # - Log standard deviation of fishery catch time-series
  param_list$catch_ln_sd = log(data_list$fleet_control$Catch_sd_prior)
  names(param_list$catch_ln_sd) <- data_list$fleet_control$Fleet_name


  # * 2.4. Comp weighting ----
  if(!is.null(data_list$fleet_control$Comp_weights)){
    param_list$comp_weights = data_list$fleet_control$Comp_weights
  }
  if(is.null(data_list$fleet_control$Comp_weights)){
    param_list$comp_weights = rep(1, nrow(data_list$fleet_control))
  }
  names(param_list$comp_weights) <- data_list$fleet_control$Fleet_name

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 3. Predation model parameters ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  # # * 3.1. Kinzery predation function parameters ----
  # param_list$logH_1 = matrix(-8.5, nrow = data_list$nspp, ncol = data_list$nspp + 1)  # Predation functional form
  # param_list$logH_1a = rep(-3, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  # param_list$logH_1b = rep(0, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  #
  # param_list$logH_2 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]
  # param_list$logH_3 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  # param_list$H_4 = matrix(1, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;
  #
  #
  # * 3.2. Suitability parameters ----
  param_list$log_gam_a = rep(0.5, data_list$nspp)  # Log predator selectivity;
  param_list$log_gam_b = rep(-.5, data_list$nspp)  # Log predator selectivity


  # * 3.3. Preference parameters ----
  param_list$log_phi = matrix(0.5, data_list$nspp, data_list$nspp)


  return(param_list)
}
