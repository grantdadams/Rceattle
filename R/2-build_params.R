#' Build parameter list from cpp file
#'
#' @description Function to read a TMB cpp file and construct parameter list object for Rceattle
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#'
#' @return a list of map arguments for each parameter
#' @export
build_params <- function(data_list) {
  # closeAllConnections()

  data_list$nspp2 = data_list$nspp + 1
  data_list$nspp_sq = data_list$nspp * data_list$nspp
  data_list$nspp_sq2 = data_list$nspp * (data_list$nspp + 1)
  param_list <- list()


  nyrs_hind <- data_list$endyr - data_list$styr + 1
  nyrs_proj <- data_list$projyr - data_list$styr + 1

  #---------------------------------------------------------------------
  # Step 1 -- Specify parameter names and dimensions used in TMB
  #---------------------------------------------------------------------

  param_list$dummy = 0  # Variable to test derived quantities given input parameters; n = [1]

  # -- 3.0. Population scalar/ sex ratio variance
  # param_list$ln_pop_scalar = matrix(0, nrow = data_list$nspp, ncol = max(data_list$nages))

  # -- 3.1. Recruitment parameters
  param_list$rec_pars = matrix(9, nrow = data_list$nspp, ncol = 3)  # col 1 = mean rec, col 2 = alpha from srr curve, col 3 = beta from srr curve
  param_list$rec_pars[,3] <- log(3) # Starting low here for beta
  param_list$rec_pars[,2] <- log(data_list$srr_prior)
  if(data_list$srr_est_mode == 3){
    param_list$rec_pars[,2] <- 3
  }

  param_list$R_ln_sd = log(as.numeric(data_list$sigma_rec_prior))  # Standard deviation of recruitment deviations; n = [1, nspp]
  param_list$rec_dev = matrix(0, nrow = data_list$nspp, ncol = nyrs_proj)  # Annual recruitment deviation; n = [nspp, nyrs_hind]

  # - Env regression parameters for recruitment
  param_list$beta_rec_pars <- matrix(0, nrow = data_list$nspp, ncol = length(data_list$srr_env_indices))


  # -- 3.2. Initial age-structure parameters
  param_list$init_dev = matrix(0, nrow = data_list$nspp, ncol = max(data_list$nages))
  for(sp in 1:data_list$nspp){
    if(data_list$initMode == 0){ # Estimate as free parameters (fill in ages above max age with -999)
      if(data_list$nages[sp] != max(data_list$nages)){
        param_list$init_dev[sp,(data_list$nages[sp]+1):max(data_list$nages)] = -999
      }
    }
    if(data_list$initMode > 0){ # Estimate as devs (fill in ages above max age - 1 with -999)
      param_list$init_dev[sp,data_list$nages[sp]:max(data_list$nages)] = -999
    }
  }

  # -- 3.3. Residual natural mortality
  m1 <- array(0, dim = c(data_list$nspp, 2, max(data_list$nages, na.rm = T))) # Set up array

  # Initialize from inputs
  for (i in 1:nrow(data_list$M1_base)) {
    sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
    sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

    # Fill in M1 array from fixed values for each sex
    if(sex == 0){ sex = c(1, 2)} # If sex = combined/both males and females, fill in both dimensions
    for(j in 1:length(sex)){
      m1[sp, sex[j], 1:max(data_list$nages, na.rm = T)] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = T)) + 2])
    }
  }
  param_list$ln_M1 <- log(m1)
  # param_list$sex_ratio_ln_sd = log(data_list$sex_ratio_sigma)


  # -- 3.4. fishing mortality parameters
  param_list$ln_mean_F = rep(-0.8, nrow(data_list$fleet_control))  # Log mean fishing mortality; n = [1, nspp]
  param_list$ln_Flimit = rep(0, data_list$nspp) # Future fishing mortality for projections for each species
  param_list$ln_Ftarget = rep(0, data_list$nspp) # Future fishing mortality for projections for each species
  param_list$ln_Finit = rep(-10, data_list$nspp)
  param_list$proj_F_prop = data_list$fleet_control$proj_F_prop  # Proportion of future fishing mortality for projections for each fleet
  param_list$F_dev = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind)  # Annual fishing mortality deviations

  # - Make ln_mean_F very low if the fleet is turned off or not a fishery
  for (i in 1:nrow(data_list$fleet_control)) {
    # Turn of F and F dev if not estimating
    if (data_list$fleet_control$Fleet_type[i] %in% c(0,2)) {
      param_list$ln_mean_F[i] <- -999
    }
  }


  # - Set Fdev for years with 0 catch to very low number
  catch_data <- data_list$catch_data
  fsh_ind <- catch_data$Fleet_code[which(catch_data$Catch == 0)]
  yr_ind <- catch_data$Year[which(catch_data$Catch == 0)] - data_list$styr + 1
  for(i in 1:length(fsh_ind)){
    param_list$F_dev[fsh_ind[i], yr_ind[i]] <- -999
  }



  # -- 3.5. Survey catchability parameters
  # Random effects version
  param_list$index_ln_q = log(data_list$fleet_control$Q_prior)   # Survey catchability; n = [sum(n_srv)]
  param_list$index_q_beta = matrix(0, nrow = nrow(data_list$fleet_control), ncol = ncol(data_list$env_data) -1) # Regression coefficients for environment-q linkage
  param_list$index_q_rho = rep(0, nrow(data_list$fleet_control)) # Rho for environment-q linkage
  # param_list$index_q_pow = rep(0, nrow(data_list$fleet_control))
  param_list$index_q_dev = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind)   # Survey catchability deviations; n = [sum(n_srv)]
  param_list$index_q_ln_sd <- log(data_list$fleet_control$Q_sd_prior)
  param_list$index_q_dev_ln_sd <- log(data_list$fleet_control$Time_varying_q_sd_prior)
  # Log standard deviation for survey selectivity random walk - used for logistic


  # -- 3.6. Selectivity parameters
  n_selectivities <- nrow(data_list$fleet_control)
  param_list$sel_coff = suppressWarnings( array(0, dim = c(n_selectivities, 2, max(1, as.numeric(c(data_list$fleet_control$Nselages) ), na.rm = T))))  # Non-parametric selectivity coef; n = [nspp, nselages]
  param_list$sel_coff_dev = suppressWarnings( array(0, dim = c(n_selectivities, 2, max(1, as.numeric(c(data_list$fleet_control$Nselages) ), na.rm = T), nyrs_hind)))  # Non-parametric selectivity coef annual deviates
  param_list$ln_sel_slp = array(0.5, dim = c(2, n_selectivities, 2))  # selectivity paramaters for logistic; n = [2, nspp, 2 sexes]
  param_list$sel_inf = array(0, dim = c(2, n_selectivities, 2))  # selectivity paramaters for logistic; n = [2, nspp, 2 sexes]
  param_list$sel_inf[1,,] <- 0
  param_list$sel_inf[2,,] <- 10

  # --- 3.5.2. Time varying selectivity parameters
  param_list$ln_sel_slp_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind))  # selectivity deviations paramaters for logistic; n = [2, nspp]
  param_list$sel_inf_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind))  # selectivity deviations paramaters for logistic; n = [2, nspp]
  param_list$sel_dev_ln_sd <- log(data_list$fleet_control$Sel_sd_prior)          # Log standard deviation for  selectivity random walk - used for logistic
  param_list$sel_curve_pen = suppressWarnings( matrix( c(data_list$fleet_control$Time_varying_sel, data_list$fleet_control$Sel_sd_prior), nrow = n_selectivities, ncol = 2)) # Non-parametric selectivity penalties


  # -- 3.7. Variance of survey and fishery time series
  param_list$index_ln_sd = log(data_list$fleet_control$Index_sd_prior)      # Log standard deviation of survey index time-series; n = [1, n_srv]
  param_list$catch_ln_sd = log(data_list$fleet_control$Catch_sd_prior)       # Log standard deviation of fishery catch time-series; n = [1, n_fsh]

  # -- 3.8. Comp weighting
  if(!is.null(data_list$fleet_control$Comp_weights)){
    param_list$comp_weights = data_list$fleet_control$Comp_weights
  }
  if(is.null(data_list$fleet_control$Comp_weights)){
    param_list$comp_weights = rep(1, nrow(data_list$fleet_control))
  }


  # # -- 3.9. Kinzery predation function parameters
  # param_list$logH_1 = matrix(-8.5, nrow = data_list$nspp, ncol = data_list$nspp2)  # Predation functional form; n = [nspp, nspp2]; # FIXME: make matrix; nspp2 = nspp + 1
  # param_list$logH_1a = rep(-3, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  # param_list$logH_1b = rep(0, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  #
  # param_list$logH_2 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]
  # param_list$logH_3 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  # param_list$H_4 = matrix(1, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;
  #
  #
  # # -- 3.9. Gamma selectivity parameters
  # param_list$log_gam_a = rep(0.5, data_list$nspp)  # Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  # param_list$log_gam_b = rep(-.5, data_list$nspp)  # Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10
  #
  #
  # # -- 3.10. Preference parameters
  # param_list$log_phi = matrix(0.5, data_list$nspp, data_list$nspp)


  return(param_list)
}
