#' Build parameter list from cpp file
#'
#' @description Function to read a TMB cpp file and construct parameter list object for Rceattle
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#' @param inits Character vector of named initial values from ADMB \code{.std} or \code{.par} files or list of previous parameter estimates from Rceattle model.
#'
#' @return a list of map arguments for each parameter
#' @export
build_params <- function(data_list, inits = NULL) {
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
  param_list$ln_pop_scalar = matrix(0, nrow = data_list$nspp, ncol = max(data_list$nages))
  param_list$ln_sex_ratio_sigma = log(data_list$sex_ratio_sigma)

  # -- 3.1. Recruitment parameters
  param_list$ln_mean_rec = rep(0, data_list$nspp)  # Mean recruitment; n = [1, nspp]
  param_list$ln_rec_sigma = log(as.numeric(data_list$sigma_rec_prior))  # Standard deviation of recruitment deviations; n = [1, nspp]
  param_list$rec_dev = matrix(0, nrow = data_list$nspp, ncol = nyrs_proj)  # Annual recruitment deviation; n = [nspp, nyrs_hind]


  # -- 3.2. Abundance parameters
  param_list$init_dev = matrix(0, nrow = data_list$nspp, ncol = max(data_list$nages)-1)  # Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this

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


  # -- 3.4. fishing mortality parameters
  param_list$ln_mean_F = rep(0, nrow(data_list$fleet_control))  # Log mean fishing mortality; n = [1, nspp]
  param_list$proj_F_prop = data_list$fleet_control$proj_F_prop  # Proportion of future fishing mortality for projections for each fleet; n = [1, nfleet]
  param_list$ln_FSPR = matrix(0, nrow = data_list$nspp, ncol = 2)                          # Future fishing mortality for projections for each species; n = [1, nfleet]
  param_list$F_dev = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind)  # Annual fishing mortality deviations; n = [nspp, nyrs_hind] # NOTE: The size of this will likely change

  # Make ln_mean_F very low if the fleet is turned off or not a fishery
  for (i in 1:nrow(data_list$fleet_control)) {
    # Turn of F and F dev if not estimating
    if (data_list$fleet_control$Fleet_type[i] %in% c(0,2)) {
      param_list$ln_mean_F[i] <- -999
    }
  }


  # Set Fdev for years with 0 catch to very low number
  fsh_biom <- data_list$fsh_biom
  fsh_ind <- fsh_biom$Fleet_code[which(fsh_biom$Catch == 0)]
  yr_ind <- fsh_biom$Year[which(fsh_biom$Catch == 0)] - data_list$styr + 1
  param_list$F_dev[fsh_ind, yr_ind] <- -999


  # -- 3.5. Survey catchability parameters
  # Random effects version
  param_list$ln_srv_q = log(data_list$fleet_control$Q_prior)   # Survey catchability; n = [sum(n_srv)]
  param_list$srv_q_pow = rep(0, nrow(data_list$fleet_control))
  param_list$ln_srv_q_dev = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind)   # Survey catchability deviations; n = [sum(n_srv)]
  #param_list$ln_srv_q_dev_re = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind)   # Survey catchability deviations; n = [sum(n_srv)]
  param_list$ln_sigma_srv_q <- log(data_list$fleet_control$Q_sd_prior)
  param_list$ln_sigma_time_varying_srv_q <- log(data_list$fleet_control$Time_varying_q_sd_prior)
  # Log standard deviation for survey selectivity random walk - used for logistic


  # -- 3.6. Selectivity parameters
  n_selectivities <- nrow(data_list$fleet_control)
  param_list$sel_coff = suppressWarnings( array(0, dim = c(n_selectivities, 2, max(1, as.numeric(c(data_list$fleet_control$Nselages) ), na.rm = T))))  # Non-parametric selectivity coef; n = [nspp, nselages]
  param_list$sel_curve_pen = suppressWarnings( matrix( c(data_list$fleet_control$Time_varying_sel, data_list$fleet_control$Sel_sd_prior), nrow = n_selectivities, ncol = 2)) # Non-parametric selectivity penalties
  param_list$ln_sel_slp = array(0.5, dim = c(2, n_selectivities, 2))  # selectivity paramaters for logistic; n = [2, nspp, 2 sexes]
  param_list$sel_inf = array(0, dim = c(2, n_selectivities, 2))  # selectivity paramaters for logistic; n = [2, nspp, 2 sexes]
  param_list$sel_inf[1,,] <- 5
  param_list$sel_inf[2,,] <- 10

  # --- 3.5.2. Time varying selectivity parameters
  param_list$ln_sel_slp_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind))  # selectivity deviations paramaters for logistic; n = [2, nspp]
  param_list$sel_inf_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind))  # selectivity deviations paramaters for logistic; n = [2, nspp]

  #param_list$ln_sel_slp_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind))  # selectivity random effect deviations paramaters for logistic; n = [2, nspp]
  #param_list$sel_inf_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind))  # selectivity random effectdeviations paramaters for logistic; n = [2, nspp]

  param_list$ln_sigma_sel <- log(data_list$fleet_control$Sel_sd_prior)          # Log standard deviation for  selectivity random walk - used for logistic


  # -- 3.7. Variance of survey and fishery time series
  param_list$ln_sigma_srv_index = log(data_list$fleet_control$Survey_sd_prior)      # Log standard deviation of survey index time-series; n = [1, n_srv]
  param_list$ln_sigma_fsh_catch = log(data_list$fleet_control$Catch_sd_prior)       # Log standard deviation of fishery catch time-series; n = [1, n_fsh]


  # -- 3.8. Kinzery predation function parameters
  # param_list$logH_1 = matrix(-8.5, nrow = data_list$nspp, ncol = data_list$nspp2)  # Predation functional form; n = [nspp, nspp2]; # FIXME: make matrix; nspp2 = nspp + 1
  # param_list$logH_1a = rep(-3, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  # param_list$logH_1b = rep(0, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  #
  # param_list$logH_2 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]
  # param_list$logH_3 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  # param_list$H_4 = matrix(1, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;


  # -- 3.9. Gamma selectivity parameters
  param_list$log_gam_a = rep(0, data_list$nspp)  # Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  param_list$log_gam_b = rep(0, data_list$nspp)  # Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10


  # -- 3.10. Preference parameters
  param_list$log_phi = matrix(0, data_list$nspp, data_list$nspp)


  # -- 3.11. Comp weighting
  if(!is.null(data_list$fleet_control$Comp_weights)){
  param_list$comp_weights = data_list$fleet_control$Comp_weights
  }
  if(is.null(data_list$fleet_control$Comp_weights)){
    param_list$comp_weights = rep(1, nrow(data_list$fleet_control))
  }


  #---------------------------------------------------------------------
  # Step 3 -- Replace inits with starting values in range
  #---------------------------------------------------------------------
  param_list$ln_mean_rec <- replace(param_list$ln_mean_rec, values = 9)
  param_list$ln_mean_F <- replace(param_list$ln_mean_F, values = -0.8)
  param_list$log_gam_a <- replace(param_list$log_gam_a, values = 0.5)
  param_list$log_gam_b <- replace(param_list$log_gam_b, values = -0.5)


  #---------------------------------------------------------------------
  # Step 4 -- Replace inits with previous parameters if desired
  #---------------------------------------------------------------------
  if(class(inits) == "list"){
    for(i in 1:length(inits)){
      param_list[[names(inits)[i]]] = inits[[names(inits)[i]]]
    }
  }


  # Start phi at 0.5
  param_list$log_phi <- replace(param_list$log_phi, values = rep(log(0.5), length(param_list$log_phi)))

  # closeAllConnections()



  return(param_list)
}
