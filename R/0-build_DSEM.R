
#' Function to fit a dynamic structural equation model related to recruitment
#'
#' @param sem Specification for time-series structural equation model structure including lagged or simultaneous effects. See Details section in \code{dsem::make_dsem_ram} for more description. All variables must be included in and named following variables in \code{env_data}. The default is assumes IID recruitment deviates. NOTE: must include \code{recdevs[spp]} for each species \code{1:nspp} (recdevs1 and recdevs2 for a 2 species model)! If no start value is provided, those model terms are not estimated.
#' @param family Character or character-vector listing the distribution used for each column of \code{env_data} used in the \code{sem}, where each element must be fixed (for no measurement error/measured exactly), normal for normal measurement error using an identity link, gamma for a gamma measurement error using a fixed CV and log-link, bernoulli for a Bernoulli measurement error using a logit-link, or poisson for a Poisson measurement error using a log-link. Default is family family="normal".
#' @param all_vars include all variables from env_data in DSEM model likelihood (estimate observation error) to allow model comparison across different SEM. Default = FALSE.
#' @param estimate_projection latent variables for projection time period are turned off. Default = FALSE.
#'
#' @description
#' The code links dynamic structural equation models to recruitment within Rceattle. The internals of \code{dsem} were copy and pasted into Rceattle. See \code{??dsem} for more description.
#'
#' @export
#'
build_DSEM <- function(sem = NULL,
                       family = "normal",
                       all_vars = FALSE,
                       estimate_projection = FALSE
){
  return(list(sem = sem, family = family, all_vars = all_vars, estimate_projection = estimate_projection))
}


#' Function to build the map and parameter objects for DSEM recruitment linkages
#'
#' @param dsem_settings dsem specifications from \code{\link{build_DSEM}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA
#' @param data_list a data_list created from \code{\link{build_dat}}.
#'
#' @export
build_dsem_objects <- function(dsem_settings = NULL, debug = FALSE, data_list = NULL){

  # Build IID sem if NULL
  if(is.null(dsem_settings$sem)){
    sem = c()
    for(sp in 1:data_list$nspp){
      sem <- c(sem, paste0("recdevs", sp, " <-> recdevs", sp, ", 0, sigmaR", sp,", 1\n")) # No space after
    }
    sem <- paste0(sem, collapse = " ")
    dsem_settings$sem <- sem
  }

  # DSEM data
  dsem_data <- data_list$env_data %>%
    # Adding NA in missing years (match assessment begining)
    dplyr::full_join(data.frame(Year=c(data_list$styr:data_list$projyr)), by = join_by(Year)) %>%
    dplyr::arrange(Year) %>%
    dplyr::select(-Year)

  # - Add column for recdev of each species
  for(sp in data_list$nspp:1){
    dsem_data <- dsem_data %>%
      dplyr::mutate(recdevs=NA) %>%
      dplyr::relocate(recdevs)
    colnames(dsem_data)[1] <- paste0("recdevs", sp)
  }

  # DSEM family
  if(length(dsem_settings$family) == 1){
    dsem_settings$family <- rep(dsem_settings$family, ncol(dsem_data))
  }

  if(length(dsem_settings$family) != ncol(dsem_data)){
    stop("Length of 'family' in 'build_DSEM' does not equal 1 or `ncol(env_data) + nspp`")
  }

  # Specify control (i.e. turn off estimation)
  DSEMcontrol <- dsem::dsem_control(use_REML = FALSE,
                                    run_model = FALSE,
                                    quiet = TRUE,
                                    getJointPrecision = TRUE,
                                    newton_loops = 0)

  # Run "dsem" to output map and parameter object
  fit_dsem = dsem::dsem(sem = dsem_settings$sem,
                        tsdata = ts(dsem_data),
                        family = dsem_settings$family,
                        control = DSEMcontrol)


  # Extract dsem map and parameter objects
  fit_dsem$tmb_inputs$map$lnsigma_j <- factor(rep(NA, length=length(fit_dsem$tmb_inputs$map$lnsigma_j))) #FIXME: Not sure why we turn this off?
  fit_dsem$tmb_inputs$parameters$lnsigma_j <- rep(log(0.1), length=length(fit_dsem$tmb_inputs$parameters$lnsigma_j)) #FIXME: Not sure why we turn this off?


  # Create mapList object
  mapList <- sapply(fit_dsem$tmb_inputs$parameters, function(x) replace(x, values = c(1:length(x))))

  # - Convert dsem map-factor to map-list
  for(i in 1:length(fit_dsem$tmb_inputs$map)){
    parname <- names(fit_dsem$tmb_inputs$map)[i]
    mapList[[parname]] <- replace(fit_dsem$tmb_inputs$parameters[[parname]], values = as.numeric(fit_dsem$tmb_inputs$map[[i]]))
  }

  # Debug mode
  if(debug){
    mapList <- sapply(mapList, function(x) replace(x, values = rep(NA, length=length(x))))
  }

  # Recruitment variance (turn off variance estimation for recdevs)
  fit_dsem$tmb_inputs$parameters$beta_z[1:data_list$nspp] <- data_list$sigma_rec_prior
  if(!data_list$random_rec){
    mapList$beta_z[1:data_list$nspp] <- NA
  }

  # Map all variables and projection
  if(!dsem_settings$all_vars){
    # DSEM parameters
    pars_off <- which(is.na(fit_dsem$sem_full$start))
    mapList$beta_z[pars_off] <- NA

    # - Latent states x_tj
    x_tj_off <- fit_dsem$sem_full %>%
      dplyr::filter(is.na(start))
    x_tj_off <- unique(c(x_tj_off$first, x_tj_off$second))
    x_tj_off <- which(colnames(dsem_data) %in% x_tj_off)
    mapList$x_tj[,x_tj_off] <- NA
    mapList$mu_j[x_tj_off] <- NA
    mapList$lnsigma_j[x_tj_off] <- NA
  }

  if(!dsem_settings$estimate_projection){
    hind_years <- data_list$styr:data_list$endyr
    all_years <- data_list$styr:data_list$projyr
    mapList$x_tj[length(hind_years+1):length(all_years),] <- NA # Turn off latent states for future
  }

  # Return
  fit_dsem$tmb_inputs$map <- sapply(mapList, function(x) factor(x))
  fit_dsem$mapList <- mapList
  fit_dsem$sem <- dsem_settings$sem
  return(fit_dsem = fit_dsem)
}
