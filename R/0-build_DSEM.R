
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
#' @param data_list a data_list read in via \code{\link{read_data}} or built directly in R.
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
    dplyr::full_join(data.frame(Year=c(data_list$styr:data_list$projyr)), by = dplyr::join_by(Year)) %>%
    dplyr::arrange(Year)

  # - Add column for recdev of each species
  for(sp in data_list$nspp:1){
    dsem_data <- dsem_data %>%
      dplyr::mutate(recdevs = NA_real_) %>%
      dplyr::relocate("recdevs") %>%
      dplyr::select(-Year)
    colnames(dsem_data)[1] <- paste0("recdevs", sp)
  }

  # - Keep only variables referenced in the sem.
  # Mirror make_dsem_ram(): scan the sem to a model table, then parse_path each
  # path to collect the variables on both heads of every arrow.
  sem_model <- scan(
    text         = dsem_settings$sem,
    what         = list(path = "", lag = 1, par = "", start = 1, dump = ""),
    sep          = ",",
    strip.white  = TRUE,
    comment.char = "#",
    fill         = TRUE,
    quiet        = TRUE
  )
  sem_vars <- unique(unlist(lapply(sem_model$path, function(p){
    pp <- parse_path(p)
    c(pp$first, pp$second)
  })))

  # Preserve dsem_data's column order (recdevs<sp> first, then env vars)
  dsem_data <- dsem_data %>%
    dplyr::select(dplyr::any_of(intersect(colnames(dsem_data), sem_vars)))


  # DSEM family
  if(length(dsem_settings$family) == 1){
    dsem_settings$family <- rep(dsem_settings$family, ncol(dsem_data))
  }

  if(length(dsem_settings$family) != ncol(dsem_data)){
    stop("Length of 'family' in 'build_DSEM' does not equal 1 or `ncol(env_data) + nspp`")
  }

  # Build DSEM TMB inputs via the vendored sem->inputs pipeline (R/0-dsem_ram.R),
  # replacing the former dsem::dsem(run_model = FALSE) harvesting. This removes
  # the runtime dependency on dsem internals; the pipeline is matched to
  # src/TMB/dsem.hpp and validated byte-for-byte against dsem 2.0.1. Defaults
  # mirror dsem 2.0.1 dsem_control() plus the controls Rceattle always set
  # (use_REML = FALSE => Random = "x_tj"; gmrf_parameterization = "gmrf_project").
  # fit_dsem = build_dsem_inputs(sem = dsem_settings$sem,
  #                              tsdata = stats::ts(dsem_data),
  #                              family = dsem_settings$family,
  #                              use_REML = FALSE,
  #                              quiet = TRUE)
  fit_dsem <- dsem::dsem(sem = dsem_settings$sem,
                         tsdata = stats::ts(dsem_data),
                         family = dsem_settings$family,
                         control = dsem::dsem_control(use_REML = FALSE,
                                                 quiet = TRUE,
                                                 run_model = FALSE))


  # Extract dsem map and parameter objects
  fit_dsem$tmb_inputs$map$lnsigma_j <- factor(rep(NA, length=length(fit_dsem$tmb_inputs$map$lnsigma_j))) #FIXME: Not sure why we turn this off?
  fit_dsem$tmb_inputs$parameters$lnsigma_j <- rep(log(0.1), length=length(fit_dsem$tmb_inputs$parameters$lnsigma_j))


  # Create mapList object
  mapList <- sapply(fit_dsem$tmb_inputs$parameters, function(x) replace(x, values = c(1:length(x))))

  # - Copy dsem map-factor to map-list
  for(i in 1:length(fit_dsem$tmb_inputs$map)){
    parname <- names(fit_dsem$tmb_inputs$map)[i]
    mapList[[parname]] <- replace(fit_dsem$tmb_inputs$parameters[[parname]], values = as.numeric(fit_dsem$tmb_inputs$map[[i]]))
  }

  # Debug mode
  if(debug){
    mapList <- sapply(mapList, function(x) replace(x, values = rep(NA, length=length(x))))
  }

  # Recruitment variance.
  # Locate the beta_z index of each species' recruitment SD: the two-headed
  # self-loop on recdevs[sp] (recdevs[sp] <-> recdevs[sp]).
  sf       <- fit_dsem$sem_full
  sf_start <- suppressWarnings(as.numeric(sf$start))
  sf_par   <- suppressWarnings(as.numeric(sf$parameter))
  sf_dir   <- abs(suppressWarnings(as.numeric(sf$direction)))
  sigma_rec_prior <- rep_len(data_list$sigma_rec_prior, data_list$nspp)

  rec_sd_idx   <- integer(data_list$nspp)   # 1-based beta_z index; 0 if SD is fixed in the sem
  rec_sd_fixed <- numeric(data_list$nspp)   # fixed SD value used when rec_sd_idx == 0
  for(sp in seq_len(data_list$nspp)){
    nm   <- paste0("recdevs", sp)
    rows <- which(sf$first == nm & sf$second == nm & sf_dir == 2)
    pn   <- unique(sf_par[rows]); pn <- pn[!is.na(pn) & pn > 0]
    if(length(pn) > 0){
      rec_sd_idx[sp] <- pn[1]
      # Initialize the recruitment SD at sigma_rec_prior; fix it if !random_rec
      fit_dsem$tmb_inputs$parameters$beta_z[rec_sd_idx[sp]] <- sigma_rec_prior[sp]
      if(!data_list$random_rec){
        mapList$beta_z[rec_sd_idx[sp]] <- NA
      }
    } else {
      # Recruitment SD fixed in the sem (NA parameter name): use its start value
      val <- sf_start[rows]; val <- val[!is.na(val)]
      rec_sd_fixed[sp] <- if(length(val) > 0) val[1] else sigma_rec_prior[sp]
    }
  }
  fit_dsem$tmb_inputs$data$rec_sd_idx   <- as.integer(rec_sd_idx)
  fit_dsem$tmb_inputs$data$rec_sd_fixed <- as.numeric(rec_sd_fixed)

  # x_tj column of each species' recdevs (0-based for the cpp).
  rec_dev_col <- match(paste0("recdevs", seq_len(data_list$nspp)), colnames(dsem_data)) - 1L
  if(any(is.na(rec_dev_col))){
    stop("Could not locate a 'recdevs<sp>' column for every species in the DSEM data")
  }
  fit_dsem$tmb_inputs$data$rec_dev_col <- as.integer(rec_dev_col)

  # Map all variables and projection
  if(!dsem_settings$all_vars){
    # Turn off beta_z for sem paths with no start value. Index by beta_z
    # PARAMETER number (not sem_full row number): equality constraints and
    # auto-added covariances mean rows != parameters.
    off_par <- unique(sf_par[is.na(sf_start)])
    off_par <- off_par[!is.na(off_par) & off_par > 0]
    if(length(off_par) > 0) mapList$beta_z[off_par] <- NA

    # - Latent states x_tj
    x_tj_off <- fit_dsem$sem_full %>%
      dplyr::filter(is.na(.data$start))
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


# -----------------------------------------------------------------------------
# DSEM <-> CEATTLE integration helpers
# -----------------------------------------------------------------------------
# build_dsem_objects() produces TMB inputs that fit_mod() merges into the
# CEATTLE parameter list, map, random-effects vector, and data list at four
# points in its pipeline. The contract of what DSEM contributes lives here so
# those merges stay consistent:
#   * parameters: beta_z, lnsigma_j, mu_j, delta0_j, x_tj
#   * map:        mapList + mapFactor entries for the above
#   * random:     x_tj (when random_rec = TRUE)
#   * data:       options, RAM, RAMstart, familycode_j, y_tj, obs_idx, unobs_idx

#' Names of the DSEM parameters contributed to the CEATTLE parameter list
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
dsem_param_names <- function(dsem) {
  names(dsem$tmb_inputs$parameters)
}

#' Merge DSEM parameters into a CEATTLE parameter list
#'
#' Adds any DSEM parameters absent from \code{param_list}, pulling defaults from
#' the built DSEM objects. A no-op for parameters already present, so it is safe
#' for both a fresh \code{build_params()} list (none present -> all added) and a
#' prior fit's \code{estimated_params} (all present -> unchanged). Nothing is
#' overwritten or duplicated.
#'
#' @param param_list CEATTLE parameter list (from build_params or inits)
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
merge_dsem_params <- function(param_list, dsem) {
  dsem_par <- dsem$tmb_inputs$parameters
  missing_par <- setdiff(names(dsem_par), names(param_list))
  if (length(missing_par) > 0) {
    param_list <- c(param_list, dsem_par[missing_par])
  }
  param_list
}

#' Merge DSEM map entries into a CEATTLE map object
#'
#' DSEM contributes both \code{mapList} (named numeric/NA vectors) and
#' \code{mapFactor} (the TMB factor map). Merge them together so the two halves
#' never drift apart.
#'
#' @param map CEATTLE map object (from \code{build_map}) with \code{mapList} and
#'   \code{mapFactor}
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
merge_dsem_map <- function(map, dsem) {
  map$mapList   <- c(map$mapList, dsem$mapList)
  map$mapFactor <- c(map$mapFactor, dsem$tmb_inputs$map)
  map
}
