#' Retrospective peels
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model. The function also evaluates retrospective forecast skill. To evaluate both retrospective bias and forecast skill, the function uses the map functionality of TMB to peel the model:
#' 1. Filters data, filters fixed inputs, and maps out time-varying parameters for the peeled years. All time-varying parameters for the peeled years are set to the terminal year of the model for that peel.
#' 2. Fits the peeled model.
#' 3. Turns off all hindcast parameters, turns on F for the peeled years, and fits to the peeled catch series to update the "forecast" dynamics given projection assumptions and observed catch from the peeled years.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param peels the number of retrospective peels to use in the calculation of rho and for model estimation
#' @param rescale TRUE/FALSE whether to subset and rescale environmental predictors for the range of peel years.
#' @param nyrs_forecast Number of forecast years to calculate Mohn's Rho in addition to terminal year
#' @param cores Number of cores to use for parallel peels. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each peel runs \code{TMB::sdreport} (standard errors).
#'   Mohn's rho uses only point estimates, so \code{FALSE} is faster with no
#'   effect on rho. Default \code{NULL} inherits the input model's setting
#'   (\code{TRUE} if it was fit with \code{getsd = TRUE}, i.e. carries an
#'   \code{sdrep}); the returned peel models then carry standard errors only
#'   when \code{getsd} is \code{TRUE}.
#'
#' @return a list of 1. list of Rceattle models and 2. vector of Mohn's rho for each species
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' retro <- retrospective(ss_run, peels = 10)
#' }
#' @export
retrospective <- function(Rceattle = NULL, peels = 5, rescale = FALSE, nyrs_forecast = 3, cores = NULL, getsd = NULL) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # Peels inherit the input model's sdreport setting unless overridden. Mohn's
  # rho reads only point estimates, so getsd = FALSE is faster and rho-neutral.
  if (is.null(getsd)) getsd <- !is.null(Rceattle$sdrep)

  # Get objects
  Rceattle$data_list$endyr_peel <- Rceattle$data_list$endyr
  data_list <- Rceattle$data_list # used by Mohn's rho block below
  endyr <- Rceattle$data_list$endyr
  styr <- Rceattle$data_list$styr
  nyrs <- length(styr:endyr)
  projyr <- Rceattle$data_list$projyr
  nyrs_proj <- projyr - styr + 1

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- peels > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-peel closure ----
  # Each peel only reads the original Rceattle, so peels are independent.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_peel <- function(i) {

    # * Get end year of peel ----
    data_list <- Rceattle$data_list
    endyr_peel <- endyr - i
    data_list$endyr_peel <- endyr_peel
    nyrs_peel <- endyr_peel - styr + 1
    peel_prj_yrs <- (endyr_peel+1):endyr
    nyrs_proj_peel <- length(peel_prj_yrs)


    # * Turn off data after endyr_peel ----
    data_list$index_data <- data_list$index_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$comp_data <- data_list$comp_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$caal_data <- data_list$caal_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$diet_data <- data_list$diet_data |>
      dplyr::filter(Year <= endyr_peel)

    peeled_catch_data <- data_list$catch_data |>
      dplyr::filter(Year > endyr_peel)
    data_list$catch_data$Catch[which(data_list$catch_data$Year > endyr_peel)] <- 0 # Set catch data in peeled years to 0 to avoid fitting


    # * Turn off fixed inputs after endyr_peel ----
    #FIXME ignores forecasted growth
    data_list$weight <- data_list$weight |>
      dplyr::filter(Year <= endyr_peel)

    data_list$emp_sel <- data_list$emp_sel |>
      dplyr::filter(Year <= endyr_peel)

    data_list$ration_data <- data_list$ration_data |>
      dplyr::filter(Year <= endyr_peel)

    # * Extend fixed inputs for "projection years"
    # - Assume weight/ration/empirical sel is same as last year of peel

    # - Weight
    #FIXME ignores forecasted growth
    proj_wt <- data_list$weight |>
      dplyr::filter(Year != 0)

    if(nrow(proj_wt) > 0){
      proj_wt <- proj_wt |>
        dplyr::group_by(Wt_index , Sex) |>
        dplyr::slice(rep(dplyr::n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
    }
    data_list$weight  <- rbind(data_list$weight, proj_wt) |>
      dplyr::arrange(Wt_index, Year)

    # - Empirical selectivity
    proj_emp_sel <- data_list$emp_sel |>
      dplyr::filter(Year != 0)

    if(nrow(proj_emp_sel) > 0){
      proj_emp_sel <- proj_emp_sel |>
        dplyr::group_by(Fleet_code, Sex) |>
        dplyr::slice(rep(n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
      data_list$emp_sel  <- rbind(data_list$emp_sel, proj_emp_sel) |>
        dplyr::arrange(Fleet_code, Year)
    }

    # - Ration data
    proj_ration_data <- data_list$ration_data |>
      dplyr::filter(Year != 0)

    if(nrow(proj_ration_data) > 0){
      proj_ration_data <- proj_ration_data |>
        dplyr::group_by(Species, Sex) |>
        dplyr::slice(rep(dplyr::n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
      data_list$ration_data  <- rbind(data_list$ration_data, proj_ration_data) |>
        dplyr::arrange(Species, Year)
    }

    # * Rescale environmental predictors ----
    if(rescale){
      data_list$env_data <- data_list$env_data |>
        dplyr::filter(Year <= endyr_peel)
      data_list$env_data[,2:ncol(data_list$env_data)]<-scale(data_list$env_data[,2:ncol(data_list$env_data)])
    }

    # * Adjust parameters ----
    #FIXME: adjust for forecasting via MVN
    inits <- Rceattle$estimated_params
    inits$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- 0
    inits$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- inits$log_M1_dev[,,,nyrs_peel]
    inits$index_q_dev[,(nyrs_peel+1):nyrs] <- inits$index_q_dev[,nyrs_peel]
    inits$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- inits$log_sel_slp_dev[,,,nyrs_peel]
    inits$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_inf_dev[,,,nyrs_peel]
    inits$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_coff_dev[,,,nyrs_peel]

    # * Adjust map size ----
    # Turn off forecasted parameters
    map <- Rceattle$map
    map$mapList$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- NA
    map$mapFactor$rec_dev <- factor(map$mapList$rec_dev)

    map$mapList$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- NA
    map$mapFactor$log_M1_dev <- factor(map$mapList$log_M1_dev)

    map$mapList$index_q_dev[,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$index_q_dev <- factor(map$mapList$index_q_dev)

    map$mapList$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$log_sel_slp_dev <- factor(map$mapList$log_sel_slp_dev)

    map$mapList$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_inf_dev <- factor(map$mapList$sel_inf_dev)

    map$mapList$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)

    # -- Map out Fdev for years with 0 catch to very low number
    zero_catch <- data_list$catch_data |>
      dplyr::filter(Year <= endyr &
                      Catch == 0) |>
      dplyr::mutate(Year = Year - styr + 1) |>
      dplyr::select(Fleet_code, Year) |>
      as.matrix()
    inits$log_F[zero_catch] <- -999
    map$mapList$log_F[zero_catch] <- NA
    map$mapFactor$log_F <- factor(map$mapList$log_F)

    # * Refit ----
    newmod <- suppressWarnings(
      suppressMessages(
        Rceattle::fit_mod(
          data_list = data_list,
          inits = inits,
          map =  map,
          bounds = NULL,
          file = NULL,
          estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
          HCR = build_hcr(HCR = data_list$HCR,
                          DynamicHCR = data_list$DynamicHCR,
                          Ftarget = data_list$Ftarget,
                          Flimit = data_list$Flimit,
                          Ptarget = data_list$Ptarget,
                          Plimit = data_list$Plimit,
                          Alpha = data_list$Alpha,
                          Pstar = data_list$Pstar,
                          Sigma = data_list$Sigma,
                          Fmult = data_list$Fmult,
                          HCRorder = data_list$HCRorder
          ),
          # suppressWarnings: legacy srr_fun = 1|3|5 / srr_indices.
          recFun = suppressWarnings(build_srr(
            srr_fun = data_list$srr_fun,
            srr_pred_fun  = data_list$srr_pred_fun,
            proj_mean_rec  = data_list$proj_mean_rec,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
            srr_hat_styr = data_list$srr_hat_styr,
            srr_hat_endyr = min(data_list$srr_hat_endyr, endyr_peel),
            srr_est_mode  = data_list$srr_est_mode,
            srr_prior  = data_list$srr_prior,
            srr_prior_sd   = data_list$srr_prior_sd,
            Bmsy_lim = data_list$Bmsy_lim,
            srr_indices = data_list$srr_indices,
            linkages = data_list$srr_linkages)),
          # suppressWarnings: legacy M1_indices may travel via data_list.
          M1Fun = suppressWarnings(build_M1(
            M1_model = data_list$M1_model,
            M1_re = data_list$M1_re,
            updateM1 = FALSE,
            M1_use_prior = data_list$M1_use_prior,
            M2_use_prior = data_list$M2_use_prior,
            M_prior = data_list$M_prior,
            M_prior_sd = data_list$M_prior_sd,
            M1_indices = data_list$M1_indices,
            linkages = data_list$M1_linkages)),
          growthFun = build_growth(fun = data_list$growth_fun,
                                   linkages = data_list$growth_linkages),
          random_rec = data_list$random_rec,
          niter = data_list$niter,
          msmMode = data_list$msmMode,
          avgnMode = data_list$avgnMode,
          suitMode = data_list$suitMode,
          suit_styr = data_list$suit_styr,
          suit_endyr = min(data_list$suit_endyr, endyr_peel),   # Update to end year if less than suit_endyr
          initMode = data_list$initMode,
          fit_control = fit_control(
            phase   = TRUE, # Phasing or else the parameters dont wanna move
            loopnum = data_list$loopnum,
            getsd   = getsd,
            verbose = 0))
      )
    )

    # Forecast ----
    peeled_pars <- newmod$estimated_params

    # - Add in peeled catch to fit to
    data_list$catch_data <- data_list$catch_data |>
      dplyr::filter(Year <= endyr_peel) |>
      rbind(peeled_catch_data) |>
      dplyr::arrange(Fleet_code, Year)

    # - Update map
    # Only new parameter we are estimating in is the log_F of the peeled years
    map <- build_map(
      data_list = data_list,
      params = newmod$estimated_params,
      debug = TRUE,
      random_rec = newmod$data_list$random_rec)
    map$mapFactor$dummy <- as.factor(NA); map$mapList$dummy <- NA

    # - Turn on F for peeled years to fit to catch (matches full model)
    peeled_pars$log_F[,(nyrs_peel+1):nyrs] <- Rceattle$estimated_params$log_F[,(nyrs_peel+1):nyrs]
    map$mapList$log_F[,(nyrs_peel+1):nyrs] <- Rceattle$map$mapList$log_F[,(nyrs_peel+1):nyrs]
    map$mapFactor$log_F <-  factor(map$mapList$log_F)

    # Adjust forecased rec_dev in new mod for bias and refit
    for(sp in 1:newmod$data_list$nspp){

      # -- where SR curve is estimated directly
      if(newmod$data_list$srr_fun == newmod$data_list$srr_pred_fun){
        rec_dev <- log(mean(newmod$quantities$R[sp,1:nyrs_peel]))  - log(newmod$quantities$R0[sp])
      }

      # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
      if(newmod$data_list$srr_fun != newmod$data_list$srr_pred_fun){
        # Already a log-scale deviation, so take the mean directly.
        rec_dev <- mean((log(newmod$quantities$R) - log(newmod$quantities$R_hat))[sp, 1:nyrs_peel])

      }

      # - Update OM with devs
      peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)] <- replace(
        peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)],
        values =  rec_dev)
    }

    newmod <- suppressMessages(
      suppressWarnings(
        Rceattle::fit_mod(
          data_list = data_list,
          inits = peeled_pars,
          map =  map,
          bounds = NULL,
          file = NULL,
          estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
          HCR = build_hcr(HCR = data_list$HCR,
                          DynamicHCR = data_list$DynamicHCR,
                          Ftarget = data_list$Ftarget,
                          Flimit = data_list$Flimit,
                          Ptarget = data_list$Ptarget,
                          Plimit = data_list$Plimit,
                          Alpha = data_list$Alpha,
                          Pstar = data_list$Pstar,
                          Sigma = data_list$Sigma,
                          Fmult = data_list$Fmult,
                          HCRorder = data_list$HCRorder
          ),
          # suppressWarnings: legacy srr_fun = 1|3|5 / srr_indices.
          recFun = suppressWarnings(build_srr(
            srr_fun = data_list$srr_fun,
            srr_pred_fun  = data_list$srr_pred_fun,
            proj_mean_rec  = data_list$proj_mean_rec,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
            srr_hat_styr = data_list$srr_hat_styr,
            srr_hat_endyr = min(data_list$srr_hat_endyr, endyr_peel),
            srr_est_mode  = data_list$srr_est_mode,
            srr_prior  = data_list$srr_prior,
            srr_prior_sd   = data_list$srr_prior_sd,
            Bmsy_lim = data_list$Bmsy_lim,
            srr_indices = data_list$srr_indices,
            linkages = data_list$srr_linkages)),
          # suppressWarnings: legacy M1_indices may travel via data_list.
          M1Fun = suppressWarnings(build_M1(
            M1_model = data_list$M1_model,
            M1_re = data_list$M1_re,
            updateM1 = FALSE,
            M1_use_prior = data_list$M1_use_prior,
            M2_use_prior = data_list$M2_use_prior,
            M_prior = data_list$M_prior,
            M_prior_sd = data_list$M_prior_sd,
            M1_indices = data_list$M1_indices,
            linkages = data_list$M1_linkages)),
          growthFun = build_growth(fun = data_list$growth_fun,
                                   linkages = data_list$growth_linkages),
          random_rec = data_list$random_rec,
          niter = data_list$niter,
          msmMode = data_list$msmMode,
          avgnMode = data_list$avgnMode,
          suitMode = data_list$suitMode,
          suit_styr = data_list$suit_styr,
          suit_endyr = min(data_list$suit_endyr, endyr_peel),   # Update to end year if less than suit_endyr
          initMode = data_list$initMode,
          fit_control = fit_control(
            phase   = TRUE, # Phasing or else the parameters dont wanna move
            loopnum = data_list$loopnum,
            getsd   = getsd,
            verbose = 0))
      )
    )

    # gc()
    #
    # map$mapFactor <- map$mapFactor[names(newmod$map$mapFactor)]
    # check <- c()
    # check_na <- c()
    # for(j in 1:length(map$mapList)){
    #   check[j] <- sum(map$mapFactor[[j]] != newmod$map$mapFactor[[j]], na.rm = TRUE)
    #   check_na[j] <- sum(is.na(map$mapFactor[[j]]) != is.na(newmod$map$mapFactor[[j]]), na.rm = TRUE)
    # }

    # Return model only if converged, else NULL (dropped post-dispatch)
    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        return(newmod)
      }
    }
    return(NULL)
  } # End run_one_peel closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch peels (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    cl <- parallel::makeCluster(min(cores, peels))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(Rceattle))
    })
    parallel::clusterExport(
      cl,
      varlist = ls(envir = environment()),
      envir = environment()
    )
    peel_results <- parallel::parLapply(cl, 1:peels, run_one_peel)
  } else {
    peel_results <- lapply(1:peels, run_one_peel)
  }

  # Drop non-converged peels and prepend the original model
  peel_results <- peel_results[!sapply(peel_results, is.null)]
  mod_list <- c(list(Rceattle), peel_results)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Calculate Mohs rho ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # * Data frame to save ----
  objects <- c("biomass", "ssb", "R", "F_spp")
  mohns <- data.frame(matrix(0, nrow = length(objects) * (nyrs_forecast+1), ncol = 3 + data_list$nspp))
  colnames(mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  # * Loop through peels ----
  for (i in 1:(length(mod_list) - 1)) {
    endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
    ind <- 1

    # * Loop output ----
    for (j in 1:length(objects)) {
      term_quantities <- mod_list[[1]]$quantities[[objects[j]]]
      retro_quantities <- mod_list[[i + 1]]$quantities[[objects[j]]]

      # * Loop forecast years ----
      for(yr in 0:nyrs_forecast){

        # If data exist for forecast (save)
        if(endyr_peel + yr <= endyr){

          # * Get full and peeled models ----
          base <- term_quantities[, nyrs_peel + yr]
          peel <- retro_quantities[, nyrs_peel + yr]
          rel_error <- ((peel - base)/base)

          # * Save and sum relative error ----
          mohns[ind, 1] <- objects[j]         # Object
          mohns[ind, 2] <- yr                 # Year
          mohns[ind, 3] <- mohns[ind, 3] + 1  # N
          mohns[ind, 4:(data_list$nspp + 3) ] <- mohns[ind, 4:(data_list$nspp + 3)] + rel_error # Relative error
        }
        ind = ind+1
      }
    }
  }

  # * Divide N ----
  mohns[, 4:(data_list$nspp + 3) ] <- mohns[, 4:(data_list$nspp + 3)]/mohns[, 3]



  # # * Beta coefficients ----
  # objects <- colnames(Rceattle$estimated_params$beta_rec_pars)
  # beta_mohns <- data.frame(matrix(0, nrow = length(objects), ncol = 3 + data_list$nspp))
  # colnames(beta_mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  #
  #   # * Loop through peels ----
  #   for (i in 1:(length(mod_list) - 1)) {
  #     endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
  #     nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
  #     ind <- 1
  #
  #     # * Loop output ----
  #     for (j in 1:length(objects)) {
  #       base <- mod_list[[1]]$estimated_params$beta_rec_pars[,j]
  #       peel <- mod_list[[i + 1]]$estimated_params$beta_rec_pars[,j]
  #       rel_error <- ((peel - base)/base)
  #
  #       # * Save and sum relative error ----
  #       beta_mohns[j, 1] <- objects[j]        # Object
  #       beta_mohns[j, 2] <- 0                 # Year
  #       beta_mohns[j, 3] <- beta_mohns[j, 3] + 1   # N
  #       beta_mohns[j, 4:(data_list$nspp + 3) ] <- beta_mohns[j, 4:(data_list$nspp + 3)] + rel_error # Relative error
  #     }
  #   }

  # * Divide N ----
  # beta_mohns[, 4:(data_list$nspp + 3) ] <- beta_mohns[, 4:(data_list$nspp + 3)]/beta_mohns[, 3]

  mod_list <- rev(mod_list)
  names(mod_list) <- paste0("Year_", (endyr - peels):endyr )

  return(list(Rceattle_list = mod_list, mohns = rbind(mohns))) #, beta_mohns)))
}




#' Jitter analysis
#'
#' @description Run's the Rceattle model at initial values that are +- N(0, sd) from the initial parameters.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param njitter the number of jitters to run
#' @param sd standard deviation for jitter (default = 0.2)
#' @param phase as in \code{\link{fit_mod}} default = FALSE
#' @param seed random number seed. Each jitter \code{i} uses \code{seed + i}
#'   so results are reproducible under both sequential and parallel execution.
#' @param cores Number of cores to use for parallel jitters. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each jitter runs \code{TMB::sdreport}. Jitter compares
#'   objectives and point estimates across starts, so \code{FALSE} is faster
#'   with no effect on that comparison. Default \code{NULL} inherits the input
#'   model's setting (\code{TRUE} only if it carries an \code{sdrep}).
#'
#' @return a list of Rceattle models
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' jitters <- jitter(ss_run, njitter = 10)
#' }
#' @export
jitter <- function(Rceattle = NULL, njitter = 50, sd = 0.2, phase = FALSE, seed = 123, cores = NULL, getsd = NULL) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # Jitters inherit the input model's sdreport setting unless overridden;
  # multimodality is judged from objectives and point estimates, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(Rceattle$sdrep)

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- njitter > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-jitter closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_jitter <- function(i) {

    set.seed(seed + i) # unique seed per jitter for reproducibility under parallel

    # * Adjust initial values ----
    inits <- Rceattle$initial_params
    mapList <- Rceattle$map$mapList
    data_list <- Rceattle$data_list

    for(j in 1:length(inits)){
      par <- names(inits)[j]
      inits[[j]] <- replace(inits[[j]],
                            values = ifelse(is.na(as.numeric(mapList[[par]])),
                                            as.numeric(inits[[j]]),
                                            as.numeric(inits[[j]]) + stats::rnorm(length(as.numeric(inits[[j]])), 0, sd))
      )
    }


    # * Refit ----
    newmod <-
      suppressMessages(
        suppressWarnings(
          Rceattle::fit_mod(
            data_list = data_list,
            inits = inits,
            map =  NULL,
            bounds = NULL,
            file = NULL,
            estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
            HCR = build_hcr(HCR = data_list$HCR,
                            DynamicHCR = data_list$DynamicHCR,
                            Ftarget = data_list$Ftarget,
                            Flimit = data_list$Flimit,
                            Ptarget = data_list$Ptarget,
                            Plimit = data_list$Plimit,
                            Alpha = data_list$Alpha,
                            Pstar = data_list$Pstar,
                            Sigma = data_list$Sigma,
                            Fmult = data_list$Fmult,
                            HCRorder = data_list$HCRorder
            ),
            # suppressWarnings: legacy srr_fun = 1|3|5 / srr_indices.
            recFun = suppressWarnings(build_srr(
              srr_fun = data_list$srr_fun,
              srr_pred_fun  = data_list$srr_pred_fun,
              proj_mean_rec  = data_list$proj_mean_rec,
              srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
              srr_hat_styr = data_list$srr_hat_styr,
              srr_hat_endyr = data_list$srr_hat_endyr,
              srr_est_mode  = data_list$srr_est_mode,
              srr_prior  = data_list$srr_prior,
              srr_prior_sd   = data_list$srr_prior_sd,
              Bmsy_lim = data_list$Bmsy_lim,
              srr_indices = data_list$srr_indices,
              linkages = data_list$srr_linkages)),
            # suppressWarnings: legacy M1_indices may travel via data_list.
            M1Fun = suppressWarnings(build_M1(
              M1_model = data_list$M1_model,
              M1_re = data_list$M1_re,
              updateM1 = FALSE,
              M1_use_prior = data_list$M1_use_prior,
              M2_use_prior = data_list$M2_use_prior,
              M_prior = data_list$M_prior,
              M_prior_sd = data_list$M_prior_sd,
              M1_indices = data_list$M1_indices,
              linkages = data_list$M1_linkages)),
            growthFun = build_growth(fun = data_list$growth_fun,
                                     linkages = data_list$growth_linkages),
            random_rec = data_list$random_rec,
            niter = data_list$niter,
            msmMode = data_list$msmMode,
            avgnMode = data_list$avgnMode,
            suitMode = data_list$suitMode,
            suit_styr = data_list$suit_styr,
            suit_endyr = min(data_list$suit_endyr, data_list$endyr),   # Update to end year if less than suit_endyr
            initMode = data_list$initMode,
            fit_control = fit_control(
              phase   = phase,
              loopnum = data_list$loopnum,
              getsd   = getsd,
              verbose = 0))
        )
      )

    # Return model only if converged, else NULL (dropped post-dispatch)
    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        return(newmod)
      }
    }
    return(NULL)
  } # End run_one_jitter closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch jitters (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    cl <- parallel::makeCluster(min(cores, njitter))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(Rceattle))
    })
    parallel::clusterExport(
      cl,
      varlist = ls(envir = environment()),
      envir = environment()
    )
    mod_list <- parallel::parLapply(cl, 1:njitter, run_one_jitter)
  } else {
    mod_list <- lapply(1:njitter, run_one_jitter)
  }

  # Drop non-converged
  mod_list <- mod_list[!sapply(mod_list, is.null)]

  # Plot ----
  jnll <- sapply(mod_list, function(x) x$opt$objective)
  # plot(x = 1:length(jnll), y = jnll)
  if (length(mod_list) > 0) {
    names(mod_list) <- paste0("Jitter_", seq_along(mod_list))
  }


  # Return ----
  return(list(Rceattle_list = mod_list, nll = jnll))
}



#' Self test simulation analysis analysis
#'
#' @description Simulates data from an Rceattle model and refits the model to the simulated data. TODO add process variation (i.e. random devs) to simulation.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param seed random number seed. Each simulation \code{i} uses \code{seed + i}
#'   so results are reproducible under both sequential and parallel execution.
#' @param nsim number of simulations
#' @param simulate passed to \code{\link{sim_mod}}. If \code{TRUE} (default),
#'   data are simulated with observation error; if \code{FALSE}, expected
#'   values from the model are used.
#' @param cores Number of cores to use for parallel simulations. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#'
#' @return a list of Rceattle models
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' sims <- self_test(ss_run, nsim = 10)
#' }
#' @export
self_test <- function(Rceattle = NULL, nsim = 50, simulate = TRUE, seed = 123, cores = NULL) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- nsim > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-simulation closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_sim <- function(i) {

    set.seed(seed + i) # unique seed per sim for reproducibility under parallel

    # * Simulate new data
    sim_data <- Rceattle::sim_mod(Rceattle, simulate = simulate)
    data_list <- sim_data

    # * Adjust initial values ----
    inits <- Rceattle$initial_params
    mapList <- Rceattle$map$mapList


    # * Refit ----
    newmod <-
      suppressMessages(
        suppressWarnings(
          Rceattle::fit_mod(
            data_list = data_list,
            inits = inits,
            map =  NULL,
            bounds = NULL,
            file = NULL,
            estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
            HCR = build_hcr(HCR = data_list$HCR,
                            DynamicHCR = data_list$DynamicHCR,
                            Ftarget = data_list$Ftarget,
                            Flimit = data_list$Flimit,
                            Ptarget = data_list$Ptarget,
                            Plimit = data_list$Plimit,
                            Alpha = data_list$Alpha,
                            Pstar = data_list$Pstar,
                            Sigma = data_list$Sigma,
                            Fmult = data_list$Fmult,
                            HCRorder = data_list$HCRorder
            ),
            # suppressWarnings: legacy srr_fun = 1|3|5 / srr_indices.
            recFun = suppressWarnings(build_srr(
              srr_fun = data_list$srr_fun,
              srr_pred_fun  = data_list$srr_pred_fun,
              proj_mean_rec  = data_list$proj_mean_rec,
              srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
              srr_hat_styr = data_list$srr_hat_styr,
              srr_hat_endyr = data_list$srr_hat_endyr,
              srr_est_mode  = data_list$srr_est_mode,
              srr_prior  = data_list$srr_prior,
              srr_prior_sd   = data_list$srr_prior_sd,
              Bmsy_lim = data_list$Bmsy_lim,
              srr_indices = data_list$srr_indices,
              linkages = data_list$srr_linkages)),
            # suppressWarnings: legacy M1_indices may travel via data_list.
            M1Fun = suppressWarnings(build_M1(
              M1_model = data_list$M1_model,
              M1_re = data_list$M1_re,
              updateM1 = FALSE,
              M1_use_prior = data_list$M1_use_prior,
              M2_use_prior = data_list$M2_use_prior,
              M_prior = data_list$M_prior,
              M_prior_sd = data_list$M_prior_sd,
              M1_indices = data_list$M1_indices,
              linkages = data_list$M1_linkages)),
            growthFun = build_growth(fun = data_list$growth_fun,
                                     linkages = data_list$growth_linkages),
            random_rec = data_list$random_rec,
            niter = data_list$niter,
            msmMode = data_list$msmMode,
            avgnMode = data_list$avgnMode,
            suitMode = data_list$suitMode,
            suit_styr = data_list$suit_styr,
            suit_endyr = min(data_list$suit_endyr, data_list$endyr),   # Update to end year if less than suit_endyr
            initMode = data_list$initMode,
            fit_control = fit_control(
              phase   = FALSE,
              loopnum = data_list$loopnum,
              getsd   = getsd,
              verbose = 0))
        )
      )

    # Return model only if converged, else NULL (dropped post-dispatch)
    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        return(newmod)
      }
    }
    return(NULL)
  } # End run_one_sim closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch sims (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    cl <- parallel::makeCluster(min(cores, nsim))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(Rceattle))
    })
    parallel::clusterExport(
      cl,
      varlist = ls(envir = environment()),
      envir = environment()
    )
    mod_list <- parallel::parLapply(cl, 1:nsim, run_one_sim)
  } else {
    mod_list <- lapply(1:nsim, run_one_sim)
  }

  # Drop non-converged
  mod_list <- mod_list[!sapply(mod_list, is.null)]

  # Plot ----
  jnll <- sapply(mod_list, function(x) x$opt$objective)
  # plot(x = 1:length(jnll), y = jnll)
  if (length(mod_list) > 0) {
    names(mod_list) <- paste0("Sim_", seq_along(mod_list))
  }


  # Return ----
  return(mod_list)
}



#' Likelihood profile across one or more parameter cells
#'
#' @description Re-fits an Rceattle model while holding selected cells of a
#'   parameter fixed at user-specified values. Supports profiling a single
#'   cell (e.g. \code{R_log_sd[species = 1]}) and arbitrary N-dimensional
#'   cross-profiles over multiple cells -- e.g. \code{log_M1[1, 1, 1]} and
#'   \code{log_M1[1, 2, 1]} jointly, to profile residual M for males against
#'   females. For each grid point the targeted cells are fixed in the TMB
#'   map and the remaining parameters are re-estimated; the result is a
#'   grid of Rceattle models for downstream NLL surfaces.
#'
#' @param fitted an Rceattle model fit using \code{\link{fit_mod}}
#' @param param Name of the parameter to profile. Two ways to specify it:
#'   \describe{
#'     \item{Raw parameter slot}{any name in
#'       \code{Rceattle$estimated_params}; tested for \code{"R_log_sd"},
#'       \code{"rec_pars"}, and \code{"log_M1"}. \code{slots} must index
#'       into the full array and \code{transform} controls the scale.}
#'     \item{Natural-scale alias}{convenience shortcut for the three
#'       documented parameters. Aliases imply \code{transform = "log"}
#'       (values are taken in natural units and log'd before being
#'       substituted) and, for \code{rec_pars}, fill in the column from
#'       the alias name so \code{slots} only needs the species index:
#'       \itemize{
#'         \item \code{"sigmaR"}, \code{"R_sd"} -> \code{R_log_sd}
#'         \item \code{"M1"} -> \code{log_M1}
#'         \item \code{"R0"} -> \code{rec_pars[, 1]}
#'         \item \code{"alpha"} -> \code{rec_pars[, 2]}
#'         \item \code{"beta"} -> \code{rec_pars[, 3]}
#'       }
#'       If \code{transform} is supplied with an alias it is ignored
#'       (with a warning).}
#'   }
#' @param slots A list whose entries are integer index vectors, one entry
#'   per cell to fix. Each entry's length must equal the number of
#'   dimensions of the resolved parameter -- 1 for vectors
#'   (\code{R_log_sd}), 2 for matrices (\code{rec_pars}), 3 for 3-D arrays
#'   (\code{log_M1}). When using the \code{"R0"}/\code{"alpha"}/\code{"beta"}
#'   aliases, supply only the species index (length 1); the column is
#'   filled in from the alias. E.g. \code{list(c(1, 2, 1))} fixes
#'   \code{log_M1[1, 2, 1]}; \code{list(c(1, 1, 1), c(1, 2, 1))} fixes both
#'   sex cells for a males-vs-females cross-profile of species 1;
#'   \code{list(1, 2)} with \code{param = "sigmaR"} cross-profiles species
#'   1 and 2. If omitted, defaults to a single species-1 slot shaped to
#'   match the resolved parameter (e.g. \code{list(1)} for
#'   \code{R_log_sd}, \code{list(c(1, 1, 1))} for \code{log_M1},
#'   \code{list(1)} for the \code{rec_pars} aliases) and emits a warning;
#'   pass \code{slots} explicitly to silence the warning. Defaulting
#'   requires \code{length(values) == 1L} (otherwise the user must
#'   explicitly say which cell each grid targets).
#' @param values A list of numeric vectors, one per entry of \code{slots}.
#'   The full grid of fits is \code{expand.grid(values)}, so a single slot
#'   gives a 1-D profile and \emph{k} slots give a \emph{k}-D cross-profile.
#' @param transform How to map user values onto the internal parameter scale
#'   before substituting them into \code{inits}. Either \code{"log"}
#'   (default), \code{"identity"}, or a unary function (e.g.
#'   \code{qlogis}). Applied element-wise to every grid value. Aliases
#'   override this with \code{"log"}.
#' @param cores Number of cores to use for parallel fits. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each grid fit runs \code{TMB::sdreport}. The profile
#'   reads only the objective (\code{nll}), so \code{FALSE} is faster with no
#'   effect on the profile. Default \code{NULL} inherits the input model's
#'   setting (\code{TRUE} only if it carries an \code{sdrep}).
#' @param ... Unused; present for consistency with the \code{stats::profile}
#'   generic.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{Rceattle_list}{list of fitted Rceattle models, one per grid
#'       row; entries for non-converged fits are \code{NULL} so positions
#'       stay aligned with \code{grid}.}
#'     \item{grid}{data frame of grid values on the user scale (before
#'       \code{transform}); one column per profiled cell, named
#'       \code{slot_1}, \code{slot_2}, ...}
#'     \item{nll}{numeric vector of joint negative log-likelihoods
#'       (\code{opt$objective}); \code{NA} where the fit did not
#'       converge.}
#'     \item{param}{the profiled parameter name (echoed).}
#'     \item{slots}{the slots list (echoed for downstream plotting).}
#'   }
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#'
#' # 1-D profile of sigmaR for species 1 (alias form -- natural scale)
#' p1 <- profile(ss_run,
#'     param  = "sigmaR",
#'     slots  = list(1),
#'     values = list(seq(0.1, 1.5, by = 0.1)))
#'
#' # Equivalent raw form (log scale -- user does the transform)
#' p1_raw <- profile(ss_run,
#'     param     = "R_log_sd",
#'     slots     = list(1),
#'     values    = list(log(seq(0.1, 1.5, by = 0.1))),
#'     transform = "identity")
#'
#' # 2-D cross-profile of M1 across species 1 and 2 (sex 1, age 1).
#' # BS2017SS is single-sex; with a multi-sex model the same form
#' # (e.g. c(1, 1, 1), c(1, 2, 1)) would cross-profile males vs females.
#' p2 <- profile(ss_run,
#'     param  = "M1",
#'     slots  = list(c(1, 1, 1), c(2, 1, 1)),
#'     values = list(seq(0.1, 0.4, length.out = 3),
#'                   seq(0.1, 0.4, length.out = 3)))
#'
#' # 1-D profile of SRR alpha for species 1 (alias drops the rec_pars column)
#' p3 <- profile(ss_run,
#'     param  = "alpha",
#'     slots  = list(1),
#'     values = list(seq(2, 80, length.out = 20)))
#' }
#' @importFrom stats profile
#' @method profile Rceattle
#' @export
profile.Rceattle <- function(fitted = NULL,
                          param = NULL,
                          slots = NULL,
                          values = NULL,
                          transform = "log",
                          cores = NULL,
                          getsd = NULL,
                          ...) {

  # -- Input validation ----
  if (!inherits(fitted, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }
  # Grid fits inherit the input model's sdreport setting unless overridden;
  # the profile reads only the objective, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(fitted$sdrep)
  if (is.null(param) || !is.character(param) || length(param) != 1L) {
    stop("`param` must be a single character string naming a parameter slot.")
  }
  if (!is.list(values) || length(values) == 0L) {
    stop("`values` must be a non-empty list of numeric grids.")
  }

  # Natural-scale aliases: each maps to a real parameter, implies log()
  # transform, and (for rec_pars aliases) fills in the column index so
  # `slots` only needs the species index.
  alias_table <- list(
    sigmaR = list(param = "R_log_sd", rec_pars_col = NA_integer_),
    R_sd   = list(param = "R_log_sd", rec_pars_col = NA_integer_),
    M1     = list(param = "log_M1",   rec_pars_col = NA_integer_),
    R0     = list(param = "rec_pars", rec_pars_col = 1L),
    alpha  = list(param = "rec_pars", rec_pars_col = 2L),
    beta   = list(param = "rec_pars", rec_pars_col = 3L)
  )

  alias_name   <- NA_character_
  rec_pars_col <- NA_integer_
  if (param %in% names(alias_table)) {
    alias_name <- param
    a <- alias_table[[alias_name]]

    # Aliases force log transform; warn if user passed something else
    if (!identical(transform, "log")) {
      warning(sprintf(
        "`param = \"%s\"` is a natural-scale alias for `%s`; ignoring the supplied `transform` (aliases imply transform = \"log\").",
        alias_name, a$param
      ))
    }
    transform    <- "log"
    rec_pars_col <- a$rec_pars_col
    param        <- a$param   # resolve to real parameter slot
  }

  if (!param %in% names(fitted$estimated_params)) {
    stop("`param` '", param, "' not found in Rceattle$estimated_params.")
  }

  par_array <- fitted$estimated_params[[param]]
  par_ndim  <- if (is.null(dim(par_array))) 1L else length(dim(par_array))

  # Default `slots` to species 1 (a single profile point shaped to match
  # the resolved parameter). For rec_pars aliases the user slot is just
  # the species index; otherwise it's a 1 for every dimension.
  if (is.null(slots)) {
    user_slot_dim <- par_ndim - if (!is.na(rec_pars_col)) 1L else 0L
    default_user_slot <- rep(1L, user_slot_dim)

    if (length(values) != 1L) {
      stop(sprintf(
        "`slots` not supplied but `values` has %d grids -- the species-1 default only covers one cell. Pass `slots` explicitly to profile multiple cells.",
        length(values)
      ))
    }

    pretty_slot <- if (length(default_user_slot) == 1L) {
      as.character(default_user_slot)
    } else {
      paste0("c(", paste(default_user_slot, collapse = ", "), ")")
    }
    warning(sprintf(
      "`slots` not supplied; defaulting to species 1 (slots = list(%s)). Pass `slots` explicitly to silence this warning.",
      pretty_slot
    ))

    slots <- list(default_user_slot)
  }

  if (!is.list(slots) || length(slots) == 0L) {
    stop("`slots` must be a non-empty list of integer index vectors.")
  }
  if (length(values) != length(slots)) {
    stop("`values` must be a list with the same length as `slots`.")
  }

  # Append rec_pars column for rec_pars aliases
  if (!is.na(rec_pars_col)) {
    for (k in seq_along(slots)) {
      if (length(slots[[k]]) != 1L) {
        stop(sprintf(
          "Under alias `\"%s\"`, slots[[%d]] should be a single species index (got length %d). The rec_pars column is filled in from the alias name.",
          alias_name, k, length(slots[[k]])
        ))
      }
      slots[[k]] <- c(as.integer(slots[[k]]), rec_pars_col)
    }
  }

  par_dim <- if (is.null(dim(par_array))) length(par_array) else dim(par_array)

  for (k in seq_along(slots)) {
    if (length(slots[[k]]) != par_ndim) {
      stop(sprintf(
        "slots[[%d]] has length %d but '%s' has %d dimension(s).",
        k, length(slots[[k]]), param, par_ndim
      ))
    }
    if (!all(is.finite(slots[[k]])) || any(slots[[k]] < 1)) {
      stop(sprintf("slots[[%d]] must be a vector of positive integers.", k))
    }
    if (any(slots[[k]] > par_dim)) {
      stop(sprintf(
        "slots[[%d]] = c(%s) is out of bounds for '%s' (dim c(%s)).",
        k,
        paste(slots[[k]], collapse = ", "),
        param,
        paste(par_dim, collapse = ", ")
      ))
    }
  }

  # Build transform fn
  trans_fun <- if (is.function(transform)) {
    transform
  } else if (identical(transform, "log")) {
    log
  } else if (identical(transform, "identity")) {
    function(x) x
  } else {
    stop("`transform` must be \"log\", \"identity\", or a function.")
  }

  # Build grid (user-scale values; transform applied at fit time)
  names(values) <- paste0("slot_", seq_along(values))
  grid <- expand.grid(values, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  ngrid <- nrow(grid)

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (mirrors jitter()/retrospective()). Respect the CRAN core limit.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- ngrid > 1L && cores > 1L

  # Generic [<-]: assign `val` into `arr` at index vector `idx`
  assign_at <- function(arr, idx, val) {
    do.call("[<-", c(list(arr), as.list(idx), list(val)))
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-grid-point closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_point <- function(i) {

    inits     <- fitted$estimated_params
    data_list <- fitted$data_list
    map_obj <- fitted$map

    # Substitute fixed values at each profiled cell
    for (k in seq_along(slots)) {
      inits[[param]] <- assign_at(inits[[param]],
                                  slots[[k]],
                                  trans_fun(grid[i, k]))
    }

    # Force profiled cells to NA
    for (k in seq_along(slots)) {
      map_obj$mapList[[param]] <- assign_at(map_obj$mapList[[param]],
                                            slots[[k]],
                                            NA)
    }
    map_obj$mapFactor <- lapply(map_obj$mapList, factor)

    newmod <-
      suppressMessages(suppressWarnings(
        Rceattle::fit_mod(
          data_list = data_list,
          inits = inits,
          map = map_obj,
          bounds = NULL,
          file = NULL,
          estimateMode = ifelse(data_list$estimateMode < 3, 1, data_list$estimateMode),
          HCR = build_hcr(HCR = data_list$HCR,
                          DynamicHCR = data_list$DynamicHCR,
                          Ftarget = data_list$Ftarget,
                          Flimit = data_list$Flimit,
                          Ptarget = data_list$Ptarget,
                          Plimit = data_list$Plimit,
                          Alpha = data_list$Alpha,
                          Pstar = data_list$Pstar,
                          Sigma = data_list$Sigma,
                          Fmult = data_list$Fmult,
                          HCRorder = data_list$HCRorder),
          # suppressWarnings: legacy srr_fun = 1|3|5 / srr_indices.
          recFun = suppressWarnings(build_srr(
            srr_fun = data_list$srr_fun,
            srr_pred_fun  = data_list$srr_pred_fun,
            proj_mean_rec  = data_list$proj_mean_rec,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
            srr_hat_styr = data_list$srr_hat_styr,
            srr_hat_endyr = data_list$srr_hat_endyr,
            srr_est_mode  = data_list$srr_est_mode,
            srr_prior  = data_list$srr_prior,
            srr_prior_sd   = data_list$srr_prior_sd,
            Bmsy_lim = data_list$Bmsy_lim,
            srr_indices = data_list$srr_indices,
            linkages = data_list$srr_linkages)),
          # suppressWarnings: legacy M1_indices may travel via data_list.
          M1Fun = suppressWarnings(build_M1(
            M1_model = data_list$M1_model,
            M1_re = data_list$M1_re,
            updateM1 = FALSE,
            M1_use_prior = data_list$M1_use_prior,
            M2_use_prior = data_list$M2_use_prior,
            M_prior = data_list$M_prior,
            M_prior_sd = data_list$M_prior_sd,
            M1_indices = data_list$M1_indices,
            linkages = data_list$M1_linkages)),
          growthFun = build_growth(fun = data_list$growth_fun,
                                   linkages = data_list$growth_linkages),
          random_rec = data_list$random_rec,
          niter = data_list$niter,
          msmMode = data_list$msmMode,
          avgnMode = data_list$avgnMode,
          suitMode = data_list$suitMode,
          suit_styr = data_list$suit_styr,
          suit_endyr = min(data_list$suit_endyr, data_list$endyr),
          initMode = data_list$initMode,
          fit_control = fit_control(
            phase   = FALSE,
            loopnum = data_list$loopnum,
            getsd   = getsd,
            verbose = 0))
      ))

    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        return(newmod)
      }
    }
    return(NULL)
  } # End run_one_point closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    cl <- parallel::makeCluster(min(cores, ngrid))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(Rceattle))
    })
    parallel::clusterExport(
      cl,
      varlist = ls(envir = environment()),
      envir = environment()
    )
    mod_list <- parallel::parLapply(cl, seq_len(ngrid), run_one_point)
  } else {
    mod_list <- lapply(seq_len(ngrid), run_one_point)
  }

  # NLL aligned with grid; NA for non-converged
  nll <- vapply(mod_list,
                function(x) if (is.null(x)) NA_real_ else x$opt$objective,
                numeric(1))

  names(mod_list) <- paste0("Fit_", seq_len(ngrid))

  return(list(
    Rceattle_list = mod_list,
    grid          = grid,
    nll           = nll,
    param         = param,
    slots         = slots
  ))
}

