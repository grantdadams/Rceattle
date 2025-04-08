#' Retrospective peels (DEPRECATED)
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param peels the number of retrospective peels to use in the calculation of rho and for model estimation
#' @param rescale rescale environmental predictors?
#' @param nyrs_forecast Number of forecast years to calculate Mohn's Rho in addition to terminal year
#'
#' @return a list of 1. list of Rceattle models and 2. vector of Mohn's rho for each species
#'
#' @examples
#' data(BS2017SS) # ?BS2017SS for more information on the data
#' data('BS2017MS') # Note: the only difference is the residual mortality is lower
#'
#' ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
#'                             inits = NULL, # Initial parameters = 0
#'                             file = NULL, # Don't save
#'                             debug = 0, # Estimate
#'                             random_rec = FALSE, # No random recruitment
#'                             msmMode = 0, # Single species mode
#'                             silent = TRUE)
#'
#' retro <- retrospective(ss_run, peels = 10)
#' @export
retrospective <- function(Rceattle = NULL, peels = NULL, rescale = FALSE, nyrs_forecast = 3) {
  if (class(Rceattle) != "Rceattle") {
    stop("Object is not of class 'Rceattle'")
  }

  library(dplyr)
  # Get objects
  mod_list <- list(Rceattle)
  endyr <- Rceattle$data_list$endyr
  styr <- Rceattle$data_list$styr
  projyr <- Rceattle$data_list$projyr
  nyrs_proj <- projyr - styr + 1

  #####################################
  # Run across retrospective peels
  #####################################
  ind <- 2
  for (i in 1:peels) {
    data_list <- Rceattle$data_list
    data_list$endyr <- endyr - i
    nyrs <- (data_list$endyr) - styr + 1


    # * Adjust data ----
    data_list$index_data <- data_list$index_data %>%
      dplyr::filter(Year <= data_list$endyr)

    data_list$catch_data <- data_list$catch_data %>%
      dplyr::filter(Year <= data_list$endyr)

    data_list$comp_data <- data_list$comp_data %>%
      dplyr::filter(Year <= data_list$endyr)

    data_list$weight <- data_list$weight %>%
      dplyr::filter(Year <= data_list$endyr)

    data_list$Pyrs <- data_list$Pyrs %>%
      dplyr::filter(Year <= data_list$endyr)

    data_list$diet_data <- data_list$diet_data %>%
      dplyr::filter(Year <= data_list$endyr)


    # * Rescale environmental predictors ----
    if(rescale){
      data_list$env_data[,2:ncol(data_list$env_data)]<-scale(data_list$env_data[,2:ncol(data_list$env_data)])
    }


    # * Adjust parameter size ----
    inits <- Rceattle$estimated_params

    inits$rec_dev[, (nyrs + 1):nyrs_proj] <- 0

    inits$ln_F <- inits$ln_F[, 1:nyrs]
    inits$index_q_dev <- inits$index_q_dev[,1:nyrs]
    inits$ln_sel_slp_dev <- inits$ln_sel_slp_dev[,,,1:nyrs]
    inits$sel_inf_dev <- inits$sel_inf_dev[,,,1:nyrs]
    inits$sel_coff_dev <- array(inits$sel_coff_dev[,,,1:nyrs], dim = c(dim(Rceattle$estimated_params$sel_coff_dev )[1:3], nyrs))

    # * Adjust map size ----
    map <- Rceattle$map

    map$mapList$rec_dev[, (nyrs + 1):nyrs_proj] <- NA
    map$mapFactor$rec_dev <- factor(map$mapList$rec_dev)

    map$mapList$ln_F <- map$mapList$ln_F[, 1:nyrs]
    map$mapFactor$ln_F <- factor(map$mapList$ln_F)

    map$mapList$index_q_dev <- map$mapList$index_q_dev[,1:nyrs]
    map$mapFactor$index_q_dev <- factor(map$mapList$index_q_dev)

    map$mapList$ln_sel_slp_dev <- map$mapList$ln_sel_slp_dev[,,,1:nyrs]
    map$mapFactor$ln_sel_slp_dev <- factor(map$mapList$ln_sel_slp_dev)

    map$mapList$sel_inf_dev <- map$mapList$sel_inf_dev[,,,1:nyrs]
    map$mapFactor$sel_inf_dev <- factor(map$mapList$sel_inf_dev)

    map$mapList$sel_coff_dev <- array(map$mapList$sel_coff_dev[,,,1:nyrs], dim = c(dim(Rceattle$estimated_params$sel_coff_dev )[1:3], nyrs))
    map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)


    # * Refit ----
    newmod <- suppressWarnings(
      Rceattle::fit_mod(
        data_list = data_list,
        inits = inits,
        map =  NULL,
        bounds = NULL,
        file = NULL,
        estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
        HCR = build_hcr(HCR = data_list$HCR, # Tier3 HCR
                        DynamicHCR = data_list$DynamicHCR,
                        FsprTarget = data_list$FsprTarget,
                        FsprLimit = data_list$FsprLimit,
                        Ptarget = data_list$Ptarget,
                        Plimit = data_list$Plimit,
                        Alpha = data_list$Alpha,
                        Pstar = data_list$Pstar,
                        Sigma = data_list$Sigma,
                        Fmult = data_list$Fmult,
                        HCRorder = data_list$HCRorder
        ),
        recFun = build_srr(srr_fun = data_list$srr_fun,
                           srr_pred_fun  = data_list$srr_pred_fun ,
                           proj_mean_rec  = data_list$proj_mean_rec ,
                           srr_meanyr = min(data_list$srr_meanyr, data_list$endyr), # Update end year if less than srr_meanyr
                           srr_hat_styr = data_list$srr_hat_styr,
                           srr_hat_endyr = data_list$srr_hat_endyr,
                           srr_est_mode  = data_list$srr_est_mode ,
                           srr_prior  = data_list$srr_prior,
                           srr_prior_sd   = data_list$srr_prior_sd,
                           Bmsy_lim = data_list$Bmsy_lim,
                           srr_env_indices = data_list$srr_env_indices),
        M1Fun =     build_M1(M1_model= data_list$M1_model,
                             updateM1 = FALSE,
                             M1_use_prior = data_list$M1_use_prior,
                             M2_use_prior = data_list$M2_use_prior,
                             M_prior = data_list$M_prior,
                             M_prior_sd = data_list$M_prior_sd),
        random_rec = data_list$random_rec,
        niter = data_list$niter,
        msmMode = data_list$msmMode,
        avgnMode = data_list$avgnMode,
        suitMode = data_list$suitMode,
        suit_styr = data_list$suit_styr,
        suit_endyr = min(data_list$suit_endyr, data_list$endyr),   # Update to end year if less than suit_endyr
        initMode = data_list$initMode,
        phase = TRUE,
        loopnum = data_list$loopnum,
        getsd = TRUE,
        verbose = 0)
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

    # Refit model If converged
    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        mod_list[[ind]] <- newmod
        ind <- ind + 1
      }
    }
  }

  #####################################
  # Calculate Mohs rho for each species
  #####################################
  objects <- c("biomass", "ssb", "R", "F_spp")

  mohns <- data.frame(matrix(0, nrow = length(objects) * (nyrs_forecast+1), ncol = 3 + data_list$nspp))
  colnames(mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  # Loop around peels that converged
  for (i in 1:(length(mod_list) - 1)) {
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr - styr + 1
    ind <- 1

    # Loop around derived quantitities
    for (j in 1:length(objects)) {
      term_quantities <- mod_list[[1]]$quantities[[objects[j]]]
      retro_quantities <- mod_list[[i + 1]]$quantities[[objects[j]]]

      # Loop around forecast
      for(yr in 0:nyrs_forecast){

        # If data exist for these years
        if(nyrs_peel + yr <= mod_list[[1]]$data_list$endyr - styr + 1){

          # Get full and peeled models
          base <- term_quantities[, nyrs_peel + yr]
          peel <- retro_quantities[, nyrs_peel + yr]
          rel_error <- ((peel - base)/base)

          # Save
          mohns[ind, 1] <- objects[j]
          mohns[ind, 2] <- yr
          mohns[ind, 3] <- mohns[ind, 3] + 1
          mohns[ind, 4:(data_list$nspp + 3) ] <- mohns[j, 4:(data_list$nspp + 3)] + rel_error
        }
        ind = ind+1
      }
    }
  }

  # Divide sum of RE by N
  mohns[, 4:(data_list$nspp + 3) ] <- mohns[, 4:(data_list$nspp + 3)]/mohns[, 3]


  return(list(Rceattle_list = rev(mod_list), mohns = mohns))
}
