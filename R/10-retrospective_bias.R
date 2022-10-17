#' Retrospective peels (DEPRECATED)
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param peels the number of retrospective peels to use in the calculation of rho and for model estimation
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
retrospective <- function(Rceattle = NULL, peels = NULL) {
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
    nyrs <- (endyr - i) - styr + 1

    # Adjust data
    data_list$srv_biom <- data_list$srv_biom %>%
      filter(Year <= data_list$endyr)

    data_list$wt <- data_list$wt %>%
      filter(Year <= data_list$endyr)

    data_list$comp_data <- data_list$comp_data %>%
      filter(Year <= data_list$endyr)

    data_list$fsh_biom <- data_list$fsh_biom %>%
      filter(Year <= data_list$endyr)

    data_list$Pyrs <- data_list$Pyrs %>%
      filter(Year <= data_list$endyr)

    # Adjust initial parameters
    inits <- Rceattle$estimated_params
    inits$rec_dev[, (nyrs + 1):nyrs_proj] <- 0
    inits$F_dev <- inits$F_dev[, 1:nyrs]

    # Adjust parameter size
    inits$ln_srv_q_dev <- inits$ln_srv_q_dev[,1:nyrs]
    inits$ln_sel_slp_dev <- inits$ln_sel_slp_dev[,,,1:nyrs]
    inits$sel_inf_dev <- inits$sel_inf_dev[,,,1:nyrs]
    inits$sel_coff_dev <- array(inits$sel_coff_dev[,,,1:nyrs], dim = c(dim(Rceattle$estimated_params$sel_coff_dev )[1:3], nyrs))

    # Refit
    newmod <- suppressWarnings(
      Rceattle::fit_mod(
        data_list = data_list,
        inits = NULL,
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
                        Sigma = data_list$Sigma
        ),
        random_rec = data_list$random_rec,
        niter = data_list$niter,
        msmMode = data_list$msmMode,
        avgnMode = data_list$avgnMode,
        minNByage = data_list$minNByage,
        suitMode = data_list$suitMode,
        phase = "default",
        meanyr = data_list$endyr, # Update end year
        updateM1 = FALSE,
        getsd = FALSE,
        verbose = 0)

    )

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
  objects <- c("biomass", "biomassSSB", "R", "F_spp")

  mohns <- data.frame(matrix(0, nrow = length(objects) + 1, ncol = 1 + data_list$nspp))
  colnames(mohns) <- c("Object", paste0("Spp/Fsh_", 1:max(c(data_list$nspp, nrow(data_list$fsh_control)))))
  mohns$Object <- c(objects, "F")

  # Loop around peels that converged
  for (i in 1:(length(mod_list) - 1)) {
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr - styr + 1

    # Loop around derived quantitities
    for (j in 1:length(objects)) {
      base <- mod_list[[1]]$quantities[[objects[j]]]
      peel <- mod_list[[i + 1]]$quantities[[objects[j]]]

      base <- base[, nyrs_peel]
      peel <- peel[, nyrs_peel]

      rel_error <- ((peel - base)/base)/peels

      mohns[j, 2:(data_list$nspp + 1) ] <- mohns[j, 2:(data_list$nspp + 1)] + rel_error
    }
  }

  return(list(Rceattle_list = rev(mod_list), mohns = mohns))
}
