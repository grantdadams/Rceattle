#' Retrospective peels
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

  # Get objects
  mod_list <- list(Rceattle)
  data_list <- Rceattle$data_list
  endyr <- data_list$endyr
  styr <- data_list$styr
  projyr <- data_list$projyr
  nyrs_proj <- projyr - styr + 1

  # Run across retrospective bits
  ind <- 2
  for (i in 1:peels) {
    data_list$endyr <- endyr - i
    nyrs <- (endyr - i) - styr + 1

    # Adjust initial parameters
    inits <- Rceattle$estimated_params
    inits$rec_dev[, (nyrs + 1):nyrs_proj] <- 0
    inits$F_dev <- inits$F_dev[, 1:nyrs]

    inits$ln_srv_q_dev <- inits$ln_srv_q_dev[, 1:nyrs]
    inits$ln_srv_q_dev_re <- inits$ln_srv_q_dev_re[, 1:nyrs]

    inits$srv_sel_inf_dev <- inits$srv_sel_inf_dev[, ,1:nyrs]
    inits$srv_sel_inf_dev_re <- inits$srv_sel_inf_dev_re[, ,1:nyrs]
    inits$srv_sel_slp_dev <- inits$srv_sel_slp_dev[, ,1:nyrs]
    inits$srv_sel_slp_dev_re <- inits$srv_sel_slp_dev_re[, ,1:nyrs]

    inits$fsh_sel_inf_dev <- inits$fsh_sel_inf_dev[, ,1:nyrs]
    inits$fsh_sel_inf_dev_re <- inits$fsh_sel_inf_dev_re[, ,1:nyrs]
    inits$fsh_sel_slp_dev <- inits$fsh_sel_slp_dev[, ,1:nyrs]
    inits$fsh_sel_slp_dev_re <- inits$fsh_sel_slp_dev_re[, ,1:nyrs]

    # Adjust map parameters
    map <- Rceattle$map
    map[[2]]$rec_dev[, (nyrs + 1):nyrs_proj] <- 0
    map[[2]]$F_dev <- map[[2]]$F_dev[, 1:nyrs]

    map[[2]]$ln_srv_q_dev <- map[[2]]$ln_srv_q_dev[, 1:nyrs]
    map[[2]]$ln_srv_q_dev_re <- map[[2]]$ln_srv_q_dev_re[, 1:nyrs]

    map[[2]]$srv_sel_inf_dev <- map[[2]]$srv_sel_inf_dev[, ,1:nyrs]
    map[[2]]$srv_sel_inf_dev_re <- map[[2]]$srv_sel_inf_dev_re[, ,1:nyrs]
    map[[2]]$srv_sel_slp_dev <- map[[2]]$srv_sel_slp_dev[, ,1:nyrs]
    map[[2]]$srv_sel_slp_dev_re <- map[[2]]$srv_sel_slp_dev_re[, ,1:nyrs]

    map[[2]]$fsh_sel_inf_dev <- map[[2]]$fsh_sel_inf_dev[, ,1:nyrs]
    map[[2]]$fsh_sel_inf_dev_re <- map[[2]]$fsh_sel_inf_dev_re[, ,1:nyrs]
    map[[2]]$fsh_sel_slp_dev <- map[[2]]$fsh_sel_slp_dev[, ,1:nyrs]
    map[[2]]$fsh_sel_slp_dev_re <- map[[2]]$fsh_sel_slp_dev_re[, ,1:nyrs]

    for (i in 1:length(map[[2]])) {
      map[[1]][[i]] <- factor(map[[2]][[i]])
    }

    # Refit
    newmod <- suppressMessages(suppressWarnings(Rceattle::fit_mod(data_list = data_list, inits = inits, file = NULL, debug = 0, map = map,
                                                                  niter = data_list$niter, random_rec = data_list$random_rec, msmMode = data_list$msmMode, suitMode = data_list$suitMode,
                                                                  avgnMode = data_list$avgnMode, silent = TRUE)))

    # Refit model If converged
    if (!is.null(newmod$opt$opt$Convergence_check)) {
      if (newmod$opt$opt$Convergence_check != "The model is definitely not converged") {
        mod_list[[ind]] <- newmod
        ind <- ind + 1
      }
    }
  }

    # Calculate Mohs rho for each species
    objects <- c("biomass", "biomassSSB", "R")

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

        # F vec
        base <- exp(mod_list[[1]]$estimated_params$F_dev + mod_list[[1]]$estimated_params$ln_mean_F)
        peel <- exp(mod_list[[1 + i]]$estimated_params$F_dev + mod_list[[1 + i]]$estimated_params$ln_mean_F)

        base <- base[, nyrs_peel]
        peel <- peel[, nyrs_peel]

        rel_error <- ((peel - base)/base)/peels

        mohns[nrow(mohns), 2:ncol(mohns) ] <- mohns[nrow(mohns), 2:ncol(mohns) ] + rel_error
    }

  return(list(Rceattle_list = rev(mod_list), mohns = mohns))
}
