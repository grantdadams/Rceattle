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
#' data("BS2017MS") # Note: the only difference is the residual mortality is lower
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
retrospective <- function( Rceattle = NULL, peels = NULL){
  if(class(Rceattle) != "Rceattle"){
    stop("Object is not of class 'Rceattle'")
  }

  # Get objects
  mod_list <- list(Rceattle)
  data_list <- Rceattle$data_list
  endyr <- data_list$endyr
  styr <- data_list$styr

  # Run across retrospective bits
  ind <- 2
  for(i in 1:peels){
    data_list$endyr <- endyr - i
    nyrs <- (endyr - i) - styr + 1

    # Adjust initial parameters
    inits <- Rceattle$estimated_params
    inits$rec_dev <- inits$rec_dev[,1:nyrs]
    inits$F_dev <- inits$F_dev[,1:nyrs]

    # Refit
    newmod <- suppressMessages(suppressWarnings( Rceattle::fit_mod(
      data_list = data_list,
      inits = inits,
      file = NULL,
      debug = data_list$debug,
      niter = data_list$niter,
      random_rec = data_list$random_rec,
      msmMode = data_list$msmMode,
      suitMode = data_list$suitMode,
      avgnMode = data_list$avgnMode,
      silent = TRUE) ))

    # Refit model
    if(newmod$opt$Convergence_check != "The model is definitely not converged" | !is.null(newmod$opt$Convergence_check)){
      mod_list[[ind]] <- newmod
      ind <- ind + 1
    }
  }

  # Calculate Mohs rho for each species
  objects <- c("biomass", "biomassSSB", "R")

  mohns <- data.frame(matrix(0, nrow = length(objects), ncol = 1 + data_list$nspp))
  colnames(mohns) <- c("Object", paste0("Species_", 1:data_list$nspp))
  mohns$Object <- objects

  # Loop around peels that converged
  for(i in 1:length(mod_list)){
    nyrs_peel <- mod_list[[i]]$data_list$endyr - styr + 1

    # Loop around derived quantitities
    for(j in 1:length(objects)){
      base <- mod_list[[1]]$quantities[[objects[j]]]
      peel <- mod_list[[i + 1]]$quantities[[objects[j]]]

      base <- base[,nyrs_peel]
      peel <- peel[,nyrs_peel]

      rel_error <- ((peel - base) / base)/peels

      mohns[,j+1] <- mohns[,j+1] + rel_error
    }
  }

  # Do it for fishing mortality
  # mohnsF <- data.frame(matrix(0, nrow = 1, ncol = 1 + data_list$nspp))
  # colnames(mohnsF) <- c("Object", paste0("Species_", 1:data_list$nspp))
  # mohnsF$Object <- "F"
  #
  # fsh_control <- data_list$fsh_control
  #
  # for(i in 1:peels){
  #
  #   # Loop around derived quantitities
  #     base <- exp(mod_list[[1]]$estimated_params$ln_mean_F + mod_list[[1]]$estimated_params$F_dev)
  #     peel <- exp(mod_list[[i+1]]$estimated_params$ln_mean_F + mod_list[[i+1]]$estimated_params$F_dev)
  #
  #     base <- base[,ncol(peel)]
  #     peel <- peel[,ncol(peel)]
  #
  #     rel_error <- ((peel - base) / base)/peels
  #   }
  # }

  return(list(Rceattle_list = rev(mod_list), mohns = mohns))
}
