#' Rerun with F = 0.
#'
#' @description
#' Function to update hindcast and set F to 0.
#' Useful for determining dynamic reference points for multi-species models under climate-change.
#'
#'
#' @param Rceattle A fitted Rceattle model object
#'
#' @export
#'
remove_F <- function(Rceattle){

  # A DSEM needs nothing special here. This sets projection F to zero and
  # rebuilds at estimateMode = 3, so nothing is re-estimated: the DSEM's fitted
  # parameters simply travel through `inits`. That only became true once
  # fit_mod() stopped dropping the dsem_* blocks out of a warm start -- before
  # that this rebuilt the model with the recruitment SD at its start value, and
  # the F = 0 projection it returned was not the fitted model's.

  # * Years for F = 0 ----
  # - don't want hindcast or it will bias suitability in Multi-species models
  # suit_endyr is per-predator; project F = 0 only after the latest suitability window.
  proj_years <- (max(Rceattle$data_list$suit_endyr)+1):Rceattle$data_list$projyr - Rceattle$data_list$styr + 1
  fdevs_cols <- 1:ncol(Rceattle$estimated_params$log_F)
  fdevs_change <- which(fdevs_cols %in% proj_years)

  # * Set F to 0 ----
  Rceattle$estimated_params$log_F[,fdevs_change] <- replace(Rceattle$estimated_params$log_F[,fdevs_change], values = -999)

  # * Update fit ----
  # Rebuild with F = 0 in the projection, reusing the model's own configuration;
  # clamp the SR-switch and suitability-end years to the hindcast terminal year.
  estMode <- Rceattle$data_list$estimateMode
  Rceattle <- .refit_like(
    data_list        = Rceattle$data_list,
    inits            = Rceattle$estimated_params,
    estimateMode     = 3,
    getsd            = TRUE,
    srr_mse_switchyr = min(Rceattle$data_list$srr_mse_switchyr, Rceattle$data_list$endyr),
    suit_endyr       = pmin(Rceattle$data_list$suit_endyr, Rceattle$data_list$endyr))

  Rceattle$data_list$estimateMode <- estMode
  return(Rceattle)

}

