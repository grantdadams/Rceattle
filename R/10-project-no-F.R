#' Rerun with F = 0.
#'
#' @description
#' Function to update hindcast and set F to 0.
#' Useful for determining dynamic reference points for multi-species models under climate-change.
#'
#'
#' @param object A fitted Rceattle model object
#' @param Rceattle deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @export
#'
remove_F <- function(object = NULL, Rceattle = NULL){
  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "remove_F")

  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).", call. = FALSE)
  }

  # * Years for F = 0 ----
  # - don't want hindcast or it will bias suitability in Multi-species models
  # suit_endyr is per-predator; project F = 0 only after the latest suitability window.
  proj_years <- (max(object$data_list$suit_endyr)+1):object$data_list$projyr - object$data_list$styr + 1
  fdevs_cols <- 1:ncol(object$estimated_params$log_F)
  fdevs_change <- which(fdevs_cols %in% proj_years)

  # * Set F to 0 ----
  object$estimated_params$log_F[,fdevs_change] <- replace(object$estimated_params$log_F[,fdevs_change], values = -999)

  # * Update fit ----
  # Rebuild with F = 0 in the projection, reusing the model's own configuration;
  # clamp the SR-switch and suitability-end years to the hindcast terminal year.
  estMode <- object$data_list$estimateMode
  object <- .refit_like(
    data_list        = object$data_list,
    inits            = object$estimated_params,
    estimateMode     = 3,
    getsd            = TRUE,
    srr_mse_switchyr = min(object$data_list$srr_mse_switchyr, object$data_list$endyr),
    suit_endyr       = pmin(object$data_list$suit_endyr, object$data_list$endyr))

  object$data_list$estimateMode <- estMode
  return(object)

}

