#' Rerun with F = 0.
#'
#' @description
#' Refits the model with fishing mortality set to 0 from `start_yr` on, keeping
#' every other parameter. The projection after `endyr` is always unfished.
#' `run_mse()` uses it for the no-fishing run (`OM_no_F`) behind the collapse
#' metrics.
#'
#' @param object A fitted Rceattle model object
#' @param start_yr First year with F = 0, from `styr` to the year after `endyr` (the default, which leaves the hindcast unchanged); under predation it must fall after the window of any predator with empirical suitability (`suitMode = 0`).
#' @param Rceattle deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @export
#'
remove_F <- function(object = NULL, start_yr = NULL, Rceattle = NULL){
  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "remove_F")

  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).", call. = FALSE)
  }

  dl <- object$data_list
  if (is.null(start_yr)) start_yr <- dl$endyr + 1
  # The projection is always unfished, so the no-F period starts by endyr + 1.
  if (!is.numeric(start_yr) || length(start_yr) != 1 || is.na(start_yr) ||
      start_yr != round(start_yr) || start_yr < dl$styr || start_yr > dl$endyr + 1) {
    stop("`start_yr` must be a single year from the first model year (", dl$styr,
         ") to the year after endyr (", dl$endyr + 1, "); the projection is always unfished.",
         call. = FALSE)
  }

  # Empirical suitability (suitMode 0) reads abundance to suit_endyr; removing F there changes it.
  if (isTRUE(.map_switch(dl$msmMode, msmMode_map, "msmMode") > 0)) {
    suit_mode <- rep_len(.map_switch(dl$suitMode, suitMode_map, "suitMode"), dl$nspp)
    suit_end  <- rep_len(pmin(dl$suit_endyr, dl$endyr), dl$nspp)[suit_mode == 0]
    if (length(suit_end) && start_yr <= max(suit_end)) {
      stop("`start_yr` (", start_yr, ") must be after the empirical suitability window ",
           "(suit_endyr ", max(suit_end), "): removing fishing inside it would re-derive ",
           "the predation suitability the model was fit with.", call. = FALSE)
    }
  }

  # * Years for F = 0 ----
  proj_years <- start_yr:dl$projyr - dl$styr + 1
  fdevs_cols <- 1:ncol(object$estimated_params$log_F)
  fdevs_change <- which(fdevs_cols %in% proj_years)

  # * Set F to 0 ----
  object$estimated_params$log_F[,fdevs_change] <- replace(object$estimated_params$log_F[,fdevs_change], values = -999)

  # * Update fit ----
  # Build-only refit (projection unfished); SR-switch and suitability years clamp to endyr.
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

