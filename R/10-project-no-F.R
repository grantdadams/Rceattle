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

  # A DSEM needs nothing special here: nothing is re-estimated below, so the
  # DSEM's fitted parameters simply travel through `inits`. That only became
  # true once fit_mod() stopped dropping the dsem_* blocks out of a warm start.

  # * Years for F = 0 ----
  # INTENTIONAL that this starts at suit_endyr, not endyr: MSVPA-derived
  # suitability is computed from the fitted dynamics up to suit_endyr, so the
  # unfished trajectory has to begin after that window or removing F would alter
  # the suitability the multispecies model was fitted with. The consequence is
  # that when suit_endyr < endyr the LATE HINDCAST YEARS are unfished too --
  # that is the intended unfished reference, not a hindcast reconstruction.
  # Measured on BS2017SS with suit_endyr = 2010: F zeroed for 2011-2017 and
  # terminal SSB 43%/187%/17% above the fitted model, which is the unfished
  # counterfactual those years are supposed to represent.
  #
  # suit_endyr is per-predator and max() takes the latest window, so no
  # predator's suitability is disturbed.
  #
  # Note this write only reaches the hindcast: log_F has nyrs_hind columns, and
  # the template reads log_F only for hindcast years. The unfished PROJECTION
  # comes from estimateMode = 3 below leaving the harvest control rule
  # unevaluated, so `forecast` stays 0 and the template forces proj_F = 0. When
  # suit_endyr == endyr the indices fall past the end of log_F and nothing is
  # written, which is correct -- there are no hindcast years to unfish.
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

