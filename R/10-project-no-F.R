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
  # PROJECTION F IS NOT A FUNCTION OF log_F. The template reads exp(log_F(flt,
  # yr)) only for yr < nyrs_hind (ceattle.cpp section 6); projection F is
  # proj_F_prop * proj_F from the harvest control rule. What actually produces
  # the unfished projection here is estimateMode = 3, which skips fit_mod()'s
  # HCR section so `forecast` stays 0 and the template forces proj_F = 0.
  #
  # This used to also write -999 into log_F for
  #   (max(suit_endyr) + 1):projyr - styr + 1
  # which are column indices PAST the end of log_F whenever the suitability
  # window runs to endyr -- so on every model where suit_endyr == endyr it
  # selected nothing and did nothing. Where suit_endyr stops EARLIER, though, the
  # same expression lands inside the hindcast and zeroed real fitted fishing
  # mortality: measured on BS2017SS with suit_endyr = 2010, it zeroed F for
  # 2011-2017 and inflated terminal hindcast SSB by 43% / 187% / 17%. The
  # returned object is used as "the fitted hindcast with the projection
  # unfished" -- including by run_mse() to build the unfished reference every
  # performance metric is normalized against -- and the truncated suitability
  # windows that trigger it are what the real GOA and hake models use. It also
  # collapsed the per-predator suit_endyr to one scalar and wrote across every
  # fleet and species. Removed: it never did anything in the intended case, and
  # what it did in the other one contradicts this function's own comment about
  # not touching the hindcast.

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

