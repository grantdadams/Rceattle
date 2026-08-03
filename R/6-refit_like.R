#' Re-fit a model reusing a source `data_list`'s structural configuration
#'
#' @description
#' Internal engine behind the repeated [fit_mod()] re-invocations in
#' [retrospective()], [jitter()], [self_test()], [profile.Rceattle()],
#' [run_mse()], and [remove_F()]. Each of those diagnostics re-fits a model that
#' shares a source model's structure but changes a few things (a data peel, a
#' jittered start, a projection year). They all rebuilt the harvest-control-rule,
#' stock-recruit, natural-mortality, and growth specifications from the source
#' `data_list` by hand -- the same ~60-line block copied at eleven call sites.
#'
#' This helper reconstructs those four `build_*()` specifications and the shared
#' process switches (`msmMode`, `suitMode`, `random_rec`, ...) from `data_list`
#' (the "...like" part), then calls [fit_mod()]. Every field that genuinely
#' varies across the callers is a named argument that defaults to the
#' `data_list` value, so each call site shows only the handful of things it
#' actually changes.
#'
#' Console suppression is intentionally left to the caller: the call sites wrap
#' this differently (`suppressMessages()`, `suppressWarnings()`, both, or
#' neither), and suppression is console-only -- it does not affect the fit.
#'
#' @param data_list Source `data_list`. Drives both the data to fit and the
#'   reconstructed HCR / recFun / M1Fun / growthFun (unless a field is overridden
#'   below).
#' @param inits Starting parameters (an `estimated_params`-shaped list).
#' @param estimateMode [fit_mod()] estimate mode for this refit. Required and
#'   passed explicitly by every caller (retro/jitter refit at 0, profile at 1,
#'   the MSE OM rebuild at 3, the input-F EM at 2, ...), so it has no default.
#' @param map Parameter map, or `NULL` to let [fit_mod()] build it.
#' @param HCR A pre-built [build_hcr()] object to use instead of reconstructing
#'   one from `data_list` (used for [run_mse()]'s input-F EM refit, whose HCR is
#'   a freshly computed average F rather than the stored rule).
#' @param phase,getsd,loopnum [fit_control()] knobs. `getsd` defaults to `FALSE`
#'   (diagnostics rarely need `sdreport`); `remove_F()` overrides it to `TRUE`.
#' @param proj_mean_rec Recruitment projection mode; defaults to the source's.
#'   [run_mse()]'s OM refit pins it to `TRUE` (project on the mean across sims).
#' @param srr_mse_switchyr,srr_hat_styr,srr_hat_endyr Stock-recruit reference
#'   fields; default to the source's values. Callers clamp `srr_mse_switchyr` /
#'   `srr_hat_endyr` to a peel or assessment year, or pin them to the pristine OM.
#' @param suit_styr,suit_endyr Suitability window; default to the source's,
#'   clamped to a peel year or pinned to the pristine OM by some callers.
#'
#' @return The fitted `Rceattle` object returned by [fit_mod()].
#' @noRd
.refit_like <- function(data_list,
                        inits,
                        estimateMode,
                        map              = NULL,
                        HCR              = NULL,
                        phase            = FALSE,
                        getsd            = FALSE,
                        loopnum          = data_list$loopnum,
                        proj_mean_rec    = data_list$proj_mean_rec,
                        srr_mse_switchyr = data_list$srr_mse_switchyr,
                        srr_hat_styr     = data_list$srr_hat_styr,
                        srr_hat_endyr    = data_list$srr_hat_endyr,
                        suit_styr        = data_list$suit_styr,
                        suit_endyr       = data_list$suit_endyr) {
  dl <- data_list

  # Reconstruct the harvest control rule from the source unless the caller
  # supplied one (the input-F EM refit builds its own from a computed average F).
  if (is.null(HCR)) {
    HCR <- build_hcr(
      HCR        = dl$HCR,
      DynamicHCR = dl$DynamicHCR,
      Ftarget    = dl$Ftarget,
      Flimit     = dl$Flimit,
      Ptarget    = dl$Ptarget,
      Plimit     = dl$Plimit,
      Alpha      = dl$Alpha,
      Pstar      = dl$Pstar,
      Sigma      = dl$Sigma,
      Fmult      = dl$Fmult,
      HCRorder   = dl$HCRorder)
  }

  fit_mod(
    data_list    = dl,
    inits        = inits,
    map          = map,
    bounds       = NULL,
    file         = NULL,
    estimateMode = estimateMode,
    HCR          = HCR,
    # suppressWarnings: legacy srr_fun = 1|3|5 / srr_indices.
    recFun = suppressWarnings(build_srr(
      srr_fun          = dl$srr_fun,
      srr_pred_fun     = dl$srr_pred_fun,
      proj_mean_rec    = proj_mean_rec,
      srr_mse_switchyr = srr_mse_switchyr,
      srr_hat_styr     = srr_hat_styr,
      srr_hat_endyr    = srr_hat_endyr,
      srr_est_mode     = dl$srr_est_mode,
      srr_prior        = dl$srr_prior,
      srr_prior_sd     = dl$srr_prior_sd,
      Bmsy_lim         = dl$Bmsy_lim,
      srr_indices      = dl$srr_indices,
      linkages         = dl$srr_linkages)),
    # suppressWarnings: legacy M1_indices may travel via data_list.
    M1Fun = suppressWarnings(build_M1(
      M1_model     = dl$M1_model,
      M1_re        = dl$M1_re,
      updateM1     = FALSE,
      M1_use_prior = dl$M1_use_prior,
      M2_use_prior = dl$M2_use_prior,
      M_prior      = dl$M_prior,
      M_prior_sd   = dl$M_prior_sd,
      M1_indices   = dl$M1_indices,
      linkages     = dl$M1_linkages)),
    growthFun = build_growth(fun      = dl$growth_fun,
                             linkages = dl$growth_linkages),
    # Reconstruct the catchability / selectivity / composition-weight linkages
    # so their beta_linkage parameters survive the refit. Without this the
    # warm-start `inits` carry linkage betas the rebuilt model lacks, and
    # MakeADFun segfaults on the length mismatch. Each build_*(NULL) is a no-op
    # for models that do not use that linkage.
    qFun       = build_catchability(linkages = dl$q_linkages),
    selFun     = build_selectivity(linkages = dl$sel_linkages),
    compFun    = build_composition(linkages = dl$comp_linkages),
    random_rec = dl$random_rec,
    random_q   = isTRUE(as.logical(dl$random_q)),
    random_sel = isTRUE(as.logical(dl$random_sel)),
    niter      = dl$niter,
    msmMode    = dl$msmMode,
    avgnMode   = dl$avgnMode,
    suitMode   = dl$suitMode,
    suit_styr  = suit_styr,
    suit_endyr = suit_endyr,
    initMode   = dl$initMode,
    fit_control = fit_control(
      phase   = phase,
      loopnum = loopnum,
      getsd   = getsd,
      verbose = 0))
}
