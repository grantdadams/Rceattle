#' Bundle the optimizer / sdreport / phasing controls for `fit_mod()`
#'
#' `fit_mod()` carries roughly a dozen optimizer- and reporting-related
#' arguments (`bias.correct`, `getsd`, `loopnum`, `newtonsteps`, ...).
#' That is a lot of surface area when the user mostly cares about
#' "what model am I fitting" rather than "how is it being fit."
#' `fit_control()` collects those knobs into a single object so calls
#' to `fit_mod()` can stay focused on the model spec:
#'
#' ```r
#' fit <- fit_mod(
#'   data_list   = BS2017SS,
#'   msmMode     = 0,
#'   fit_control = fit_control(loopnum = 1, getsd = FALSE)
#' )
#' ```
#'
#' Pass the result via the `fit_control` argument to [fit_mod()]. When
#' supplied, the values in the `fit_control` object override the
#' corresponding individual arguments to `fit_mod()`. Individual
#' arguments are kept for backward compatibility.
#'
#' @param bias.correct logical. If `TRUE`, applies bias correction via
#'   [TMB::sdreport()]. Default `FALSE`.
#' @param getsd logical. If `TRUE`, run [TMB::sdreport()] after
#'   optimization. Default `TRUE`.
#' @param getJointPrecision logical. Return the full Hessian of fixed
#'   and random effects. Default `TRUE` (matches `fit_mod()` default).
#' @param getReportCovariance logical. Return the variance-covariance
#'   of `ADREPORT` variables. Default `FALSE`.
#' @param projection_uncertainty logical. If `TRUE`, accounts for hindcast
#'   parameter uncertainty in projections when using an HCR (refits with all
#'   hindcast and biological-reference-point parameters turned on). Default
#'   `FALSE` for speed.
#' @param osa logical. If `TRUE`, assemble the full one-step-ahead (OSA)
#'   observation data during the fit so [osa_residuals()] can be computed for
#'   composition / caal / diet as well as the index/catch series. Default `FALSE`:
#'   the fitted objective is identical either way, but skipping the (sizeable)
#'   composition OSA metadata keeps fits fast -- important for
#'   simulation testing such as [run_mse()], [self_test()], [jitter()], and
#'   [retrospective()]. Set `TRUE` when you intend to call [osa_residuals()] on
#'   composition data.
#' @param comp_offset Numeric or `NULL`. Small proportion offset added to the
#'   observed and predicted age/length composition and conditional-age-at-length
#'   (caal) bins before the multinomial / Dirichlet-multinomial likelihood (to
#'   avoid `log(0)` for empty bins). Stored on `data_list` so the fitted
#'   likelihood and the OSA observation vector use the same offset, and so
#'   internal re-fits inherit it. `NULL` (the default) inherits
#'   `data_list$comp_offset` if set, otherwise `1e-5`; a number overrides it.
#'   Does not apply to diet (stomach-content) compositions.
#' @param bias_adjust_obs logical with default TRUE. Whether to
#'   apply a bias adjustment (mean-sd^2/2) to lognormal data
#'   likelihoods
#' @param bias_adjust_proc logical with default TRUE. Whether to
#'   apply a bias adjustment (mean-sd^2/2) to lognormal process
#'   likelihoods
#' @param use_gradient logical. Use the analytic gradient during
#'   phasing. Default `TRUE`.
#' @param rel_tol Numeric tolerance used to flag discontinuous
#'   likelihood warnings (compares the TMB and `nlminb` objectives).
#'   Default `1`.
#' @param loopnum Integer. Number of times to re-start optimization
#'   (`loopnum = 3` sometimes achieves a lower final gradient than
#'   `loopnum = 1`). Default `5`.
#' @param newtonsteps Integer. Number of extra Newton steps to take
#'   after optimization (alternative to `loopnum`). Default `0`.
#' @param phase `TRUE`/`FALSE` or a list. If `FALSE`, the model is not
#'   phased. If `TRUE`, default phasing is used. Can also accept a
#'   list of parameter object names with corresponding phase. Default
#'   `FALSE`.
#' @param TMBfilename Optional character. Path (without `.cpp`) to an
#'   alternate TMB template for development. Default `NULL` (use the
#'   bundled `ceattle_v01_11`).
#' @param verbose `0` = silent, `1` = print updates of model fit, `2` =
#'   print updates of model fit and TMB estimation progress. Default
#'   `1`.
#' @param nlminb_control A list of control parameters passed to
#'   [stats::nlminb()]. See `?nlminb`. Default
#'   `list(eval.max = 1e9, iter.max = 1e9, trace = 0)`.
#'
#' @return A list of class `"Rceattle_fit_control"`.
#' @export
#'
#' @examples
#' # Quick-and-dirty fit: skip sdreport, single optimizer pass
#' ctl <- fit_control(getsd = FALSE, loopnum = 1)
#'
#' # Production fit with bias correction and joint precision
#' ctl <- fit_control(bias.correct = TRUE, getJointPrecision = TRUE)
fit_control <- function(
  bias.correct        = FALSE,
  getsd               = TRUE,
  getJointPrecision   = TRUE,
  getReportCovariance = FALSE,
  projection_uncertainty = FALSE,
  osa                 = FALSE,
  comp_offset         = NULL,
  bias_adjust_obs     = TRUE,
  bias_adjust_proc    = TRUE,
  use_gradient        = TRUE,
  rel_tol             = 1,
  loopnum             = 5,
  newtonsteps         = 0,
  phase               = FALSE,
  TMBfilename         = NULL,
  verbose             = 1,
  nlminb_control      = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)
) {
  ans <- list(
    bias.correct        = bias.correct,
    getsd               = getsd,
    getJointPrecision   = getJointPrecision,
    getReportCovariance = getReportCovariance,
    projection_uncertainty = projection_uncertainty,
    osa                 = osa,
    comp_offset         = comp_offset,
    bias_adjust_obs     = bias_adjust_obs,
    bias_adjust_proc    = bias_adjust_proc,
    use_gradient        = use_gradient,
    rel_tol             = rel_tol,
    loopnum             = loopnum,
    newtonsteps         = newtonsteps,
    phase               = phase,
    TMBfilename         = TMBfilename,
    verbose             = verbose,
    nlminb_control      = nlminb_control
  )
  class(ans) <- c("Rceattle_fit_control", "list")
  ans
}
