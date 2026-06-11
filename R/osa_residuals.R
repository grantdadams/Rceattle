#' One-step-ahead (OSA) residuals for an Rceattle model
#'
#' @description
#' Computes one-step-ahead (OSA) residuals -- also called forecast or quantile
#' residuals (Thygesen et al. 2017) -- for a fitted [Rceattle] model via
#' [TMB::oneStepPredict()]. Unlike Pearson residuals, OSA residuals are
#' distributed iid standard normal under a correctly specified model even when
#' observations are correlated (through composition bins) or when the model
#' contains random effects, so they support objective goodness-of-fit testing
#' (Trijoulet et al. 2023; Stewart and Monnahan 2025).
#'
#' These are *internal* OSA residuals: the residualization is integrated into
#' the assessment via TMB, so it also accounts for correlation induced by the
#' model's random effects across years -- the gold standard relative to the
#' *external* `compResidual` approach (Stewart and Monnahan 2025).
#'
#' OSA residuals are computed *post hoc* and are expensive (TMB re-optimizes the
#' random effects as each observation is added), so they are not produced during
#' [fit_mod()]. The model must have been optimized with `estimateMode < 3`.
#'
#' Supported observation types are the aggregate `"catch"` and `"index"` series,
#' the `"comp"` (age/length composition) and `"caal"` (conditional age-at-length)
#' compositions, and `"diet"` (predator stomach-content composition, for
#' multispecies models with estimated suitability). Diet is opt-in (not in the
#' default `types`) because it applies only to multispecies models and can be
#' expensive.
#'
#' For composition data the multivariate multinomial / Dirichlet-multinomial is
#' decomposed into a sequence of univariate conditional residuals (binomial /
#' beta-binomial; Trijoulet et al. 2023). The final bin of each composition is
#' fixed by the sum-to-N constraint and so has no residual (returned as `NA`).
#' Composition OSA uses an internal model rebuild with unweighted, proper
#' densities (the `osa_mode` switch); fleets fit with the AFSC `MultinomialAFSC`
#' pseudo-likelihood are residualized under the full multinomial.
#'
#' @param fit A fitted object of class `Rceattle` (from [fit_mod()]).
#' @param source Character vector of observation sources to residualize: any of
#'   `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, or `"all"`. Defaults to
#'   the four non-diet sources (`diet` is opt-in because it applies only to
#'   multispecies models and can be expensive); pass `"all"` to include `diet`.
#'   Sources with no observations in the model are silently skipped. Mirrors the
#'   `source` argument of [residuals.Rceattle()] and [plot.rceattle_osa()].
#' @param method Passed to [TMB::oneStepPredict()]. Defaults to
#'   `"oneStepGaussianOffMode"` (the WHAM/SAM default), appropriate for the
#'   lognormal aggregate series.
#' @param discrete Logical, passed to [TMB::oneStepPredict()]. Default `FALSE`.
#'   Discrete-distribution OSA residuals are randomized quantile residuals
#'   (Dunn and Smyth 1996) and so are stochastic; set `seed` for reproducibility.
#' @param seed Random seed passed to [TMB::oneStepPredict()] for reproducibility
#'   of randomized-quantile residuals. Default `123`.
#' @param trace Logical; print [TMB::oneStepPredict()] progress. Default `FALSE`.
#' @param ... Further arguments passed to [TMB::oneStepPredict()].
#'
#' @return A data frame of class `rceattle_osa` with one row per residualized
#'   observation and columns `type`, `fleet`, `fleet_name`, `species`, `sex`,
#'   `year`, `age_or_length`, `length` (the conditioning length bin for caal;
#'   `NA` otherwise), `index_label` (`"age"`/`"length"`/`NA`), `observed`,
#'   `predicted`, `sd`, and `residual`. For aggregate series `observed` and
#'   `predicted` are on the model (log) scale; for compositions they are bin
#'   counts. Carries `method` and `seed` attributes, and (when composition types
#'   are present) a `"pearson"` attribute holding the matching Pearson residuals
#'   so [plot.rceattle_osa()] can show both. Summarize it with
#'   [osa_diagnostics()] and plot it with [plot.rceattle_osa()].
#'
#' @references
#' Thygesen, U.H., et al. 2017. Validation of ecological state space models using
#'   the Laplace approximation. Environ. Ecol. Stat. 24:317-339.
#'
#' Trijoulet, V., et al. 2023. Model validation for compositional data in stock
#'   assessment models. Fish. Res. 257:106487.
#'
#' Stewart, I.J., and Monnahan, C.C. 2025. Diagnosing common sources of lack of
#'   fit to composition data using one-step-ahead residuals. Can. J. Fish. Aquat.
#'   Sci. 82:1-13.
#'
#' @seealso [osa_diagnostics()], [plot.rceattle_osa()], [process_residuals()]
#' @export
osa_residuals <- function(fit,
                          source  = c("index", "catch", "comp", "caal"),
                          method  = "oneStepGaussianOffMode",
                          discrete = FALSE,
                          seed    = 123,
                          trace   = FALSE,
                          ...) {

  # ---- Validate the fit ----
  if (!inherits(fit, "Rceattle")) {
    stop("'fit' must be a fitted Rceattle model (from fit_mod()).")
  }
  if (is.null(fit$obj)) {
    stop("'fit' has no TMB object ($obj); OSA residuals require the fitted ",
         "model object.")
  }
  if (is.null(fit$obs_ctl)) {
    stop("'fit' has no OSA observation map ($obs_ctl). Re-fit with a current ",
         "version of Rceattle that builds the OSA observation vector.")
  }
  em <- fit$data_list$estimateMode
  if (!is.null(em) && em >= 3) {
    stop("OSA residuals require a model optimized with estimateMode < 3 ",
         "(the returned objective for estimateMode >= 3 is a debug placeholder ",
         "that oneStepPredict cannot differentiate).")
  }
  if (is.null(fit$obj$env$last.par.best)) {
    stop("'fit' does not appear to have been optimized; OSA residuals require ",
         "a converged fit.")
  }

  # "diet" is supported but opt-in: it applies only to multispecies models with
  # estimated suitability and can be expensive, so it is not in the default set.
  # "all" is a synonym for every source including diet.
  valid_sources <- c("index", "catch", "comp", "caal", "diet", "all")
  source <- match.arg(source, choices = valid_sources, several.ok = TRUE)
  if ("all" %in% source) source <- c("index", "catch", "comp", "caal", "diet")

  # Composition OSA data (comp / caal / diet) is only assembled when the model
  # was fit with fit_control(osa = TRUE); the default fast path builds just the
  # aggregate index/catch entries (see build_osa_data()). Give a clear message
  # rather than silently returning fewer residuals than requested.
  comp_sources <- intersect(source, c("comp", "caal", "diet"))
  if (!isTRUE(fit$osa) && length(comp_sources) > 0) {
    msg <- paste0(
      "OSA residuals for source(s) ", paste(comp_sources, collapse = ", "),
      " require fitting with fit_control(osa = TRUE); this model was fit with ",
      "osa = FALSE (the fast default used for simulation testing). Refit with ",
      "fit_mod(..., fit_control = fit_control(osa = TRUE)) to enable them.")
    if (all(source %in% comp_sources)) stop(msg) else warning(msg)
  }

  # ---- Select the obsvec positions to residualize ----
  obs_ctl <- fit$obs_ctl
  sel <- obs_ctl[obs_ctl$type %in% source & !obs_ctl$is_last_bin, , drop = FALSE]
  if (nrow(sel) == 0) {
    stop("No observations of source(s) ", paste(source, collapse = ", "),
         " are available for OSA residuals in this model.")
  }

  # Order observations chronologically (then by source and bin) so the
  # conditional 'one-step-ahead' sequence is in time order, as recommended by
  # Stewart and Monnahan (2025). The obsvec storage order is independent of this;
  # subset defines the conditioning order.
  sel <- sel[order(sel$year, match(sel$type, source), sel$bin_index,
                   na.last = TRUE), , drop = FALSE]

  # oneStepPredict()'s 'subset' uses 1-based R indices into obsvec; obs_pos is
  # the 0-based TMB position.
  subset_pos <- sel$obs_pos + 1L

  # Rebuild the model in OSA mode (osa_mode = 1) at the fitted parameters. This
  # is required so the composition (comp/caal) likelihoods read their counts
  # from obsvec and use the unweighted, proper density that oneStepPredict needs.
  # It leaves the aggregate catch/index likelihood unchanged, so the same
  # rebuilt object serves all observation types.
  obj_osa <- .osa_build_obj(fit)

  # ---- Compute the OSA residuals ----
  osa <- TMB::oneStepPredict(
    obj                 = obj_osa,
    observation.name    = "obsvec",
    data.term.indicator = "keep",
    method              = method,
    subset              = subset_pos,
    discrete            = discrete,
    seed                = seed,
    trace               = trace,
    ...)

  # oneStepPredict returns rows in 'subset' order (see SAM::residuals.sam), so
  # they align with the rows of 'sel'. The 'residual' column is always present;
  # 'observation'/'mean'/'sd' depend on the method, so pull them defensively.
  get_col <- function(df, nm) {
    if (!is.null(df[[nm]])) df[[nm]] else rep(NA_real_, nrow(df))
  }
  index_label <- c("age", "length")[sel$comp_type + 1L]   # NA for aggregates
  fc <- fit$data_list$fleet_control                        # fleet code -> name
  fleet_name <- if (!is.null(fc)) {
    fc$Fleet_name[match(sel$fleet_code, fc$Fleet_code)]
  } else NA_character_

  out <- data.frame(
    type          = sel$type,
    fleet         = sel$fleet_code,
    fleet_name    = fleet_name,
    species       = sel$species,
    sex           = sel$sex,
    year          = sel$year,
    age_or_length = sel$age_or_length,
    length        = sel$length,
    index_label   = index_label,
    observed      = get_col(osa, "observation"),
    predicted     = get_col(osa, "mean"),
    sd            = get_col(osa, "sd"),
    residual      = osa$residual,
    stringsAsFactors = FALSE)

  rownames(out) <- NULL
  class(out) <- c("rceattle_osa", "data.frame")
  attr(out, "method") <- method
  attr(out, "seed")   <- seed

  # Attach the matching Pearson residuals for composition sources so the
  # plot() method can show OSA and Pearson bubbles side by side.
  comp_types <- intersect(unique(out$type), c("comp", "caal"))
  if (length(comp_types) > 0) {
    attr(out, "pearson") <- tryCatch(
      stats::residuals(fit, type = "pearson", source = comp_types),
      error = function(e) NULL)
  }

  n_bad <- sum(!is.finite(out$residual))
  if (n_bad > 0) {
    warning(n_bad, " of ", nrow(out), " OSA residual(s) are non-finite. This ",
            "usually indicates a poorly converged fit or very sparse / ",
            "degenerate compositions (common for conditional age-at-length); ",
            "check model convergence before interpreting the residuals.")
  }
  out
}


#' Rebuild a fitted Rceattle TMB object in OSA mode at the fitted parameters
#'
#' @description
#' Returns a TMB ADFun object equivalent to `fit$obj` but with `osa_mode = 1`,
#' built at the fitted parameter values and the same map / random-effect
#' structure, ready for [TMB::oneStepPredict()]. In OSA mode the composition
#' likelihoods read their counts from `obsvec` and use unweighted proper
#' densities; the aggregate catch/index likelihood is identical in both modes,
#' so this single object serves every observation type.
#'
#' @param fit A fitted `Rceattle` object.
#' @return A TMB ADFun object with `osa_mode = 1`.
#' @keywords internal
.osa_build_obj <- function(fit) {
  obj <- fit$obj
  if (!is.null(obj$env$data$osa_mode) && obj$env$data$osa_mode == 1L) {
    return(obj)
  }
  data2 <- obj$env$data
  # obj$env$data is already sanitized (stored as double with a 'check.passed'
  # attribute). Overwrite osa_mode as a double (DATA_INTEGER reads it via
  # CppAD::Integer) and drop 'check.passed' so MakeADFun re-sanitizes cleanly.
  data2$osa_mode <- 1
  attr(data2, "check.passed") <- NULL
  random_names <- if (length(obj$env$random)) {
    unique(names(obj$env$par)[obj$env$random])
  } else {
    NULL
  }
  # parList() returns the current (fitted) parameters as a list structured to
  # match the map, which is the reliable way to rebuild at the fitted values.
  obj2 <- TMB::MakeADFun(
    data       = data2,
    parameters = obj$env$parList(),
    map        = obj$env$map,
    random     = random_names,
    DLL        = fit$TMBfilename,
    silent     = TRUE)
  obj2$fn(obj2$par)   # evaluate once so last.par is populated
  obj2
}


#' Statistical diagnostics for OSA residuals
#'
#' @description
#' Computes the Stewart and Monnahan (2025) statistical diagnostics for a set of
#' OSA residuals: the standard deviation of the normalized residuals (SDNR) and
#' the lower/upper tail statistics, each with the 95% interval expected under
#' the standard normal null hypothesis (so departures can be judged objectively
#' rather than by eye). Computed per data source (type x fleet) and overall.
#'
#' Under a correctly specified model OSA residuals are already iid standard
#' normal, so the SDNR is simply their sample standard deviation. Its null
#' interval follows the chi-square result for the sample standard deviation of
#' `n` standard normals (Francis 2014); the tail-statistic null intervals are
#' obtained by simulation.
#'
#' @param osa An `rceattle_osa` object from [osa_residuals()], or a data frame
#'   with `residual` and (optionally) `type`/`fleet` columns.
#' @param nsim Number of simulations for the tail-statistic null intervals.
#'   Default 10000.
#' @param probs Lower/upper tail probabilities. Default `c(0.025, 0.975)`.
#' @param seed Seed for the tail-interval simulation (reproducibility).
#'   Default 123.
#'
#' @return A data frame with one row per data source plus an `"all"` row, with
#'   columns: `source`, `type`, `fleet`, `n`, `sdnr`, `sdnr_lo`, `sdnr_hi`,
#'   `lower`, `lower_lo`, `lower_hi`, `upper`, `upper_lo`, `upper_hi`, and the
#'   logical flags `sdnr_ok`, `lower_ok`, `upper_ok` (TRUE when the statistic is
#'   inside its null interval).
#'
#' @references
#' Francis, R.I.C.C. 2014. Replacing the multinomial in stock assessment models:
#'   a first step. Fish. Res. 151:70-84.
#'
#' Stewart, I.J., and Monnahan, C.C. 2025. Can. J. Fish. Aquat. Sci. 82:1-13.
#'
#' @seealso [osa_residuals()]
#' @export
osa_diagnostics <- function(osa, nsim = 10000, probs = c(0.025, 0.975),
                            seed = 123) {

  if (is.null(osa$residual)) {
    stop("'osa' must have a 'residual' column (e.g. output of osa_residuals()).")
  }

  has_groups <- !is.null(osa$type) && !is.null(osa$fleet)
  groups <- if (has_groups) {
    split(osa, list(osa$type, osa$fleet), drop = TRUE)
  } else {
    list(all = osa)
  }

  rows <- lapply(names(groups), function(g) {
    grp <- groups[[g]]
    stat <- .osa_sdnr_tails(grp$residual, nsim = nsim, probs = probs, seed = seed)
    data.frame(
      source = if (has_groups) paste(grp$type[1], "fleet", grp$fleet[1]) else "all",
      type   = if (has_groups) grp$type[1] else NA_character_,
      fleet  = if (has_groups) grp$fleet[1] else NA_integer_,
      stat,
      stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)

  # Overall row across all residuals.
  overall <- .osa_sdnr_tails(osa$residual, nsim = nsim, probs = probs, seed = seed)
  out <- rbind(out, data.frame(source = "all", type = NA_character_,
                               fleet = NA_integer_, overall,
                               stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


#' SDNR and tail statistics with standard-normal null intervals
#'
#' @param resid Numeric vector of residuals (assumed standard normal under H0).
#' @param nsim,probs,seed See [osa_diagnostics()].
#' @return A one-row data frame of statistics and their null intervals.
#' @keywords internal
.osa_sdnr_tails <- function(resid, nsim = 10000, probs = c(0.025, 0.975),
                            seed = 123) {
  resid <- resid[is.finite(resid)]
  n <- length(resid)
  if (n < 2) {
    return(data.frame(n = n, sdnr = NA_real_, sdnr_lo = NA_real_,
                      sdnr_hi = NA_real_, lower = NA_real_, lower_lo = NA_real_,
                      lower_hi = NA_real_, upper = NA_real_, upper_lo = NA_real_,
                      upper_hi = NA_real_, sdnr_ok = NA, lower_ok = NA,
                      upper_ok = NA))
  }

  # SDNR and its chi-square null interval (Francis 2014).
  sdnr    <- stats::sd(resid)
  df      <- n - 1
  sdnr_lo <- sqrt(stats::qchisq(probs[1], df) / df)
  sdnr_hi <- sqrt(stats::qchisq(probs[2], df) / df)

  # Observed lower/upper tail statistics.
  lower <- stats::quantile(resid, probs[1], names = FALSE)
  upper <- stats::quantile(resid, probs[2], names = FALSE)

  # Null intervals for the tail statistics by simulation: draw nsim sets of n
  # standard normals, take each tail quantile, then summarize across draws.
  withr_seed <- function(s, expr) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
    }
    set.seed(s)
    expr
  }
  sim <- withr_seed(seed, matrix(stats::rnorm(n * nsim), nrow = nsim, ncol = n))
  sim_lower <- apply(sim, 1L, stats::quantile, probs = probs[1], names = FALSE)
  sim_upper <- apply(sim, 1L, stats::quantile, probs = probs[2], names = FALSE)
  lower_int <- stats::quantile(sim_lower, probs, names = FALSE)
  upper_int <- stats::quantile(sim_upper, probs, names = FALSE)

  data.frame(
    n        = n,
    sdnr     = sdnr,    sdnr_lo  = sdnr_lo,     sdnr_hi  = sdnr_hi,
    lower    = lower,   lower_lo = lower_int[1], lower_hi = lower_int[2],
    upper    = upper,   upper_lo = upper_int[1], upper_hi = upper_int[2],
    sdnr_ok  = sdnr  >= sdnr_lo     & sdnr  <= sdnr_hi,
    lower_ok = lower >= lower_int[1] & lower <= lower_int[2],
    upper_ok = upper >= upper_int[1] & upper <= upper_int[2])
}


#' @export
print.rceattle_osa <- function(x, ...) {
  cat("Rceattle one-step-ahead (OSA) residuals\n")
  cat("  method:", attr(x, "method"), " seed:", attr(x, "seed"), "\n")
  cat("  ", nrow(x), " residuals across ",
      length(unique(paste(x$type, x$fleet))), " data source(s)\n", sep = "")
  print(utils::head(as.data.frame(x), ...))
  invisible(x)
}


#' @export
summary.rceattle_osa <- function(object, ...) {
  osa_diagnostics(object, ...)
}
