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
#' @param discrete Logical; whether to treat *composition* (comp / caal / diet)
#'   observations as discrete. Default `FALSE` (continuous, matching how CEATTLE
#'   fits the composition likelihood with effective-sample-size-scaled counts).
#'   When `TRUE`, composition residuals are randomized quantile residuals (Dunn
#'   and Smyth 1996) and so are stochastic; set `seed` for reproducibility. The
#'   aggregate index/catch series are always continuous (lognormal); the
#'   [TMB::oneStepPredict()] call is split by observation type so `discrete` is
#'   applied correctly per type (the discrete group uses the generic CDF-based
#'   method rather than `method`).
#' @param parallel Logical; compute the per-observation OSA loop in parallel via
#'   [parallel::mclapply()]. Default `TRUE`. This is the main speedup for models
#'   with random effects, where each observation triggers a Laplace
#'   re-evaluation -- it gives a near-linear speedup across cores (set
#'   `options(mc.cores = )` to choose how many; forking falls back to serial on
#'   Windows). Only the continuous group is parallelized; the discrete
#'   (randomized-quantile) path always runs serially so it stays reproducible
#'   given `seed`.
#' @param seed Random seed passed to [TMB::oneStepPredict()] for reproducibility
#'   of randomized-quantile residuals. Default `123`.
#' @param trace Logical; print [TMB::oneStepPredict()] progress. Default `FALSE`.
#' @param ... Further arguments passed to [TMB::oneStepPredict()].
#'
#' @return A data frame of class `rceattle_osa` with one row per residualized
#'   observation and columns `source` (the data source: index/catch/comp/caal/
#'   diet), `fleet`, `fleet_name`, `species`, `sex`, `year`, `age_length_bin`
#'   (the age or length bin index), `length` (the conditioning length bin for
#'   caal; `NA` otherwise), `index_label` (`"age"`/`"length"`/`NA`), `observed`,
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
                          source   = c("index", "catch", "comp", "caal"),
                          method   = "oneStepGaussianOffMode",
                          discrete = FALSE,
                          parallel = TRUE,
                          seed     = 123,
                          trace    = FALSE,
                          ...) {

  # ---- Validate the fit ----
  if (!inherits(fit, "Rceattle")) {
    stop("'fit' must be a fitted Rceattle model (from fit_mod()).")
  }
  if (is.null(fit$obj)) {
    stop("'fit' has no TMB object ($obj); OSA residuals require the fitted ",
         "model object.")
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

  # Build the full OSA observation data (comp / caal / diet segments) on demand.
  # This works from any fit and no longer requires fitting with
  # fit_control(osa = TRUE): build_osa_data() reads only the *_ctl / *_obs
  # arrays the template already carries, and `obs_ctl` maps each obsvec position
  # back to its source. The same regenerated data is reused by .osa_build_obj().
  osa_dat <- build_osa_data(fit$obj$env$data, build_osa = TRUE)
  obs_ctl <- osa_dat$obs_ctl

  # ---- Select the obsvec positions to residualize ----
  sel <- obs_ctl[obs_ctl$source %in% source & !obs_ctl$is_last_bin, , drop = FALSE]
  if (nrow(sel) == 0) {
    stop("No observations of source(s) ", paste(source, collapse = ", "),
         " are available for OSA residuals in this model.")
  }

  # ---- CONDITIONING ORDER (this is the one-step-ahead sequence) ----------------
  # `subset` order = the order in which oneStepPredict() conditions each
  # observation on the previously-residualized ones. It does not change the joint
  # likelihood, but for compositions it changes the per-bin residual values (the
  # within-multinomial conditional binomials are sequenced in this order).
  #
  # NOTE: for a fixed-effects fit (random_rec = FALSE, no random effects) the
  # observations are independent given the parameters, so the OSA residuals are
  # invariant to this subset order
  #
  # Optional -- chronological: year first, then data source, then bin. Puts the
  # conditional sequence in time order. 
  # sel <- sel[order(sel$year, match(sel$source, source), sel$bin_index,
  #                 na.last = TRUE), , drop = FALSE]

  # WHAM-style -- type-blocked: source first (aggregate index/catch, then comp,
  # then caal), then year, then bin. This reproduces
  # WHAM's make_osa_residuals() conditioning (aggregate -> comp -> CAAL,
  # each conditional on the previous types).
  sel <- sel[order(match(sel$source, source), sel$year, sel$bin_index,
                    na.last = TRUE), , drop = FALSE]

  # Rebuild the model in OSA mode (osa = TRUE) at the fitted parameters. This
  # is required so the composition (comp/caal) likelihoods read their counts
  # from obsvec and use the unweighted, proper density that oneStepPredict needs.
  # It leaves the aggregate catch/index likelihood unchanged, so the same
  # rebuilt object serves all observation types.
  obj_osa <- .osa_build_obj(fit, osa_dat)

  # ---- Compute the OSA residuals ----
  # oneStepPredict() applies a single `discrete` setting per call, but the
  # observation types differ: the aggregate index/catch series are continuous
  # (lognormal), while composition (comp/caal/diet) residuals can be treated as
  # discrete via `discrete = TRUE`. Because composition data is normalized, default
  # is continuous. Residualize each discrete group in its own call so the
  # setting is correct per type; with the default `discrete = FALSE` every type
  # is continuous and this is a single call. `subset` uses 1-based indices into
  # obsvec (obs_pos is 0-based); '.row' maps each result back to its 'sel' row.
  get_col <- function(df, nm) {
    if (!is.null(df[[nm]])) df[[nm]] else rep(NA_real_, nrow(df))
  }
  is_comp      <- sel$source %in% c("comp", "caal", "diet")
  sel_discrete <- ifelse(is_comp, isTRUE(discrete), FALSE)

  # oneStepPredict(parallel = TRUE) calls TMB::openmp(), which can only resolve
  # the active model when a single TMB DLL is loaded; with several loaded (e.g.
  # Rceattle alongside WHAM, or the full test suite) it errors with "Multiple TMB
  # models loaded". Probe openmp() once and silently fall back to serial when it
  # is unavailable, so parallel = TRUE stays a safe default everywhere.
  parallel_ok <- isTRUE(parallel) &&
    !inherits(try(TMB::openmp(), silent = TRUE), "try-error")
  if (isTRUE(parallel) && !parallel_ok) {
    message("osa_residuals(): parallel unavailable (multiple TMB models loaded?); ",
            "computing one-step-ahead residuals serially.")
  }

  .run_osp <- function(rows, dsc) {
    # The Gaussian methods are continuous-only, so discrete groups use the
    # generic (CDF-based) method, which supports randomized quantile residuals.
    # Parallelize only the continuous group; the discrete path uses the seeded
    # RNG, so keep it serial to stay bit-reproducible across runs.
    res <- TMB::oneStepPredict(
      obj                 = obj_osa,
      observation.name    = "obsvec",
      data.term.indicator = "keep",
      method              = if (dsc) "oneStepGeneric" else method,
      subset              = sel$obs_pos[rows] + 1L,
      discrete            = dsc,
      parallel            = parallel_ok && !dsc,
      seed                = seed,
      trace               = trace,
      ...)
    data.frame(.row = rows, observed = get_col(res, "observation"),
               predicted = get_col(res, "mean"), sd = get_col(res, "sd"),
               residual = res$residual)
  }
  osa <- do.call(rbind, lapply(unique(sel_discrete),
                  function(dsc) .run_osp(which(sel_discrete == dsc), dsc)))
  osa <- osa[order(osa$.row), , drop = FALSE]   # restore chronological 'sel' order
  index_label <- c("age", "length")[sel$comp_type + 1L]   # NA for aggregates
  fc <- fit$data_list$fleet_control                        # fleet code -> name
  fleet_name <- if (!is.null(fc)) {
    fc$Fleet_name[match(sel$fleet_code, fc$Fleet_code)]
  } else NA_character_

  out <- data.frame(
    source         = sel$source,
    fleet          = sel$fleet_code,
    fleet_name     = fleet_name,
    species        = sel$species,
    sex            = sel$sex,
    year           = sel$year,
    age_length_bin = sel$age_length_bin,
    length         = sel$length,
    index_label    = index_label,
    observed      = osa$observed,
    predicted     = osa$predicted,
    sd            = osa$sd,
    residual      = osa$residual,
    stringsAsFactors = FALSE)

  rownames(out) <- NULL
  class(out) <- c("rceattle_osa", "data.frame")
  attr(out, "method") <- method
  attr(out, "seed")   <- seed
  # Per-species bin counts, so plot() can split joint-sex (Sex == 3) composition
  # bins onto a single age/length axis (males are stored as bins nbin+1 .. 2*nbin).
  attr(out, "nages")    <- fit$data_list$nages
  attr(out, "nlengths") <- fit$data_list$nlengths

  # Attach the matching Pearson residuals for composition sources so the
  # plot() method can show OSA and Pearson bubbles side by side.
  comp_types <- intersect(unique(out$source), c("comp", "caal"))
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
#' @param osa_dat Optional pre-built OSA observation data (the list returned by
#'   [build_osa_data()] with `build_osa = TRUE`) to reuse instead of rebuilding
#'   it. `NULL` (the default) rebuilds it from `fit$obj$env$data`.
#' @return A TMB ADFun object with `osa_mode = 1`.
#' @keywords internal
.osa_build_obj <- function(fit, osa_dat = NULL) {
  obj <- fit$obj
  if (!is.null(obj$env$data$osa_mode) && obj$env$data$osa_mode == 1L) {
    return(obj)
  }
  # Regenerate the full OSA observation vector (comp / CAAL / diet segments) on
  # demand, so residuals no longer require fitting with fit_control(osa = TRUE).
  # build_osa_data() reads only the *_ctl / *_obs arrays the template already
  # carries in obj$env$data, so the result is identical to an osa = TRUE fit.
  # `osa_dat` lets the caller pass a pre-built copy to avoid recomputing it.
  data2 <- if (is.null(osa_dat)) build_osa_data(obj$env$data, build_osa = TRUE) else osa_dat
  data2$obs_ctl <- NULL   # R-side metadata table, not a TMB DATA input
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
  # parList(last.par.best) returns the best (MLE) parameters as a list structured
  # to match the map -- the reliable way to rebuild at the fitted values. Use
  # last.par.best (not the parList() default last.par) so OSA residuals are
  # computed at the optimum even when getsd = FALSE (no sdreport re-eval).
  obj2 <- TMB::MakeADFun(
    data       = data2,
    parameters = obj$env$parList(par = obj$env$last.par.best),
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

  has_groups <- !is.null(osa$source) && !is.null(osa$fleet)
  groups <- if (has_groups) {
    split(osa, list(osa$source, osa$fleet), drop = TRUE)
  } else {
    list(all = osa)
  }

  rows <- lapply(names(groups), function(g) {
    grp <- groups[[g]]
    stat <- .osa_sdnr_tails(grp$residual, nsim = nsim, probs = probs, seed = seed)
    data.frame(
      group  = if (has_groups) paste(grp$source[1], "fleet", grp$fleet[1]) else "all",
      source = if (has_groups) grp$source[1] else NA_character_,
      fleet  = if (has_groups) grp$fleet[1] else NA_integer_,
      stat,
      stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)

  # Overall row across all residuals.
  overall <- .osa_sdnr_tails(osa$residual, nsim = nsim, probs = probs, seed = seed)
  out <- rbind(out, data.frame(group = "all", source = NA_character_,
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
      length(unique(paste(x$source, x$fleet))), " data source(s)\n", sep = "")
  print(utils::head(as.data.frame(x), ...))
  invisible(x)
}


#' @export
summary.rceattle_osa <- function(object, ...) {
  osa_diagnostics(object, ...)
}
