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
#' Survey-index OSA residuals are supported for every index likelihood family
#' (`Index_distribution`). Lognormal IID (`"Lognormal"`) residualizes on the log
#' scale, and the natural-scale `"Normal"` and `"TruncatedNormal"` on the natural
#' scale. The correlated covariance families (`"MVN"` / `"MVNORM"`) are whitened
#' by the lower Cholesky of the fleet's survey covariance Sigma = L L', so the
#' residuals are the multivariate-Gaussian one-step-ahead innovations
#' L^-1 (obs - q*pred) -- the closed form [TMB::oneStepPredict()] reproduces for a
#' Gaussian block.
#'
#' `"TruncatedNormal"` rows are residualized in their own
#' [TMB::oneStepPredict()] call, with `method = "oneStepGeneric"` and a range
#' starting at zero, whatever `method` is passed. Its density differs from
#' `"Normal"` only by `log Phi(mu/sd)`, which is a function of the prediction and
#' not of the observation, so a Gaussian method -- which reads the curvature of
#' the density in the observation -- cannot see the truncation at all and returns
#' the untruncated `(obs - mu)/sd`. Integrating over the family's own support
#' instead gives the truncated CDF
#' `F(x) = [Phi((x - mu)/sd) - Phi(-mu/sd)] / Phi(mu/sd)`,
#' so `qnorm(F(x))` is standard normal by the probability integral transform
#' however hard the truncation bites. The upper limit is finite rather than
#' `Inf` -- ten standard deviations past the largest fitted index in the group,
#' which leaves under 1e-23 of the mass outside while keeping the Laplace inner
#' problem in a region it can solve. That group also runs with
#' `splineApprox = FALSE`, because the spline shortcut integrates over whatever
#' range its profile slice covered. The range is a
#' property of the family, so it cannot be shared with the other fleets:
#' `"Normal"` is genuinely untruncated, and a lognormal fleet's stored observation
#' is `log(obs)`, which is negative for a small index.
#'
#' The size of the correction is the truncated mass: on a fleet predicting 100
#' with an absolute sd of 150 (a quarter of the density below zero), an
#' observation exactly at the prediction has an untruncated residual of 0 and a
#' truncated one of -0.44.
#'
#' Three consequences of that group being residualized separately, worth knowing
#' before reading the output:
#'
#' * **On a model with random effects the exact integration can fail.** It
#'   evaluates the Laplace marginal at arbitrary values of the observation, and
#'   the inner problem does not always converge across the whole support. Rather
#'   than return `NA` for those rows -- which would shrink the sample
#'   [osa_diagnostics()] passes verdict on, without saying so -- the fleet is
#'   recomputed under TMB's spline approximation and a warning says the residuals
#'   for it are approximate. Fixed-effect models are unaffected.
#' * **`sd` is `NA` and `predicted` means something different for this group.**
#'   `oneStepGeneric` returns neither a standard deviation nor the fitted value;
#'   `predicted` is the truncated conditional mean `E[x | x > 0]`, which sits
#'   above the fitted index (163.8 against a fitted 100.0 on the fixture above),
#'   not the fitted index itself. The residual is unaffected.
#' * **Adding a `"TruncatedNormal"` fleet moves the other fleets' residuals on a
#'   random-effects model.** Each [TMB::oneStepPredict()] call marks the rows it
#'   is not residualizing as unconditional, which zeroes their data terms and so
#'   changes the conditional distribution of the latent states. Both orderings
#'   are valid probability-integral-transform sequences, so no value is wrong,
#'   but catch and lognormal-index residuals will not match those from an
#'   otherwise identical model with no truncated fleet. Fixed-effect models are
#'   again unaffected, since their observations are independent.
#'
#' @param object A fitted object of class `Rceattle` (from [fit_mod()]).
#' @param fit deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @param source Character vector of observation sources to residualize: any of
#'   `"ecov"`, `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, or `"all"`.
#'   Defaults to the five non-diet sources (`diet` is opt-in because it applies
#'   only to multispecies models and can be expensive); pass `"all"` to include
#'   `diet`. `"ecov"` is the state-space covariate (QAR1 `observe=` term),
#'   residualized first against its own series, as in WHAM's
#'   `make_osa_residuals()`.
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
#'   \code{\link[parallel]{mclapply}}. Default `TRUE`. This is the main speedup for models
#'   with random effects, where each observation triggers a Laplace
#'   re-evaluation -- it gives a near-linear speedup across cores (set
#'   `options(mc.cores = )` to choose how many; forking falls back to serial on
#'   Windows). Some models -- heavy random-effect structures such as a
#'   time-varying catchability -- abort the forked worker instead of returning;
#'   the loop then recomputes serially, after rebuilding, and prints the worker's
#'   own "irrecoverable exception" message, which comes from C and cannot be
#'   suppressed. That message does not mean the call failed. Pass `FALSE` to skip
#'   the attempt. Only the continuous group is parallelized; the discrete
#'   (randomized-quantile) path always runs serially so it stays reproducible
#'   given `seed`.
#'
#'   **A DSEM fit always computes serially, whatever you pass.** Fitting one
#'   loads the `dsem` DLL alongside `ceattle`, and TMB will not choose between
#'   two loaded models ("Multiple TMB models loaded. Failed to guess DLL name."),
#'   so the parallel attempt is skipped and a message says so. The residuals are
#'   unaffected; only the run time is.
#' @param seed Random seed passed to [TMB::oneStepPredict()] for reproducibility
#'   of randomized-quantile residuals. Default `123`.
#' @param trace Logical; print [TMB::oneStepPredict()] progress. Default `FALSE`.
#' @param ... Further arguments passed to [TMB::oneStepPredict()].
#'
#' @return A data frame of class `rceattle_osa` with one row per residualized
#'   observation and columns `source` (the data source: index/catch/comp/caal/
#'   diet), `fleet`, `fleet_name`, `species`, `sex`, `year`, `age_length_bin`
#'   (the age or length bin the value stands for), `accumulated` (`TRUE` where
#'   tail accumulation folded neighbouring bins into this one, so it covers a
#'   range of ages rather than the one named), `length` (the conditioning length
#'   bin for caal; `NA` otherwise), `index_label` (`"age"`/`"length"`/`NA`), `observed`,
#'   `predicted`, `sd`, and `residual`. For aggregate series `observed` and
#'   `predicted` are on the residualization scale -- log for lognormal catch/index,
#'   natural scale for a `"Normal"` or `"TruncatedNormal"` index, and the
#'   whitened (`L^-1`) scale for an
#'   `"MVN"`/`"MVNORM"` index; for compositions they are bin counts. Carries
#'   `method` and `seed` attributes -- `method` is the string that was passed,
#'   except on a model with a `"TruncatedNormal"` index fleet, where it is the
#'   named vector `c(default = <method>, TruncatedNormal = "oneStepGeneric")`
#'   because that family is residualized with its own method (see above) -- and
#'   (when composition types
#'   are present) a `"pearson"` attribute holding the matching Pearson residuals
#'   so [plot.rceattle_osa()] can show both. The attribute uses this data
#'   frame's column names rather than the data-sheet names
#'   [residuals.Rceattle()] returns, so the two halves of one object read alike.
#'   Note the shared names do not mean a shared scale: in the attribute
#'   `observed` and `predicted` are proportions summing to one within a
#'   fleet-year, with the sample size in `sample_size`, because composition
#'   Pearson residuals are defined on proportions -- not the bin counts the
#'   columns above carry. Do not compare the two directly.
#'   Both describe the bins the likelihood fit, so a fleet with tail
#'   accumulation reports the folded window in each -- with one asymmetry: the
#'   one-step-ahead decomposition drops each group's last bin (it is fixed by
#'   sum-to-N), and under an *old*-tail accumulation that dropped bin is the
#'   upper accumulated one. Such a fleet therefore shows its upper boundary bin
#'   in the Pearson residuals and not in the OSA residuals. Summarize it with
#'   [osa_diagnostics()] and plot it with [plot.rceattle_osa()].
#'
#' @section Negative composition `predicted` values:
#' A composition `predicted` is an expected bin count and cannot truly be
#' negative, but it goes slightly negative where a bin holds almost no fish.
#' Composition observations enter as counts, `(proportion + comp_offset) * N`,
#' and [TMB::oneStepPredict()]'s conditional mean is a numerical step away from
#' the observation, which overshoots below zero when the count is near it.
#' The function warns, naming the count and the years.
#'
#' It is the *bin's* count that drives this, not the composition's sample size:
#' on EBS pollock the negative rows have a median observed count of 0.05 against
#' 4.9 for the rest, while their sample sizes span the same 1 to 821 as
#' everything else (69 of 353 occur above a sample size of 100). A rare age in a
#' well-sampled year does it as readily as a poorly sampled year, so the warning
#' is not by itself evidence of thin data.
#'
#' The values are reported rather than clamped, because a negative expected
#' count is the signal that the bin is too sparse for the decomposition to
#' describe and clamping would hide it. Treat `predicted` on those rows as
#' uninformative; their `residual` standardises the observation against that
#' same mean and is therefore biased positive.
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
#' @examples
#' \dontrun{
#' data(BS2017SS)
#' fit <- fit_mod(BS2017SS, estimateMode = "Hindcast")
#' osa <- osa_residuals(fit, source = c("index", "comp"))
#' plot(osa)
#' }
#' @export
osa_residuals <- function(object = NULL,
                          source   = c("ecov", "index", "catch", "comp", "caal"),
                          method   = "oneStepGaussianOffMode",
                          discrete = FALSE,
                          parallel = TRUE,
                          seed     = 123,
                          trace    = FALSE,
                          ..., fit = NULL) {
  # `fit` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(fit))
    object <- .rce_deprecated_arg(fit, !missing(object), "fit", "object", "osa_residuals")


  # ---- Validate the fit ----
  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).")
  }
  # A DSEM needs nothing special here: these are one-step-ahead residuals for
  # the DATA, and TMB::oneStepPredict() integrates the random effects out
  # whatever their structure -- the GMRF is just another random block, reached
  # through obj$env$random like any other. Unlike process_residuals(), nothing
  # here standardizes a deviation by an assumed per-year prior, so there is no
  # IID assumption to violate.
  if (is.null(object$obj)) {
    stop("`object` has no TMB object ($obj); OSA residuals require the fitted ",
         "model object.")
  }
  em <- object$data_list$estimateMode
  if (!is.null(em) && em >= 3) {
    stop("OSA residuals require a model optimized with estimateMode < 3 ",
         "(the returned objective for estimateMode >= 3 is a debug placeholder ",
         "that oneStepPredict cannot differentiate).")
  }
  if (is.null(object$obj$env$last.par.best)) {
    stop("`object` does not appear to have been optimized; OSA residuals require ",
         "a converged fit.")
  }

  # "diet" is supported but opt-in: it applies only to multispecies models with
  # estimated suitability and can be expensive, so it is not in the default set.
  # "all" is a synonym for every source including diet.
  valid_sources <- c("ecov", "index", "catch", "comp", "caal", "diet", "all")
  source <- match.arg(source, choices = valid_sources, several.ok = TRUE)
  if ("all" %in% source) source <- c("ecov", "index", "catch", "comp", "caal", "diet")

  # Build the full OSA observation data (comp / caal / diet segments) on demand.
  # This works from any fit and no longer requires fitting with
  # fit_control(osa = TRUE): build_osa_data() reads only the *_ctl / *_obs
  # arrays the model already carries, and `obs_ctl` maps each obsvec position
  # back to its source. The same regenerated data is reused by .osa_build_obj().
  osa_dat <- build_osa_data(object$obj$env$data, build_osa = TRUE)
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
  # `fleet_code` is included as a tie-break so the sequence is fully determined:
  # within a given (source, year) -- or (year, source) below -- several fleets can
  # supply observations (e.g. multiple surveys), and without an explicit key their
  # order would fall to the incidental row order of obs_ctl. Under random effects
  # that within-year fleet order is not cosmetic -- it shifts individual residuals
  # (though not the overall N(0,1) properties / validity conclusion; Trijoulet et
  # al. 2023). Ascending fleet_code also matches WHAM's within-stage ordering
  # (index fleets before the fishery), so it does not disturb the WHAM cross-check.
  #
  # Optional -- chronological: year first, then data source, then fleet, then bin.
  # Puts the conditional sequence in time order.
  # sel <- sel[order(sel$year, match(sel$source, source), sel$fleet_code,
  #                 sel$bin_index, na.last = TRUE), , drop = FALSE]

  # Order by source, then year, fleet, bin: the covariate (ecov) is residualized
  # first and standalone, then index/catch -> comp -> CAAL each conditional on the
  # earlier types.
  #
  # This is the order that agrees with WHAM, MEASURED rather than assumed. WHAM's
  # own obsvec is year-major (its `ind` cycles logindex, logcatch, indexpal,
  # catchpal, indexcaal within each year), so the chronological alternative below
  # looks like the closer match -- and is not. Against WHAM's cached residuals on
  # the GOAcod growth comparison, source-major gives survey-lencomp r = 0.99995
  # (median |difference| 3.6e-5, 2.4% of bins over 0.01) and CAAL r = 0.99998;
  # switching to year-major degrades every one of them -- survey lencomp to
  # r = 0.99764 with 34% of bins over 0.01, catch from r = 0.909 to 0.718. Do not
  # "fix" this to match WHAM's obsvec layout; it is already the better match.
  sel <- sel[order(match(sel$source, source), sel$year, sel$fleet_code,
                    sel$bin_index, na.last = TRUE), , drop = FALSE]

  # Rebuild the model in OSA mode (osa = TRUE) at the fitted parameters. This
  # is required so the composition (comp/caal) likelihoods read their counts
  # from obsvec and use the unweighted, proper density that oneStepPredict needs.
  # It leaves the aggregate catch/index likelihood unchanged, so the same
  # rebuilt object serves all observation types.
  obj_osa <- .osa_build_obj(object, osa_dat)

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

  # `Index_distribution = "TruncatedNormal"` is residualized on its own, with
  # oneStepGeneric over a (0, Inf) range. That integrates the density over the
  # range and normalizes by the integral, giving the truncated CDF
  #   F(x) = [Phi((x-mu)/sd) - Phi(-mu/sd)] / Phi(mu/sd)
  # so qnorm(F(x)) is standard normal. The Gaussian methods, including the
  # default, cannot see the truncation -- it enters the density only through
  # log Phi(mu/sd), a function of the prediction and not the observation -- and
  # would return the untruncated residual wherever truncation carries real mass.
  #
  # The range belongs to the FAMILY, so only these rows may have it: "Normal"
  # is genuinely untruncated, and a lognormal fleet's obsvec entry is log(obs),
  # which is negative for a small index. Keyed off index_ll_type, the same
  # vector build_osa_data() used to decide whether the entry holds obs or
  # log(obs), so the two cannot disagree.
  fc_osa   <- object$data_list$fleet_control
  ill_osa  <- object$obj$env$data$index_ll_type
  sel_trunc <- rep(FALSE, nrow(sel))
  if (!is.null(ill_osa) && length(ill_osa)) {
    sel_trunc <- sel$source == "index" & ill_osa[sel$fleet_code] %in% 4L
    sel_trunc[is.na(sel_trunc)] <- FALSE
  }
  sel_group <- paste0(sel_discrete, "|", sel_trunc)

  # Say so. `method` is overridden for these rows whatever the caller passed, and
  # a silent override is exactly the kind of thing that makes a Q-Q plot hard to
  # account for later. Also flag the cost: the exact integration is several times
  # slower per observation than a Gaussian approximation.
  if (any(sel_trunc)) {
    trunc_fleets <- unique(as.character(
      fc_osa$Fleet_name[match(sel$fleet_code[sel_trunc], fc_osa$Fleet_code)]))
    message("osa_residuals(): fleet(s) ", paste(trunc_fleets, collapse = ", "),
            " use Index_distribution = \"TruncatedNormal\", whose truncation a ",
            "Gaussian method cannot see. Those ", sum(sel_trunc), " observation(s) ",
            "are residualized with method = \"oneStepGeneric\" over the ",
            "truncated support instead of \"", method,
            "\" -- exact, and slower. See ?osa_residuals.")
  }

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

  # Upper limit for the truncated group's exact integration. `range = c(0, Inf)`
  # is mathematically what the support is, but integrate()'s infinite-limit
  # substitution probes observations far out in the tail, and on a random-effects
  # model the Laplace inner problem does not converge there: the integrand comes
  # back non-finite, integrate() errors, and TMB's try() writes NA. On a 15-year
  # fixture with random recruitment that lost 5 of 15 residuals -- silently, and
  # under a warning blaming the fit's convergence.
  #
  # Bound it instead. The upper limit covers every fitted index in the group by
  # ten of its own standard deviations, and every observation being residualized,
  # so the mass discarded is below 1e-23 while the integrand stays in a region the
  # inner problem can solve. oneStepGeneric normalizes over the range it is given,
  # so what matters is that the range contains the observation and effectively all
  # the density -- not that it reaches infinity.
  .trunc_range <- function(rows) {
    dr <- sel$data_row[rows]
    mu <- suppressWarnings(as.numeric(object$quantities$index_hat[dr]))
    sg <- suppressWarnings(as.numeric(.observation_sd(object$quantities, "index")[dr]))
    ob <- suppressWarnings(as.numeric(osa_dat$obsvec[sel$obs_pos[rows] + 1L]))
    hi <- suppressWarnings(max(c(mu + 10 * sg, ob * 1.5), na.rm = TRUE))
    if (!is.finite(hi) || hi <= 0) hi <- Inf   # nothing usable; fall back
    c(0, hi)
  }

  .run_osp <- function(rows, dsc, trunc = FALSE, spline = FALSE) {
    # The Gaussian methods are continuous-only, so discrete groups use the
    # generic (CDF-based) method, which supports randomized quantile residuals.
    # Parallelize only the continuous group; the discrete path uses the seeded
    # RNG, so keep it serial to stay bit-reproducible across runs.
    osp <- function(par) {
      args <- list(
        obj                 = obj_osa,
        observation.name    = "obsvec",
        data.term.indicator = "keep",
        method              = if (dsc || trunc) "oneStepGeneric" else method,
        # Only the truncated family restricts the support; every other group
        # keeps TMB's default (-Inf, Inf).
        range               = if (trunc) .trunc_range(rows) else c(-Inf, Inf),
        subset              = sel$obs_pos[rows] + 1L,
        discrete            = dsc,
        parallel            = par,
        seed                = seed,
        trace               = trace)
      dots <- list(...)
      if (length(dots) && (is.null(names(dots)) || !all(nzchar(names(dots))))) {
        stop("osa_residuals(): all arguments passed through `...` to ",
             "oneStepPredict() must be named; an unnamed one would be matched ",
             "positionally against an argument this function already sets.",
             call. = FALSE)
      }
      clash <- intersect(names(dots), names(args))
      if (length(clash)) {
        warning("osa_residuals(): ignoring ", paste(clash, collapse = ", "),
                " passed through `...` -- ", if (length(clash) > 1) "these are"
                else "this is", " set per observation group by this function ",
                "(see ?osa_residuals). Use the `method` argument to choose a ",
                "method.", call. = FALSE)
      }
      args <- c(args, dots[setdiff(names(dots), names(args))])
      # `range` only binds the integration limits when oneStepGeneric integrates
      # the density itself. Under its default splineApprox = TRUE it instead
      # splines a tmbprofile() slice and integrates over the range that slice
      # happened to cover, which does not respect the (0, Inf) support: on a
      # fixture at mu = x = 100, sd = 150 that returns -0.790 where the exact
      # truncated transform is -0.437. So the exact path is the one that makes
      # this group worth splitting out. `spline` lets the caller below retry with
      # the approximation when the exact path cannot finish -- see .run_osp().
      if (trunc) args$splineApprox <- isTRUE(spline)
      do.call(TMB::oneStepPredict, args)
    }

    # oneStepPredict(parallel = TRUE) forks via mclapply. A worker that dies --
    # some model/observation combinations abort the child rather than return --
    # leaves an error object where a gradient should be, which surfaces from deep
    # inside TMB as "non-numeric argument to mathematical function". Retry
    # serially instead: slower, but it returns residuals rather than a message
    # that points nowhere near the cause.
    #
    # The retry must build a FRESH object. Re-entering the one the failed attempt
    # used ends the R session outright ("An irrecoverable exception occurred")
    # instead of recovering, which is worse than the failure it is handling.
    #
    # What is established, from GOA pollock 2025 (164 random effects, 109 of them
    # from an ar1 and a rw catchability linkage): forked workers abort at any
    # width above one, and mclapply normally absorbs that -- a direct
    # oneStepPredict() over the same 90 index observations returned all of them
    # at 1, 2, 4, 8 and 11 cores despite aborts at every width above one.
    # Sometimes the absorption fails and the error surfaces here; that case did
    # not reproduce on demand, so it looks load-dependent rather than a property
    # of the model. Either way the parent must not reuse the object afterwards.
    # Serial (`parallel = FALSE`) never fails.
    want_par <- parallel_ok && !dsc
    res <- if (!want_par) osp(FALSE) else
      tryCatch(osp(TRUE), error = function(e) {
        message("osa_residuals(): the parallel one-step-ahead loop failed (",
                conditionMessage(e), "); recomputing serially.")
        obj_osa <<- .osa_build_obj(object, osa_dat, force = TRUE)
        osp(FALSE)
      })
    # oneStepGeneric returns no `observation` column, so the generic groups (the
    # discrete compositions and the truncated index) would otherwise report NA
    # for a value that is simply the obsvec entry that was residualized.
    obs_col <- get_col(res, "observation")
    if (all(is.na(obs_col))) {
      obs_col <- as.numeric(osa_dat$obsvec[sel$obs_pos[rows] + 1L])
    }
    data.frame(.row = rows, observed = obs_col,
               predicted = get_col(res, "mean"), sd = get_col(res, "sd"),
               residual = res$residual)
  }
  osa <- do.call(rbind, lapply(unique(sel_group), function(g) {
    rows  <- which(sel_group == g)
    trunc <- sel_trunc[rows[1]]
    res   <- .run_osp(rows, dsc = sel_discrete[rows[1]], trunc = trunc)

    # The exact integration evaluates the Laplace marginal at arbitrary values of
    # the observation, and on a random-effects model the inner Newton problem
    # does not always converge there: integrate() errors and TMB writes NA. That
    # is worse than an approximate residual -- the rows vanish, and
    # osa_diagnostics() then passes verdict on whatever survived without saying
    # how many it lost. Retry the whole group under TMB's spline approximation,
    # which is robust but does NOT respect the (0, Inf) support, so the result is
    # approximate in the same direction the Gaussian methods are. Retrying the
    # whole group rather than the failed rows keeps one method per fleet.
    if (trunc && any(!is.finite(res$residual))) {
      n_bad <- sum(!is.finite(res$residual))
      warning("osa_residuals(): exact integration of the TruncatedNormal ",
              "density failed for ", n_bad, " of ", nrow(res), " observation(s) ",
              "-- the Laplace inner problem does not converge across the whole ",
              "support, which happens on models with random effects. The fleet ",
              "was recomputed with TMB's spline approximation, so ITS RESIDUALS ",
              "ARE APPROXIMATE and do not carry the truncation exactly. Treat ",
              "them as indicative; the other fleets are unaffected.",
              call. = FALSE)
      res <- .run_osp(rows, dsc = sel_discrete[rows[1]], trunc = TRUE,
                      spline = TRUE)
      attr(res, "approx") <- TRUE
    }
    res
  }))
  osa <- osa[order(osa$.row), , drop = FALSE]   # restore chronological 'sel' order
  index_label <- c("age", "length")[sel$comp_type + 1L]   # NA for aggregates
  fc <- object$data_list$fleet_control                        # fleet code -> name
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
    accumulated    = if (is.null(sel$accumulated)) FALSE else sel$accumulated,
    length         = sel$length,
    index_label    = index_label,
    observed      = osa$observed,
    predicted     = osa$predicted,
    sd            = osa$sd,
    residual      = osa$residual,
    stringsAsFactors = FALSE)

  rownames(out) <- NULL
  class(out) <- c("rceattle_osa", "data.frame")
  # Record what was actually used, not what was asked for: the truncated index
  # family is residualized with its own method regardless of `method`, so a
  # single string would misdescribe those rows. Stays a plain string when nothing
  # was overridden.
  attr(out, "method") <- if (any(sel_trunc)) {
    c(default = method, TruncatedNormal = "oneStepGeneric")
  } else method
  attr(out, "seed")   <- seed
  # Per-species bin counts, so plot() can split joint-sex (Sex == 3) composition
  # bins onto a single age/length axis (males are stored as bins nbin+1 .. 2*nbin).
  attr(out, "nages")    <- object$data_list$nages
  attr(out, "nlengths") <- object$data_list$nlengths

  # Attach the matching Pearson residuals for composition sources so the
  # plot() method can show OSA and Pearson bubbles side by side.
  comp_types <- intersect(unique(out$source), c("comp", "caal"))
  if (length(comp_types) > 0) {
    pear <- tryCatch(
      stats::residuals(object, type = "pearson", source = comp_types),
      error = function(e) NULL)
    # residuals() is a general-purpose method and names its columns in the
    # data-sheet style (Fleet_code, Year, Observed, ...); this data frame names
    # them in the style of the object it is attached to. Carrying both
    # conventions on one object means a reader has to know which half they are
    # holding, so rename to match here.
    if (!is.null(pear)) {
      nm <- c(Source = "source", Fleet_code = "fleet", Fleet_name = "fleet_name",
              Species = "species", Sex = "sex", Year = "year",
              Age0_Length1 = "index_label", Bin = "age_length_bin",
              Length = "length", Sample_size = "sample_size",
              Accumulated = "accumulated",
              Observed = "observed", Fitted = "predicted",
              Residual = "residual")
      hit <- names(pear) %in% names(nm)
      names(pear)[hit] <- unname(nm[names(pear)[hit]])

      # Same form as the OSA frame, not a second encoding under a second name.
      if (!is.null(pear$index_label)) {
        pear$index_label <- ifelse(pear$source == "caal", "age",
                                   ifelse(!is.na(pear$index_label) &
                                            pear$index_label == 1,
                                          "length", "age"))
      }
    }
    attr(out, "pearson") <- pear
  }

  n_bad <- sum(!is.finite(out$residual))
  if (n_bad > 0) {
    warning(n_bad, " of ", nrow(out), " OSA residual(s) are non-finite. This ",
            "usually indicates a poorly converged fit or very sparse / ",
            "degenerate compositions (common for conditional age-at-length); ",
            "check model convergence before interpreting the residuals.")
  }

  # A composition `predicted` is an expected bin count, so it cannot be
  # negative. It goes slightly negative where the bin holds almost no fish:
  # composition observations enter as counts, (proportion + comp_offset) * N,
  # and oneStepPredict()'s conditional mean is a numerical step from the
  # observation, which overshoots below zero when the count is near it.
  #
  # It is the BIN's count that matters, not the composition's sample size. On
  # EBS pollock the negative rows have a median observed count of 0.05 against
  # 4.9 for the rest, while their sample sizes span the same 1 to 821 as
  # everything else -- 69 of 353 occur above a sample size of 100. A rare age in
  # a well-sampled year does this just as readily as a poorly sampled year.
  #
  # Reported rather than clamped: the negative value is the signal that the bin
  # is too sparse for the decomposition to describe, and clamping hides it.
  is_comp <- out$source %in% c("comp", "caal", "diet")
  bad <- is_comp & is.finite(out$predicted) & out$predicted < 0
  if (any(bad)) {
    yrs <- sort(unique(out$year[bad]))
    warning(sum(bad), " composition `predicted` value(s) are negative, in ",
            "year(s) ", paste(utils::head(yrs, 5), collapse = ", "),
            if (length(yrs) > 5) ", ..." else "",
            ". An expected count cannot be negative: those bins hold too few ",
            "fish for the one-step-ahead decomposition to describe, so treat ",
            "`predicted` there as uninformative and their residuals as biased ",
            "positive. See ?osa_residuals.", call. = FALSE)
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
#' @param force Build a new object even when `fit` is already in OSA mode, where
#'   this otherwise returns `fit$obj` itself. The retry after a failed parallel
#'   one-step-ahead loop needs a genuinely new one.
#' @return A TMB ADFun object with `osa_mode = 1`.
#' @keywords internal
.osa_build_obj <- function(fit, osa_dat = NULL, force = FALSE) {
  obj <- fit$obj
  # A model already fitted in OSA mode can use its own object -- except when the
  # caller needs a genuinely new one. The retry after a failed parallel loop does:
  # handing back `fit$obj` there would reuse the object the failure touched,
  # which is the thing that ends the R session.
  if (!force && !is.null(obj$env$data$osa_mode) &&
      obj$env$data$osa_mode == 1L) {
    return(obj)
  }
  # Regenerate the full OSA observation vector (comp / CAAL / diet segments) on
  # demand, so residuals no longer require fitting with fit_control(osa = TRUE).
  # build_osa_data() reads only the *_ctl / *_obs arrays the model already
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
#' @return A data frame (class `"rceattle_osa_diagnostics"`, so it prints as a
#'   compact severity-tagged summary; every column is still there and `$` works
#'   as before) with one row per data source plus an `"all"` row, with columns:
#'   `group` (the `"<source> fleet <n>"` label), `source`, `fleet`, `n`, `sdnr`,
#'   `sdnr_lo`, `sdnr_hi`, `lower`, `lower_lo`, `lower_hi`, `upper`, `upper_lo`,
#'   `upper_hi`, and the logical flags `sdnr_ok`, `lower_ok`, `upper_ok` (TRUE
#'   when the statistic is inside its null interval). On the `"all"` row
#'   `source` and `fleet` are `NA`.
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
  # Still a data frame -- every column and `$` access is unchanged. The class
  # only adds a print method, so the sixteen columns stop wrapping across three
  # screen-widths with no verdict; see print.rceattle_osa_diagnostics().
  class(out) <- c("rceattle_osa_diagnostics", "data.frame")
  out
}


#' @export
print.rceattle_osa_diagnostics <- function(x, ...) {
  df <- as.data.frame(x)
  per <- df[df$group != "all" | !is.na(df$source), , drop = FALSE]
  all_row <- df[df$group == "all" & is.na(df$source), , drop = FALSE]

  # SDNR is the headline statistic, so it carries WARN; a tail outside its null
  # interval with an acceptable SDNR is a NOTE. Both intervals are simulated
  # under the standard-normal null, so "outside" already means "further than
  # chance", and neither is a FAIL: these are diagnostics on fit, not a broken
  # model.
  sev <- rep("OK", nrow(per))
  sev[!is.na(per$lower_ok) & !per$lower_ok] <- "NOTE"
  sev[!is.na(per$upper_ok) & !per$upper_ok] <- "NOTE"
  sev[!is.na(per$sdnr_ok)  & !per$sdnr_ok]  <- "WARN"

  n_bad <- sum(sev != "OK")
  .rce_diag_header(
    "OSA diagnostics", .rce_worst(sev),
    paste0(.rce_n_of(n_bad, nrow(per)),
           " source(s) outside a null interval",
           if (nrow(all_row)) paste0("; overall SDNR ",
                                     formatC(all_row$sdnr[1], format = "f", digits = 2),
                                     " (", formatC(all_row$sdnr_lo[1], format = "f", digits = 2),
                                     "-", formatC(all_row$sdnr_hi[1], format = "f", digits = 2),
                                     ")") else ""))

  # Worst first, so a long fleet list opens with what needs attention.
  ord <- order(match(sev, .CONV_SEVERITY), decreasing = TRUE)
  per <- per[ord, , drop = FALSE]; sev <- sev[ord]
  per$.tag <- .rce_sev_tag(sev)
  per$.sdnr_int <- sprintf("%.2f-%.2f", per$sdnr_lo, per$sdnr_hi)

  .rce_diag_table(per, c(" " = ".tag", "source" = "source", "fleet" = "fleet",
                         "n" = "n", "sdnr" = "sdnr", "null" = ".sdnr_int"))
  cat("  tails and their intervals are in the data frame",
      "(lower/upper, *_lo, *_hi)\n")
  invisible(x)
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
  # `method` is a named vector when a family overrode the caller's choice.
  m <- attr(x, "method")
  m_txt <- if (length(m) > 1L && !is.null(names(m))) {
    paste(paste0(names(m), " = ", m), collapse = ", ")
  } else as.character(m)
  cat("  method:", m_txt, " seed:", attr(x, "seed"), "\n")
  cat("  ", nrow(x), " residuals across ",
      length(unique(paste(x$source, x$fleet))), " data source(s)\n", sep = "")
  print(utils::head(as.data.frame(x), ...))
  invisible(x)
}


#' @export
summary.rceattle_osa <- function(object, ...) {
  osa_diagnostics(object, ...)
}
