#' Process residuals for an Rceattle model's random-effect processes
#'
#' @description
#' One-sample process residuals for the model's random-effect (or penalized)
#' process deviations, in the style of SAM's `procres()` (Nielsen and Berg
#' 2014). They validate the *process* model -- whether the deviations behave like
#' their assumed iid normal process -- a complement to the observation-based
#' [osa_residuals()].
#'
#' Each set of deviations carries a Gaussian process prior in the model. The
#' posterior *mode* of the deviations is shrunk toward that prior, so -- following
#' SAM -- a single draw is taken from the joint posterior of the deviations (from
#' the joint precision when they are random effects, or the fixed-effect
#' covariance when they are penalized fixed effects) and standardized by the
#' process standard deviation. Under a correctly specified process these are
#' approximately iid N(0, 1). Because a random draw is used the residuals are
#' stochastic; set `seed` for reproducibility.
#'
#' Supported processes and their deviations:
#' \describe{
#'   \item{`"recruitment"`}{`rec_dev`, prior `N(R_sd^2/2, R_sd^2)` per species.}
#'   \item{`"initial"`}{`init_dev`, prior `N(R_sd^2/2, R_sd^2)` per species.}
#'   \item{`"catchability"`}{`index_q_dev`, prior `N(0, q_dev_sd^2)` per index.}
#' }
#' `"all"` returns every supported process present in the fit. Selectivity and
#' natural-mortality deviations (which use random-walk / 2D-AR1 priors) are not
#' yet supported. For processes with iid priors the standardization is exact; for
#' AR1-configured catchability it standardizes by the marginal SD.
#'
#' @param fit A fitted `Rceattle` model. The targeted deviations must be
#'   estimated -- as random effects (e.g. `random_rec = TRUE`) or as penalized
#'   fixed effects -- with a usable covariance.
#' @param process One of `"recruitment"`, `"initial"`, `"catchability"`, or
#'   `"all"`.
#' @param seed Seed for the posterior draw. Default 123.
#'
#' @return A data frame of class `rceattle_osa` (so it can be passed to
#'   [osa_diagnostics()] and [plot.rceattle_osa()]) with one row per deviation:
#'   columns `type` (the process), `fleet`, `species`, `year`, `age_or_length`,
#'   and `residual`.
#'
#' @references
#' Nielsen, A., and Berg, C.W. 2014. Estimation of time-varying selectivity in
#'   stock assessments using state-space models. Fish. Res. 158:96-101.
#'
#' @seealso [osa_residuals()], [osa_diagnostics()]
#' @export
process_residuals <- function(fit,
                              process = c("recruitment", "initial",
                                          "catchability", "all"),
                              seed = 123) {
  if (!inherits(fit, "Rceattle")) {
    stop("'fit' must be a fitted Rceattle model (from fit_mod()).")
  }
  if (is.null(fit$obj)) stop("'fit' has no TMB object ($obj).")
  process <- match.arg(process)

  specs <- list(recruitment  = "rec_dev",
                initial      = "init_dev",
                catchability = "index_q_dev")
  todo <- if (process == "all") names(specs) else process

  out <- lapply(todo, function(p) {
    res <- tryCatch(.process_residual_one(fit, p, specs[[p]], seed),
                    error = function(e) {
                      if (process == "all") NULL else stop(e)
                    })
    res
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) {
    stop("No supported process deviations were found in this fit.")
  }
  rownames(out) <- NULL
  class(out) <- c("rceattle_osa", "data.frame")
  attr(out, "method") <- "procres"
  attr(out, "seed")   <- seed
  out
}


#' Process residuals for a single deviation parameter
#' @keywords internal
.process_residual_one <- function(fit, process, par_name, seed) {
  obj <- fit$obj

  # ---- Posterior mode + covariance of the estimated deviations ----
  # Random effects: use the joint precision (the deviation block of its inverse
  # is the marginal posterior covariance). Penalized fixed effects: use the
  # fixed-effect covariance (inverse Hessian) directly.
  jp <- fit$sdrep$jointPrecision
  if (!is.null(jp) && par_name %in% rownames(jp)) {
    idx   <- which(rownames(jp) == par_name)
    lp    <- obj$env$last.par.best
    mode  <- lp[names(lp) == par_name]
    Sigma <- as.matrix(Matrix::solve(jp))[idx, idx, drop = FALSE]
  } else {
    sdr <- fit$sdrep
    if (is.null(sdr) || is.null(sdr$cov.fixed)) sdr <- TMB::sdreport(obj)
    cf <- sdr$cov.fixed
    if (is.null(cf) || !(par_name %in% colnames(cf))) {
      stop("'", par_name, "' is not estimated in this fit (no usable ",
           "covariance). Refit with the matching option (e.g. random_rec for ",
           "recruitment) so the deviations are estimated.")
    }
    cols  <- which(colnames(cf) == par_name)
    mode  <- sdr$par.fixed[names(sdr$par.fixed) == par_name]
    Sigma <- cf[cols, cols, drop = FALSE]
  }
  n <- length(mode)
  if (n == 0) stop("No estimated '", par_name, "' deviations.")

  # ---- Prior mean / SD and labels for the estimated elements ----
  spec <- .process_prior_spec(fit, process, par_name, n)
  if (length(spec$sd) != n) {
    stop("Could not align the prior for '", par_name, "' (", length(spec$sd),
         " priors vs ", n, " estimated deviations).")
  }

  # ---- One posterior draw, standardized by the prior ----
  L <- t(chol(Sigma))
  set.seed(seed)
  draw  <- as.numeric(mode + L %*% stats::rnorm(n))
  resid <- (draw - spec$mean) / spec$sd

  data.frame(
    type          = process,
    fleet         = spec$fleet,
    species       = spec$species,
    sex           = NA_integer_,
    year          = spec$year,
    age_or_length = spec$age,
    length        = NA_integer_,
    index_label   = NA_character_,
    observed      = NA_real_,
    predicted     = NA_real_,
    sd            = NA_real_,
    residual      = resid,
    stringsAsFactors = FALSE)
}


#' Prior mean/SD and labels for the estimated elements of a deviation parameter
#'
#' Maps the estimated (non-fixed) elements of `par_name` -- in the column-major
#' order used by TMB and the covariance matrices -- to their species/fleet,
#' year/age labels, and the process prior mean and SD.
#' @keywords internal
.process_prior_spec <- function(fit, process, par_name, n) {
  d    <- fit$data_list
  obj  <- fit$obj
  full <- obj$env$parList()[[par_name]]          # full parameter array
  mapf <- fit$map$mapFactor[[par_name]]
  est  <- if (is.null(mapf)) seq_along(full) else which(!is.na(as.integer(mapf)))
  if (length(est) != n) {
    # Fall back to the first n (e.g. when projection columns are simply absent
    # from the estimated set in a consistent leading block).
    est <- seq_len(n)
  }
  ai   <- arrayInd(est, dim(full))               # one row per estimated element
  styr <- d$styr

  if (process %in% c("recruitment", "initial")) {
    R_sd <- exp(as.numeric(obj$env$parList()$R_log_sd))   # per species
    sp   <- ai[, 1]
    list(species = sp, fleet = sp,
         year = if (process == "recruitment") styr + ai[, 2] - 1L else NA_integer_,
         age  = if (process == "initial") ai[, 2] + 1L else NA_integer_,
         mean = R_sd[sp]^2 / 2, sd = R_sd[sp])
  } else if (process == "catchability") {
    q_sd <- exp(as.numeric(obj$env$parList()$index_q_dev_log_sd))  # per fleet
    flt  <- ai[, 1]
    list(species = NA_integer_, fleet = flt,
         year = styr + ai[, 2] - 1L, age = NA_integer_,
         mean = rep(0, n), sd = q_sd[flt])
  } else {
    stop("Process '", process, "' is not supported.")
  }
}
