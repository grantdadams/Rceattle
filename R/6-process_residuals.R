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
#'   \item{`"recruitment"`}{`rec_dev`, prior `N(-bias_adjust_proc * R_sd^2/2, R_sd^2)` per species.}
#'   \item{`"initial"`}{`init_dev`, prior `N(-bias_adjust_proc * R_sd^2/2, R_sd^2)` per species.}
#'   \item{`"catchability"`}{`index_q_dev`, prior `N(0, q_dev_sd^2)` per index.}
#' }
#' `"all"` returns every supported process present in the fit. Selectivity and
#' natural-mortality deviations (which use random-walk / 2D-AR1 priors) are not
#' yet supported. Recruitment and initial-abundance residuals are exact (their
#' priors are iid). Catchability residuals are exact only for the iid deviate
#' prior (`Time_varying_q = 1` or `2`); for a random-walk or AR1 catchability
#' prior the marginal-SD standardization ignores the prior correlation, so those
#' residuals are approximate and a warning is emitted.
#'
#' @section On a DSEM fit:
#' A DSEM makes the recruitment deviations the latent states of a Gaussian
#' Markov random field, so there is no per-year normal prior to divide by.
#' `"recruitment"` instead whitens the states with the process precision:
#' conditioning on the states the model was given -- the observed covariates --
#' the estimated states have precision \eqn{Q_{EE}}, and
#' \eqn{r = U (x_E - m)} with \eqn{Q_{EE} = U^{\top} U} has covariance
#' \eqn{I}. So the claim is the same as on the standard path, and exact in the
#' same sense. One residual is returned per hindcast year per species; a
#' covariate's own latent states are not residualized here.
#'
#' `"initial"` standardizes by the recruitment SD the density actually used,
#' which under a DSEM comes from the SEM's variance path rather than from
#' `R_log_sd` (mapped out, and a placeholder). `"catchability"` is untouched by
#' a DSEM.
#'
#' @param object A fitted `Rceattle` model. The targeted deviations must be estimated -- as random effects (e.g. `random_rec = TRUE`) or as penalized fixed effects -- with a usable covariance.
#' @param fit deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @param process One of `"recruitment"`, `"initial"`, `"catchability"`, or
#'   `"all"`.
#' @param seed Seed for the posterior draw. Default 123.
#'
#' @return A data frame of class `rceattle_osa` (so it can be passed to
#'   [osa_diagnostics()] and [plot.rceattle_osa()]) with one row per deviation:
#'   columns `source` (the process), `fleet`, `species`, `year`,
#'   `age_length_bin`, and `residual`.
#'
#' @references
#' Nielsen, A., and Berg, C.W. 2014. Estimation of time-varying selectivity in
#'   stock assessments using state-space models. Fish. Res. 158:96-101.
#'
#' @seealso [osa_residuals()], [osa_diagnostics()]
#' @export
process_residuals <- function(object = NULL,
                              process = c("recruitment", "initial",
                                          "catchability", "all"),
                              seed = 123, fit = NULL) {
  # `fit` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(fit))
    object <- .rce_deprecated_arg(fit, !missing(object), "fit", "object", "process_residuals")

  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).")
  }
  process <- match.arg(process)

  if (is.null(object$obj)) stop("`object` has no TMB object ($obj).")

  # Under a DSEM the recruitment deviations are the latent states of a GMRF, so
  # `rec_dev` is mapped out and there is no per-year normal prior to standardize
  # by. The residual is defined all the same: whiten the states with the
  # process precision instead of dividing by a scalar SD (.process_residual_dsem).
  # "initial" and "catchability" take the ordinary path -- init_dev is still
  # estimated, and a DSEM does not touch the catchability deviates.
  dsem_rec <- .has_dsem(object)

  specs <- list(recruitment  = "rec_dev",
                initial      = "init_dev",
                catchability = "index_q_dev")
  todo <- if (process == "all") names(specs) else process

  # Seed once here -- saving and restoring the caller's RNG -- and let the
  # processes draw from a single advancing stream. Seeding inside the per-process
  # loop instead would make equal-length deviation blocks (e.g. recruitment vs
  # initial) draw identical standard normals (spurious cross-process correlation)
  # and would clobber the caller's global .Random.seed.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  out <- lapply(todo, function(p) {
    one <- if (dsem_rec && p == "recruitment") {
      function() .process_residual_dsem(object)
    } else {
      function() .process_residual_one(object, p, specs[[p]])
    }
    tryCatch(one(), error = function(e) if (process == "all") NULL else stop(e))
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
.process_residual_one <- function(fit, process, par_name) {
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
  # The RNG is seeded once by the caller (process_residuals()); drawing from that
  # advancing stream keeps each process's normals distinct and reproducible.
  L <- t(chol(Sigma))
  draw  <- as.numeric(mode + L %*% stats::rnorm(n))
  resid <- (draw - spec$mean) / spec$sd

  data.frame(
    source         = process,
    fleet          = spec$fleet,
    species        = spec$species,
    sex            = NA_integer_,
    year           = spec$year,
    age_length_bin = spec$age,
    length         = NA_integer_,
    index_label    = NA_character_,
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
  # Estimated (non-fixed) elements, in column-major order -- the same order as
  # the covariance matrix rows. If this does not match the covariance dimension
  # the element-to-(species/year/age) mapping would be guesswork, so error
  # rather than mislabel.
  est <- if (is.null(mapf)) seq_along(full) else which(!is.na(as.integer(mapf)))
  if (length(est) != n) {
    stop("Could not map the ", n, " estimated '", par_name, "' deviations to ",
         "model elements (found ", length(est), " estimated positions). ",
         "Process residuals are not supported for this parameter's ",
         "configuration.")
  }
  ai <- arrayInd(est, dim(full))                 # one row per estimated element
  styr <- d$styr

  if (process %in% c("recruitment", "initial")) {
    # The recruitment SD the DENSITY used. Under a DSEM the template sets it
    # from the SEM's variance path and maps R_log_sd out, so exp(R_log_sd) is a
    # placeholder: on BS2017SS it reads 0.7071 against the 0.9888 / 1.3321 /
    # 0.7976 actually in force. Take the reported R_sd wherever it is available
    # and fall back to the parameter only if it is not.
    R_sd <- suppressWarnings(as.numeric(fit$quantities$R_sd))
    if (length(R_sd) != d$nspp || anyNA(R_sd))
      R_sd <- exp(as.numeric(obj$env$parList()$R_log_sd))   # per species
    # rec_dev / init_dev prior mean is the lognormal bias correction
    # -bias_adjust_proc * sigma^2/2 (matches ceattle.cpp slots 9 & 10).
    ba   <- obj$env$data$bias_adjust_proc
    if (is.null(ba)) ba <- 1
    sp   <- ai[, 1]
    list(species = sp, fleet = sp,
         year = if (process == "recruitment") styr + ai[, 2] - 1L else NA_integer_,
         age  = if (process == "initial") ai[, 2] + 1L else NA_integer_,
         mean = -as.numeric(ba) * R_sd[sp]^2 / 2, sd = R_sd[sp])
  } else if (process == "catchability") {
    q_sd <- exp(as.numeric(obj$env$parList()$index_q_dev_log_sd))  # per fleet
    flt  <- ai[, 1]
    # The (mean 0, marginal SD) standardization below is exact only for the iid
    # catchability deviate prior (Time_varying_q = 1 or 2). Warn when an involved
    # fleet uses a correlated prior -- random walk (index_varying_q == 4) or AR1
    # (est_index_q == 6) -- because ignoring the prior correlation makes those
    # residuals approximate.
    ivq  <- obj$env$data$index_varying_q
    eqd  <- obj$env$data$est_index_q
    uflt <- unique(flt)
    corr <- rep(FALSE, length(uflt))
    if (!is.null(ivq)) corr <- corr | (ivq[uflt] %in% 4L)
    if (!is.null(eqd)) corr <- corr | (eqd[uflt] %in% 6L)
    if (any(corr, na.rm = TRUE)) {
      warning("process = 'catchability': index fleet(s) ",
              paste(uflt[which(corr)], collapse = ", "),
              " use a correlated catchability deviate prior (random walk or AR1); ",
              "their residuals standardize by the marginal SD and are therefore ",
              "approximate. They are exact only for the iid prior ",
              "(Time_varying_q = 1 or 2).", call. = FALSE)
    }
    list(species = NA_integer_, fleet = flt,
         year = styr + ai[, 2] - 1L, age = NA_integer_,
         mean = rep(0, n), sd = q_sd[flt])
  } else {
    stop("Process '", process, "' is not supported.")
  }
}


#' Process residuals for the recruitment deviations of a DSEM
#'
#' @description
#' The DSEM makes the recruitment deviations the latent states of a GMRF, so
#' they are not independent and cannot be standardized by a per-year SD. The
#' generalization is the same one SAM's `procres()` makes for a scalar prior,
#' written for a correlated process: take one draw from the joint posterior of
#' the latent states and whiten it with the process precision.
#'
#' Conditioning on the states the model was GIVEN -- the observed covariates,
#' which are fixed in the map -- the estimated states have precision
#' \eqn{Q_{EE}} and mean \eqn{\mu_E - Q_{EE}^{-1} Q_{EF} (x_F - \mu_F)}. With
#' \eqn{Q_{EE} = U^{\top} U} from a Cholesky factor, \eqn{r = U (x_E - m)} has
#' covariance \eqn{U Q_{EE}^{-1} U^{\top} = I}, so under a correctly specified
#' SEM the residuals are approximately iid N(0, 1) -- the same claim the scalar
#' path makes, and exact in the same sense.
#'
#' @param fit a fitted `Rceattle` model carrying a DSEM.
#' @return A data.frame in the `rceattle_osa` schema, one row per estimated
#'   latent state.
#' @keywords internal
#' @noRd
.process_residual_dsem <- function(fit) {
  q  <- fit$quantities
  Q  <- q$dsem_Q
  if (is.null(Q) || !length(Q)) {
    stop("This fit reports no DSEM precision (`dsem_Q`), so its recruitment ",
         "deviations cannot be whitened. Refit with the current version.",
         call. = FALSE)
  }
  Q  <- as.matrix(Q)
  X  <- as.matrix(fit$estimated_params$dsem_x_tj)
  mu <- as.matrix(q$dsem_xhat_tj) + as.matrix(q$dsem_delta_tj)
  n_t <- nrow(X); n_j <- ncol(X); n_k <- n_t * n_j
  if (nrow(Q) != n_k) {
    stop("The reported DSEM precision is ", nrow(Q), "x", ncol(Q),
         " but the latent field is ", n_t, "x", n_j, "; they cannot be ",
         "aligned. Usually this means a variable in the sem has no exogenous ",
         "variance -- check_dsem_spec() names it.", call. = FALSE)
  }

  # The RECRUITMENT-deviation columns only, conditional on every other column at
  # its fitted value. Selecting by the map instead would also return one residual
  # per unobserved covariate year -- dsem does not map an observed covariate cell
  # off, so the map is not the latent/observed boundary here -- and those rows
  # would be labelled "recruitment" and pooled by osa_diagnostics(). A covariate
  # state is a different diagnostic.
  rc <- fit$dsem$tmb_inputs$data$rec_dev_col + 1L   # 0-based in the C++
  if (is.null(rc) || !length(rc)) {
    stop("This fit does not record which latent columns carry the recruitment ",
         "deviations, so its process residuals cannot be computed. Refit ",
         "before residualizing.", call. = FALSE)
  }
  # HINDCAST years only. build_DSEM(estimate_projection = TRUE) -- which
  # sample_rec() requires, so it is not an exotic setting -- runs the latent
  # field to projyr, and those states are scored by no data. Whitening them
  # returns draws from the prior, which are N(0, 1) by construction, and
  # osa_diagnostics() pools them into the SDNR: on BS2017SS with a 5-year
  # projection that is 5 of every 44 rows per species, pulling the SDNR toward 1
  # for a reason that has nothing to do with fit. They stay in the CONDITIONING
  # set, the same treatment an unobserved covariate cell already gets here.
  n_hind <- fit$data_list$endyr - fit$data_list$styr + 1L
  mp  <- fit$map$mapList$dsem_x_tj
  est <- if (is.null(mp)) rep(TRUE, n_k) else !is.na(as.numeric(mp))
  in_rec <- rep(FALSE, n_k)
  in_rec[unlist(lapply(rc, function(j)
    (j - 1L) * n_t + seq_len(min(n_t, n_hind))))] <- TRUE
  E   <- which(est & in_rec); F_ <- setdiff(seq_len(n_k), E)
  if (!length(E)) stop("Every recruitment-deviation state is fixed; nothing to ",
                       "residualize.", call. = FALSE)

  # Conditional mean of the estimated states given the ones the model was given.
  m <- as.numeric(mu)[E]
  if (length(F_)) {
    m <- m - solve(Q[E, E, drop = FALSE],
                   Q[E, F_, drop = FALSE] %*%
                     (as.numeric(X)[F_] - as.numeric(mu)[F_]))
  }

  # One draw from the joint posterior of the states, as the scalar path does.
  # The `dsem_x_tj` rows of the joint precision are the MAP-ESTIMATED cells in
  # column-major order, which is a superset of the recruitment cells, so the
  # draw is taken over that block and then subset to E.
  jp   <- fit$sdrep$jointPrecision
  keep <- match(E, which(est))
  if (!is.null(jp) && "dsem_x_tj" %in% rownames(jp)) {
    idx   <- which(rownames(jp) == "dsem_x_tj")
    lp    <- fit$obj$env$last.par.best
    mode  <- lp[names(lp) == "dsem_x_tj"]
    if (length(mode) != sum(est)) {
      stop("The DSEM latent block has ", length(mode), " entries in the ",
           "parameter vector but the map says ", sum(est), "; they cannot be ",
           "aligned.", call. = FALSE)
    }
    Sigma <- as.matrix(Matrix::solve(jp))[idx, idx, drop = FALSE]
    draw  <- as.numeric(mode + t(chol(Sigma)) %*% stats::rnorm(length(mode)))[keep]
  } else {
    # No joint precision (getsd = FALSE): the posterior mode is the best
    # available point, so the residuals are the mode's rather than a draw's.
    # Say so -- shrunk residuals look better than they are.
    warning("This fit carries no joint precision, so the DSEM process ",
            "residuals use the posterior MODE rather than a draw. Mode ",
            "residuals are shrunk toward the process and look better behaved ",
            "than they are. Refit with fit_control(getsd = TRUE) for the ",
            "one-sample residuals this function is meant to return.",
            call. = FALSE)
    draw <- as.numeric(X)[E]
  }

  resid <- as.numeric(chol(Q[E, E, drop = FALSE]) %*% (draw - m))

  # Label each cell. k = (j - 1) * n_t + t, column-major over x_tj.
  t_of <- ((E - 1L) %% n_t) + 1L
  j_of <- ((E - 1L) %/% n_t) + 1L
  rc   <- fit$dsem$tmb_inputs$data$rec_dev_col + 1L   # 0-based in the C++
  sp   <- match(j_of, rc)                             # NA for a covariate column
  vnm  <- colnames(X)
  data.frame(
    source         = "recruitment",
    fleet          = sp,
    species        = sp,
    sex            = NA_integer_,
    year           = fit$data_list$styr + t_of - 1L,
    age_length_bin = NA_integer_,
    length         = NA_integer_,
    index_label    = if (is.null(vnm)) NA_character_ else vnm[j_of],
    observed       = NA_real_,
    predicted      = NA_real_,
    sd             = NA_real_,
    residual       = resid,
    stringsAsFactors = FALSE)
}
