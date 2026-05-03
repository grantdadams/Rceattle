#' Fit a TMB object with optional TMBhelper enhancements
#'
#' Internal wrapper used by [fit_mod()]. Delegates to
#' \code{TMBhelper::fit_tmb} when the suggested package is installed,
#' and otherwise performs an equivalent core fit using \code{stats::nlminb}
#' followed by \code{TMB::sdreport}. The fallback path covers `loopnum`
#' restarts, optional `newtonsteps`, and the `getsd` /
#' `getJointPrecision` / `getReportCovariance` reporting options.
#'
#' @keywords internal
#' @noRd
.fit_tmb <- function(obj, lower = -Inf, upper = Inf,
                     control = list(eval.max = 1e9, iter.max = 1e9, trace = 0),
                     loopnum = 1, newtonsteps = 0,
                     getsd = TRUE, bias.correct = FALSE,
                     getJointPrecision = FALSE, getReportCovariance = FALSE,
                     quiet = TRUE) {

  if (requireNamespace("TMBhelper", quietly = TRUE)) {
    return(
      TMBhelper::fit_tmb(
        obj                 = obj,
        fn                  = obj$fn,
        gr                  = obj$gr,
        startpar            = obj$par,
        lower               = lower,
        upper               = upper,
        loopnum             = loopnum,
        newtonsteps         = newtonsteps,
        getsd               = getsd,
        control             = control,
        bias.correct        = bias.correct,
        bias.correct.control = list(sd = getsd),
        getJointPrecision   = getJointPrecision,
        getReportCovariance = getReportCovariance,
        quiet               = quiet
      )
    )
  }

  # Fallback: plain nlminb + sdreport
  start_par <- obj$par
  best <- NULL
  for (i in seq_len(max(1, loopnum))) {
    init <- if (is.null(best)) start_par else best$par
    fit_i <- tryCatch(
      stats::nlminb(start = init, objective = obj$fn, gradient = obj$gr,
                    lower = lower, upper = upper, control = control),
      error = function(e) NULL
    )
    if (!is.null(fit_i) && (is.null(best) || fit_i$objective < best$objective)) {
      best <- fit_i
    }
  }
  if (is.null(best)) {
    stop("Optimization failed in .fit_tmb fallback (install TMBhelper for richer diagnostics).")
  }

  # Optional Newton steps to refine the gradient
  if (newtonsteps > 0) {
    par <- best$par
    for (i in seq_len(newtonsteps)) {
      g <- obj$gr(par)
      h <- tryCatch(stats::optimHess(par, fn = obj$fn, gr = obj$gr),
                    error = function(e) NULL)
      if (is.null(h)) break
      step <- tryCatch(solve(h, as.numeric(g)), error = function(e) NULL)
      if (is.null(step)) break
      par <- par - step
    }
    best$par <- par
    best$objective <- obj$fn(par)
  }

  best$diagnostics <- data.frame(
    Param = names(best$par),
    starting_value = start_par,
    Lower = if (length(lower) == 1) rep(lower, length(start_par)) else lower,
    MLE = best$par,
    Upper = if (length(upper) == 1) rep(upper, length(start_par)) else upper,
    final_gradient = as.numeric(obj$gr(best$par))
  )
  best$max_gradient <- max(abs(best$diagnostics$final_gradient))
  best$Convergence_check <- if (best$max_gradient < 1e-4) {
    "There is no evidence that the model is not converged"
  } else {
    "The maximum gradient is high; the model may not be converged"
  }
  best$number_of_coefficients <- c(
    Total = length(start_par),
    Fixed_effects = length(start_par) - length(obj$env$random),
    Random_effects = length(obj$env$random)
  )

  if (isTRUE(getsd)) {
    best$SD <- tryCatch(
      TMB::sdreport(obj,
                    bias.correct        = bias.correct,
                    bias.correct.control = list(sd = getsd),
                    getJointPrecision   = getJointPrecision,
                    getReportCovariance = getReportCovariance),
      error = function(e) NULL
    )
  }

  best
}


#' Check parameter estimability (TMBhelper::check_estimability fallback)
#'
#' @keywords internal
#' @noRd
.check_estimability <- function(obj) {
  if (requireNamespace("TMBhelper", quietly = TRUE)) {
    return(TMBhelper::check_estimability(obj))
  }
  # Minimal in-package fallback: test which fixed effects move the gradient
  par <- obj$env$last.par.best
  fixed <- setdiff(seq_along(par), obj$env$random)
  g <- as.numeric(obj$gr(par[fixed]))
  bad <- which(abs(g) > 1)
  list(
    BadParams = data.frame(
      Param = names(par[fixed]),
      MLE = unname(par[fixed]),
      Param_check = ifelse(seq_along(fixed) %in% bad, 2L, 1L),
      Final_gradient = g
    ),
    WhichBad = bad
  )
}
