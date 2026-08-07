# -----------------------------------------------------------------------------
# Convergence diagnostics
# -----------------------------------------------------------------------------
# Each check emits one record -- list(id, tier, severity, message, data),
# severity in c("OK","NOTE","WARN","FAIL") -- so the battery is a single object
# (fit$convergence) rather than scattered messages. fit_mod() attaches the
# result and surfaces non-OK checks via message() (never errors).
# -----------------------------------------------------------------------------

# Severity ordering (low -> high). Used to roll individual records up into an
# overall status (the worst severity present).
.CONV_SEVERITY <- c("OK", "NOTE", "WARN", "FAIL")

# Build one check record.
.conv_record <- function(id, tier, severity, message, data = NULL) {
  severity <- match.arg(severity, .CONV_SEVERITY)
  list(id = id, tier = tier, severity = severity,
       message = message, data = data)
}

# Worst severity across a list of records ("OK" when empty).
.conv_overall <- function(checks) {
  if (length(checks) == 0) return("OK")
  sev <- vapply(checks, function(x) x$severity, character(1))
  .CONV_SEVERITY[max(match(sev, .CONV_SEVERITY))]
}

# Snapshot of the hindcast optimizer result, captured in fit_mod() before the
# projection re-optimization overwrites `opt`. `obj` lets us recompute the
# gradient and MLE vector directly -- essential in the hard-failure case, where
# opt$SD/$max_gradient/$diagnostics are often absent. `getsd` flags whether an
# sdreport was requested.
#
# `bounds` (the full build_bounds() list, $lower/$upper), `mapFactor`
# (map$mapFactor) and `random_vars` are used to assemble the per-parameter
# (mle, lower, upper) table.
.capture_opt_convergence <- function(opt, obj = NULL, bounds = NULL,
                                     mapFactor = NULL, random_vars = NULL,
                                     getsd = NA) {
  if (is.null(opt) && is.null(obj)) return(NULL)

  # Fixed-effect MLE vector (named), recomputed from the object.
  par_fixed <- NULL
  if (!is.null(obj)) {
    par_fixed <- tryCatch({
      par   <- obj$env$last.par.best
      ran   <- obj$env$random
      fixed <- if (length(ran)) setdiff(seq_along(par), ran) else seq_along(par)
      stats::setNames(as.numeric(par[fixed]), names(par)[fixed])
    }, error = function(e) NULL)
  }

  diag  <- opt$diagnostics
  worst <- NULL
  gg <- NULL
  if (!is.null(diag) && !is.null(diag$final_gradient)) {
    gg <- diag$final_gradient
    i <- which.max(abs(diag$final_gradient))
    if (length(i) == 1L) {
      worst <- list(param = as.character(diag$Param[i]),
                    gradient = as.numeric(diag$final_gradient[i]))
    }
  }
  mg <- opt$max_gradient
  if (is.null(mg) || !is.finite(mg)) {
    mg <- if (!is.null(worst)) abs(worst$gradient) else NA_real_
  }

  # Fallback: recompute the marginal gradient (and the worst parameter) from the
  # hindcast object when opt did not provide it.
  if ((is.na(mg) || is.null(worst)) && !is.null(par_fixed)) {
    gg <- tryCatch(stats::setNames(as.numeric(obj$gr(par_fixed)),
                                   names(par_fixed)),
                   error = function(e) NULL)
    if (!is.null(gg) && length(gg) > 0 && any(is.finite(gg))) {
      i <- which.max(abs(gg))
      worst <- list(param = names(gg)[i], gradient = unname(gg[i]))
      if (is.na(mg)) mg <- max(abs(gg))
    }
  }

  # Per-parameter (mle, lower, upper) for the bounds check. Build all three from
  # a single shared shape so they cannot drift out of alignment.
  bnd_par <- bnd_lo <- bnd_hi <- NULL
  if (!is.null(bounds) && !is.null(bounds$lower) && !is.null(bounds$upper) &&
      !is.null(mapFactor) && !is.null(obj)) {
    # parList(x = par[lfixed()], par = last.par): pass last.par.best as `par`
    # (the full vector) so `x` defaults to its fixed-effects subset. Passing it
    # positionally as `x` triggers `par[lfixed()] <- x` with a length mismatch.
    full <- tryCatch(obj$env$parList(par = obj$env$last.par.best),
                     error = function(e) NULL)
    if (!is.null(full)) {
      pv <- lv <- uv <- nv <- list()
      for (nm in names(mapFactor)) {
        if (nm %in% random_vars) next                 # random effects are unbounded
        if (is.null(full[[nm]]) || is.null(bounds$lower[[nm]])) next
        mf <- as.numeric(unlist(mapFactor[[nm]]))
        pb <- as.numeric(unlist(full[[nm]]))
        lb <- as.numeric(unlist(bounds$lower[[nm]]))
        ub <- as.numeric(unlist(bounds$upper[[nm]]))
        # one row per estimated coefficient (drop NA-mapped/fixed and duplicates)
        keep <- which(!is.na(mf) & !duplicated(mf))
        if (length(keep) == 0L) next
        pv[[nm]] <- pb[keep]; lv[[nm]] <- lb[keep]; uv[[nm]] <- ub[keep]
        nv[[nm]] <- rep(nm, length(keep))
      }
      if (length(pv)) {
        bnd_par <- stats::setNames(unlist(pv, use.names = FALSE),
                                   unlist(nv, use.names = FALSE))
        bnd_lo  <- unlist(lv, use.names = FALSE)
        bnd_hi  <- unlist(uv, use.names = FALSE)
      }
    }
  }

  list(
    par                = if (!is.null(bnd_par)) bnd_par else par_fixed,
    gradient           = gg,
    lower              = bnd_lo,
    upper              = bnd_hi,
    sd_requested       = isTRUE(getsd),
    sd_present         = !is.null(opt$SD),
    max_gradient       = as.numeric(mg),
    worst              = worst,
    nlminb_convergence = opt$convergence,
    message            = opt$message,
    convergence_check  = opt$Convergence_check,
    pdHess             = if (!is.null(opt$SD)) isTRUE(opt$SD$pdHess) else NA
  )
}


# --- checks ------------------------------------------------------------------

# Optimizer convergence: max |gradient| (+ the parameter carrying it) and
# Hessian positive-definiteness. Reads the hindcast snapshot so the result is
# not clobbered by the projection re-optimization.
.check_optimizer <- function(object) {
  ch <- object$.conv_hindcast
  out <- list()
  if (is.null(ch)) return(out)

  mg <- ch$max_gradient
  if (!is.null(mg) && is.finite(mg)) {
    sev <- if (mg > 1) "FAIL" else if (mg > 1e-3) "WARN" else "OK"
    worst_txt <- if (!is.null(ch$worst)) {
      sprintf(" (largest on '%s')", ch$worst$param)
    } else ""
    out$max_gradient <- .conv_record(
      "max_gradient", "fit", sev,
      sprintf("Maximum absolute marginal gradient = %.3g%s.", mg, worst_txt),
      ch)
  }

  if (!is.null(ch$pdHess) && !is.na(ch$pdHess)) {
    out$pdHess <- .conv_record(
      "pdHess", "fit",
      if (isTRUE(ch$pdHess)) "OK" else "FAIL",
      if (isTRUE(ch$pdHess)) "Hessian is positive definite."
      else "Hessian is not positive definite; standard errors are unavailable.")
  }
  out
}

# Hessian eigen-identifiability. Reads the fixed-effect covariance
# (sdrep$cov.fixed); its eigenvalues are 1/(Hessian eigenvalues), so the
# condition number kappa = lambda_max/lambda_min of cov.fixed equals the Hessian
# condition number, and the eigenvector of the LARGEST covariance eigenvalue is
# the least-identified parameter direction. Complements (does not duplicate)
# TMBhelper::check_estimability: that gives a per-parameter verdict; this gives
# the severity (kappa) and the offending linear *combination*.
.check_hessian_eigen <- function(object) {
  out <- list()
  cov <- tryCatch(object$sdrep$cov.fixed, error = function(e) NULL)
  if (is.null(cov) || !is.matrix(cov) || nrow(cov) < 2L) return(out)
  ev <- tryCatch(eigen(cov, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) return(out)
  vals <- ev$values
  pos  <- vals[is.finite(vals) & vals > 0]
  if (length(pos) < 2L) return(out)

  kappa <- max(pos) / min(pos)               # = condition number of the Hessian
  v  <- ev$vectors[, which.max(vals)]         # flattest Hessian direction (unit norm)
  nm <- rownames(cov)
  if (is.null(nm) || length(nm) != length(v)) nm <- names(object$sdrep$par.fixed)
  if (is.null(nm) || length(nm) != length(v)) nm <- paste0("p", seq_along(v))

  # The eigenvector is unit-norm, so v^2 partitions the direction across
  # coefficients. Aggregate by parameter block (a vector parameter shares one
  # name across its coefficients): each block's share is the sum of its squared
  # loadings, and the shares sum to 1. This names the block(s) the flat direction
  # lives in even when it is spread diffusely over many coefficients -- the case
  # a single-coefficient threshold missed, leaving the message blank.
  share <- tapply(v^2, factor(nm, levels = unique(nm)), sum)
  share <- sort(share, decreasing = TRUE)
  cum   <- cumsum(as.numeric(share))
  ntop  <- which(cum >= 0.90)[1]                          # blocks explaining >=90%
  if (is.na(ntop)) ntop <- length(share)
  ntop  <- max(1L, min(ntop, 5L))                         # always name >=1, cap at 5
  combo <- paste(sprintf("%s (%.0f%%)", names(share)[seq_len(ntop)],
                         100 * as.numeric(share)[seq_len(ntop)]),
                 collapse = " + ")
  top   <- data.frame(param = names(share), share = round(as.numeric(share), 3))

  sev <- if (kappa > 1e10) "FAIL" else if (kappa > 1e6) "WARN" else "OK"
  msg <- sprintf("Hessian condition number = %.2g.%s", kappa,
    if (sev != "OK")
      sprintf(" Least-identified direction loads on: %s.", combo) else "")
  out$hessian_conditioning <- .conv_record(
    "hessian_conditioning", "fit", sev, msg,
    list(condition_number = kappa, loadings = top))
  out
}

# sdreport failed: requested but did not return (Hessian not invertible). A
# strong non-convergence signal even when no gradient is available.
.check_sdreport_failed <- function(object) {
  ch <- object$.conv_hindcast
  if (is.null(ch) || !isTRUE(ch$sd_requested) || isTRUE(ch$sd_present)) {
    return(list())
  }
  list(sdreport_failed = .conv_record(
    "sdreport_failed", "fit", "FAIL",
    "sdreport failed (Hessian not invertible); standard errors unavailable."))
}

# Parameters at a configured bound. Optimization is unbounded in fit_mod(), so a
# parameter at/beyond its build_bounds() range means the MLE hit the edge of the
# plausible range -- often unidentified or mis-scaled.
.check_bounds <- function(object) {
  ch <- object$.conv_hindcast
  if (is.null(ch) || is.null(ch$par) || is.null(ch$lower) || is.null(ch$upper)) {
    return(list())
  }
  par <- ch$par; lo <- ch$lower; hi <- ch$upper
  rng <- hi - lo
  tol <- pmax(1e-6, 1e-3 * ifelse(is.finite(rng) & rng > 0, rng, 1))
  at_lo <- is.finite(lo) & par <= lo + tol
  at_hi <- is.finite(hi) & par >= hi - tol
  hit <- which((at_lo | at_hi) & par > -900)   # skip -999 sentinels (e.g. log_F)
  if (length(hit) == 0) return(list())
  tab <- data.frame(param = names(par)[hit], mle = signif(par[hit], 4),
                    lower = signif(lo[hit], 4), upper = signif(hi[hit], 4),
                    bound = ifelse(at_lo[hit], "lower", "upper"))
  list(parameters_on_bounds = .conv_record(
    "parameters_on_bounds", "fit", "WARN",
    sprintf("%d parameter(s) at a configured bound: %s.",
            length(hit), paste(unique(names(par)[hit]), collapse = ", ")),
    tab))
}

# Phasing log: phases that ended with a high gradient localize which parameter
# block breaks. Reads the per-phase log captured from TMBphase().
.check_phasing <- function(object) {
  pl <- object$.conv_phase
  if (is.null(pl) || length(pl) == 0) return(list())
  df <- tryCatch(do.call(rbind, lapply(pl, as.data.frame)),
                 error = function(e) NULL)
  if (is.null(df) || is.null(df$max_grad)) return(list())
  bad <- df[is.finite(df$max_grad) & df$max_grad > 1, , drop = FALSE]
  if (nrow(bad) == 0) return(list())
  list(phasing = .conv_record(
    "phasing", "fit", "WARN",
    sprintf("Phase(s) ended with max|grad| > 1: %s.",
            paste(sprintf("phase %s (%.2g)", bad$phase, bad$max_grad),
                  collapse = ", ")),
    df))
}

# Surface TMBhelper::check_estimability (run by fit_mod() and stored as
# object$identified when sdreport fails). Binary per-parameter verdict; pairs
# with .check_hessian_eigen()'s continuous view.
.check_estimability_record <- function(object) {
  out <- list()
  id <- object$identified
  if (is.null(id) || is.character(id)) return(out)
  bad <- tryCatch(id$WhichBad, error = function(e) NULL)
  if (is.null(bad)) {
    pc  <- tryCatch(id$BadParams$Param_check, error = function(e) NULL)
    bad <- if (!is.null(pc)) which(pc == 2L) else NULL
  }
  if (!is.null(bad) && length(bad) > 0) {
    nms <- tryCatch(unique(as.character(id$BadParams$Param[bad])),
                    error = function(e) NULL)
    out$estimability <- .conv_record(
      "estimability", "fit", "FAIL",
      sprintf("check_estimability flagged %d non-identifiable fixed parameter(s)%s.",
              length(bad),
              if (!is.null(nms)) paste0(" (", paste(nms, collapse = ", "), ")")
              else ""),
      id)
  }
  out
}

# Stock-recruit curve sanity. Beverton-Holt steepness is derived as
# h = alpha * SPR0 / (4 + alpha * SPR0), and alpha carries no lower bound, so the
# optimizer can reach alpha < 1/SPR0 -- below the replacement line. Steepness
# then falls under 0.2, the stock cannot replace itself, and the implied unfished
# recruitment (alpha - 1/SPR0)/beta is negative, which carries into the initial
# age structure.
.check_stock_recruit <- function(object) {
  dl <- object[["data_list"]]
  q  <- object[["quantities"]]
  if (is.null(dl) || is.null(q)) return(list())
  # Keyed on srr_fun, the curve the hindcast is actually fit with. Under the
  # Ianelli configuration (srr_fun = 0, srr_pred_fun > 0) recruitment is annual
  # deviates about mean R and the reported steepness / R0 describe that, not the
  # penalty curve, so there is nothing here to check.
  srr_fun <- suppressWarnings(as.integer(dl[["srr_fun"]])[1])
  if (is.na(srr_fun) || srr_fun < 2) return(list())

  # Both are [nspp, nyrs]; the stock-recruit curve is summarised by its first
  # year, so read column 1.
  col1 <- function(x) {
    if (is.null(x)) return(numeric(0))
    if (is.matrix(x)) as.numeric(x[, 1]) else as.numeric(x)
  }
  h  <- suppressWarnings(col1(q[["steepness"]]))
  r0 <- suppressWarnings(col1(q[["R0"]]))
  if (!length(h) && !length(r0)) return(list())

  bad_r0 <- which(is.finite(r0) & r0 <= 0)
  bad_h  <- if (srr_fun %in% c(2L, 3L)) which(is.finite(h) & h < 0.2) else integer(0)
  if (!length(bad_r0) && !length(bad_h)) {
    return(list(stock_recruit = .conv_record(
      "stock_recruit", "fit", "OK",
      sprintf("Stock-recruit curve is well posed (steepness %s).",
              paste(sprintf("%.3f", h), collapse = ", ")),
      list(steepness = h, R0 = r0))))
  }

  spp <- dl[["spnames"]]
  nm  <- function(i) if (length(spp) >= max(i)) paste(spp[i], collapse = ", ") else paste(i, collapse = ", ")
  msg <- character(0)
  if (length(bad_h)) {
    msg <- c(msg, sprintf(
      "Beverton-Holt steepness below 0.2 for %s (h = %s): alpha is under the replacement line 1/SPR0, so the stock cannot replace itself.",
      nm(bad_h), paste(sprintf("%.3f", h[bad_h]), collapse = ", ")))
  }
  if (length(bad_r0)) {
    msg <- c(msg, sprintf(
      "Implied unfished recruitment R0 is non-positive for %s (R0 = %s).",
      nm(bad_r0), paste(signif(r0[bad_r0], 4), collapse = ", ")))
  }
  msg <- c(msg, "Refit from stock-recruit starting values on the right scale -- see build_srr(srr_alpha_init =, srr_beta_init =).")

  list(stock_recruit = .conv_record(
    "stock_recruit", "fit", "FAIL", paste(msg, collapse = " "),
    list(steepness = h, R0 = r0)))
}



#' Convergence diagnostics for a fitted Rceattle model
#'
#' Runs the post-fit convergence battery and returns a single structured
#' object. Each check yields a record with a common schema
#' (\code{id}, \code{tier}, \code{severity}, \code{message}, \code{data});
#' \code{severity} is one of \code{"OK"}, \code{"NOTE"}, \code{"WARN"},
#' \code{"FAIL"}. The object's \code{status} is the worst severity present.
#'
#' \code{fit_mod()} runs this automatically and attaches the result as
#' \code{fit$convergence}; call \code{convergence_diagnostics()} directly to
#' re-run it on any fit. Checks cover the optimizer gradient, Hessian
#' positive-definiteness and conditioning, parameters on bounds, phasing, and
#' parameter estimability.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return An object of class \code{"Rceattle_convergence"}: a list with
#'   \code{status} (overall worst severity) and \code{checks} (named list of
#'   records).
#' @export
convergence_diagnostics <- function(object, ...) {
  checks <- c(
    .check_phasing(object),
    .check_optimizer(object),
    .check_sdreport_failed(object),
    .check_hessian_eigen(object),
    .check_bounds(object),
    .check_estimability_record(object),
    .check_stock_recruit(object)
  )
  structure(
    list(status = .conv_overall(checks), checks = checks),
    class = "Rceattle_convergence"
  )
}


#' @export
print.Rceattle_convergence <- function(x, all = FALSE, ...) {
  cat("<Rceattle convergence>  status:", x$status, "\n")
  if (length(x$checks) == 0) {
    cat("  (no checks run)\n")
    return(invisible(x))
  }
  show <- if (isTRUE(all)) x$checks else
    Filter(function(c) c$severity != "OK", x$checks)
  if (length(show) == 0) {
    cat("  all checks OK\n")
    return(invisible(x))
  }
  for (c in show) {
    cat(sprintf("  [%-4s] %-16s %s\n", c$severity, c$id, c$message))
  }
  invisible(x)
}
