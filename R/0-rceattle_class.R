#' Print method for fitted Rceattle models
#'
#' Provides a compact summary so that auto-printing inside R Markdown /
#' knitr / RStudio does not recurse into the (very deep) data and TMB
#' objects stored on the fit. Only structural metadata, convergence
#' status, and headline derived quantities are printed.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return The input `x`, invisibly.
#' @export
print.Rceattle <- function(x, ...) {
  dat <- x$data_list
  cat("<Rceattle model>\n")
  cat("  Species   :", paste(dat$spnames, collapse = ", "), "\n")
  cat("  Years     :", dat$styr, "-", dat$endyr,
      if (!is.null(dat$projyr)) paste0(" (projyr ", dat$projyr, ")") else "",
      "\n")
  cat("  msmMode   :", dat$msmMode, "\n")
  cat("  HCR       :", if (is.null(dat$HCR)) "(none)" else dat$HCR, "\n")
  cat("  initMode  :", if (is.null(dat$initMode)) NA else dat$initMode, "\n")

  if (!is.null(x$opt) && !is.null(x$opt$objective)) {
    cat("  -log L    :", signif(x$opt$objective, 6), "\n")
  }
  if (!is.null(x$opt$max_gradient)) {
    cat("  max |grad|:", signif(x$opt$max_gradient, 4), "\n")
  } else if (!is.null(x$obj)) {
    g <- tryCatch(max(abs(as.numeric(x$obj$gr()))), error = function(e) NA_real_)
    if (is.finite(g)) cat("  max |grad|:", signif(g, 4), "\n")
  }
  if (!is.null(x$run_time)) {
    cat("  Run time  :", format(x$run_time), "\n")
  }
  invisible(x)
}


#' Compact summary method for Rceattle fits
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @export
summary.Rceattle <- function(object, ...) {
  print(object, ...)
  invisible(object)
}


# Save and restore graphics par() in the caller frame.
#
# Plotting functions in Rceattle modify graphics state with par(...) but
# do not own the user's device. Calling .save_par() at the top of each
# plot function snapshots par() and registers an on.exit() handler in
# the caller so the user's device is restored even on error.
#
# Implementation note: do.call("on.exit", ..., envir = parent.frame())
# attaches the on.exit hook to the caller, not to .save_par() itself.
.save_par <- function() {
  oldpar <- graphics::par(no.readonly = TRUE)
  do.call(
    "on.exit",
    list(substitute(graphics::par(oldpar)), add = TRUE),
    envir = parent.frame()
  )
  invisible(oldpar)
}
