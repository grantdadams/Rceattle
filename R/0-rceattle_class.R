#' Print method for fitted Rceattle models
#'
#' Provides a compact summary so that auto-printing inside R Markdown /
#' knitr / RStudio does not recurse into the (very deep) data and TMB
#' objects stored on the fit. Only structural metadata, convergence
#' status, headline derived quantities, and the package / TMB-DLL
#' versions used to produce the fit are printed.
#'
#' For operational use, the package version line is meant to make it
#' obvious which version of `Rceattle` produced an archived fit so that
#' results can be reproduced even if `master` has moved on. Tag a
#' release (`devtools::install_github("grantdadams/Rceattle@vX.Y.Z")`)
#' and the same version string will reappear here on a fresh run.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return The input `x`, invisibly.
#' @export
print.Rceattle <- function(x, ...) {
  dat <- x$data_list
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("Rceattle")),
    error = function(e) NA_character_
  )
  tmb_dll <- if (!is.null(x$TMBfilename)) x$TMBfilename else "ceattle_v01_11"

  cat("<Rceattle model>\n")
  cat("  Rceattle  :", pkg_ver, "  TMB DLL:", tmb_dll, "\n")
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


#' Plot method for fitted Rceattle models
#'
#' Thin S3 dispatcher around the package's existing `plot_*()` functions
#' so that `plot(fit)` works the way users expect. Pick the panel with
#' `what`; everything in `...` is forwarded to the underlying function.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param what Character. One of `"biomass"` (default), `"ssb"`,
#'   `"recruitment"`, `"depletion"`, `"index"`, `"catch"`,
#'   `"selectivity"`, `"mortality"`, or `"data"`.
#' @param ... Passed to the underlying plotting function.
#'
#' @return Invisibly returns `NULL`. Called for the side effect of
#'   producing a plot.
#' @export
plot.Rceattle <- function(x, what = "biomass", ...) {
  what <- match.arg(
    what,
    c("biomass", "ssb", "recruitment", "depletion",
      "index", "catch", "selectivity", "mortality", "data")
  )
  switch(
    what,
    biomass      = plot_biomass(x, ...),
    ssb          = plot_ssb(x, ...),
    recruitment  = plot_recruitment(x, ...),
    depletion    = plot_depletion(x, ...),
    index        = plot_index(x, ...),
    catch        = plot_catch(x, ...),
    selectivity  = plot_selectivity(x, ...),
    mortality    = plot_mortality(x, ...),
    data         = plot_data(x, ...)
  )
  invisible(NULL)
}


#' Extract estimated parameters from an Rceattle fit
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return A named numeric vector of estimated fixed-effect parameters
#'   (the optimizer's `par`). Returns `NULL` if the model was not
#'   estimated.
#' @export
coef.Rceattle <- function(object, ...) {
  if (is.null(object$opt) || is.null(object$opt$par)) return(NULL)
  object$opt$par
}


#' Variance-covariance matrix for an Rceattle fit
#'
#' Returns the fixed-effect covariance matrix produced by
#' [TMB::sdreport()]. Random-effect covariance is not returned here —
#' use `object$sdrep` for the full report.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return A numeric matrix, or `NULL` if `sdreport` was not run
#'   (i.e. the fit was produced with `getsd = FALSE`).
#' @export
vcov.Rceattle <- function(object, ...) {
  sdr <- object$sdrep
  if (is.null(sdr) || is.null(sdr$cov.fixed)) return(NULL)
  sdr$cov.fixed
}


#' Log-likelihood of an Rceattle fit
#'
#' Returns the joint negative-objective from `nlminb`, sign-flipped, as
#' a `"logLik"` object. `df` is the number of estimated fixed-effect
#' parameters. The `nobs` attribute is intentionally omitted: counting
#' "observations" in a stock assessment likelihood (with composition
#' cells, indices, catches, and priors) is not well-defined, so
#' [stats::AIC()] works (uses `df`) while [stats::BIC()] does not.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return An object of class `"logLik"`, or `NULL` if the model was
#'   not estimated.
#' @export
logLik.Rceattle <- function(object, ...) {
  if (is.null(object$opt) || is.null(object$opt$objective)) return(NULL)
  ll <- -as.numeric(object$opt$objective)
  attr(ll, "df") <- length(object$opt$par)
  class(ll) <- "logLik"
  ll
}


#' Observed-vs-fitted residuals from an Rceattle fit
#'
#' Returns a long-format data frame of residuals across one or more of
#' the four fitted data sources: `"index"` (survey indices),
#' `"catch"` (fishery catches), `"comp"` (age- or length-composition
#' proportions), and `"caal"` (conditional age-at-length proportions).
#'
#' For `"index"` and `"catch"`, the `Residual` column is on the log
#' scale by default (matching the lognormal observation likelihood) and
#' can be switched to the natural scale via `scale = "natural"`. For
#' `"comp"` and `"caal"`, residuals are Pearson residuals on the
#' fitted proportions:
#' \deqn{r = (p - \hat p)/\sqrt{\hat p (1 - \hat p)/N}}
#' where N is the input sample size. Composition rows are returned in
#' long form: one row per (observation, age/length bin).
#'
#' Composition rows carry the `Age0_Length1` flag from `comp_data`
#' (`0` for age comps, `1` for length comps) so age and length comps
#' can be filtered apart. CAAL rows carry both the conditioning
#' `Length` and the age `Bin`.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param type One or more of `"index"`, `"catch"`, `"comp"`, `"caal"`,
#'   or `"all"` (default `"index"`).
#' @param scale `"log"` (default) or `"natural"`. Only affects
#'   `"index"` and `"catch"` residuals.
#' @param ... Currently unused.
#'
#' @return A `data.frame` with columns `Source`, `Fleet_code`,
#'   `Fleet_name`, `Species`, `Sex`, `Year`, `Length`, `Bin`,
#'   `Age0_Length1`, `Sample_size`, `Observed`, `Fitted`, `Residual`.
#'   `Sex`, `Length`, `Bin`, `Age0_Length1`, and `Sample_size` are
#'   `NA` where they do not apply (e.g. for index/catch rows).
#' @export
residuals.Rceattle <- function(object, type = "index", scale = "log", ...) {
  valid <- c("index", "catch", "comp", "caal", "all")
  type  <- match.arg(type, valid, several.ok = TRUE)
  if ("all" %in% type) type <- c("index", "catch", "comp", "caal")
  scale <- match.arg(scale, c("log", "natural"))

  q <- object$quantities
  d <- object$data_list

  empty_row <- function(n) {
    data.frame(
      Source       = character(n),
      Fleet_code   = integer(n),
      Fleet_name   = character(n),
      Species      = integer(n),
      Sex          = rep(NA_integer_, n),
      Year         = integer(n),
      Length       = rep(NA_real_, n),
      Bin          = rep(NA_real_, n),
      Age0_Length1 = rep(NA_integer_, n),
      Sample_size  = rep(NA_real_, n),
      Observed     = numeric(n),
      Fitted       = numeric(n),
      Residual     = numeric(n),
      stringsAsFactors = FALSE
    )
  }

  out <- list()

  if ("index" %in% type &&
      !is.null(d$index_data) && !is.null(q$index_hat)) {
    idx <- d$index_data
    obs <- idx$Observation
    hat <- as.numeric(q$index_hat)
    res <- if (scale == "log") log(obs) - log(hat) else obs - hat
    df <- empty_row(nrow(idx))
    df$Source     <- "index"
    df$Fleet_code <- idx$Fleet_code
    df$Fleet_name <- idx$Fleet_name
    df$Species    <- idx$Species
    df$Year       <- idx$Year
    df$Observed   <- obs
    df$Fitted     <- hat
    df$Residual   <- res
    out$index <- df
  }

  if ("catch" %in% type &&
      !is.null(d$catch_data) && !is.null(q$catch_hat)) {
    ctc <- d$catch_data
    obs <- ctc$Catch
    hat <- as.numeric(q$catch_hat)
    res <- if (scale == "log") log(obs) - log(hat) else obs - hat
    df <- empty_row(nrow(ctc))
    df$Source     <- "catch"
    df$Fleet_code <- ctc$Fleet_code
    df$Fleet_name <- ctc$Fleet_name
    df$Species    <- ctc$Species
    df$Year       <- ctc$Year
    df$Observed   <- obs
    df$Fitted     <- hat
    df$Residual   <- res
    out$catch <- df
  }

  if ("comp" %in% type &&
      !is.null(d$comp_data) && !is.null(q$comp_hat) &&
      nrow(d$comp_data) > 0) {
    cd <- d$comp_data
    bin_cols <- grep("^Comp_", colnames(cd), value = TRUE)
    obs_mat <- as.matrix(cd[, bin_cols, drop = FALSE])
    row_tot <- rowSums(obs_mat, na.rm = TRUE)
    row_tot[row_tot == 0] <- NA_real_
    obs_prop <- obs_mat / row_tot
    hat_prop <- as.matrix(q$comp_hat)
    bin_vals <- suppressWarnings(as.numeric(sub("^Comp_", "", bin_cols)))
    n_obs <- nrow(cd); n_bin <- length(bin_cols)
    a0l1 <- if (!is.null(cd$Age0_Length1)) cd$Age0_Length1 else NA_integer_
    df <- empty_row(n_obs * n_bin)
    df$Source       <- "comp"
    df$Fleet_code   <- rep(cd$Fleet_code, times = n_bin)
    df$Fleet_name   <- rep(cd$Fleet_name, times = n_bin)
    df$Species      <- rep(cd$Species,    times = n_bin)
    df$Sex          <- rep(if (!is.null(cd$Sex)) cd$Sex else NA_integer_,
                           times = n_bin)
    df$Year         <- rep(cd$Year,       times = n_bin)
    df$Bin          <- rep(bin_vals, each = n_obs)
    df$Age0_Length1 <- rep(a0l1, times = n_bin)
    df$Sample_size  <- rep(cd$Sample_size, times = n_bin)
    df$Observed     <- as.numeric(obs_prop)
    df$Fitted       <- as.numeric(hat_prop)
    df$Residual     <- (df$Observed - df$Fitted) /
      sqrt(df$Fitted * (1 - df$Fitted) / df$Sample_size)
    df <- df[!is.na(df$Observed) & !is.na(df$Fitted), , drop = FALSE]
    out$comp <- df
  }

  if ("caal" %in% type &&
      !is.null(d$caal_data) && !is.null(q$caal_hat) &&
      nrow(d$caal_data) > 0) {
    cd <- d$caal_data
    bin_cols <- grep("^CAAL_", colnames(cd), value = TRUE)
    obs_mat <- as.matrix(cd[, bin_cols, drop = FALSE])
    row_tot <- rowSums(obs_mat, na.rm = TRUE)
    row_tot[row_tot == 0] <- NA_real_
    obs_prop <- obs_mat / row_tot
    hat_prop <- as.matrix(q$caal_hat)
    if (ncol(hat_prop) >= length(bin_cols)) {
      hat_prop <- hat_prop[, seq_along(bin_cols), drop = FALSE]
    }
    bin_vals <- suppressWarnings(as.numeric(sub("^CAAL_", "", bin_cols)))
    n_obs <- nrow(cd); n_bin <- length(bin_cols)
    df <- empty_row(n_obs * n_bin)
    df$Source      <- "caal"
    df$Fleet_code  <- rep(cd$Fleet_code, times = n_bin)
    df$Fleet_name  <- rep(cd$Fleet_name, times = n_bin)
    df$Species     <- rep(cd$Species,    times = n_bin)
    df$Sex         <- rep(if (!is.null(cd$Sex)) cd$Sex else NA_integer_,
                          times = n_bin)
    df$Year        <- rep(cd$Year,       times = n_bin)
    df$Length      <- rep(if (!is.null(cd$Length)) cd$Length else NA_real_,
                          times = n_bin)
    df$Bin         <- rep(bin_vals, each = n_obs)
    df$Sample_size <- rep(cd$Sample_size, times = n_bin)
    df$Observed    <- as.numeric(obs_prop)
    df$Fitted      <- as.numeric(hat_prop)
    df$Residual    <- (df$Observed - df$Fitted) /
      sqrt(df$Fitted * (1 - df$Fitted) / df$Sample_size)
    df <- df[!is.na(df$Observed) & !is.na(df$Fitted), , drop = FALSE]
    out$caal <- df
  }

  if (length(out) == 0) return(empty_row(0))
  do.call(rbind, out)
}
