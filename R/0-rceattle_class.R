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

  cat("<Rceattle model>\n")
  cat("  Rceattle  :", pkg_ver, "\n")

  # Model specification as an indented spec tree (dimensions -> fleets ->
  # processes -> linkages -> config), shared with print.Rceattle_data(). Guarded
  # so a rendering edge case can never break auto-printing a fitted model.
  tree <- tryCatch(.rce_spec_tree(dat),
                   error = function(e) paste0("  (spec tree unavailable: ",
                                              conditionMessage(e), ")"))
  cat(paste(tree, collapse = "\n"), "\n", sep = "")

  cat("  fit\n")
  # Alias the integer switch codes to their string forms so the fit block reads the
  # same regardless of the underlying code (matching the spec tree above).
  cat("  \u251c\u2500 msmMode  :", if (is.null(dat$msmMode)) NA else .rce_alias_show(dat$msmMode, msmMode_map), "\n")
  cat("  \u251c\u2500 HCR      :", if (is.null(dat$HCR)) "(none)" else .rce_alias_show(dat$HCR, hcr_map), "\n")
  cat("  \u251c\u2500 initMode :", if (is.null(dat$initMode)) NA else .rce_alias_show(dat$initMode, initMode_map), "\n")

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
  if (!is.null(x$convergence)) {
    cat("  Converged :", x$convergence$status,
        if (x$convergence$status != "OK")
          "(see fit$convergence)" else "", "\n")
  }
  invisible(x)
}


#' Compact summary method for Rceattle fits
#'
#' Prints the compact model overview (via [print.Rceattle()]). When the fit
#' carries a dynamic structural equation model (DSEM) on recruitment deviations,
#' also prints and returns the estimated SEM path coefficients — with standard
#' errors, z-values and Wald p-values when an `sdreport` is available — plus the
#' per-species recruitment SD. Adapted from `summary.dsem`; see
#' `?dsem::summary.dsem` for the underlying parameterization.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Passed to [print.Rceattle()].
#'
#' @return Invisibly, either:
#' \itemize{
#'   \item a list with elements \code{coefficients} and \code{recruitment_sd},
#'     when the fit contains a DSEM. \code{coefficients} is a data.frame of path
#'     coefficients with columns \code{path}, \code{lag}, \code{name},
#'     \code{start}, \code{parameter}, \code{first}, \code{second},
#'     \code{direction}, \code{Estimate}, \code{Std_Error}, \code{z_value} and
#'     \code{p_value}. The last three are \code{NA} when the fit carries no
#'     \code{sdreport}, but the columns are always present so that tables from
#'     several fits can be row-bound. \code{recruitment_sd} has one row per
#'     species with columns \code{Species}, \code{R_sd}, \code{Estimated}
#'     (whether it was estimated rather than fixed) and \code{Std_Error}
#'     (\code{NA} for a fixed SD); or
#'   \item the input \code{object}, when no DSEM is attached.
#' }
#'
#' @section Reading the variance rows:
#' Two-headed (\code{<->}) paths are variance parameters, and their reported
#' statistics follow `summary.dsem` rather than being rewritten here. Two
#' consequences are worth knowing when quoting these tables:
#' \itemize{
#'   \item Their sign is not identified — the model uses \eqn{|beta_z|} (the C++
#'     takes \code{sqrt(square(.))}), so a variance row's \code{Estimate} may
#'     appear negative while \code{recruitment_sd$R_sd} reports the same
#'     quantity as positive. Compare magnitudes, not signs.
#'   \item \code{z_value} / \code{p_value} on a variance are a test against a
#'     boundary of the parameter space, so the usual Wald interpretation does not
#'     apply; they are typically extreme regardless of evidence. Do not read a
#'     small p-value on a \code{<->} row as support for that variance term.
#' }
#'
#' @export
summary.Rceattle <- function(object, ...) {
  if (!inherits(object, "Rceattle")) {
    stop("Input is not an Rceattle model.")
  }

  print(object, ...)

  # No DSEM attached: behave exactly as the non-DSEM summary always has.
  if (is.null(object$dsem) || is.null(object$dsem$sem_full)) {
    return(invisible(object))
  }

  # sem_full is a data.frame from dsem::dsem(); as.data.frame() keeps it one if a
  # future dsem returns a matrix, whose columns would otherwise all coerce to
  # character and silently change every downstream column's type.
  model  <- as.data.frame(object$dsem$sem_full, stringsAsFactors = FALSE)
  beta_z <- object$estimated_params$beta_z
  par_idx <- as.numeric(model[, "parameter"])

  # sem_full$parameter is a 1-based index into beta_z, or 0 for a path fixed at
  # its start value. Refuse to guess if beta_z is missing or too short: without
  # it every estimate silently falls back to its start value, which looks like a
  # converged table of zeros rather than an error.
  if (is.null(beta_z)) {
    stop("The fit has no estimated_params$beta_z, so DSEM path coefficients ",
         "cannot be reported. Was this fit produced with a `dsem` argument?")
  }
  if (length(beta_z) < max(par_idx, 0, na.rm = TRUE)) {
    stop("estimated_params$beta_z has ", length(beta_z), " entries but the SEM ",
         "references parameter index ", max(par_idx, na.rm = TRUE), ".")
  }

  # Prepending NA lets index 0 select it; the start value then replaces it.
  # Keyed on `parameter == 0` rather than is.na(Estimate): is.na(NaN) is TRUE, so
  # testing the estimate would rewrite a path whose MLE came back NaN (diverged
  # fit, non-PD Hessian) into its start value and report it as if it were fixed.
  coefs <- data.frame(model, Estimate = c(NA, beta_z)[par_idx + 1])
  fixed_path <- !is.na(par_idx) & par_idx == 0
  coefs$Estimate[fixed_path] <- as.numeric(model[, "start"])[fixed_path]

  # Std_Error / z_value / p_value are always present, NA when no sdreport was
  # run. Making them conditional changes the column count between fits, and the
  # consumer scripts rbind() several models' tables together.
  coefs$Std_Error <- NA_real_
  if (!is.null(object$sdrep)) {
    SE <- as.list(object$sdrep, report = FALSE, what = "Std. Error")
    if (!is.null(SE$beta_z)) coefs$Std_Error <- c(NA, SE$beta_z)[par_idx + 1]
  }
  coefs$z_value <- coefs$Estimate / coefs$Std_Error
  coefs$p_value <- stats::pnorm(-abs(coefs$z_value)) * 2

  # Drop paths whose beta_z entry is mapped off. The map is indexed by beta_z
  # entry and coefs by sem row, so the two are bridged by par_idx rather than by
  # position: sem rows outnumber beta_z entries whenever a path is fixed at its
  # start value (parameter == 0) or several paths share one parameter, both of
  # which are ordinary. A fixed path owns no beta_z entry and is always kept --
  # it cannot be mapped off.
  bmap <- object$map$mapList$beta_z
  if (!is.null(bmap)) {
    if (length(bmap) >= max(par_idx, 0, na.rm = TRUE)) {
      estimated_path <- !is.na(par_idx) & par_idx >= 1
      keep <- rep(TRUE, nrow(coefs))
      keep[estimated_path] <- !is.na(bmap[par_idx[estimated_path]])
      coefs <- coefs[keep, , drop = FALSE]
    } else {
      warning("The DSEM map has ", length(bmap), " beta_z entries but the SEM ",
              "references parameter index ", max(par_idx, na.rm = TRUE),
              "; reporting all paths unfiltered.", call. = FALSE)
    }
  }

  if (nrow(coefs) > 0) {
    cat("\n<DSEM path coefficients>\n")
    print(coefs, row.names = FALSE)
  }

  rec_sd <- .rceattle_rec_sd(object)
  if (!is.null(rec_sd)) {
    cat("\n<Recruitment SD (R_sd)>\n")
    rec_sd_print <- rec_sd
    num <- vapply(rec_sd_print, is.numeric, logical(1))
    rec_sd_print[num] <- lapply(rec_sd_print[num], round, 4)
    print(rec_sd_print, row.names = FALSE)
  }

  invisible(list(coefficients = coefs, recruitment_sd = rec_sd))
}


# Per-species recruitment SD (R_sd) for summary.Rceattle(): the value the model
# used, whether it was estimated or fixed, and its Std_Error when an sdreport is
# available. Returns NULL if R_sd cannot be located.
#
# R_sd(sp) = |beta_z(rec_sd_idx(sp))| when estimated, else the fixed sem value;
# see build_dsem_objects() and the DSEM block in ceattle.cpp.
.rceattle_rec_sd <- function(object) {
  d <- object$data_list
  nspp <- d$nspp
  if (is.null(nspp) || nspp < 1) return(NULL)
  spnames <- if (!is.null(d$spnames)) d$spnames else as.character(seq_len(nspp))

  # Prefer the ADREPORTed R_sd (carries an SE), then quantities, then a report().
  rsd <- rep(NA_real_, nspp)
  se  <- rep(NA_real_, nspp)
  sdr <- object$sdrep
  if (!is.null(sdr) && !is.null(sdr$value)) {
    i <- which(names(sdr$value) == "R_sd")
    if (length(i) == nspp) {
      rsd <- as.numeric(sdr$value[i])
      se  <- as.numeric(sdr$sd[i])
    }
  }
  if (all(is.na(rsd)) && !is.null(object$quantities$R_sd)) {
    rsd <- as.numeric(object$quantities$R_sd)[seq_len(nspp)]
  }
  if (all(is.na(rsd)) && !is.null(object$obj)) {
    rsd <- tryCatch(as.numeric(object$obj$report()$R_sd)[seq_len(nspp)],
                    error = function(e) rep(NA_real_, nspp))
  }
  if (all(is.na(rsd))) return(NULL)

  # Estimated per species when that species' recruitment-SD beta_z entry is
  # mapped on; fall back to the model-level random_rec flag.
  estimated <- rep(isTRUE(as.logical(d$random_rec)), nspp)
  idx  <- object$dsem$tmb_inputs$data$rec_sd_idx
  bmap <- object$map$mapList$beta_z
  if (!is.null(idx) && !is.null(bmap) && length(idx) == nspp) {
    estimated <- vapply(seq_len(nspp), function(sp) {
      isTRUE(idx[sp] >= 1) && !is.na(bmap[idx[sp]])
    }, logical(1))
  }

  # Fixed SDs carry no estimation uncertainty. Mask first, then attach the column
  # unconditionally: testing any(!is.na(se)) beforehand made the column's very
  # presence depend on the fit, which destabilizes downstream rbind()s.
  se[!estimated] <- NA_real_
  data.frame(Species = spnames, R_sd = rsd, Estimated = estimated,
             Std_Error = se, stringsAsFactors = FALSE)
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
#' [TMB::sdreport()]. Random-effect covariance is not returned here --
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


#' Residuals from an Rceattle fit
#'
#' Returns a long-format data frame of residuals, following the convention of
#' [stats::residuals.glm()] where `type` selects the *kind* of residual. By
#' default residuals are returned for every applicable data source -- survey
#' indices, fishery catches, age/length composition (`comp`), and conditional
#' age-at-length (`caal`); use `source` to restrict to particular ones.
#'
#' Residual kinds (`type`):
#' \describe{
#'   \item{`"response"`}{Observed minus fitted. For `index` / `catch` this is on
#'     the log scale by default (matching the lognormal likelihood; set
#'     `scale = "natural"` for the arithmetic difference); for `comp` / `caal`
#'     it is the difference in proportions, observed minus fitted.}
#'   \item{`"pearson"`}{Standardized residuals. For `index` / `catch`,
#'     \eqn{(\log o - (\log\hat{o} - b\,\sigma^2/2))/\sigma} using the model's
#'     realized observation log-SD \eqn{\sigma} and the observation
#'     bias-adjustment flag \eqn{b} (`bias_adjust_obs`, default 1); for
#'     `comp` / `caal`,
#'     \eqn{(p - \hat{p})/\sqrt{\hat{p}(1 - \hat{p})/N}} with input sample size
#'     N.}
#'   \item{`"osa"`}{One-step-ahead residuals via [osa_residuals()], which builds
#'     the composition observation data on demand from any fit.}
#'   \item{`"process"`}{Process residuals via [process_residuals()] for the
#'     model's random-effect deviations; `source` does not apply.}
#' }
#'
#' Composition rows are returned in long form (one row per observation x
#' age/length bin) and carry the `Age0_Length1` flag from `comp_data` (`0` age,
#' `1` length); CAAL rows carry both the conditioning `Length` and the age `Bin`.
#'
#' Where a fleet uses tail accumulation (`Comp_accum_young` /
#' `Comp_accum_old`), composition residuals describe the bins the likelihood
#' actually fit: the tails are folded into their boundary bin and the bins
#' outside the window are not reported, because the model never fit them
#' separately. `Bin` names the age or length each residual belongs to and
#' `Accumulated` marks the boundary bins, which stand for a range rather than
#' the single bin they are named for. Fleets without accumulation are unchanged.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param type Residual kind: one of `"response"` (default), `"pearson"`,
#'   `"osa"`, or `"process"`.
#' @param source Data source(s) to include: any of `"index"`, `"catch"`,
#'   `"comp"`, `"caal"`, or `"all"` (default). `"diet"` (predator
#'   stomach-content composition) is also accepted but, because it uses a
#'   predator/prey schema rather than the fleet/bin layout, must be requested on
#'   its own. Ignored when `type = "process"`.
#' @param scale `"log"` (default) or `"natural"`. Only affects `"response"`
#'   residuals for `index` / `catch`.
#' @param species Optional species code(s) to include (matched against the
#'   `Species` column). Default `NULL` keeps all species. Mirrors the `species`
#'   argument of [plot.rceattle_osa()].
#' @param ... Passed to [osa_residuals()] (e.g. `method`, `seed`) when
#'   `type = "osa"`, or to [process_residuals()] when `type = "process"`.
#'
#' @return A `data.frame` with columns `Source`, `Fleet_code`, `Fleet_name`,
#'   `Species`, `Sex`, `Year`, `Length`, `Bin`, `Age0_Length1`, `Sample_size`,
#'   `Accumulated`, `Observed`, `Fitted`, `Residual`. Columns are `NA` where they
#'   do not apply; `Accumulated` is `FALSE` except on a composition bin that
#'   absorbed a folded tail.
#' @export
residuals.Rceattle <- function(object, type = "response", source = "all",
                               scale = "log", species = NULL, ...) {
  # Following stats::residuals.glm(), `type` selects the residual *kind* and
  # `source` the data source(s).
  type   <- match.arg(type, c("response", "pearson", "osa", "process"))
  source <- match.arg(source, c("all", "index", "catch", "comp", "caal", "diet"),
                      several.ok = TRUE)
  # TODO(review): osa path -- osa_residuals("all") includes diet, but this "all"
  # expansion (forwarded when type = "osa") omits it, and source = c("all",
  # "diet") silently drops diet too. Add "diet" here or document the asymmetry.
  if ("all" %in% source) source <- c("index", "catch", "comp", "caal")
  scale  <- match.arg(scale, c("log", "natural"))

  # Optional species filter, applied to whichever data frame is returned. Keep
  # rows with NA species (e.g. catchability process residuals, which belong to a
  # fleet, not a species) so a species filter does not silently drop them.
  .sp_filter <- function(df) {
    if (is.null(species)) df
    else df[df$Species %in% species | is.na(df$Species), , drop = FALSE]
  }

  # Diet (predator stomach-content) composition uses a predator/prey schema that
  # does not fit the fleet/bin layout of the other sources, so it is returned on
  # its own. (For type = "osa" it flows through to osa_residuals() below.)
  if ("diet" %in% source && type %in% c("response", "pearson")) {
    if (!all(source == "diet")) {
      stop("source = 'diet' uses a predator/prey schema and must be requested ",
           "on its own (not combined with index/catch/comp/caal).")
    }
    dd <- object$data_list$diet_data
    if (is.null(dd) || nrow(dd) == 0 || is.null(object$quantities$diet_hat)) {
      stop("This model has no fitted diet composition (need msmMode > 0 with ",
           "diet_data).")
    }
    obs <- dd$Stomach_proportion_by_weight
    hat <- object$quantities$diet_hat[, 2]
    res <- if (type == "pearson") .pearson_proportion(obs, hat, dd$Sample_size) else obs - hat
    return(.sp_filter(data.frame(
      Source = "diet", Species = dd$Pred, Pred_sex = dd$Pred_sex,
      Prey = dd$Prey, Prey_sex = dd$Prey_sex, Pred_age = dd$Pred_age,
      Prey_age = dd$Prey_age, Year = dd$Year, Sample_size = dd$Sample_size,
      Observed = obs, Fitted = hat, Residual = res, stringsAsFactors = FALSE)))
  }

  # Process residuals (recruitment / random-effect deviations) are computed
  # separately and reshaped into the common schema; `source` does not apply.
  if (type == "process") {
    pr <- process_residuals(object, ...)
    return(.sp_filter(data.frame(
      Source       = pr$source,
      Fleet_code   = NA_integer_,
      Fleet_name   = NA_character_,
      Species      = pr$species,
      Sex          = pr$sex,
      Year         = pr$year,
      Length       = rep(NA_real_, nrow(pr)),
      Bin          = pr$age_length_bin,
      Age0_Length1 = rep(NA_integer_, nrow(pr)),
      Sample_size  = rep(NA_real_, nrow(pr)),
      Observed     = pr$observed,
      Fitted       = pr$predicted,
      Residual     = pr$residual,
      stringsAsFactors = FALSE)))
  }

  # One-step-ahead residuals are computed separately (oneStepPredict) and
  # reshaped into the common schema. They are standard-normal, so `scale` does
  # not apply. `source` selects the observation sources; extra args (method,
  # seed) flow through `...` to osa_residuals().
  if (type == "osa") {
    osa <- osa_residuals(object, source = source, ...)
    fc  <- object$data_list$fleet_control
    return(.sp_filter(data.frame(
      Source       = osa$source,
      Fleet_code   = osa$fleet,
      Fleet_name   = fc$Fleet_name[match(osa$fleet, fc$Fleet_code)],
      Species      = osa$species,
      Sex          = osa$sex,
      Year         = osa$year,
      Length       = osa$length,
      Bin          = osa$age_length_bin,
      Age0_Length1 = ifelse(osa$index_label == "length", 1L,
                            ifelse(osa$index_label == "age", 0L, NA_integer_)),
      Sample_size  = rep(NA_real_, nrow(osa)),
      Observed     = osa$observed,
      Fitted       = osa$predicted,
      Residual     = osa$residual,
      stringsAsFactors = FALSE)))
  }

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
      # TRUE where tail accumulation folded neighbouring bins into this one, so
      # `Bin` names a range rather than the single age or length it appears to.
      Accumulated  = rep(FALSE, n),
      Observed     = numeric(n),
      Fitted       = numeric(n),
      Residual     = numeric(n),
      stringsAsFactors = FALSE
    )
  }

  # Observation bias-adjustment flag (default 1): the index/catch lognormal
  # likelihood fits to mean log(hat) - bias_adjust_obs * sigma^2/2, so the Pearson
  # residual must subtract the same offset.
  ba_obs <- d$bias_adjust_obs
  if (is.null(ba_obs) && !is.null(object$obj)) ba_obs <- object$obj$env$data$bias_adjust_obs
  if (is.null(ba_obs)) ba_obs <- 1
  ba_obs <- as.numeric(ba_obs)[1]

  out <- list()

  if ("index" %in% source &&
      !is.null(d$index_data) && !is.null(q$index_hat)) {
    idx <- d$index_data
    obs <- idx$Observation
    hat <- as.numeric(q$index_hat)
    if (type == "pearson") {
      sigma <- if (!is.null(q$log_index_sd)) as.numeric(q$log_index_sd) else idx$Log_sd
      res   <- (log(obs) - (log(hat) - ba_obs * sigma^2 / 2)) / sigma
    } else {
      res <- if (scale == "log") log(obs) - log(hat) else obs - hat
    }
    df <- empty_row(nrow(idx))
    df$Source     <- "index"
    df$Fleet_code <- idx$Fleet_code
    df$Fleet_name <- idx$Fleet_name
    df$Species    <- idx$Species
    df$Year       <- idx$Year
    df$Observed   <- obs
    df$Fitted     <- hat
    df$Residual   <- res
    # Keep only genuine observations: drop projection / NA rows and non-positive
    # values, which carry no (lognormal) residual.
    # TODO(review): held-out rows (Year <= 0) with a positive observation still
    # pass this filter and appear with a negative Year, though the C++ likelihood
    # excludes them (Year>0 & Year<=endyr & flt_type>0). Gate on that same rule;
    # the "drop projection rows" comment above is not fully accurate for them.
    df <- df[is.finite(df$Observed) & is.finite(df$Fitted) &
             df$Observed > 0 & df$Fitted > 0, , drop = FALSE]
    out$index <- df
  }

  if ("catch" %in% source &&
      !is.null(d$catch_data) && !is.null(q$catch_hat)) {
    ctc <- d$catch_data
    obs <- ctc$Catch
    hat <- as.numeric(q$catch_hat)
    if (type == "pearson") {
      sigma <- if (!is.null(q$log_catch_sd)) as.numeric(q$log_catch_sd) else ctc$Log_sd
      res   <- (log(obs) - (log(hat) - ba_obs * sigma^2 / 2)) / sigma
    } else {
      res <- if (scale == "log") log(obs) - log(hat) else obs - hat
    }
    df <- empty_row(nrow(ctc))
    df$Source     <- "catch"
    df$Fleet_code <- ctc$Fleet_code
    df$Fleet_name <- ctc$Fleet_name
    df$Species    <- ctc$Species
    df$Year       <- ctc$Year
    df$Observed   <- obs
    df$Fitted     <- hat
    df$Residual   <- res
    # Keep only genuine observations: drop projection / NA rows and non-positive
    # values, which carry no (lognormal) residual.
    # TODO(review): held-out rows (Year <= 0) with a positive observation still
    # pass this filter and appear with a negative Year, though the C++ likelihood
    # excludes them (Year>0 & Year<=endyr & flt_type>0). Gate on that same rule;
    # the "drop projection rows" comment above is not fully accurate for them.
    df <- df[is.finite(df$Observed) & is.finite(df$Fitted) &
             df$Observed > 0 & df$Fitted > 0, , drop = FALSE]
    out$catch <- df
  }

  if ("comp" %in% source &&
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

    # Fold the tails onto the bins the likelihood fit (see .fold_comp_props()).
    # NULL means no fleet accumulates, so the rectangular path below is used and
    # every existing model is unchanged.
    folded <- .fold_comp_props(d, cd, obs_prop, hat_prop, a0l1, n_bin)

    if (is.null(folded)) {
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
    } else {
      df <- empty_row(length(folded$obs))
      df$Source       <- "comp"
      df$Fleet_code   <- cd$Fleet_code[folded$row]
      df$Fleet_name   <- cd$Fleet_name[folded$row]
      df$Species      <- cd$Species[folded$row]
      df$Sex          <- if (!is.null(cd$Sex)) cd$Sex[folded$row] else NA_integer_
      df$Year         <- cd$Year[folded$row]
      df$Bin          <- folded$bin
      df$Age0_Length1 <- a0l1[folded$row]
      df$Sample_size  <- cd$Sample_size[folded$row]
      df$Accumulated  <- folded$acc
      df$Observed     <- folded$obs
      df$Fitted       <- folded$hat
    }
    df$Residual     <- if (type == "pearson")
      .pearson_proportion(df$Observed, df$Fitted, df$Sample_size)
    else df$Observed - df$Fitted
    # Drop zero-padded phantom bins (ragged multispecies comps have
    # Observed == Fitted == 0), which otherwise give 0/0 = NaN Pearson residuals.
    df <- df[!is.na(df$Observed) & !is.na(df$Fitted) &
             !(df$Observed == 0 & df$Fitted == 0), , drop = FALSE]
    out$comp <- df
  }

  if ("caal" %in% source &&
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
    df$Residual    <- if (type == "pearson")
      .pearson_proportion(df$Observed, df$Fitted, df$Sample_size)
    else df$Observed - df$Fitted
    # Drop zero-padded phantom bins (ragged multispecies comps have
    # Observed == Fitted == 0), which otherwise give 0/0 = NaN Pearson residuals.
    df <- df[!is.na(df$Observed) & !is.na(df$Fitted) &
             !(df$Observed == 0 & df$Fitted == 0), , drop = FALSE]
    out$caal <- df
  }

  if (length(out) == 0) return(empty_row(0))
  .sp_filter(do.call(rbind, out))
}


#' Pearson residual for an observed proportion
#'
#' \eqn{(p - \hat p)/\sqrt{\hat p (1 - \hat p)/N}} -- the multinomial-bin Pearson
#' residual shared by the composition, conditional-age-at-length, and diet
#' branches of [residuals.Rceattle()].
#' @param observed,fitted Observed and fitted proportions.
#' @param n Input sample size.
#' @keywords internal
.pearson_proportion <- function(observed, fitted, n) {
  (observed - fitted) / sqrt(fitted * (1 - fitted) / n)
}


#' Fold composition proportions onto the bins the likelihood actually fit
#'
#' Tail accumulation (`Comp_accum_young` / `Comp_accum_old`) folds a
#' composition's tails into a boundary bin and fits only `[yng, old]`, per fleet
#' and per sex block. A residual taken on the unfolded row therefore describes
#' bins the model never fit separately.
#'
#' Reuses [build_osa_data()]'s fold, so the Pearson and OSA residuals of a fleet
#' cover the same bins with the same labels.
#'
#' @param d The model's `data_list`.
#' @param cd Its `comp_data`.
#' @param obs_prop,hat_prop Observed / predicted proportions, `[row, bin]`.
#' @param a0l1 Per-row `Age0_Length1` (0 = age bins, 1 = length bins).
#' @param n_bin Number of bin columns, including any joint-sex doubling.
#' @return `NULL` when no fleet accumulates (the caller keeps its vectorized
#'   path); otherwise a list of equal-length vectors `row` (index into `cd`),
#'   `bin` (the ordinal each value stands for), `obs`, `hat`, and `acc` (whether
#'   the bin absorbed a folded tail).
#' @keywords internal
#' @noRd
.fold_comp_props <- function(d, cd, obs_prop, hat_prop, a0l1, n_bin) {
  acc_y <- d$comp_accum_young
  acc_o <- d$comp_accum_old
  if (is.null(acc_y)) acc_y <- d$fleet_control$Comp_accum_young
  if (is.null(acc_o)) acc_o <- d$fleet_control$Comp_accum_old
  if (is.null(acc_y) && is.null(acc_o)) return(NULL)

  nages    <- d$nages
  nlengths <- d$nlengths
  fleets   <- cd$Fleet_code
  sexes    <- if (!is.null(cd$Sex)) cd$Sex else rep(NA_integer_, nrow(cd))

  # Per-row window, clamped exactly as build_osa_data() and the cpp do.
  win <- lapply(seq_len(nrow(cd)), function(r) {
    flt <- fleets[r]
    sp  <- cd$Species[r]
    nb  <- if (isTRUE(a0l1[r] == 1)) nlengths[sp] else nages[sp]
    if (is.na(nb) || nb < 1L) nb <- n_bin
    yng <- if (!is.null(acc_y)) acc_y[flt] else 1L
    old <- if (!is.null(acc_o)) acc_o[flt] else 0L
    if (is.na(yng) || yng < 1L) yng <- 1L
    if (is.na(old) || old < 1L || old > nb) old <- nb
    if (yng > old) yng <- old
    list(nb = nb, yng = yng, old = old,
         blk = if (isTRUE(sexes[r] == 3)) 2L else 1L)
  })
  if (!any(vapply(win, function(w) w$yng > 1L || w$old < w$nb, logical(1)))) {
    return(NULL)   # nothing folds: leave the caller's fast path alone
  }

  parts <- lapply(seq_len(nrow(cd)), function(r) {
    w <- win[[r]]
    n_use <- w$blk * w$nb
    # A row declaring more bins than it has columns is inconsistent, and
    # data_check() is where that is judged. Report its raw bins rather than
    # index past the end.
    if (n_use > n_bin) {
      o <- as.numeric(obs_prop[r, seq_len(n_bin)])
      h <- as.numeric(hat_prop[r, seq_len(n_bin)])
      return(list(row = rep(r, n_bin), bin = seq_len(n_bin),
                  obs = o, hat = h, acc = rep(FALSE, n_bin)))
    }
    o <- as.numeric(obs_prop[r, seq_len(n_use)])
    h <- as.numeric(hat_prop[r, seq_len(n_use)])
    does_fold <- w$yng > 1L || w$old < w$nb
    if (does_fold) {
      o <- .fold_comp_bins(o, w$nb, w$blk, w$yng, w$old)
      h <- .fold_comp_bins(h, w$nb, w$blk, w$yng, w$old)
    }
    keep <- seq.int(w$yng, w$old)
    bin  <- rep((seq_len(w$blk) - 1L) * w$nb, each = length(keep)) +
      rep(keep, times = w$blk)
    # Only the boundary bins absorbed a tail, so only they cover more than the
    # age they are named for.
    acc1 <- rep(FALSE, length(keep))
    if (does_fold) {
      if (w$yng > 1L)     acc1[1] <- TRUE
      if (w$old < w$nb)   acc1[length(acc1)] <- TRUE
    }
    list(row = rep(r, length(o)), bin = bin, obs = o, hat = h,
         acc = rep(acc1, times = w$blk))
  })

  list(row = unlist(lapply(parts, `[[`, "row"), use.names = FALSE),
       bin = unlist(lapply(parts, `[[`, "bin"), use.names = FALSE),
       obs = unlist(lapply(parts, `[[`, "obs"), use.names = FALSE),
       hat = unlist(lapply(parts, `[[`, "hat"), use.names = FALSE),
       acc = unlist(lapply(parts, `[[`, "acc"), use.names = FALSE))
}


# Catalogue of derived quantities that as.data.frame.Rceattle() can extract,
# keyed by REPORT name, with their shape: "sy" = matrix(nspp, nyrs);
# "ssay" = array(nspp, max_sex, max_age, nyrs). Quantities that are also
# ADREPORT'd in the TMB template are flagged so we can pull standard errors
# out of sdrep$value / sdrep$sd by name.
.RCEATTLE_QUANTITIES <- list(
  biomass             = list(shape = "sy",   adreport = TRUE),
  ssb                 = list(shape = "sy",   adreport = TRUE),
  R                   = list(shape = "sy",   adreport = TRUE),
  biomass_depletion   = list(shape = "sy",   adreport = FALSE),
  ssb_depletion       = list(shape = "sy",   adreport = FALSE),
  B0                  = list(shape = "sy",   adreport = FALSE),
  SB0                 = list(shape = "sy",   adreport = FALSE),
  DynamicB0           = list(shape = "sy",   adreport = FALSE),
  DynamicSB0          = list(shape = "sy",   adreport = FALSE),
  DynamicSBF          = list(shape = "sy",   adreport = FALSE),
  exploitable_biomass = list(shape = "sy",   adreport = FALSE),
  F_spp               = list(shape = "sy",   adreport = FALSE),
  N_at_age            = list(shape = "ssay", adreport = FALSE),
  biomass_at_age      = list(shape = "ssay", adreport = FALSE),
  Z_at_age            = list(shape = "ssay", adreport = FALSE),
  M_at_age            = list(shape = "ssay", adreport = FALSE),
  M1_at_age           = list(shape = "ssay", adreport = FALSE),
  M2_at_age           = list(shape = "ssay", adreport = FALSE),
  F_at_age            = list(shape = "ssay", adreport = FALSE),
  consumption_at_age  = list(shape = "ssay", adreport = FALSE),
  B_eaten_as_prey     = list(shape = "ssay", adreport = FALSE),
  NByage0             = list(shape = "ssay", adreport = FALSE),
  NByageF             = list(shape = "ssay", adreport = FALSE)
)


#' Tidy long-format derived quantities from an Rceattle fit
#'
#' Returns the model's derived population quantities in long form so that
#' custom plots and post-processing don't have to walk the deeply nested
#' `quantities` list or inherit the dimnames decisions in
#' [rename_output()]. Two shapes are flattened into one tidy frame:
#' species-by-year quantities (e.g. `biomass`, `ssb`, `R`, `F_spp`) and
#' species-by-sex-by-age-by-year quantities (e.g. `N_at_age`,
#' `biomass_at_age`, `M_at_age`). For the species-level shapes, `sex`
#' and `age` are returned as `NA`. Cells of the 4D arrays that are
#' padded out to `max(nsex)` / `max(nages)` for species with fewer
#' sexes or ages are dropped.
#'
#' Standard errors (`se`) and confidence intervals (`lwr`, `upr`) are
#' populated from the TMB `sdreport` for any quantity that was
#' `ADREPORT`'d (currently `biomass`, `ssb`, `R`); other quantities and
#' fits produced with `getsd = FALSE` get `NA` for `se` / `lwr` / `upr`.
#' Set `ci_level` to widen or narrow the band.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param row.names,optional Ignored; present for the [as.data.frame()] generic.
#' @param which Character vector of quantity names to extract, or
#'   `"all"` for every quantity with a known shape. Defaults to a common
#'   population-level summary. See `names(Rceattle:::.RCEATTLE_QUANTITIES)`
#'   for the full list.
#' @param ci_level Confidence level for `lwr` and `upr`. Default `0.95`.
#' @param ... Currently unused.
#'
#' @return A `data.frame` with columns `year`, `species`, `sex`, `age`,
#'   `quantity`, `value`, `se`, `lwr`, `upr`. `species` is the character
#'   species name from `data_list$spnames`. Rows are sorted in the order
#'   `which` was given.
#' @export
as.data.frame.Rceattle <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   which = c("biomass", "ssb", "R",
                                             "biomass_depletion",
                                             "ssb_depletion", "F_spp"),
                                   ci_level = 0.95,
                                   ...) {
  known <- names(.RCEATTLE_QUANTITIES)
  if (length(which) == 1 && identical(which, "all")) which <- known
  bad <- setdiff(which, known)
  if (length(bad)) {
    stop("Unknown quantity name(s): ", paste(bad, collapse = ", "),
         ".\nKnown: ", paste(known, collapse = ", "))
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).")
  }

  d <- x$data_list
  q <- x$quantities
  yrs   <- d$styr:d$projyr
  nyrs  <- length(yrs)
  nspp  <- d$nspp
  nages <- d$nages
  nsex  <- d$nsex
  minage  <- if (!is.null(d$minage)) d$minage else rep(1L, nspp)
  spnames <- if (!is.null(d$spnames)) d$spnames else as.character(seq_len(nspp))

  z <- stats::qnorm(0.5 + ci_level / 2)

  # sdreport vectors are stored in the same column-major order TMB used to
  # stamp the matrix/array, so reshaping `sd` against the report's dim
  # lines up cell-for-cell. A fit run with getsd=FALSE has no $sdrep.
  sdrep <- x$sdrep
  sd_lookup <- function(name, dims) {
    if (is.null(sdrep) || is.null(sdrep$value)) return(NULL)
    idx <- which(names(sdrep$value) == name)
    if (!length(idx) || length(idx) != prod(dims)) return(NULL)
    array(sdrep$sd[idx], dim = dims)
  }

  empty <- data.frame(
    year     = integer(),
    species  = character(),
    sex      = integer(),
    age      = integer(),
    quantity = character(),
    value    = numeric(),
    se       = numeric(),
    lwr      = numeric(),
    upr      = numeric(),
    stringsAsFactors = FALSE
  )

  out <- vector("list", length(which))
  names(out) <- which

  for (qn in which) {
    spec <- .RCEATTLE_QUANTITIES[[qn]]
    arr  <- q[[qn]]
    if (is.null(arr)) next

    if (spec$shape == "sy") {
      mat <- as.matrix(arr)
      if (length(dim(mat)) != 2 ||
          dim(mat)[1] != nspp || dim(mat)[2] != nyrs) next
      sd_mat <- if (spec$adreport) sd_lookup(qn, dim(mat)) else NULL
      grid <- expand.grid(species_idx = seq_len(nspp),
                          year_idx    = seq_len(nyrs),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
      cell <- cbind(grid$species_idx, grid$year_idx)
      val  <- mat[cell]
      sdv  <- if (!is.null(sd_mat)) sd_mat[cell] else NA_real_
      out[[qn]] <- data.frame(
        year     = yrs[grid$year_idx],
        species  = spnames[grid$species_idx],
        sex      = NA_integer_,
        age      = NA_integer_,
        quantity = qn,
        value    = as.numeric(val),
        se       = as.numeric(sdv),
        lwr      = as.numeric(val - z * sdv),
        upr      = as.numeric(val + z * sdv),
        stringsAsFactors = FALSE
      )
    } else if (spec$shape == "ssay") {
      a4 <- arr
      if (length(dim(a4)) != 4 || dim(a4)[1] != nspp || dim(a4)[4] != nyrs) next
      max_sex <- dim(a4)[2]
      max_age <- dim(a4)[3]
      sd4 <- if (spec$adreport) sd_lookup(qn, dim(a4)) else NULL
      grid <- expand.grid(species_idx = seq_len(nspp),
                          sex_idx     = seq_len(max_sex),
                          age_idx     = seq_len(max_age),
                          year_idx    = seq_len(nyrs),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
      keep <- grid$sex_idx <= nsex[grid$species_idx] &
              grid$age_idx <= nages[grid$species_idx]
      grid <- grid[keep, , drop = FALSE]
      cell <- cbind(grid$species_idx, grid$sex_idx,
                    grid$age_idx, grid$year_idx)
      val  <- a4[cell]
      sdv  <- if (!is.null(sd4)) sd4[cell] else NA_real_
      out[[qn]] <- data.frame(
        year     = yrs[grid$year_idx],
        species  = spnames[grid$species_idx],
        sex      = as.integer(grid$sex_idx),
        age      = as.integer(grid$age_idx + minage[grid$species_idx] - 1L),
        quantity = qn,
        value    = as.numeric(val),
        se       = as.numeric(sdv),
        lwr      = as.numeric(val - z * sdv),
        upr      = as.numeric(val + z * sdv),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(empty)
  do.call(rbind, c(unname(out), list(make.row.names = FALSE)))
}
