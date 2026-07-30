# =============================================================================
# Fit diagnostics: survey-index and fishery-catch fits, and index residuals.
# All plotters build a tidy data frame, render it with the shared ggplot theme
# (theme_classic) + colourblind palette, and return the ggplot object.
# =============================================================================

# Assemble a tidy observed-vs-predicted data frame for a fleet "fit" plot.
# `kind` selects the index (survey) or catch (fishery) data/quantities. Returns
# one row per (model, fleet, year) with the observation, its lognormal 95%
# band, and the model prediction.
.fleet_fit_df <- function(Rceattle, kind = c("index", "catch"),
                          species = NULL, incl_proj = FALSE, maxyr = NULL,
                          model_names = NULL) {
  kind <- match.arg(kind)
  cfg <- switch(kind,
    index = list(data = "index_data", obs = "Observation",
                 hat = "index_hat", sd = "log_index_sd", type = "Survey"),
    catch = list(data = "catch_data", obs = "Catch",
                 hat = "catch_hat", sd = "log_catch_sd", type = "Fishery"))

  model_names_use <- .model_labels(Rceattle, model_names)
  fc <- Rceattle[[1]]$data_list$fleet_control
  codes <- fc$Fleet_code[fc$Fleet_type == cfg$type]
  if (is.null(species)) species <- seq_len(Rceattle[[1]]$data_list$nspp)

  out <- list()
  for (i in seq_along(Rceattle)) {
    dat <- Rceattle[[i]]$data_list[[cfg$data]]
    dat$.obs <- dat[[cfg$obs]]
    dat$.sd  <- Rceattle[[i]]$quantities[[cfg$sd]]
    dat$.hat <- Rceattle[[i]]$quantities[[cfg$hat]]
    # Year < 0 flags a prediction-only row (excluded from the likelihood; see
    # ceattle.cpp "Likelihood (yr > 0) vs prediction (yr < 0)"). Drop
    # them so they are not drawn as fitted observations.
    keep <- dat$Species %in% species & dat$Fleet_code %in% codes & dat$Year > 0
    if (!incl_proj) keep <- keep & dat$Year <= Rceattle[[i]]$data_list$endyr
    if (!is.null(maxyr)) keep <- keep & dat$Year <= maxyr
    dat <- dat[keep, , drop = FALSE]
    if (nrow(dat) == 0) next
    pos <- dat$.obs > 0
    upr <- lwr <- rep(NA_real_, nrow(dat))
    upr[pos] <- stats::qlnorm(0.975, log(dat$.obs[pos]), dat$.sd[pos])
    lwr[pos] <- stats::qlnorm(0.025, log(dat$.obs[pos]), dat$.sd[pos])
    out[[length(out) + 1L]] <- data.frame(
      Model       = model_names_use[i],
      Fleet       = as.character(fc$Fleet_name[match(dat$Fleet_code,
                                                     fc$Fleet_code)]),
      Year        = abs(dat$Year),
      Observation = dat$.obs,
      Lower95     = lwr,
      Upper95     = upr,
      Predicted   = dat$.hat,
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0L) {
    stop("No data to plot: the selected species/fleet(s) have no ",
         "likelihood observations (all rows are projections or excluded).",
         call. = FALSE)
  }
  df <- do.call(rbind, out)
  df$Model <- factor(df$Model, levels = unique(model_names_use))
  df$Fleet <- factor(df$Fleet,
                     levels = unique(fc$Fleet_name[match(codes, fc$Fleet_code)]))
  df
}

# Shared renderer for index/catch fit plots.
.render_fleet_fit <- function(df, ylab, error = TRUE, single_model = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Year, colour = .data$Model,
                                        fill = .data$Model))
  if (error) {
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Lower95,
                                          ymax = .data$Upper95),
                             width = 0, colour = "grey60") +
      ggplot2::geom_point(ggplot2::aes(y = .data$Observation),
                          shape = 21, fill = "white", colour = "grey30")
  }
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$Predicted), linewidth = 1)
  p <- .rceattle_scale(
    p +
      ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
      ggplot2::labs(x = "Year", y = ylab) +
      .rceattle_theme())
  if (single_model) p <- p + ggplot2::guides(colour = "none", fill = "none")
  p
}


#' Survey index fits
#'
#' Plots fitted survey CPUE indices: observed points with lognormal 95%
#' intervals and the model-predicted index, faceted by survey fleet.
#'
#' @param Rceattle A single [fit_mod()] object or a list of them (overlaid).
#' @param file Optional file stem; if given the figure is written to
#'   `<file>_survey_indices.png`.
#' @param model_names Legend labels for the models.
#' @param species Species (indices) to include. Default all.
#' @param incl_proj Include projection years.
#' @param width,height Saved figure size (inches).
#' @param error Draw observed points and error bars.
#' @param log Plot the index on the log scale.
#' @param line_col,right_adj,top_adj,single.plots Deprecated base-graphics
#'   arguments, retained for back-compatibility and ignored.
#' @return A `ggplot` object.
#' @export
plot_index <- function(Rceattle,
                       file = NULL,
                       model_names = NULL,
                       species = NULL,
                       incl_proj = FALSE,
                       width = 7,
                       height = 6.5,
                       error = TRUE,
                       log = FALSE,
                       line_col = NULL,
                       right_adj = 0,
                       top_adj = 0.05,
                       single.plots = FALSE) {
  Rceattle <- .as_model_list(Rceattle)
  df <- .fleet_fit_df(Rceattle, kind = "index", species = species,
                      incl_proj = incl_proj, model_names = model_names)
  ylab <- "Index"
  if (log) {
    for (col in c("Observation", "Lower95", "Upper95", "Predicted")) {
      df[[col]] <- log(df[[col]])
    }
    ylab <- "log(Index)"
  }
  p <- .render_fleet_fit(df, ylab = ylab, error = error,
                         single_model = nlevels(df$Model) < 2L)
  .save_ggplot(p, file = file, suffix = "survey_indices",
               width = width, height = height)
}


#' Survey index fits on the log scale (deprecated)
#'
#' @description Deprecated. Use `plot_index(..., log = TRUE)`, which draws the
#'   survey index on the log scale.
#' @param ... Passed through to [plot_index()].
#' @return A `ggplot` object.
#' @keywords internal
#' @export
plot_logindex <- function(...) {
  .Deprecated("plot_index")
  plot_index(..., log = TRUE)
}


#' Fishery catch fits
#'
#' Plots fitted fishery catch: observed points with lognormal 95% intervals and
#' the model-predicted catch, faceted by fishery fleet. For MSE objects the
#' projection period is summarised with 50% / 95% ribbons across simulations.
#'
#' @inheritParams plot_index
#' @param maxyr Last year to plot.
#' @param mse Is `Rceattle` an MSE object (operating models)?
#' @param fleets Fishery fleets to include (indices into the sorted fleet
#'   codes). Default all.
#' @param alpha,lwd,ymax Deprecated base-graphics arguments, ignored.
#' @return A `ggplot` object.
#' @export
plot_catch <- function(Rceattle,
                       file = NULL,
                       model_names = NULL,
                       species = NULL,
                       incl_proj = FALSE,
                       width = 7,
                       height = 6.5,
                       error = TRUE,
                       maxyr = NULL,
                       mse = FALSE,
                       fleets = NULL,
                       line_col = NULL,
                       right_adj = 0,
                       top_adj = 1.2,
                       single.plots = FALSE,
                       alpha = 0.4,
                       lwd = 2,
                       ymax = NULL) {
  if (mse) {
    Rceattle <- .as_model_list(Rceattle, mse = TRUE, OM = TRUE)
    incl_proj <- TRUE
  } else {
    Rceattle <- .as_model_list(Rceattle)
  }

  df <- .fleet_fit_df(Rceattle, kind = "catch", incl_proj = incl_proj,
                      maxyr = maxyr, model_names = model_names)

  if (!is.null(fleets)) {
    lev <- levels(df$Fleet)
    df <- df[df$Fleet %in% lev[fleets], , drop = FALSE]
    df$Fleet <- droplevels(df$Fleet)
  }

  if (mse) {
    # Summarise predicted catch across operating models per fleet/year
    agg <- stats::aggregate(Predicted ~ Fleet + Year, data = df,
      FUN = function(v) c(mean = mean(v),
                          l95 = stats::quantile(v, 0.025, names = FALSE),
                          u95 = stats::quantile(v, 0.975, names = FALSE),
                          l50 = stats::quantile(v, 0.25, names = FALSE),
                          u50 = stats::quantile(v, 0.75, names = FALSE)))
    mse_df <- data.frame(Fleet = agg$Fleet, Year = agg$Year,
                         agg$Predicted)
    p <- ggplot2::ggplot(mse_df, ggplot2::aes(x = .data$Year)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$l95, ymax = .data$u95),
                           alpha = 0.25, fill = "grey40") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$l50, ymax = .data$u50),
                           alpha = 0.4, fill = "grey40") +
      ggplot2::geom_line(ggplot2::aes(y = .data$mean), linewidth = 1) +
      ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Catch") +
      .rceattle_theme()
    return(.save_ggplot(p, file = file, suffix = "fishery_catch",
                        width = width, height = height))
  }

  p <- .render_fleet_fit(df, ylab = "Catch", error = error,
                         single_model = nlevels(df$Model) < 2L)
  .save_ggplot(p, file = file, suffix = "fishery_catch",
               width = width, height = height)
}


#' Survey index residuals
#'
#' Plots log residuals `log(predicted) - log(observed)` of the survey index by
#' year, faceted by survey fleet (`residual_type = "pearson"`, the default), or
#' one-step-ahead (OSA) residual diagnostics for a single fit
#' (`residual_type = "osa"`, via [osa_residuals()] / [plot.rceattle_osa()]).
#'
#' @inheritParams plot_index
#' @param residual_type `"pearson"` (log index residuals) or `"osa"`.
#' @return A `ggplot` object.
#' @export
plot_indexresidual <- function(Rceattle,
                               file = NULL,
                               model_names = NULL,
                               species = NULL,
                               width = 7,
                               height = 6.5,
                               residual_type = c("pearson", "osa"),
                               line_col = NULL,
                               right_adj = 0,
                               top_adj = 0.05,
                               incl_proj = FALSE,
                               single.plots = FALSE) {
  residual_type <- match.arg(residual_type)
  if (residual_type == "osa") {
    fit <- if (inherits(Rceattle, "Rceattle")) Rceattle else Rceattle[[1]]
    osa <- osa_residuals(fit, source = c("index", "catch"))
    return(plot(osa))
  }

  Rceattle <- .as_model_list(Rceattle)
  model_names_use <- .model_labels(Rceattle, model_names)
  fc <- Rceattle[[1]]$data_list$fleet_control
  if (is.null(species)) species <- seq_len(Rceattle[[1]]$data_list$nspp)

  out <- list()
  for (i in seq_along(Rceattle)) {
    dat <- Rceattle[[i]]$data_list$index_data
    dat$.hat <- Rceattle[[i]]$quantities$index_hat
    # Drop prediction-only rows (Year < 0); they are not in the likelihood and
    # have no residual to plot.
    keep <- dat$Species %in% species & dat$Year > 0
    dat <- dat[keep, , drop = FALSE]
    if (nrow(dat) == 0) next
    out[[length(out) + 1L]] <- data.frame(
      Model    = model_names_use[i],
      Fleet    = as.character(fc$Fleet_name[match(dat$Fleet_code,
                                                  fc$Fleet_code)]),
      Year     = abs(dat$Year),
      Residual = log(dat$.hat) - log(dat$Observation),
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0L) {
    stop("No index residuals to plot: the selected species have no ",
         "likelihood observations (all rows are projections or excluded).",
         call. = FALSE)
  }
  df <- do.call(rbind, out)
  df$Model <- factor(df$Model, levels = unique(model_names_use))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Year, y = .data$Residual,
                                        colour = .data$Model)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$Year, yend = 0),
                          position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6))
  p <- .rceattle_scale(
    p +
      ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "log(index) residual") +
      .rceattle_theme(),
    aesthetics = "colour")
  if (nlevels(df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")
  .save_ggplot(p, file = file, suffix = "index_residuals",
               width = width, height = height)
}
