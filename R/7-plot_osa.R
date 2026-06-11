#' Plot one-step-ahead (OSA) residual diagnostics
#'
#' @description
#' Diagnostic plots for an `rceattle_osa` object (from [osa_residuals()]),
#' following the recommendations of Stewart and Monnahan (2025) and the styling
#' of the NOAA-AFSC `afscOSA` package. Under a correctly specified model OSA
#' residuals are iid standard normal, so the headline diagnostic is statistical
#' (the Q-Q plot with its standard-normal null envelope and annotated SDNR /
#' tail statistics) rather than the residual values themselves.
#'
#' The panels drawn adapt to the residual types present, all faceted by data
#' source: a Q-Q panel is always produced; composition residuals (comp / caal)
#' add a signed-bubble panel by year and bin, while aggregate (catch / index)
#' residuals add a residual-versus-year panel.
#'
#' @param x An `rceattle_osa` object from [osa_residuals()].
#' @param which Which panels to draw: any of `"qq"`, `"bubble"`, and
#'   `"resid_year"`. Defaults to the Q-Q panel plus a bubble panel when
#'   composition residuals are present, or a residual-vs-year panel otherwise.
#' @param ... Unused.
#'
#' @return Invisibly, the assembled `ggplot`/`cowplot` object. Called for its
#'   side effect of drawing the plot.
#'
#' @references Stewart, I.J., and Monnahan, C.C. 2025. Can. J. Fish. Aquat. Sci.
#'   82:1-13.
#' @seealso [osa_residuals()], [osa_diagnostics()]
#' @export
plot.rceattle_osa <- function(x, which = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot OSA residuals.")
  }
  x <- x[is.finite(x$residual), , drop = FALSE]
  x$source <- .osa_source_label(x)

  has_comp <- any(x$type %in% c("comp", "caal"))
  if (is.null(which)) which <- c("qq", if (has_comp) "bubble" else "resid_year")
  which <- match.arg(which, c("qq", "bubble", "resid_year"), several.ok = TRUE)

  panels <- list()
  if ("qq" %in% which)         panels$qq <- .osa_qqplot(x)
  if ("bubble" %in% which)     panels$bubble <- .osa_bubble_plot(x)
  if ("resid_year" %in% which) panels$resid_year <- .osa_resid_year_plot(x)

  if (length(panels) == 1L) {
    print(panels[[1]])
    return(invisible(panels[[1]]))
  }

  if (requireNamespace("cowplot", quietly = TRUE)) {
    g <- cowplot::plot_grid(plotlist = panels, ncol = 1L,
                            align = "v", axis = "lr")
    print(g)
    invisible(g)
  } else {
    for (p in panels) print(p)
    invisible(panels)
  }
}


#' Q-Q plot of OSA residuals with standard-normal null envelope
#'
#' @param osa An `rceattle_osa` data frame with a `source` column.
#' @param nsim,seed Passed to the null-envelope / tail-statistic simulation.
#' @return A `ggplot` object.
#' @keywords internal
.osa_qqplot <- function(osa, nsim = 10000, seed = 123) {
  # Per-source quantile points and SDNR/tail annotation.
  srcs <- split(osa, osa$source, drop = TRUE)
  qq <- do.call(rbind, lapply(srcs, function(g) {
    n <- nrow(g)
    data.frame(source     = g$source[1],
               theoretical = stats::qnorm(stats::ppoints(n)),
               sample      = sort(g$residual),
               stringsAsFactors = FALSE)
  }))

  ann <- do.call(rbind, lapply(srcs, function(g) {
    s <- .osa_sdnr_tails(g$residual, nsim = nsim, seed = seed)
    data.frame(source = g$source[1],
               label  = sprintf("SDNR=%.2f\n(%.2f-%.2f)",
                                s$sdnr, s$sdnr_lo, s$sdnr_hi),
               stringsAsFactors = FALSE)
  }))

  ggplot2::ggplot(qq, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    ggplot2::geom_point(colour = "#2c7fb8", alpha = 0.7) +
    ggplot2::geom_text(data = ann,
                       ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
                       hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ source, nrow = 1L) +
    ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                  title = "OSA residual Q-Q") +
    ggplot2::theme_bw(base_size = 10)
}


#' Residual-versus-year plot for OSA residuals
#'
#' @param osa An `rceattle_osa` data frame with a `source` column.
#' @return A `ggplot` object.
#' @keywords internal
.osa_resid_year_plot <- function(osa) {
  osa$sign <- ifelse(osa$residual >= 0, "positive", "negative")
  ggplot2::ggplot(osa, ggplot2::aes(x = .data$year, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_hline(yintercept = c(-2, 2), colour = "grey80",
                        linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$sign), alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c(positive = "#d7301f",
                                            negative = "#2c7fb8"),
                                 guide = "none") +
    ggplot2::facet_wrap(~ source, nrow = 1L) +
    ggplot2::labs(x = "Year", y = "OSA residual",
                  title = "OSA residual by year") +
    ggplot2::theme_bw(base_size = 10)
}


#' Bubble plot of composition OSA residuals (afscOSA styling)
#'
#' @param osa An `rceattle_osa` data frame with a `source` column. Bubbles are
#'   placed at (year, age/length bin); red = positive, blue = negative; size and
#'   transparency scale with the absolute residual; outliers (`|resid| > 3`) are
#'   drawn as triangles. Mirrors the NOAA-AFSC `afscOSA` package.
#' @return A `ggplot` object.
#' @keywords internal
.osa_bubble_plot <- function(osa) {
  osa$sign  <- ifelse(osa$residual >= 0, "positive", "negative")
  osa$shape <- ifelse(abs(osa$residual) > 3, "outlier", "normal")
  ggplot2::ggplot(osa, ggplot2::aes(x = .data$year, y = .data$age_or_length)) +
    ggplot2::geom_point(ggplot2::aes(size = abs(.data$residual),
                                     alpha = abs(.data$residual),
                                     colour = .data$sign,
                                     shape = .data$shape)) +
    ggplot2::scale_colour_manual(values = c(positive = "#d7301f",
                                            negative = "#2c7fb8"),
                                 guide = "none") +
    ggplot2::scale_shape_manual(values = c(normal = 16L, outlier = 17L),
                                guide = "none") +
    ggplot2::scale_size_continuous(range = c(0.5, 5), guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
    ggplot2::facet_wrap(~ source, nrow = 1L) +
    ggplot2::labs(x = "Year", y = "Age / length bin",
                  title = "OSA residual bubbles") +
    ggplot2::theme_bw(base_size = 10)
}


#' Build a human-readable data-source label for OSA residual rows
#' @param osa An `rceattle_osa` data frame.
#' @return Character vector of labels (one per row).
#' @keywords internal
.osa_source_label <- function(osa) {
  lab <- ifelse(is.na(osa$index_label) | osa$index_label == "",
                osa$type,
                paste0(osa$type, " (", osa$index_label, ")"))
  paste0(lab, " - fleet ", osa$fleet)
}
