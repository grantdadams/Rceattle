#' Plot one-step-ahead (OSA) residual diagnostics
#'
#' @description
#' Diagnostic plots for an `rceattle_osa` object (from [osa_residuals()]),
#' following the recommendations of Stewart and Monnahan (2025) and the styling
#' of the NOAA-AFSC `afscOSA` package. Under a correctly specified model OSA
#' residuals are iid standard normal, so the headline diagnostic is statistical
#' (the Q-Q plot with its annotated SDNR / tail statistics).
#'
#' Up to two separate figures are drawn, depending on which data are present:
#' \enumerate{
#'   \item **Aggregate** (`index` / `catch`): a Q-Q panel faceted by data source.
#'     These series have no age/length bin, so no bubble plots are drawn.
#'   \item **Composition** (`comp` / `caal`): a Q-Q panel, a signed OSA-residual
#'     bubble panel, and a signed Pearson-residual bubble panel (the Pearson
#'     residuals carried on the `rceattle_osa` object). By default age-based bins
#'     (age composition and conditional age-at-length) are shown in the left
#'     column and length-based bins in the right column, each with its own bin
#'     axis; set `combine = FALSE` to draw the age and length composition as two
#'     separate figures instead (useful when a model has many fleets/species).
#' }
#' Panel headers use the fleet name from `fleet_control`. Process residuals
#' (from [process_residuals()]) are drawn as a Q-Q panel plus a
#' residual-by-year panel.
#'
#' @param x An `rceattle_osa` object from [osa_residuals()] or
#'   [process_residuals()].
#' @param source Data source(s) to plot: any of `"index"`, `"catch"`, `"comp"`,
#'   `"caal"`, `"diet"`, or `"all"` (default). Mirrors the `source` argument of
#'   [residuals.Rceattle()]; filters which figures are produced.
#' @param species Optional species code(s) to include (matched against the
#'   `species` column). Default `NULL` keeps all species.
#' @param combine Logical. When `TRUE` (default), age and length composition
#'   share one figure (age in the left column, length in the right). When
#'   `FALSE`, they are drawn as separate `composition_age` / `composition_length`
#'   figures.
#' @param ... Unused.
#'
#' @return Invisibly, a named list of the assembled `ggplot` / `cowplot`
#'   objects (some of `aggregate`, `composition` / `composition_age` /
#'   `composition_length`, and `process`). Called for its side effect of drawing
#'   the plot(s).
#'
#' @references Stewart, I.J., and Monnahan, C.C. 2025. Can. J. Fish. Aquat. Sci.
#'   82:1-13.
#' @seealso [osa_residuals()], [osa_diagnostics()]
#' @export
plot.rceattle_osa <- function(x, source = "all", species = NULL,
                              combine = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot OSA residuals.")
  }
  pearson  <- attr(x, "pearson")
  nages    <- attr(x, "nages")      # per-species, for joint-sex bin rebasing
  nlengths <- attr(x, "nlengths")

  # ---- Subset by data source and species (like residuals.Rceattle()) ----
  # TODO(review): process-residual objects (from process_residuals()) carry
  # source names like "recruitment"/"initial"/"catchability", which are not in
  # valid_src, so plot(x, source = "recruitment") errors in match.arg(); only
  # source = "all" works for those. Allow process names when x is a procres.
  valid_src  <- c("index", "catch", "comp", "caal", "diet")
  keep_types <- if (identical(source, "all")) unique(x$source) else
    match.arg(source, valid_src, several.ok = TRUE)
  x <- x[x$source %in% keep_types, , drop = FALSE]
  if (!is.null(species)) x <- x[x$species %in% species, , drop = FALSE]
  if (!is.null(pearson)) {
    pearson <- pearson[pearson$source %in% keep_types, , drop = FALSE]
    if (!is.null(species)) {
      pearson <- pearson[pearson$species %in% species, , drop = FALSE]
    }
  }

  x <- x[is.finite(x$residual), , drop = FALSE]
  if (nrow(x) == 0) {
    warning("No finite residuals to plot for the requested source / species.")
    return(invisible(NULL))
  }

  agg  <- x[x$source %in% c("index", "catch"), , drop = FALSE]
  comp <- x[x$source %in% c("comp", "caal", "diet"), , drop = FALSE]
  proc <- x[!x$source %in% valid_src, , drop = FALSE]

  plots <- list()

  # Aggregate (index / catch): Q-Q only -- no age/length bin to bubble.
  if (nrow(agg) > 0) {
    agg$source <- .osa_source_label(agg)
    plots$aggregate <- .osa_qqplot(agg)
  }

  # Composition (comp / caal): Q-Q + OSA bubbles + Pearson bubbles. One combined
  # figure (age | length columns) by default, or separate figures per bin type.
  if (nrow(comp) > 0) {
    plots <- c(plots, .osa_composition_figures(comp, pearson, combine = combine,
                                               nages = nages, nlengths = nlengths))
  }

  # Process residuals: Q-Q + residual-by-year.
  if (nrow(proc) > 0) {
    proc$source <- paste0(proc$source,
                          ifelse(is.na(proc$species), "",
                                 paste0(" - sp ", proc$species)))
    plots$process <- .osa_stack(list(.osa_qqplot(proc),
                                     .osa_resid_year_plot(proc)))
  }

  for (g in plots) print(g)
  invisible(plots)
}


#' Composition OSA figure(s): Q-Q + OSA bubbles + Pearson bubbles (age / length)
#'
#' @param comp The composition rows of an `rceattle_osa` object.
#' @param pearson Matching Pearson residuals (the `"pearson"` attribute of the
#'   `rceattle_osa` object), or `NULL`.
#' @param combine When `TRUE`, return a single `composition` figure with age
#'   bins in the left column and length bins in the right; when `FALSE`, return
#'   separate `composition_age` and `composition_length` figures.
#' @param nages,nlengths Per-species bin counts (the `rceattle_osa` `"nages"` /
#'   `"nlengths"` attributes), used to split joint-sex (Sex == 3) bins onto a
#'   single age/length axis, matching [plot_comp()].
#' @return A named list of `cowplot`/`ggplot` objects.
#' @keywords internal
.osa_composition_figures <- function(comp, pearson = NULL, combine = TRUE,
                                     nages = NULL, nlengths = NULL) {
  comp$source <- .osa_source_label(comp)
  comp$.side  <- .osa_bin_side(comp$index_label)
  comp        <- .osa_jointsex(comp, nages, nlengths)

  # The attached Pearson residuals already carry the OSA column names
  # (osa_residuals() renames them once, where the attribute is attached), so
  # only the derived index label is added here.
  pear <- NULL
  if (!is.null(pearson) && nrow(pearson) > 0) {
    pear <- pearson
    pear$index_label <- ifelse(pear$source == "caal", "age",
                               ifelse(!is.na(pear$index_label_code) &
                                        pear$index_label_code == 1,
                                      "length", "age"))
    pear <- pear[is.finite(pear$residual), , drop = FALSE]
    pear$source <- .osa_source_label(pear)
    pear$.side  <- .osa_bin_side(pear$index_label)
    pear        <- .osa_jointsex(pear, nages, nlengths)
  }

  age_fig <- .osa_comp_side(comp, pear, "age")
  len_fig <- .osa_comp_side(comp, pear, "length")

  if (!combine) {
    out <- list()
    if (!is.null(age_fig)) out$composition_age    <- age_fig
    if (!is.null(len_fig)) out$composition_length <- len_fig
    return(out)
  }

  cols <- Filter(Negate(is.null), list(age_fig, len_fig))
  if (length(cols) == 0L) return(list())
  fig <- if (length(cols) == 1L) {
    cols[[1]]
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(plotlist = cols, ncol = length(cols))
  } else {
    cols[[1]]   # cowplot unavailable: fall back to the first column
  }
  list(composition = fig)
}


#' One bin-side (age or length) composition column: Q-Q + OSA + Pearson bubbles
#'
#' @param comp Composition rows with `source` and `.side` columns.
#' @param pear Reshaped Pearson rows with `source` and `.side` columns, or `NULL`.
#' @param side `"age"` or `"length"`.
#' @return A stacked `cowplot`/`ggplot` object, or `NULL` if no rows on that side.
#' @keywords internal
.osa_comp_side <- function(comp, pear, side) {
  cs <- comp[comp$.side == side, , drop = FALSE]
  if (nrow(cs) == 0) return(NULL)
  ylab   <- if (side == "age") "Age bin" else "Length bin"
  panels <- list(
    .osa_qqplot(cs),
    .osa_bubble_plot(cs, ylab = ylab, title = "OSA residuals"))
  if (!is.null(pear)) {
    ps <- pear[pear$.side == side, , drop = FALSE]
    if (nrow(ps) > 0) {
      panels[[length(panels) + 1L]] <-
        .osa_bubble_plot(ps, ylab = ylab, title = "Pearson residuals")
    }
  }
  .osa_stack(panels)
}


#' Stack ggplot panels vertically (cowplot if available)
#' @param panels A list of ggplot objects.
#' @keywords internal
.osa_stack <- function(panels) {
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 1L) return(panels[[1]])
  if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(plotlist = panels, ncol = 1L, align = "v", axis = "lr")
  } else {
    for (p in panels) print(p)
    panels[[1]]
  }
}


#' Map an `index_label` to the bin side ("age" left, "length" right)
#' @param index_label Character vector (`"age"`, `"length"`, or `NA`).
#' @keywords internal
.osa_bin_side <- function(index_label) {
  ifelse(!is.na(index_label) & index_label == "length", "length", "age")
}


#' Split joint-sex (Sex == 3) composition bins onto a single age/length axis
#'
#' Joint-sex compositions stack females in bins `1..nbin` and males in bins
#' `nbin+1..2*nbin` (where `nbin` is `nages` or `nlengths` for the species).
#' This re-bases the male bins to `1..nbin` and tags the source label by sex so
#' males and females face the same bin axis -- matching [plot_comp()]. Rows with
#' Sex != 3 (single-sex or combined) are returned unchanged.
#' @param df A data frame with `species`, `sex`, `index_label`, `age_length_bin`,
#'   and `source` columns.
#' @param nages,nlengths Per-species bin counts (or `NULL` to skip the split).
#' @keywords internal
.osa_jointsex <- function(df, nages, nlengths) {
  if (is.null(df) || nrow(df) == 0 || is.null(nages) ||
      is.null(df$sex) || is.null(df$species)) {
    return(df)
  }
  bins_per_sex <- ifelse(!is.na(df$index_label) & df$index_label == "length",
                         nlengths[df$species], nages[df$species])
  joint <- !is.na(df$sex) & df$sex == 3 & !is.na(bins_per_sex)
  male  <- joint & df$age_length_bin > bins_per_sex
  df$age_length_bin[male] <- df$age_length_bin[male] - bins_per_sex[male]
  df$source <- paste0(df$source,
                      ifelse(joint, ifelse(male, " - male", " - female"), ""))
  df
}


#' Q-Q plot of OSA residuals with standard-normal null envelope
#'
#' @param osa An `rceattle_osa` data frame with a `source` column.
#' @param nsim,seed Passed to the SDNR / tail-statistic annotation.
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


#' Residual-versus-year plot for OSA / process residuals
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
    ggplot2::labs(x = "Year", y = "Residual",
                  title = "Residual by year") +
    ggplot2::theme_bw(base_size = 10)
}


#' Bubble plot of composition residuals (afscOSA styling)
#'
#' @param osa A data frame with `source`, `year`, `age_length_bin`, and
#'   `residual` columns. Bubbles are placed at (year, age/length bin); red =
#'   positive, blue = negative; size and transparency scale with the absolute
#'   residual; outliers (`|resid| > 3`) are drawn as triangles.
#' @param ylab Y-axis label (e.g. `"Age bin"` or `"Length bin"`).
#' @param title Panel title.
#' @return A `ggplot` object.
#' @keywords internal
.osa_bubble_plot <- function(osa, ylab = "Bin", title = "OSA residuals") {
  osa$sign  <- ifelse(osa$residual >= 0, "positive", "negative")
  osa$shape <- ifelse(abs(osa$residual) > 3, "outlier", "normal")
  ggplot2::ggplot(osa, ggplot2::aes(x = .data$year, y = .data$age_length_bin)) +
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
    ggplot2::labs(x = "Year", y = ylab, title = title) +
    ggplot2::theme_bw(base_size = 10)
}


#' Build a two-row data-source label for residual panel headers
#'
#' The first row is the fleet name from `fleet_control` (falling back to the
#' fleet code); the second row is the data type, tagged with whether composition
#' bins are ages or lengths. The two rows are separated by a newline so they
#' render as a two-line facet strip.
#' @param osa A data frame with `source`, `fleet`, and (optionally) `fleet_name`
#'   and `index_label` columns.
#' @return Character vector of labels (one per row).
#' @keywords internal
.osa_source_label <- function(osa) {
  flt <- if (!is.null(osa$fleet_name)) {
    ifelse(is.na(osa$fleet_name), paste("fleet", osa$fleet), osa$fleet_name)
  } else {
    paste("fleet", osa$fleet)
  }
  lab <- if (!is.null(osa$index_label)) {
    ifelse(is.na(osa$index_label) | osa$index_label == "",
           osa$source, paste0(osa$source, " (", osa$index_label, ")"))
  } else {
    osa$source
  }
  paste0(flt, "\n", lab)   # row 1: fleet name, row 2: data type
}
