# =============================================================================
# Shared helpers for the Rceattle ggplot2 plotting functions
# -----------------------------------------------------------------------------
# The exported plot_*() functions follow one standardized contract:
#   1. coerce the `Rceattle` argument to a list of fits      (.as_model_list)
#   2. assemble a tidy data frame of exactly what is plotted (per function)
#   3. render it with the shared theme + colourblind palette  (.rceattle_theme,
#      .rceattle_scale)
#   4. optionally write the figure to `file`                  (.save_ggplot)
#   5. RETURN the ggplot object (so plots are composable and the plotted data
#      is testable as `p$data`).
#
# Keeping the data assembly separate from rendering is deliberate: tests assert
# the assembled data equals the model's source quantities, so a ggplot plotter
# is verified to draw the *correct numbers*, not merely to render.
# =============================================================================

#' Coerce the `Rceattle` plotting argument to a list of fits
#'
#' Accepts a single `Rceattle` fit, a list of fits (multi-model overlay), or an
#' MSE object, and returns a plain list of `Rceattle` fits.
#'
#' @param Rceattle A fit, list of fits, or (when `mse = TRUE`) MSE object.
#' @param mse,OM When `mse = TRUE`, pull the operating model (`OM = TRUE`) or
#'   the terminal estimation model from each MSE element.
#' @return A list of `Rceattle` fits.
#' @keywords internal
.as_model_list <- function(Rceattle, mse = FALSE, OM = TRUE) {
  if (mse) {
    if (OM) {
      return(lapply(Rceattle, function(x) x$OM))
    }
    return(lapply(Rceattle, function(x) x$EM[[length(x$EM)]]))
  }
  if (inherits(Rceattle, "Rceattle")) {
    return(list(Rceattle))
  }
  Rceattle
}

#' Model display names for a list of fits
#'
#' Returns `model_names` if supplied, otherwise `"Model 1"`, `"Model 2"`, ... so
#' the colour/legend mapping always has labels.
#' @keywords internal
.model_labels <- function(model_list, model_names = NULL) {
  if (!is.null(model_names)) {
    return(rep(model_names, length.out = length(model_list)))
  }
  if (length(model_list) == 1L) {
    return("Model")
  }
  paste("Model", seq_along(model_list))
}

#' Standard Rceattle ggplot theme
#'
#' Clean, legible, and colorblind-friendly: `theme_classic` with bold facet
#' strips, a light panel border, and the legend along the bottom.
#' @param base_size Base font size.
#' @keywords internal
.rceattle_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size),
      panel.border     = ggplot2::element_rect(colour = "black", fill = NA,
                                               linewidth = 0.5),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank()
    )
}

#' Okabe-Ito colorblind-safe qualitative palette
#'
#' Eight fixed hues in a fixed order, reordered so the strongest-contrast,
#' best-separated hues come first (2-4 model overlays are the common case) and
#' the low-contrast yellow comes last. Unlike a sequential ramp (viridis), this
#' encodes *identity* (which model), not magnitude: the colours carry no implied
#' ordering and every hue holds a legible contrast against a white panel.
#'
#' Reference: Okabe & Ito (2008), "Color Universal Design".
#' @keywords internal
.okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#56B4E9", "#E69F00", "#000000", "#F0E442")

#' Discrete palette function: Okabe-Ito for up to 8 series, else interpolate
#'
#' Mirrors the arity contract of `scale_*_viridis_d` (accepts any `n`) so the
#' scale never errors on a many-model overlay; beyond 8 series it interpolates
#' the Okabe-Ito anchors (rare enough that CVD-optimality is not guaranteed).
#' @keywords internal
.okabe_ito_pal <- function(n) {
  if (n <= length(.okabe_ito)) return(.okabe_ito[seq_len(n)])
  grDevices::colorRampPalette(.okabe_ito)(n)
}

#' Normalize a user colour specification to something ggplot2 accepts
#'
#' `line_col` is a base-graphics argument, so users pass base-graphics colours:
#' `line_col = 1` and `line_col = c(1, 2, 4)` are as common in real scripts as
#' `"black"` / `"#0072B2"`. ggplot2's manual scales take names and hex only, so
#' an integer has to be resolved through the palette first -- otherwise
#' `line_col = 1` errors as an invalid colour instead of drawing black.
#'
#' @param x Colours as a character vector, or integers indexing
#'   [grDevices::palette()], or `NULL`.
#' @return A character vector of colours, or `NULL` for `NULL` input.
#' @keywords internal
#' @noRd
.as_colour <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.numeric(x)) {
    pal <- grDevices::palette()
    # Base recycles the palette by (i - 1) %% n, and treats 0 as "transparent".
    out <- ifelse(x == 0, NA_character_, pal[(as.integer(x) - 1L) %% length(pal) + 1L])
    return(out)
  }
  as.character(x)
}


#' Add the standard colorblind-safe scales to a plot
#'
#' Discrete aesthetics (series identity, e.g. model) use the Okabe-Ito
#' qualitative palette; continuous aesthetics (ordered magnitude, e.g. year)
#' use viridis. Both are colorblind-safe.
#'
#' `line_col` overrides the default. It supplies the palette for whichever
#' variable the plot already maps to colour -- which is *not* always the model:
#' see the shared argument documentation (`?rceattle-plot-args`). On a
#' continuous colour aesthetic (a year fan) there is no set of levels to name,
#' so the colours become the ramp anchors instead: one colour draws the whole
#' fan in that colour, two or more interpolate between them in order.
#'
#' @param p A ggplot.
#' @param discrete Use the discrete (Okabe-Ito) or continuous (viridis) scale.
#' @param aesthetics Which of `"colour"`/`"fill"` to add.
#' @param line_col User-supplied colours, or `NULL` for the package default.
#' @param levels Levels of the mapped discrete variable, in plotting order.
#'   Colours are recycled over these, so the caller must pass the levels of the
#'   variable actually mapped to colour.
#' @keywords internal
.rceattle_scale <- function(p, discrete = TRUE,
                            aesthetics = c("colour", "fill"),
                            line_col = NULL, levels = NULL) {
  line_col <- .as_colour(line_col)

  # Build the scale for one aesthetic. Kept as a closure so "colour" and "fill"
  # cannot drift apart -- a user's colours must reach both or the ribbon and its
  # line disagree.
  scale_for <- function(aes_name) {
    if (is.null(line_col)) {
      return(if (discrete)
               ggplot2::discrete_scale(aes_name, palette = .okabe_ito_pal)
             else if (aes_name == "colour") ggplot2::scale_colour_viridis_c()
             else ggplot2::scale_fill_viridis_c())
    }
    if (discrete) {
      # Without levels there is nothing to name the colours by; fall back to a
      # palette function, which ggplot2 applies in level order anyway.
      if (is.null(levels)) {
        return(ggplot2::discrete_scale(
          aes_name,
          palette = function(n) rep(line_col, length.out = n)))
      }
      vals <- stats::setNames(rep(line_col, length.out = length(levels)), levels)
      return(if (aes_name == "colour") ggplot2::scale_colour_manual(values = vals)
             else ggplot2::scale_fill_manual(values = vals))
    }
    # Continuous: the colours are ramp anchors. A single colour gives a flat
    # ramp, i.e. "draw it all in this colour", which is what a user asking for
    # `line_col = "black"` means.
    if (length(line_col) == 1L) {
      return(if (aes_name == "colour")
               ggplot2::scale_colour_gradient(low = line_col, high = line_col)
             else ggplot2::scale_fill_gradient(low = line_col, high = line_col))
    }
    if (aes_name == "colour") ggplot2::scale_colour_gradientn(colours = line_col)
    else ggplot2::scale_fill_gradientn(colours = line_col)
  }

  for (a in c("colour", "fill")) {
    if (a %in% aesthetics) p <- p + scale_for(a)
  }
  p
}

#' Save a ggplot to `file` if a file stem is supplied; always return the plot
#'
#' @param p A ggplot.
#' @param file File stem (no extension), or `NULL` to skip saving.
#' @param suffix Appended to `file` to form the PNG name.
#' @param width,height Figure size in inches.
#' @keywords internal
.save_ggplot <- function(p, file = NULL, suffix = "plot",
                         width = 7, height = 6.5) {
  if (!is.null(file)) {
    ggplot2::ggsave(paste0(file, "_", suffix, ".png"), plot = p,
                    width = width, height = height, dpi = 300)
  }
  p
}


#' Facet by species, but only when there is more than one
#'
#' With a single species the strip label just names the only thing on the plot.
#' Faceting is what draws that strip, so the single-species case omits the facet
#' entirely rather than blanking the strip and leaving its gap.
#'
#' Returns `NULL` when there is nothing to facet, which ggplot2 accepts as a
#' no-op inside a `+` chain -- so call sites read the same either way.
#'
#' @param df The plot data frame; its `Species` column supplies the count.
#' @param ... Passed to [ggplot2::facet_wrap()] (e.g. `scales`, `ncol`).
#' @return A `facet_wrap` layer, or `NULL` for a single-species model.
#' @keywords internal
#' @noRd
.facet_species <- function(df, ...) {
  n <- if (is.null(df$Species)) 0L else length(unique(stats::na.omit(df$Species)))
  if (n > 1L) ggplot2::facet_wrap(~ Species, ...) else NULL
}


# =============================================================================
# Standardized argument handling
# -----------------------------------------------------------------------------
# One resolver per shared argument, so every plot_*() reads `species`, `lwd`,
# `lty`, `line_col`, `minyr`/`maxyr`, `incl_proj`, `incl_mean` and `add_ci` the
# same way. Before these existed, plot_timeseries() was the only plotter that
# honoured any of them and the rest declared them as formals and ignored them.
# =============================================================================

#' Resolve the `species` / `spnames` pair to indices and labels
#'
#' `species` selects species by **index**; `spnames` supplies their **labels**.
#' The two were historically confused -- some plotters read `species` as names
#' and had no `spnames` at all -- so this accepts every spelling that used to
#' work and normalizes it:
#'
#' * `NULL` or `"all"` -- every species, in model order.
#' * numeric -- indices, validated against `1:nspp`.
#' * logical -- a mask over `1:nspp`.
#' * character matching `spnames` -- selection by name.
#' * character matching *nothing* -- back-compatibility: these are labels, not a
#'   selection, which is how `plot_selectivity()` and `plot_maturity()` used the
#'   argument. Emits a message pointing at `spnames` and keeps every species.
#' * character matching *some* names -- a selection with typos. Warns, naming
#'   the unmatched entries, and keeps the matches.
#'
#' The returned `index` preserves the order the caller asked for, so
#' `species = c(3, 1)` yields facets in that order rather than model order.
#'
#' @param models A list of `Rceattle` fits (from [.as_model_list()]).
#' @param species Species selection; see above.
#' @param spnames Species labels, length `nspp`. `NULL` takes the model's own.
#' @return `list(index = <integer>, spnames = <character, nspp>,
#'   labels = <character, length(index)>)`.
#' @keywords internal
#' @noRd
.resolve_species <- function(models, species = NULL, spnames = NULL) {
  dl   <- models[[1]]$data_list
  nspp <- dl$nspp
  if (is.null(nspp) || is.na(nspp)) nspp <- length(dl$spnames)

  # Default labels come from the model; a model without them gets positional
  # ones so the legend is never blank.
  default_names <- dl$spnames
  if (is.null(default_names) || length(default_names) == 0L) {
    default_names <- paste("Species", seq_len(nspp))
  }
  user_named <- !is.null(spnames)
  if (!user_named) spnames <- default_names
  spnames <- rep(as.character(spnames), length.out = nspp)

  idx <- NULL
  if (is.null(species)) {
    idx <- seq_len(nspp)
  } else if (is.character(species)) {
    if (length(species) == 1L && identical(tolower(species), "all")) {
      idx <- seq_len(nspp)
    } else {
      hit <- match(species, spnames)
      if (all(!is.na(hit))) {
        idx <- hit
      } else if (all(is.na(hit))) {
        # Nothing matched: the old "species = species names" spelling. Keep every
        # species and take the strings as labels, which is what they were.
        if (!user_named) spnames <- rep(as.character(species), length.out = nspp)
        idx <- seq_len(nspp)
        message("`species` matched no species name, so it is being read as ",
                "display labels. Pass labels as `spnames` and a selection ",
                "(indices or names) as `species`.")
      } else {
        warning("`species` did not match: ",
                paste(species[is.na(hit)], collapse = ", "),
                ". Plotting the ", sum(!is.na(hit)), " that did.", call. = FALSE)
        idx <- hit[!is.na(hit)]
      }
    }
  } else if (is.logical(species)) {
    idx <- which(rep(species, length.out = nspp))
  } else if (is.numeric(species)) {
    idx <- as.integer(species)
  } else {
    stop("`species` must be species indices, names, a logical mask, or \"all\".",
         call. = FALSE)
  }

  if (length(idx) == 0L) {
    stop("`species` selected no species.", call. = FALSE)
  }
  bad <- idx[is.na(idx) | idx < 1L | idx > nspp]
  if (length(bad)) {
    stop("`species` out of range: ", paste(bad, collapse = ", "),
         ". This model has ", nspp, " species (valid indices 1:", nspp, ").",
         call. = FALSE)
  }

  list(index = idx, spnames = spnames, labels = spnames[idx])
}


#' Line width / type / transparency parameters for a `geom_line()`
#'
#' Generalizes the per-model `lwd` / `lty` handling that only `plot_timeseries()`
#' had. A value that does not vary is passed as a **fixed** geom parameter, so
#' the uniform default adds no spurious legend; a value that varies is **mapped**
#' to the keying variable and given a manual scale.
#'
#' `lwd` keeps the base-graphics convention, where the package default of 3
#' renders as a standard-weight ggplot line. That is the `lwd / 3` below, and it
#' is load-bearing: every plotter converted to this helper previously hardcoded
#' `linewidth = 1`, so `lwd = 3` must keep producing exactly that or every
#' existing default figure changes.
#'
#' @param lwd Line width(s), base-graphics scale (default 3 == linewidth 1).
#' @param lty Line type(s).
#' @param alpha Fixed transparency, or `NULL` to leave it unset.
#' @param lwd_by,lty_by Name of the data column keying each aesthetic when it
#'   varies, or `NULL` if the plot has nothing to key it to.
#' @param lwd_levels,lty_levels Levels of the corresponding keying variable, in
#'   plotting order.
#' @return `list(mapping = <aes or NULL>, args = <fixed geom params>,
#'   scales = <list of scale objects>)`.
#' @keywords internal
#' @noRd
.rce_line_params <- function(lwd = 3, lty = 1, alpha = NULL,
                             lwd_by = NULL, lwd_levels = NULL,
                             lty_by = NULL, lty_levels = NULL) {
  args   <- list()
  scales <- list()

  lwd <- if (length(lwd)) lwd else 3
  lty <- if (length(lty)) lty else 1

  # An aesthetic can only be mapped if there is a variable to map it to AND
  # levels to name the values by. Otherwise the first value is used for all,
  # which is what the caller gets when they vary an aesthetic on a plot that has
  # no corresponding series.
  map_lwd <- length(unique(lwd)) > 1L && !is.null(lwd_by) && length(lwd_levels)
  map_lty <- length(unique(lty)) > 1L && !is.null(lty_by) && length(lty_levels)

  if (length(unique(lwd)) > 1L && !map_lwd) {
    warning("`lwd` varies but this plot has no per-series line width to map it ",
            "to; using lwd = ", lwd[1], " throughout.", call. = FALSE)
  }
  if (length(unique(lty)) > 1L && !map_lty) {
    warning("`lty` varies but this plot has no per-series line type to map it ",
            "to; using lty = ", lty[1], " throughout.", call. = FALSE)
  }

  if (map_lwd) {
    scales <- c(scales, list(ggplot2::scale_linewidth_manual(
      values = stats::setNames(
        rep(lwd, length.out = length(lwd_levels)) / 3, lwd_levels))))
  } else {
    args$linewidth <- lwd[1] / 3
  }
  if (map_lty) {
    scales <- c(scales, list(ggplot2::scale_linetype_manual(
      values = stats::setNames(
        rep(lty, length.out = length(lty_levels)), lty_levels))))
  } else {
    args$linetype <- lty[1]
  }
  if (!is.null(alpha)) args$alpha <- alpha

  # Built explicitly rather than by deleting elements, so the quosures capture
  # the column names from this frame.
  mapping <- if (map_lwd && map_lty) {
    ggplot2::aes(linewidth = .data[[lwd_by]], linetype = .data[[lty_by]])
  } else if (map_lwd) {
    ggplot2::aes(linewidth = .data[[lwd_by]])
  } else if (map_lty) {
    ggplot2::aes(linetype = .data[[lty_by]])
  } else {
    NULL
  }

  list(mapping = mapping, args = args, scales = scales)
}


#' Add a `geom_line()` built from [.rce_line_params()], plus its scales
#'
#' The three-line incantation every converted plotter would otherwise repeat.
#' @param p A ggplot.
#' @param lp The list returned by [.rce_line_params()].
#' @param ... Further arguments to [ggplot2::geom_line()].
#' @keywords internal
#' @noRd
.rce_add_line <- function(p, lp, ...) {
  p <- p + do.call(ggplot2::geom_line,
                   c(list(mapping = lp$mapping), lp$args, list(...)))
  for (s in lp$scales) p <- p + s
  p
}


#' Year-axis limits, for plots with Year on the x axis
#'
#' Clips rather than filters, so the returned plot's `$data` still holds the
#' whole series -- the accuracy tests read `p$data` and compare it against the
#' model's quantities, which a filtered frame would silently narrow.
#'
#' Use [.rce_year_filter()] instead where Year is a grouping or colour variable
#' rather than the x axis; there, dropping rows is the point.
#'
#' @param minyr,maxyr First / last year to show, or `NULL` for no limit.
#' @return A `coord_cartesian` layer, or `NULL` when neither limit is set.
#' @keywords internal
#' @noRd
.rce_year_limits <- function(minyr = NULL, maxyr = NULL) {
  if (is.null(minyr) && is.null(maxyr)) return(NULL)
  ggplot2::coord_cartesian(xlim = c(minyr, maxyr))
}


#' Restrict a plot frame to a year range, for plots where Year is not the x axis
#'
#' Thins a year fan (selectivity, mortality-at-age) to the requested window.
#' @param df A plot data frame.
#' @param minyr,maxyr First / last year to keep, or `NULL` for no limit.
#' @param year Name of the year column.
#' @keywords internal
#' @noRd
.rce_year_filter <- function(df, minyr = NULL, maxyr = NULL, year = "Year") {
  if (is.null(df[[year]])) return(df)
  keep <- rep(TRUE, nrow(df))
  if (!is.null(minyr)) keep <- keep & df[[year]] >= minyr
  if (!is.null(maxyr)) keep <- keep & df[[year]] <= maxyr
  if (!any(keep)) {
    stop("No years left to plot between minyr = ", minyr, " and maxyr = ",
         maxyr, ".", call. = FALSE)
  }
  df[keep, , drop = FALSE]
}


#' Dashed divider at the last hindcast year
#'
#' Marks where the projection starts. Consolidates the identical `geom_vline()`
#' that six plotters carried inline.
#'
#' @param models A list of `Rceattle` fits.
#' @param incl_proj Whether projection years are being drawn.
#' @return A `geom_vline` layer, or `NULL`.
#' @keywords internal
#' @noRd
.rce_proj_divider <- function(models, incl_proj = FALSE) {
  if (!isTRUE(incl_proj)) return(NULL)
  endyr <- models[[length(models)]]$data_list$endyr
  if (is.null(endyr) || is.na(endyr)) return(NULL)
  ggplot2::geom_vline(xintercept = endyr, linetype = 2, colour = "grey50")
}


#' Horizontal long-term mean line(s)
#'
#' `incl_mean` was documented as "include time series mean as horizontal line"
#' on four plotters and implemented on none. The mean is taken over **hindcast
#' years only**: a mean that folded in the projection would not be a historical
#' reference, and would move when the projection horizon changed.
#'
#' @param df The plot data frame.
#' @param incl_mean Draw it?
#' @param by Columns defining a group (a facet and/or a series), e.g.
#'   `c("Species", "Model")`.
#' @param value Name of the value column.
#' @param hind_endyr Last hindcast year; rows after it are excluded.
#' @param year Name of the year column.
#' @return A `geom_hline` layer, or `NULL`.
#' @keywords internal
#' @noRd
.rce_mean_line <- function(df, incl_mean = FALSE, by = NULL, value = "value",
                           hind_endyr = NULL, year = "Year") {
  if (!isTRUE(incl_mean)) return(NULL)
  hind <- df
  if (!is.null(hind_endyr) && !is.null(hind[[year]])) {
    hind <- hind[hind[[year]] <= hind_endyr, , drop = FALSE]
  }
  hind <- hind[!is.na(hind[[value]]), , drop = FALSE]
  if (nrow(hind) == 0L) return(NULL)

  by <- by[by %in% names(hind)]
  agg <- if (length(by)) {
    stats::aggregate(hind[[value]], by = hind[by], FUN = mean)
  } else {
    data.frame(x = mean(hind[[value]]))
  }
  names(agg)[ncol(agg)] <- ".mean"

  # Colour follows Model when the plot separates models, so the mean line sits
  # under the series it belongs to rather than in an unrelated colour.
  mapping <- if ("Model" %in% by) {
    ggplot2::aes(yintercept = .data$.mean, colour = .data$Model)
  } else {
    ggplot2::aes(yintercept = .data$.mean)
  }
  args <- list(data = agg, mapping = mapping, inherit.aes = FALSE,
               linetype = 3)
  if (!("Model" %in% by)) args$colour <- "grey40"
  do.call(ggplot2::geom_hline, args)
}


#' Standard errors for an ADREPORTed series, by name
#'
#' Factored out of `plot_timeseries()` so any plotter drawing a confidence
#' interval reads `sdrep` the same way.
#'
#' @param model An `Rceattle` fit.
#' @param name Name of the ADREPORTed series.
#' @param n_need How many values the caller needs.
#' @return The first `n_need` standard errors, or `NULL` when the fit does not
#'   carry them (no `sdreport`, or the series is `REPORT()`-only).
#' @keywords internal
#' @noRd
.rce_series_sd <- function(model, name, n_need) {
  sdrep <- model$sdrep
  if (is.null(sdrep) || is.null(sdrep$value)) return(NULL)
  rows <- which(names(sdrep$value) == name)
  if (length(rows) < n_need) return(NULL)
  sdrep$sd[rows][seq_len(n_need)]
}


#' Decline an `add_ci` request that the model cannot support
#'
#' Some plotted quantities are products or sums of ADREPORTed series (ration is
#' consumption x biomass; the M2 proportions are ratios), so their standard
#' error needs a covariance the fit does not carry by default. Silently ignoring
#' `add_ci = TRUE` is what let these arguments look functional for so long, so
#' say it once instead.
#'
#' @param add_ci The user's request.
#' @param quantity Name of the plotted quantity, for the message.
#' @param reason Why it cannot be drawn.
#' @return `FALSE`, always -- so the caller can write `add_ci <- .rce_no_ci(...)`.
#' @keywords internal
#' @noRd
.rce_no_ci <- function(add_ci, quantity, reason) {
  if (isTRUE(add_ci)) {
    warning("`add_ci = TRUE` is not available for ", quantity, ": ", reason,
            ". Drawing no confidence interval.", call. = FALSE)
  }
  FALSE
}


#' Shared plotting arguments
#'
#' The `plot_*()` functions take a common vocabulary, documented here once and
#' inherited by each of them.
#'
#' @section How `line_col` and `lty` are applied:
#'
#' They supply values for whichever **discrete variable the plot already
#' encodes with that aesthetic**, matched in level order -- which is not always
#' the model. In [plot_b_eaten_prop()] colour separates predators; in
#' [plot_ration()] line type separates sexes. Check what the figure's legend
#' names before assuming a colour vector maps to your models.
#'
#' Where colour encodes a continuous variable -- the year fan in
#' [plot_selectivity()] and [plot_mortality()] -- `line_col` supplies the ramp
#' anchors instead: one colour draws the fan in that colour, several interpolate
#' between them.
#'
#' @param Rceattle A single [fit_mod()] object or a list of them (overlaid).
#' @param file Optional file stem; the figure is written to
#'   `<file>_<suffix>.png` if given.
#' @param model_names Legend labels for the models.
#' @param species Species to include, as indices (`c(1, 3)`), names, a logical
#'   mask, or `"all"`. Default `NULL` plots every species. Species **labels**
#'   belong in `spnames`; a character `species` that matches no species name is
#'   read as labels for back-compatibility.
#' @param spnames Species labels, length `nspp`. Default: the model's own.
#' @param line_col Line colours; names, hex, or base-graphics integers. `NULL`
#'   uses the colorblind-safe Okabe-Ito palette. See above for what it maps to.
#' @param lwd Line width on the base-graphics scale: the default `3` renders as
#'   a standard-weight ggplot line. A vector varies it across series.
#' @param lty Line type. A vector varies it across series.
#' @param alpha Transparency of confidence ribbons and shaded areas.
#' @param add_ci Add a 95% confidence interval. Only available where the
#'   plotted quantity carries standard errors; warns and draws none otherwise.
#' @param minyr,maxyr First / last year to plot.
#' @param incl_proj Include the projection years, with a dashed divider at the
#'   last hindcast year.
#' @param incl_mean Add a horizontal line at the hindcast mean of each series.
#' @param width,height Saved figure size in inches.
#' @param right_adj,top_adj,mod_cex,legend.pos,single.plots,save,theta
#'   Base-graphics arguments from before the ggplot2 rewrite. Accepted so
#'   existing scripts keep running, and ignored. The returned object is a
#'   ggplot, so margins, legend placement and text size are set on it directly
#'   (e.g. `p + ggplot2::theme(legend.position = "right")`).
#' @return A `ggplot` object.
#' @name rceattle-plot-args
#' @keywords internal
NULL
