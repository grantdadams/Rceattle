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
#' `"black"` / `"#0072B2"`. R's colour parser already reads the *string* `"1"` as
#' a palette index, so a single small integer happens to survive; resolving it
#' here makes that explicit and fixes the cases that do not, namely indices past
#' the end of the palette (base wraps, ggplot2 does not).
#'
#' Invalid colours fail here, where the argument can be named, rather than deep
#' inside `ggplot_build()` -- or, worse, silently: an `NA` colour renders as a
#' missing line, so a bad `line_col` would otherwise produce a blank panel with
#' no message.
#'
#' @param x Colours as a character vector, or integers indexing
#'   [grDevices::palette()], or `NULL`.
#' @return A character vector of colours, or `NULL` for `NULL`/empty input.
#' @keywords internal
#' @noRd
.as_colour <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NULL)
  if (is.logical(x) && all(is.na(x))) {
    stop("`line_col` is NA. Pass colours, palette indices, or NULL for the ",
         "default palette.", call. = FALSE)
  }
  if (is.numeric(x)) {
    if (anyNA(x) || any(x < 0)) {
      stop("`line_col` must be non-negative palette indices; got ",
           paste(utils::head(x, 5), collapse = ", "), ".", call. = FALSE)
    }
    pal <- grDevices::palette()
    # Base wraps the palette by (i - 1) %% n and reads 0 as the background,
    # i.e. transparent.
    return(ifelse(x == 0, NA_character_,
                  pal[(as.integer(x) - 1L) %% length(pal) + 1L]))
  }
  if (!is.character(x) && !is.factor(x)) {
    stop("`line_col` must be colour names, hex codes, or palette indices; got ",
         class(x)[1], ".", call. = FALSE)
  }
  x <- as.character(x)
  # Reject unrenderable colours now. grDevices::col2rgb() is the same parser
  # ggplot2 eventually uses, so anything it accepts will draw.
  bad <- vapply(x, function(one) {
    if (is.na(one)) return(FALSE)   # NA is a legitimate "draw nothing"
    inherits(try(grDevices::col2rgb(one), silent = TRUE), "try-error")
  }, logical(1))
  if (any(bad)) {
    stop("`line_col` is not a valid colour: ",
         paste(unique(x[bad]), collapse = ", "), ".", call. = FALSE)
  }
  x
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
#'   Recycled over the mapped variable's levels in plotting order.
#' @keywords internal
.rceattle_scale <- function(p, discrete = TRUE,
                            aesthetics = c("colour", "fill"),
                            line_col = NULL) {
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
      # Applied in level order by a palette function rather than by naming the
      # values. A named scale_*_manual() silently renders every series grey50
      # when the names do not match the data's levels, so the caller would have
      # to know the factor's levels exactly to use it safely -- and the
      # documented contract is level order anyway.
      return(ggplot2::discrete_scale(
        aes_name,
        palette = function(n) rep(line_col, length.out = n)))
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
  if (!length(nspp) || nspp < 1L) {
    # Blame the model, not `species` -- otherwise a malformed fit reports
    # "`species` selected no species" for an argument the caller never passed.
    stop("This model reports no species (`data_list$nspp` is ",
         format(dl$nspp), ").", call. = FALSE)
  }

  # Default labels come from the model; a model without them gets positional
  # ones so the legend is never blank.
  default_names <- dl$spnames
  if (is.null(default_names) || length(default_names) != nspp) {
    default_names <- paste("Species", seq_len(nspp))
  }
  user_named <- !is.null(spnames)
  if (!user_named) spnames <- default_names
  # Recycling a short `spnames` would label species 3 with species 1's name --
  # a plausible-looking wrong answer on a figure that goes into an assessment
  # document. The un-helped behaviour (an NA strip) is at least visibly broken,
  # so refuse instead of quietly inventing a label.
  spnames <- as.character(spnames)
  if (length(spnames) != nspp) {
    stop("`spnames` has ", length(spnames), " name(s) but the model has ",
         nspp, " species. Supply one name per species.", call. = FALSE)
  }

  idx <- NULL
  if (is.null(species)) {
    idx <- seq_len(nspp)
  } else if (is.factor(species)) {
    # A factor of species names is a plausible input (e.g. carried out of a
    # data frame); treat it as the character vector it prints as.
    return(.resolve_species(models, as.character(species), spnames))
  } else if (is.character(species)) {
    if (length(species) == 1L && identical(tolower(species), "all")) {
      idx <- seq_len(nspp)
    } else {
      # Match against the model's own names as well as any relabelling, so
      # renaming species for display does not make them unselectable.
      hit <- match(species, spnames)
      hit[is.na(hit)] <- match(species[is.na(hit)], default_names)
      if (all(!is.na(hit))) {
        idx <- hit
      } else if (all(is.na(hit)) && length(species) == nspp && !user_named) {
        # The legacy spelling: `species` held one display label per species,
        # which is how plot_selectivity() and plot_maturity() used it. Only
        # this exact shape -- one string per species, no separate `spnames` --
        # is read that way. Anything else is a selection, and a selection that
        # matches nothing is a typo, not a relabelling: silently plotting every
        # species under the user's misspelled names is the worst outcome
        # available, so it is an error.
        spnames <- as.character(species)
        idx <- seq_len(nspp)
        message("`species` matched no species name and supplies one string per ",
                "species, so it is being read as display labels. Pass labels ",
                "as `spnames` and a selection as `species`.")
      } else if (all(is.na(hit))) {
        stop("`species` matched no species name: ",
             paste(species, collapse = ", "), ". This model has: ",
             paste(default_names, collapse = ", "), ".", call. = FALSE)
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
#'   varies, or `NULL` if the plot has nothing to key it to. Values are applied
#'   to that column's levels in plotting order.
#' @return `list(mapping = <aes or NULL>, args = <fixed geom params>,
#'   scales = <list of scale objects>)`.
#' @keywords internal
#' @noRd
.rce_line_params <- function(lwd = 3, lty = 1, alpha = NULL,
                             lwd_by = NULL, lty_by = NULL) {
  args   <- list()
  scales <- list()

  # An NA or negative width does not draw: the line silently disappears, or the
  # device errors deep inside the render. Catch it where the argument has a name.
  if (is.null(lwd) || !length(lwd)) lwd <- 3
  if (is.null(lty) || !length(lty)) lty <- 1
  if (!is.numeric(lwd) || anyNA(lwd) || any(lwd < 0)) {
    stop("`lwd` must be non-negative numbers; got ",
         paste(utils::head(format(lwd), 5), collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(lty)) stop("`lty` must not be NA.", call. = FALSE)

  # An aesthetic can only be mapped if the plot has a variable to key it to.
  # Otherwise the first value is used throughout, which is what a caller gets
  # for varying an aesthetic on a plot with no corresponding series.
  map_lwd <- length(unique(lwd)) > 1L && !is.null(lwd_by)
  map_lty <- length(unique(lty)) > 1L && !is.null(lty_by)

  if (length(unique(lwd)) > 1L && !map_lwd) {
    warning("`lwd` varies but this plot draws one line width; using lwd = ",
            lwd[1], " throughout.", call. = FALSE)
  }
  if (length(unique(lty)) > 1L && !map_lty) {
    warning("`lty` varies but this plot draws one line type; using lty = ",
            lty[1], " throughout.", call. = FALSE)
  }

  # Values are applied in level order by a palette function rather than by
  # naming them: a named scale_*_manual() whose names miss the data's levels
  # renders every line blank, and the caller cannot be expected to know the
  # factor's levels.
  if (map_lwd) {
    scales <- c(scales, list(ggplot2::discrete_scale(
      "linewidth", palette = function(n) rep(lwd, length.out = n) / 3)))
  } else {
    args$linewidth <- lwd[1] / 3
  }
  if (map_lty) {
    scales <- c(scales, list(ggplot2::discrete_scale(
      "linetype", palette = function(n) rep(lty, length.out = n))))
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
#' One-sided limits are the common case (`plot_ssb(fit, minyr = 1990)`), so the
#' unset end is passed as `NA` -- `coord_cartesian()` reads that as "take it from
#' the data". Passing `NULL` through `c()` instead drops the element and yields a
#' length-1 `xlim`, which ggplot2 rejects outright.
#'
#' @param minyr,maxyr First / last year to show, or `NULL` for no limit.
#' @return A `coord_cartesian` layer, or `NULL` when neither limit is set.
#' @keywords internal
#' @noRd
.rce_year_limits <- function(minyr = NULL, maxyr = NULL) {
  if (is.null(minyr) && is.null(maxyr)) return(NULL)
  ggplot2::coord_cartesian(
    xlim = c(if (is.null(minyr)) NA_real_ else as.numeric(minyr),
             if (is.null(maxyr)) NA_real_ else as.numeric(maxyr)))
}


#' Restrict a plot frame to a year range, for plots where Year is not the x axis
#'
#' Thins a year fan (selectivity, mortality-at-age) to the requested window.
#'
#' Unlike [.rce_year_limits()], which clips the axis and keeps every row, this
#' drops rows -- that is the point where Year is the colour or grouping variable.
#' Emptying the frame entirely is an error rather than a blank panel, because a
#' figure with nothing in it is never what the caller wanted; a window that only
#' empties *some* panels just loses those rows, matching `plot_timeseries()`'s
#' convention of clipping a `maxyr` past the end of a model's data.
#'
#' @param df A plot data frame.
#' @param minyr,maxyr First / last year to keep, or `NULL` for no limit.
#' @param year Name of the year column.
#' @keywords internal
#' @noRd
.rce_year_filter <- function(df, minyr = NULL, maxyr = NULL, year = "Year") {
  if (is.null(df[[year]])) return(df)
  if (!is.null(minyr) && (length(minyr) != 1L || is.na(minyr))) {
    stop("`minyr` must be a single year.", call. = FALSE)
  }
  if (!is.null(maxyr) && (length(maxyr) != 1L || is.na(maxyr))) {
    stop("`maxyr` must be a single year.", call. = FALSE)
  }
  keep <- rep(TRUE, nrow(df))
  if (!is.null(minyr)) keep <- keep & df[[year]] >= minyr
  if (!is.null(maxyr)) keep <- keep & df[[year]] <= maxyr
  # A comparison against an NA year gives NA, which subsets to a fabricated
  # all-NA row rather than dropping it.
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) {
    stop("No years left to plot between minyr = ",
         if (is.null(minyr)) "NULL" else minyr, " and maxyr = ",
         if (is.null(maxyr)) "NULL" else maxyr, ".", call. = FALSE)
  }
  df[keep, , drop = FALSE]
}


#' Dashed divider at the last hindcast year
#'
#' Marks where the projection starts. Consolidates the identical `geom_vline()`
#' that six plotters carried inline.
#'
#' Takes the **latest** hindcast year across the models rather than the last
#' model's. On a retrospective peel list the models end in different years, and
#' keying on whichever happens to be last puts the divider mid-hindcast and
#' labels real data as projection.
#'
#' @param models A list of `Rceattle` fits.
#' @param incl_proj Whether projection years are being drawn.
#' @return A `geom_vline` layer, or `NULL`.
#' @keywords internal
#' @noRd
.rce_proj_divider <- function(models, incl_proj = FALSE) {
  if (!isTRUE(incl_proj)) return(NULL)
  endyr <- vapply(models, function(m) {
    e <- m$data_list$endyr
    if (is.null(e) || !length(e)) NA_real_ else as.numeric(e[1])
  }, numeric(1))
  if (all(is.na(endyr))) return(NULL)
  ggplot2::geom_vline(xintercept = max(endyr, na.rm = TRUE),
                      linetype = 2, colour = "grey50")
}


#' Horizontal long-term mean line(s)
#'
#' `incl_mean` was documented as "include time series mean as horizontal line"
#' on four plotters and implemented on none. The mean is taken over **hindcast
#' years only**: a mean that folded in the projection would not be a historical
#' reference, and would move when the projection horizon changed. `hind_endyr`
#' is therefore required, not optional -- a caller that forgot it would silently
#' get the projection-contaminated mean this exists to avoid.
#'
#' The mean line takes a colour only when the caller says colour encodes the
#' model. Inferring that from `by` is wrong: this helper's own callers map
#' colour to predator or to species, and adding a `colour = Model` layer to such
#' a plot trains model names into that scale's legend.
#'
#' @param df The plot data frame.
#' @param incl_mean Draw it?
#' @param by Columns defining a group (a facet and/or a series), e.g.
#'   `c("Species", "Model")`.
#' @param value Name of the value column.
#' @param hind_endyr Last hindcast year; rows after it are excluded. Pass `NA`
#'   to average every row deliberately.
#' @param year Name of the year column.
#' @param colour_by Column the plot maps to colour, or `NULL` if the mean line
#'   should be drawn in a neutral grey.
#' @return A `geom_hline` layer, or `NULL`.
#' @keywords internal
#' @noRd
.rce_mean_line <- function(df, incl_mean = FALSE, by = NULL, value = "value",
                           hind_endyr, year = "Year", colour_by = NULL) {
  if (!isTRUE(incl_mean)) return(NULL)
  if (missing(hind_endyr)) {
    stop("`hind_endyr` is required: the mean is over hindcast years only.",
         call. = FALSE)
  }
  hind <- df
  if (!is.null(hind_endyr) && !is.na(hind_endyr) && !is.null(hind[[year]])) {
    hind <- hind[!is.na(hind[[year]]) & hind[[year]] <= hind_endyr, ,
                 drop = FALSE]
  }
  hind <- hind[!is.na(hind[[value]]), , drop = FALSE]
  if (nrow(hind) == 0L) return(NULL)

  # A grouping column must survive into the layer's data for faceting to place
  # the line, so drop any the frame does not carry.
  by <- by[by %in% names(hind)]
  by <- setdiff(by, ".mean")
  agg <- if (length(by)) {
    stats::aggregate(hind[[value]], by = hind[by], FUN = mean)
  } else {
    data.frame(x = mean(hind[[value]]))
  }
  names(agg)[ncol(agg)] <- ".mean"

  colour_by <- if (!is.null(colour_by) && colour_by %in% names(agg)) colour_by
  args <- list(data = agg, inherit.aes = FALSE, linetype = 3)
  if (is.null(colour_by)) {
    args$mapping <- ggplot2::aes(yintercept = .data$.mean)
    args$colour  <- "grey40"
  } else {
    args$mapping <- ggplot2::aes(yintercept = .data$.mean,
                                 colour = .data[[colour_by]])
  }
  do.call(ggplot2::geom_hline, args)
}


#' Standard errors for an ADREPORTed series, by name
#'
#' Factored out of `plot_timeseries()` so any plotter drawing a confidence
#' interval reads `sdrep` the same way.
#'
#' The block length must match `n_total` **exactly**. `sdrep$value` holds the
#' whole flattened series, so taking the first `n_need` of a block of unverified
#' length silently returns the standard errors of different cells -- an interval
#' that is wrong rather than absent.
#'
#' Taking a leading slice is legitimate only when the caller knows the block's
#' full shape, which is why `n_total` has to be stated rather than inferred. The
#' species-by-year series are flattened column-major with the hindcast years
#' first, so the first `nspp * nyrs_hindcast` values are exactly the hindcast --
#' that is the one prefix any caller here needs.
#'
#' @param model An `Rceattle` fit.
#' @param name Name of the ADREPORTed series.
#' @param n_need How many leading values the caller wants.
#' @param n_total Length the whole block must have. Defaults to `n_need`, i.e.
#'   no slicing.
#' @return `n_need` standard errors, or `NULL` when the fit does not carry a
#'   block of exactly `n_total` (no `sdreport`, or a `REPORT()`-only series).
#' @keywords internal
#' @noRd
.rce_series_sd <- function(model, name, n_need, n_total = n_need) {
  sdrep <- model$sdrep
  if (is.null(sdrep) || is.null(sdrep$value)) return(NULL)
  rows <- which(names(sdrep$value) == name)
  if (length(rows) != n_total || n_need > n_total) return(NULL)
  sdrep$sd[rows][seq_len(n_need)]
}


#' Can this quantity carry a confidence interval at all?
#'
#' Reads the same `.RCEATTLE_QUANTITIES` registry that `as.data.frame.Rceattle()`
#' uses, so the figure and the table agree on which series have standard errors.
#' A quantity absent from the registry is a derived combination (ration is
#' consumption x biomass; the M2 proportions are ratios) and has none.
#'
#' @param quantity REPORT name of the plotted quantity.
#' @keywords internal
#' @noRd
.rce_quantity_adreport <- function(quantity) {
  isTRUE(.RCEATTLE_QUANTITIES[[quantity]]$adreport)
}


#' Decline an `add_ci` request that the model cannot support
#'
#' Silently ignoring `add_ci = TRUE` is what let these arguments look functional
#' for so long, so say it once instead. Kept quiet for a quantity the model was
#' never going to report -- `F_spp` and the other `REPORT()`-only series would
#' otherwise warn on every fit, which is the gate `plot_timeseries()` already
#' applies.
#'
#' @param add_ci The user's request.
#' @param quantity Name of the plotted quantity, for the message.
#' @param reason Why it cannot be drawn.
#' @param warn Emit the warning? Pass `FALSE` where the quantity is
#'   `REPORT()`-only by design.
#' @return `FALSE`, always -- so the caller can write `add_ci <- .rce_no_ci(...)`.
#' @keywords internal
#' @noRd
.rce_no_ci <- function(add_ci, quantity, reason, warn = TRUE) {
  if (isTRUE(add_ci) && isTRUE(warn)) {
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
#' @param right_adj Ignored. Base-graphics leftover: the figure widened its
#'   right margin to fit the legend. Set margins on the returned ggplot instead.
#' @param top_adj Ignored. Base-graphics leftover; see `right_adj`.
#' @param mod_cex Ignored. Base-graphics leftover: legend text size. Use
#'   `p + ggplot2::theme(legend.text = ggplot2::element_text(size = ...))`.
#' @param legend.pos Ignored. Base-graphics leftover: legend placement. Use
#'   `p + ggplot2::theme(legend.position = "right")`.
#' @param single.plots Ignored. Base-graphics leftover: one device per panel.
#'   The ggplot facets instead.
#' @param theta Ignored. Base-graphics leftover: viewing angle of a `persp`
#'   surface, which is no longer drawn.
#' @param ymax Ignored. Base-graphics leftover: y-axis maximum. Use
#'   `p + ggplot2::coord_cartesian(ylim = c(0, ymax))`.
#' @param cex Ignored. Base-graphics leftover: point expansion.
#' @return A `ggplot` object.
#' @name rceattle-plot-args
#' @keywords internal
NULL
