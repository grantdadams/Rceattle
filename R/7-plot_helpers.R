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
#' Clean, legible, and colourblind-friendly: `theme_classic` with bold facet
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

#' Okabe-Ito colourblind-safe qualitative palette
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

#' Add the standard colourblind-safe scales to a plot
#'
#' Discrete aesthetics (series identity, e.g. model) use the Okabe-Ito
#' qualitative palette; continuous aesthetics (ordered magnitude, e.g. year)
#' use viridis. Both are colourblind-safe.
#'
#' @param p A ggplot.
#' @param discrete Use the discrete (Okabe-Ito) or continuous (viridis) scale.
#' @param aesthetics Which of `"colour"`/`"fill"` to add.
#' @keywords internal
.rceattle_scale <- function(p, discrete = TRUE,
                            aesthetics = c("colour", "fill")) {
  if ("colour" %in% aesthetics) {
    p <- p + (if (discrete)
                ggplot2::discrete_scale("colour", palette = .okabe_ito_pal)
              else ggplot2::scale_colour_viridis_c())
  }
  if ("fill" %in% aesthetics) {
    p <- p + (if (discrete)
                ggplot2::discrete_scale("fill", palette = .okabe_ito_pal)
              else ggplot2::scale_fill_viridis_c())
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
