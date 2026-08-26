# =============================================================================
# Likelihood-profile figures
# -----------------------------------------------------------------------------
# The component profile is a conflict diagnostic, not a precision one. The total
# alone answers "how well is this parameter determined"; the components answer
# "do the data sources agree about where it should be", and a smooth quadratic
# total can hide two of them pulling in opposite directions.
#
# Layout follows r4ss::SSplotProfile(), which is what most readers of a US
# federal assessment document will have seen: change in negative log-likelihood
# on the y axis, each series re-zeroed at its own minimum, small components
# filtered out, and the total drawn heavier than the rest.
# =============================================================================

#' Build the plot frame for one or more profiles
#'
#' Concatenates [profile_components()] over a list of profiles, adding the
#' `Model` column the facets key on, and checks that they profile the same
#' thing -- overlaying a profile of M on a profile of sigmaR would draw two
#' unrelated curves on one axis.
#'
#' @param profiles List of `"Rceattle_profile"` objects.
#' @param labels Model labels, one per profile.
#' @param weighted,relative,minfraction Passed to [profile_components()].
#' @return A data frame with a `Model` column and an `x` column holding the
#'   profiled value.
#' @keywords internal
#' @noRd
.profile_plot_frame <- function(profiles, labels, weighted, relative,
                                minfraction) {
  params <- vapply(profiles, function(p) as.character(p$param), character(1))
  if (length(unique(params)) > 1L) {
    stop("These profiles are over different parameters (",
         paste(unique(params), collapse = ", "),
         "), so they do not share an x axis.", call. = FALSE)
  }
  ncell <- vapply(profiles, function(p) ncol(p$grid), integer(1))
  if (any(ncell != 1L)) {
    stop("plot_profile() draws a 1-D profile; profile ",
         which(ncell != 1L)[1], " is a cross-profile over ",
         ncell[which(ncell != 1L)[1]],
         " cells. Read the surface from `$grid` and `$nll` instead.",
         call. = FALSE)
  }

  # `param` is the internal slot, which is the same for a natural-scale and a
  # log-scale run of the same parameter: profile(param = "M1") and
  # profile(param = "log_M1") both report "log_M1" while their grids hold M and
  # log M. Overlaid, one curve would sit at 0.3 and the other at -1.2 on an axis
  # labelled from the first.
  aliases <- vapply(profiles, function(p) {
    a <- p$alias
    if (is.null(a) || is.na(a)) NA_character_ else as.character(a)
  }, character(1))
  if (length(unique(aliases)) > 1L) {
    stop("These profiles are over the same parameter but not on the same ",
         "scale: ", paste(unique(ifelse(is.na(aliases), params, aliases)),
                          collapse = " and "),
         ". Run them the same way before overlaying them.", call. = FALSE)
  }

  # Different cells of the same parameter -- species 1 against species 2 -- is a
  # legitimate comparison, but only the first one's cell reaches the axis label.
  cells <- vapply(profiles, function(p) paste(p$slots[[1]], collapse = ", "),
                  character(1))
  if (length(unique(cells)) > 1L) {
    warning("These profiles are over different cells of ", params[1], " (",
            paste(unique(cells), collapse = "; "),
            "), so the shared x axis is labelled from the first. Pass `xlab` ",
            "to say what it is.", call. = FALSE)
  }

  frames <- lapply(seq_along(profiles), function(i) {
    d <- profile_components(profiles[[i]], weighted = weighted,
                            relative = relative, minfraction = minfraction)
    # The grid column is named slot_1 by profile(); the plot wants one name
    # whatever the profile.
    d$x <- d[[1]]
    d$Model <- labels[i]
    d[, c("Model", "x", "fit", "component", "unit", "axis", "series", "value")]
  })

  out <- do.call(rbind, frames)
  # rbind() of factors with different levels drops to the union in the wrong
  # order, so the series ordering is rebuilt across all models at once.
  span <- vapply(split(out$value, as.character(out$series)),
                 function(v) suppressWarnings(diff(range(v, na.rm = TRUE))),
                 numeric(1))
  span[!is.finite(span)] <- 0
  ord <- names(sort(span, decreasing = TRUE))
  if ("Total" %in% ord) ord <- c("Total", setdiff(ord, "Total"))
  out$series <- factor(as.character(out$series), levels = ord)
  out$Model  <- factor(out$Model, levels = unique(labels))
  out
}


#' Default x-axis label for a profile
#'
#' Names the parameter in the units the grid is in: the alias when one was
#' used, since `param` has by then been resolved to the log-scale slot the
#' model estimates, and the profiled cell in brackets so a multi-species or
#' multi-sex profile says which one it is.
#'
#' @param prof An `"Rceattle_profile"`.
#' @keywords internal
#' @noRd
.profile_xlab <- function(prof) {
  nm <- prof$alias
  if (is.null(nm) || is.na(nm) || !nzchar(nm)) nm <- prof$param
  slot <- prof$slots[[1]]
  # Under an alias profile() has already appended the rec_pars column, which is
  # part of the alias name rather than something to show.
  if (!is.null(prof$alias) && !is.na(prof$alias) &&
      prof$param == "rec_pars" && length(slot) == 2L) {
    slot <- slot[1]
  }
  paste0(nm, "[", paste(slot, collapse = ", "), "]")
}


#' Plot the likelihood components across a profile
#'
#' @description Draws each negative log-likelihood component against the
#'   profiled parameter, with the total overlaid, in the style of
#'   `r4ss::SSplotProfile()`. Where the curves disagree about which value of the
#'   parameter they prefer, the data sources are in conflict -- which is what
#'   the total on its own cannot show.
#'
#' @details
#' Read the figure by where each curve bottoms out, not by how deep it is. Every
#' series is re-zeroed at its own minimum (`relative = "own"`) and a point marks
#' that minimum, so the spread of the points along the x axis *is* the
#' disagreement: components whose points sit together support the same value,
#' and a component whose point sits far from the total's is pulling against the
#' rest. `relative = "minimum"` instead re-zeroes every series at the total's
#' minimum, which shows what each component gives up by moving away from the
#' fitted value.
#'
#' The total is drawn in black, heavier than the components, and is kept out of
#' the colour legend so the palette separates only the components being
#' compared. Components whose change over the grid is under `minfraction` of the
#' total's are dropped: they are flat on this scale and would only crowd the
#' legend. Non-converged grid points leave a gap in the curve.
#'
#' `line_col`, `lwd` and `lty` apply to the **components**. The total is drawn
#' solid and black at 1.6 times `lwd`, so that it stays the reference line
#' whatever the components are styled as. Series are ordered, in the legend and
#' in the palette, by how much each moves over the grid.
#'
#' Under `random_rec = TRUE` the total is the Laplace-approximated marginal
#' likelihood while the components are the inner joint negative
#' log-likelihood, so they will not sum; [profile_components()] says so when
#' they differ. The shapes are still comparable.
#'
#' @param Rceattle_profile A single `"Rceattle_profile"` from
#'   [profile.Rceattle()], or a list of them to compare models in facets (e.g.
#'   the same profile run on two model configurations).
#' @param weighted,relative,minfraction Passed to [profile_components()].
#'   `minfraction` defaults to `0.01` here, as in `r4ss::SSplotProfile()`.
#' @param add_cutoff Draw a horizontal line at `cutoff`. Off by default: the
#'   cutoff is a statement about the total, and drawing it across the
#'   components invites reading it as one about them.
#' @param cutoff Height of that line. Default `1.92`, the 95%
#'   profile-likelihood cutoff for one parameter, \eqn{\chi^2_1(0.95)/2}.
#' @param xlab X-axis label. Default names the profiled parameter and cell in
#'   the units the grid is in.
#' @param ylab Y-axis label. Default follows `relative`.
#' @inheritParams rceattle-plot-args
#'
#' @return A `ggplot` object. Its `$data` is the plotted frame, so the numbers
#'   behind the figure can be read back off it.
#'
#' @seealso [profile_components()] for the same numbers as a data frame,
#'   [profile.Rceattle()] to run the profile, [print.Rceattle_profile()] for
#'   whether the grid brackets the minimum.
#'
#' @references
#' Taylor, I.G., et al. (2021) Beyond visualizing catch-at-age models: Lessons
#' learned from the r4ss package about software to support stock assessments.
#' Fisheries Research 239: 105924.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE, msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#'
#' # Profile M for species 1 and see which data sources disagree about it
#' prof <- profile(ss_run, param = "M1", slots = list(c(1, 1, 1)),
#'                 values = list(seq(0.2, 0.5, by = 0.05)))
#'
#' plot_profile(prof, xlab = "M at age 1")
#'
#' # Two model configurations side by side
#' plot_profile(list(prof, prof), model_names = c("Base", "Alternative"))
#' }
#' @export
plot_profile <- function(Rceattle_profile,
                         model_names = NULL,
                         weighted = TRUE,
                         relative = c("own", "minimum", "none"),
                         minfraction = 0.01,
                         add_cutoff = FALSE,
                         cutoff = 1.92,
                         xlab = NULL,
                         ylab = NULL,
                         line_col = NULL,
                         lwd = 3,
                         lty = 1,
                         file = NULL,
                         width = 7,
                         height = 5) {

  relative <- match.arg(relative)
  profiles <- if (inherits(Rceattle_profile, "Rceattle_profile")) {
    list(Rceattle_profile)
  } else {
    Rceattle_profile
  }
  if (!length(profiles) ||
      !all(vapply(profiles, inherits, logical(1), "Rceattle_profile"))) {
    stop("`Rceattle_profile` must be an \"Rceattle_profile\" from profile(), ",
         "or a list of them.", call. = FALSE)
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a single number.", call. = FALSE)
  }
  labels <- .model_labels(profiles, model_names)
  df <- .profile_plot_frame(profiles, labels, weighted, relative, minfraction)

  comps <- df[df$series != "Total", , drop = FALSE]
  total <- df[df$series == "Total", , drop = FALSE]
  n_comp <- length(unique(droplevels(comps$series)))
  if (!n_comp) {
    # The figure exists to compare components, so drawing the total alone is a
    # blank answer to the question asked. Say which argument emptied it.
    warning("Every likelihood component moves less than `minfraction` = ",
            minfraction, " of the total over this grid, so only the total is ",
            "drawn. Lower `minfraction` to see them.", call. = FALSE)
  }
  .rce_check_line_col(line_col, n_comp, "likelihood components")

  # The point marks each series' minimum: with relative = "own" that is where
  # the component sits at zero, and the spread of those points along x is the
  # disagreement the figure exists to show.
  at_min <- function(d) {
    marks <- do.call(rbind, lapply(split(d, list(d$Model, d$series), drop = TRUE),
                                   function(g) {
                                     if (all(is.na(g$value))) return(NULL)
                                     g[which.min(g$value), , drop = FALSE]
                                   }))
    # An empty layer must still be a data frame: geom_point(data = NULL) means
    # "inherit the plot's data", so returning NULL would draw a point at every
    # grid point instead of none.
    if (is.null(marks)) d[0, , drop = FALSE] else marks
  }

  if (is.null(xlab)) xlab <- .profile_xlab(profiles[[1]])
  if (is.null(ylab)) {
    ylab <- if (relative == "none") {
      "Negative log-likelihood"
    } else {
      "Change in negative log-likelihood"
    }
  }

  lp <- .rce_line_params(lwd = lwd, lty = lty, lty_by = "series",
                         lty_n = n_comp)

  p <- ggplot2::ggplot(comps, ggplot2::aes(x = .data$x, y = .data$value,
                                           colour = .data$series))
  if (relative != "none") {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey70",
                                 linewidth = 0.4)
  }
  if (add_cutoff) {
    p <- p + ggplot2::geom_hline(yintercept = cutoff, colour = "grey70",
                                 linetype = 2, linewidth = 0.4)
  }
  p <- .rce_add_line(p, lp, na.rm = TRUE)
  p <- p + ggplot2::geom_point(data = at_min(comps), size = 2, na.rm = TRUE)

  # The total is its own black layer rather than another level of `series`: it
  # is the reference the components are read against, not one of the things
  # being compared, and leaving it out of the palette keeps the component
  # colours stable as `minfraction` changes what is shown.
  p <- p +
    ggplot2::geom_line(data = total, colour = "black",
                       linewidth = lwd[1] / 3 * 1.6,
                       inherit.aes = FALSE, na.rm = TRUE,
                       mapping = ggplot2::aes(x = .data$x, y = .data$value)) +
    ggplot2::geom_point(data = at_min(total), colour = "black", size = 2.5,
                        inherit.aes = FALSE, na.rm = TRUE,
                        mapping = ggplot2::aes(x = .data$x, y = .data$value)) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL) +
    .rceattle_theme()

  # A fleet-by-component legend runs to long labels ("EIT_Pollock: Composition
  # data") and there are usually more than a handful, so the single bottom row
  # ggplot2 would draw runs off the page. Wrap to as many columns as fit a
  # ~90-character row.
  lab_width <- max(nchar(levels(droplevels(comps$series))), 0L)
  p <- p + ggplot2::guides(colour = ggplot2::guide_legend(
    ncol = max(1L, min(4L, floor(90 / (lab_width + 6)))), byrow = TRUE))

  p <- .rceattle_scale(p, discrete = TRUE, aesthetics = "colour",
                       line_col = line_col)

  if (nlevels(df$Model) > 1L) p <- p + ggplot2::facet_wrap(~ Model)

  .save_ggplot(p, file, "profile", width = width, height = height)
}


#' Plot method for a likelihood profile
#'
#' @description Shorthand for [plot_profile()]: draws the likelihood components
#'   against the profiled parameter, with the total overlaid.
#'
#' @param x An `"Rceattle_profile"` object from [profile.Rceattle()].
#' @param ... Passed to [plot_profile()].
#' @return A `ggplot` object.
#' @method plot Rceattle_profile
#' @export
plot.Rceattle_profile <- function(x, ...) {
  plot_profile(x, ...)
}
