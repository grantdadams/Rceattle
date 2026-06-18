#' Plot stock recruit function
#'
#' @description Plots the stock-recruit relationship estimated by Rceattle:
#'   spawning stock biomass (x) against recruitment (y) as points, with the
#'   fitted stock-recruit curve (mean recruitment, Beverton-Holt, or Ricker)
#'   overlaid, faceted by species. A 95% normal data ellipse of the SSB-R
#'   cloud is added when `add_ci = TRUE`.
#'
#' @param Rceattle A single [fit_mod()] object or a list of them (overlaid).
#' @param file Optional file stem; the figure is written to
#'   `<file>_stock_recruit.png` if given.
#' @param model_names Legend labels for the models.
#' @param species Species (indices) to include. Default all.
#' @param spnames Species names for the facet labels.
#' @param incl_proj Currently unused (kept for back-compatibility).
#' @param width,height Saved figure size (inches).
#' @param add_ci Add a 95% normal data ellipse of the SSB-R points.
#' @param line_col,lwd,lty,plot_env,mod_cex Deprecated base-graphics
#'   arguments, retained for back-compatibility and ignored.
#'
#' @return A `ggplot` object.
#' @export
plot_stock_recruit <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           width = 7,
           height = 6.5,
           species = NULL,
           spnames = NULL,
           lwd = 3,
           lty = 1,
           incl_proj = FALSE,
           plot_env = FALSE,
           mod_cex = 1,
           add_ci = TRUE) {

    models <- .as_model_list(Rceattle)
    model_names_use <- .model_labels(models, model_names)
    nspp <- models[[1]]$data_list$nspp
    if (is.null(spnames)) spnames <- models[[1]]$data_list$spnames
    if (is.null(species)) species <- 1:nspp

    pts <- list()
    crv <- list()
    for (mod in seq_along(models)) {
      dl       <- models[[mod]]$data_list
      nyrshind <- dl$endyr - dl$styr + 1L
      srr_pred <- dl$srr_pred_fun
      # SSB (million mt) and recruitment (millions) over the hindcast.
      ssb <- models[[mod]]$quantities$ssb[, seq_len(nyrshind), drop = FALSE] / 1e6
      R   <- models[[mod]]$quantities$R[,   seq_len(nyrshind), drop = FALSE] / 1e6
      rp  <- models[[mod]]$estimated_params$rec_pars
      for (sp in species) {
        pts[[length(pts) + 1L]] <- data.frame(
          Model = model_names_use[mod], Species = spnames[sp],
          SSB = as.numeric(ssb[sp, ]), R = as.numeric(R[sp, ]),
          stringsAsFactors = FALSE)

        # Fitted stock-recruit curve over the observed SSB range. Same forms
        # (and parameter transforms) as the original base-graphics version.
        xmax <- max(ssb[sp, ], na.rm = TRUE)
        xg <- seq(0, xmax, length.out = 100)
        yg <- if (srr_pred == 0) {
          rep(exp(rp[sp, 1]), length(xg))                                # mean rec
        } else if (srr_pred %in% c(2, 3)) {                             # Beverton-Holt
          exp(rp[sp, 2]) * xg / (1 + exp(rp[sp, 3]) * xg * 1e6)
        } else if (srr_pred %in% c(4, 5)) {                             # Ricker
          exp(rp[sp, 2]) * xg * exp(-exp(rp[sp, 3]) * xg)
        } else {
          rep(NA_real_, length(xg))
        }
        crv[[length(crv) + 1L]] <- data.frame(
          Model = model_names_use[mod], Species = spnames[sp],
          SSB = xg, R = yg, stringsAsFactors = FALSE)
      }
    }
    pts_df <- do.call(rbind, pts)
    crv_df <- do.call(rbind, crv)
    pts_df$Species <- factor(pts_df$Species, levels = spnames[species])
    crv_df$Species <- factor(crv_df$Species, levels = spnames[species])

    p <- ggplot2::ggplot(
      pts_df, ggplot2::aes(x = .data$SSB, y = .data$R, colour = .data$Model)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_line(data = crv_df, linewidth = 1) +
      ggplot2::facet_wrap(~ Species, scales = "free") +
      ggplot2::labs(x = "Spawning stock biomass (million mt)",
                    y = "Recruitment (millions)")
    if (add_ci) {
      p <- p + ggplot2::stat_ellipse(type = "norm", linetype = 2, na.rm = TRUE)
    }
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (nlevels(pts_df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")

    .save_ggplot(p, file = file, suffix = "stock_recruit",
                 width = width, height = height)
  }


#' Make a transparent colour
#'
#' From <https://www.dataanalytics.org.uk/make-transparent-colors-in-r/>.
#'
#' @param color color name
#' @param percent percent transparency
#' @param name an optional name for the color
#' @keywords internal
t_col <- function(color, percent = 50, name = NULL) {

  ## Get RGB values for named color
  rgb.val <- grDevices::col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- grDevices::rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                          maxColorValue = 255,
                          alpha = (100 - percent) * 255 / 100,
                          names = name)

  ## Save the color
  invisible(t.col)
}
