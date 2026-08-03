#' Plot fishery selectivity and maturity
#'
#' @description Overlays terminal-year fishery selectivity-at-age against input
#'   maturity-at-age, faceted by species. Useful for debugging SPR-based
#'   reference points (the two curves drive spawning-potential ratio).
#'
#' @param Rceattle A single [fit_mod()] object or a list of them (the first is
#'   used).
#' @param file Optional file stem; the figure is written to
#'   `<file>_selectivity_vs_maturity.png` if given.
#' @param model_names,line_col,lwd Deprecated base-graphics arguments, retained
#'   for back-compatibility and ignored.
#' @param species Species names for the facet labels. Default: model species
#'   names.
#' @param width,height Saved figure size (inches).
#'
#' @return A `ggplot` object.
#' @export
plot_selectivity_vs_maturity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           width = 7,
           height = 6.5,
           species = NULL,
           lwd = 3) {

    Rceattle <- .as_model_list(Rceattle)
    fc      <- Rceattle[[1]]$data_list$fleet_control
    nages   <- Rceattle[[1]]$data_list$nages
    minage  <- Rceattle[[1]]$data_list$minage
    nspp    <- Rceattle[[1]]$data_list$nspp
    spnames <- if (is.null(species)) Rceattle[[1]]$data_list$spnames else species
    nyrshind <- Rceattle[[1]]$data_list$endyr - Rceattle[[1]]$data_list$styr + 1L

    sel <- Rceattle[[1]]$quantities$sel_at_age           # [fleet, sex, age, yr]
    mat <- Rceattle[[1]]$data_list$maturity[, -1, drop = FALSE]  # [species, age]

    df_list <- list()
    for (sp in seq_len(nspp)) {
      ages <- seq_len(nages[sp]) - 1 + minage[sp]
      # terminal-year selectivity for each fishery fleet of this species
      fleets <- fc$Fleet_code[fc$Species == sp & fc$Fleet_type == "Fishery"]
      for (flt in fleets) {
        df_list[[length(df_list) + 1L]] <- data.frame(
          Species = spnames[sp],
          Age     = ages,
          Series  = as.character(fc$Fleet_name[fc$Fleet_code == flt]),
          value   = as.numeric(sel[flt, 1, seq_len(nages[sp]), nyrshind]),
          stringsAsFactors = FALSE)
      }
      # input maturity-at-age
      df_list[[length(df_list) + 1L]] <- data.frame(
        Species = spnames[sp],
        Age     = ages,
        Series  = "Maturity",
        value   = as.numeric(mat[sp, seq_len(nages[sp])]),
        stringsAsFactors = FALSE)
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Species <- factor(plot_df$Species, levels = spnames)

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$Age, y = .data$value,
                                      colour = .data$Series)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~ Species) +
      ggplot2::labs(x = "Age", y = "Selectivity / Maturity")
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")

    .save_ggplot(p, file = file, suffix = "selectivity_vs_maturity",
                 width = width, height = height)
  }
