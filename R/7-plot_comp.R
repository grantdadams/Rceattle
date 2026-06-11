#' Plot composition fits and residuals
#'
#' @description
#' Diagnostic plots for age / length composition data from a fitted [Rceattle]
#' model, drawn with ggplot2 for a consistent look with [plot.rceattle_osa()].
#' With the default `residual_type = "pearson"` three figures are produced per
#' fleet x composition type:
#' \itemize{
#'   \item **Pearson residual bubbles** by year and bin, faceted by fleet (and,
#'     for joint-sex data, by sex); red = positive, blue = negative, sized by
#'     magnitude. The Pearson residual is
#'     \eqn{(p - \hat p)/\sqrt{\hat p (1 - \hat p)/N}}, the same form used by
#'     [residuals.Rceattle()].
#'   \item **Annual composition** -- observed (shaded area) vs fitted (line)
#'     proportion at age / length, one panel per year. Joint-sex data are
#'     mirrored (females up, males down).
#'   \item **Aggregated composition** -- the same, summed over (hindcast) years.
#' }
#' The shaded area and fitted line span only the observed bins (they do not
#' extend past the first/last bin), and bins with zero observed proportion are
#' retained (only `NA` bins are dropped), so the curves are not interpolated
#' across gaps.
#'
#' @param Rceattle A single fitted `Rceattle` model.
#' @param file Optional filename stem; when supplied, each figure is also saved
#'   as a PNG.
#' @param model_names Unused (kept for back-compatibility).
#' @param species Optional species code(s) to plot (matched against the
#'   `Species` column); `NULL` (default) plots all. Mirrors the `species`
#'   argument of [residuals.Rceattle()] and [plot.rceattle_osa()].
#' @param cex Unused (kept for back-compatibility).
#' @param lwd Width of the fitted-composition line. Default `3`.
#' @param right_adj Unused (kept for back-compatibility).
#' @param residual_type `"pearson"` (default) for the ggplot2 Pearson-residual
#'   and composition-fit figures drawn here, or `"osa"` to instead draw the
#'   one-step-ahead residual diagnostics via [osa_residuals()] and
#'   [plot.rceattle_osa()] -- a Q-Q plot (with SDNR / tail annotation) alongside
#'   signed OSA- and Pearson-residual bubbles. The `"osa"` path requires a model
#'   fit with `fit_control(osa = TRUE)`.
#'
#' @return Invisibly, a named list of the `ggplot` objects. Called for its side
#'   effect of drawing (and optionally saving) the figures.
#' @importFrom rlang .data
#' @export
plot_comp <- function(Rceattle, file = NULL, model_names = NULL, species = NULL,
                      cex = 3, lwd = 3, right_adj = 0,
                      residual_type = c("pearson", "osa")) {

  if (!inherits(Rceattle, "Rceattle")) {
    stop("Please only use one Rceattle model")
  }
  residual_type <- match.arg(residual_type)

  # One-step-ahead residual diagnostics: delegate to the rceattle_osa plot.
  if (residual_type == "osa") {
    osa <- osa_residuals(Rceattle, source = "comp")
    return(invisible(plot(osa, species = species)))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot composition data.")
  }

  long <- .comp_resid_long(Rceattle, species = species)
  if (is.null(long) || nrow(long) == 0) {
    warning("No composition data to plot for the requested species.")
    return(invisible(NULL))
  }

  sex_cols  <- c(combined = "black", female = "#d7301f", male = "#2c7fb8")
  sign_cols <- c(positive = "#d7301f", negative = "#2c7fb8")
  save_png  <- function(p, tag, w = 9, h = 6.5) {
    if (!is.null(file)) {
      ggplot2::ggsave(paste0(file, "_comp_", tag, ".png"), p,
                      width = w, height = h, units = "in", dpi = 300)
    }
  }
  plots <- list()

  # ---- Pearson residual bubbles (faceted by fleet x type [x sex]) ----
  pear <- long[is.finite(long$pearson), , drop = FALSE]
  if (nrow(pear) > 0) {
    pear$sign <- ifelse(pear$pearson >= 0, "positive", "negative")
    g <- ggplot2::ggplot(pear, ggplot2::aes(.data$Year, .data$bin)) +
      ggplot2::geom_point(ggplot2::aes(size = abs(.data$pearson),
                                       colour = .data$sign), alpha = 0.8) +
      ggplot2::scale_colour_manual(values = sign_cols, guide = "none") +
      ggplot2::scale_size_continuous(range = c(0.5, 6), name = "|Pearson|") +
      ggplot2::facet_wrap(~ source, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Age / length bin",
                    title = "Composition Pearson residuals") +
      ggplot2::theme_bw(base_size = 10)
    plots$pearson <- g; print(g); save_png(g, "pearson")
  }

  # ---- Composition fit: observed (area) vs fitted (line) ----
  # Mirror joint-sex (females up, males down). geom_area spans only the observed
  # bins, so it cannot extend past the data.
  long$y_obs <- ifelse(long$sex_grp == "male", -long$obs, long$obs)
  long$y_hat <- ifelse(long$sex_grp == "male", -long$hat, long$hat)

  # Observed proportion as a shaded area, fitted as a line; joint-sex males are
  # mirrored below zero (y_obs / y_hat). geom_area spans only the observed bins.
  comp_area_plot <- function(df, ylab, title) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$bin)) +
      ggplot2::geom_area(ggplot2::aes(y = .data$y_obs, group = .data$sex_grp),
                         position = "identity", fill = "grey80") +
      ggplot2::geom_line(ggplot2::aes(y = .data$y_hat, group = .data$sex_grp,
                                      colour = .data$sex_grp),
                         linewidth = lwd / 3) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
      ggplot2::scale_colour_manual(values = sex_cols, guide = "none") +
      ggplot2::labs(x = df$bin_lab[1], y = ylab, title = title) +
      ggplot2::theme_bw(base_size = 9)
  }

  for (nm in unique(long$panel)) {
    d <- long[long$panel == nm, , drop = FALSE]

    # Annual (one panel per year)
    g_a <- comp_area_plot(d, "Proportion", nm) + ggplot2::facet_wrap(~ Year)
    plots[[paste0("annual_", nm)]] <- g_a; print(g_a)
    save_png(g_a, paste0("annual_fleet", d$Fleet[1], "_", d$type_lab[1]))

    # Aggregated across hindcast years
    agg <- .comp_aggregate(d)
    if (nrow(agg) > 0) {
      g_g <- comp_area_plot(agg, "Mean proportion", paste(nm, "(aggregated)"))
      plots[[paste0("aggregated_", nm)]] <- g_g; print(g_g)
      save_png(g_g, paste0("aggregated_fleet", d$Fleet[1], "_", d$type_lab[1]))
    }
  }

  invisible(plots)
}


#' Tidy long-format composition observed / fitted / Pearson table
#'
#' Reuses [residuals.Rceattle()] (`type = "pearson"`, `source = "comp"`) for the
#' observed and fitted proportions and the Pearson residual -- the single source
#' of truth -- then adds the plotting-only columns: joint-sex bins re-based onto
#' a single age/length axis (males stored as bins `nbin+1 .. 2*nbin` are mapped
#' to `1 .. nbin` and tagged `sex_grp = "male"`) and the facet labels. Zero
#' observed proportions are kept; only `NA` observed/fitted are dropped (by
#' `residuals()`).
#' @param Rceattle A fitted `Rceattle` model.
#' @param species Optional species code(s) to keep.
#' @return A data frame, or `NULL` when there are no composition residuals.
#' @keywords internal
.comp_resid_long <- function(Rceattle, species = NULL) {
  r <- stats::residuals(Rceattle, type = "pearson", source = "comp",
                        species = species)
  if (is.null(r) || nrow(r) == 0) return(NULL)
  nages    <- Rceattle$data_list$nages
  nlengths <- Rceattle$data_list$nlengths

  long <- data.frame(
    Fleet      = r$Fleet_code,
    Fleet_name = r$Fleet_name,
    Species    = r$Species,
    Sex        = r$Sex,
    comp_type  = r$Age0_Length1,
    Year       = r$Year,
    N          = r$Sample_size,
    bin        = r$Bin,
    obs        = r$Observed,      # observed proportion
    hat        = r$Fitted,        # fitted proportion
    pearson    = r$Residual,      # (obs - hat) / sqrt(hat (1 - hat) / N)
    stringsAsFactors = FALSE)

  # Joint-sex (Sex == 3) stacks females in bins 1..nbin and males in
  # nbin+1..2*nbin; re-base the male bins onto the same 1..nbin axis and label
  # each row's sex so the two can be drawn together (males mirrored below 0).
  bins_per_sex <- ifelse(long$comp_type == 1, nlengths[long$Species],
                         nages[long$Species])
  is_male <- long$Sex == 3 & long$bin > bins_per_sex
  long$bin[is_male] <- long$bin[is_male] - bins_per_sex[is_male]
  long$sex_grp <- ifelse(long$Sex == 1, "female",
                  ifelse(long$Sex == 2, "male",
                  ifelse(long$Sex == 3, ifelse(is_male, "male", "female"),
                         "combined")))

  long$type_lab <- ifelse(long$comp_type == 1, "length", "age")
  long$bin_lab  <- ifelse(long$comp_type == 1, "Length bin", "Age bin")
  long$panel    <- paste0(long$Fleet_name, " - ", long$type_lab, " comp")
  long$source   <- paste0(long$Fleet_name, "\n", long$type_lab, " comp",
                          ifelse(long$Sex == 3, paste0(" - ", long$sex_grp), ""))
  long
}


#' Aggregate composition proportions across hindcast years
#'
#' Sums observed and fitted proportions over years (for one fleet x type panel)
#' and renormalizes to the mean composition; joint-sex groups keep their shared
#' normalization (females + males sum to 1).
#' @param d One panel's rows from [.comp_resid_long()].
#' @return A data frame with `bin`, `sex_grp`, `obs`, `hat`, `y_obs`, `y_hat`.
#' @keywords internal
.comp_aggregate <- function(d) {
  d <- d[d$Year > 0, , drop = FALSE]              # hindcast only
  if (nrow(d) == 0) return(d[, c("bin", "sex_grp"), drop = FALSE])
  agg <- stats::aggregate(cbind(obs, hat) ~ bin + sex_grp, data = d, FUN = sum)
  agg$obs <- agg$obs / sum(agg$obs)               # joint normalization (mean prop)
  agg$hat <- agg$hat / sum(agg$hat)
  agg$bin_lab <- d$bin_lab[1]
  agg$y_obs <- ifelse(agg$sex_grp == "male", -agg$obs, agg$obs)
  agg$y_hat <- ifelse(agg$sex_grp == "male", -agg$hat, agg$hat)
  agg
}




#' Plot diet composition fits
#'
#' @description
#' If year == 0, diet data are averaged from suit_styr to suit_endyr
#' If prey_age >= 0 diet data are diet proportion of prey-at-age in predator-at-age
#' If prey_age < 0 diet data are diet proportion of prey-spp in predator-at-age (sum across prey ages)
#' If prey_age < 0 and pred_age < 0, diet data are mean diet proportion of prey-spp in predator-spp (sum across prey ages and take mean across predator ages)
#' If prey_age < 0 and pred_age < -500, diet data are weighted mean diet proportion of prey-spp in predator-spp (sum across prey ages and take weighted mean across predator ages)
#'
#'
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param file name of a file to identified the files exported by the function.
#' @param species Species names for legend
#'
#' @return Returns and saves a figure
#' @importFrom rlang .data
#' @importFrom ggplot2 theme element_blank aes ggplot geom_point theme_classic ylim ylab xlim xlab ggtitle ggsave
#' @export
plot_diet_comp <-
  function(Rceattle,
           file = NULL,
           species = NULL) {

    # Make sure we are using only one model
    if(!inherits(Rceattle, "Rceattle")){
      stop("Please only use one Rceattle model")
    }
    data_list <- Rceattle$data_list

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }

    # Get colors
    colvec=c("red", "blue", "black")

    # * Extract data objects ----
    # Diet observed proportions, fitted, and the Pearson residual all come from
    # residuals(type = "pearson", source = "diet") -- the single source of truth
    # for residual computation -- then aliased to the column names this plotting
    # loop expects.
    comp_data <- stats::residuals(Rceattle, type = "pearson", source = "diet")
    comp_data$Pred <- comp_data$Species
    comp_data$Stomach_proportion_by_weight <- comp_data$Observed
    comp_data$Est <- comp_data$Fitted
    comp_data$pearson <- comp_data$Residual

    # Loop around predators ----
    for(pred in 1:data_list$nspp) {
      for(pred_sex_idx in 1:data_list$nsex[pred]){
        for(prey in 1:data_list$nspp) {

          # * Get sex for legend ----
          if(data_list$nsex[pred] > 1){
            # Fixed: sex variable was undefined in your original snippet
            pred_legend <- paste("Pred-", species[pred], ifelse(pred_sex_idx == 1, "female", "male"))
            current_pred_sex = pred_sex_idx - 1
          } else{
            pred_legend <- paste("Pred-", species[pred])
            current_pred_sex = 0
          }

          # * Extract comps ----
          comp_tmp <- comp_data |>
            dplyr::filter(.data$Pred == pred & .data$Pred_sex == current_pred_sex & .data$Prey == prey)

          yrs <- sort(unique(comp_tmp$Year))
          nyrs <- length(yrs)

          # * Plot annual comps ----
          for(yr in 1:nyrs){
            comp_tmp_yr <- comp_tmp |>
              dplyr::filter(.data$Year == yrs[yr] ) |>
              dplyr::mutate(Prey_age = ifelse(.data$Prey_sex == 2, -.data$Prey_age, .data$Prey_age))

            plot_obs <- comp_tmp_yr |>
              dplyr::filter(.data$Stomach_proportion_by_weight > 0) |>
              ggplot2::ggplot(ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$Stomach_proportion_by_weight)) +
              ggplot2::geom_point(alpha = 1) +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Observed diet: year", yrs[yr])) +
              ggplot2::theme(legend.position = c(0.25, 0.7),
                             legend.title = ggplot2::element_blank())

            plot_est <- comp_tmp_yr |>
              dplyr::filter(.data$Stomach_proportion_by_weight > 0) |>
              ggplot2::ggplot(ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$Est)) +
              ggplot2::geom_point(alpha = 1) +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Estimated diet: year", yrs[yr])) +
              ggplot2::theme(legend.position = "none")

            plot_pear <- comp_tmp_yr |>
              dplyr::filter(.data$Stomach_proportion_by_weight > 0) |>
              ggplot2::ggplot(ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age,
                                           size = abs(.data$pearson),
                                           color = .data$pearson < 0)) +
              ggplot2::geom_point(alpha = 1) +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Pearson residual: year", yrs[yr])) +
              ggplot2::theme(legend.position = c(0.25, 0.8))

            p1 <- cowplot::plot_grid(plot_obs, plot_est, plot_pear, nrow = 1)
            print(p1)

            if (!is.null(file)) {
              filename <- paste0(file, "_aggregated_diet_comps_year", yr,"_", pred_legend, "_prey_", species[prey],".png")
              ggplot2::ggsave(filename = filename, plot = p1, width = 10, height = 6.5, units = "in", dpi = 300)
            }
          }
        }
      }
    }
  }
