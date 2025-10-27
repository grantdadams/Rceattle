#' @title match_diet_preds (Corrected for Matrix)
#' @description Binds the estimated diet proportions from quantities$diet_hat
#' to the observed diet data from data_list$diet_data.
match_diet_preds <- function(data_list, quantities) {

  if (is.null(data_list$diet_data) || nrow(data_list$diet_data) == 0) {
    message("No diet data found in data_list$diet_data")
    return(NULL)
  }

  if (is.null(quantities$diet_hat)) {
    message("No estimated diet found in quantities$diet_hat")
    return(NULL)
  }

  plot_data <- data_list$diet_data
  n_rows_data <- nrow(plot_data)

  # --- NEW: Check if diet_hat is a matrix ---
  est_diet_values <- NULL

  if (is.matrix(quantities$diet_hat)) {
    if (nrow(quantities$diet_hat) == n_rows_data && ncol(quantities$diet_hat) == 2) {
      # This is the 170x2 matrix case, select column 2
      est_diet_values <- quantities$diet_hat[, 2]
    } else {
      # Matrix dimensions don't match
      message(paste(
        "Warning: 'diet_hat' is a matrix but dimensions [",
        nrow(quantities$diet_hat), ",", ncol(quantities$diet_hat),
        "] do not match data rows (", n_rows_data, ")."
      ))
      return(NULL)
    }
  } else if (is.vector(quantities$diet_hat)) {
    # This is the vector case (like before)
    if (length(quantities$diet_hat) == n_rows_data) {
      est_diet_values <- quantities$diet_hat
    }
  }

  # Check for mismatches
  if (is.null(est_diet_values) || length(est_diet_values) != n_rows_data) {
    message(paste(
      "Warning: Mismatch between diet data rows (", n_rows_data,
      ") and processed diet estimates (", length(est_diet_values), ")."
    ))
    return(NULL)
  }
  # --- END NEW BLOCK ---

  # Bind the estimated diet as a new column
  plot_data$Est_diet <- est_diet_values

  return(plot_data)
}


plot_diet_comp2 <- function(Rceattle, file = NULL, species = NULL, add_ci = TRUE) {

  # 1. SETUP & DATA PREPARATION ----
  if(!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for this function.")
  if(!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required for this function.")
  if(!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required for this function.")
  if(!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot is required for this function.")
  if (!inherits(Rceattle, "Rceattle")) stop("Input 'Rceattle' must be a single object of class Rceattle.")
  if (is.null(species)) species <- Rceattle$data_list$spnames

  # Use the helper function
  plot_data <- match_diet_preds(
    data_list = Rceattle$data_list,
    quantities = Rceattle$quantities
  )
  if (is.null(plot_data) || nrow(plot_data) == 0) { message("No diet data to plot."); return(invisible(NULL)) }

  # Get row count for SD matching
  n_rows_data <- nrow(plot_data)

  # <--- NEW: ADD STANDARD DEVIATION & CIs ---
  if (add_ci) {
    if (is.null(Rceattle$sdrep)) {
      message("add_ci = TRUE but Rceattle$sdrep is NULL. Skipping CIs.")
      add_ci <- FALSE
    } else {
      sd_indices <- which(names(Rceattle$sdrep$value) == "diet_hat")

      if (length(sd_indices) > 0) {
        all_sd_values <- Rceattle$sdrep$sd[sd_indices]

        # --- NEW: Check length and select correct SDs ---
        sd_values_to_use <- NULL

        if (length(all_sd_values) == n_rows_data) {
          # Case 1: 170 SDs match 170 rows (diet_hat was a vector)
          sd_values_to_use <- all_sd_values
        } else if (length(all_sd_values) == n_rows_data * 2) {
          # Case 2: 340 SDs for 170 rows (diet_hat was 170x2 matrix)
          # We want the second half, corresponding to the second column
          sd_values_to_use <- all_sd_values[(n_rows_data + 1):(n_rows_data * 2)]
          message("Using second half of 'diet_hat' standard deviations (340 -> 170).")
        }
        # --- END NEW BLOCK ---

        if (!is.null(sd_values_to_use)) {
          plot_data$Est_sd <- sd_values_to_use

          plot_data <- plot_data %>%
            dplyr::mutate(
              Est_Lower = pmax(0, Est_diet - 1.96 * Est_sd),
              Est_Upper = pmin(1, Est_diet + 1.96 * Est_sd)
            )
        } else {
          message(paste(
            "Warning: Length mismatch. SD values (", length(all_sd_values),
            ") do not match diet rows (", n_rows_data, "). Skipping CIs."
          ))
          add_ci <- FALSE
        }
      } else {
        message("Warning: Could not find 'diet_hat' in Rceattle$sdrep$value. Skipping CIs.")
        add_ci <- FALSE
      }
    }
  }

  # Rename columns using names from diet_data
  plot_data <- plot_data %>%
    dplyr::rename(Observed = Stomach_proportion_by_weight, Est = Est_diet) %>%
    dplyr::mutate(
      Est_clipped = pmin(0.9999, pmax(0.0001, Est)),
      Pearson = (Observed - Est) / sqrt((Est_clipped * (1 - Est_clipped)) / Sample_size),
      AbsPearson = abs(Pearson)
    )

  # 2. PLOTTING LOGIC ----
  plot_list <- list()

  for (pred_ind in 1:Rceattle$data_list$nspp) {
    for (prey_ind in 1:Rceattle$data_list$nspp) {

      subset_data <- plot_data %>% dplyr::filter(Pred == pred_ind, Prey == prey_ind)
      if(nrow(subset_data) == 0) next

      is_prey_age_agg <- any(subset_data$Prey_age < 0)
      is_pred_age_agg <- any(subset_data$Pred_age < 0)

      pred_legend <- paste("Predator:", species[pred_ind])

      # --- PATHWAY FOR PREY-AGE AGGREGATED (Line plot) ---
      if(is_prey_age_agg && !is_pred_age_agg) {

        message(paste("Generating line plot for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))

        p <- ggplot2::ggplot(subset_data, ggplot2::aes(x = Pred_age)) +

          {
            if(add_ci && "Est_Lower" %in% names(subset_data)) {
              ggplot2::geom_ribbon(
                ggplot2::aes(ymin = Est_Lower, ymax = Est_Upper, fill = "Est"),
                alpha = 0.2
              )
            }
          } +

          ggplot2::geom_line(ggplot2::aes(y = Est, color = "Est", linetype = "Est"), linewidth = 1) +
          ggplot2::geom_point(ggplot2::aes(y = Est, color = "Est"), size = 2.5, alpha = 0.7) +

          ggplot2::geom_line(ggplot2::aes(y = Observed, color = "Observed", linetype = "Observed"), linewidth = 1) +
          ggplot2::geom_point(ggplot2::aes(y = Observed, color = "Observed"), size = 2.5) +

          ggplot2::facet_wrap(~ Year, scales = "free_y", labeller = ggplot2::labeller(Year = ~paste("Year:", .))) +
          ggplot2::facet_grid(~ Pred_sex, scales = "free_y", labeller = ggplot2::labeller(Pred_sex = ~paste("Sex:", .))) +
          ggplot2::scale_color_manual(name = "Source", values = c("Observed" = "black", "Est" = "darkred")) +
          ggplot2::scale_linetype_manual(name = "Source", values = c("Observed" = "dashed", "Est" = "solid")) +
          ggplot2::scale_fill_manual(name = "Source", values = c("Est" = "darkred"), guide = "none") +
          ggplot2::labs(x = "Predator Age", y = "Diet Proportion", title = paste("Diet of", species[pred_ind], "on", species[prey_ind])) +
          ggplot2::theme_bw()

        print(p)
        plot_list[[length(plot_list) + 1]] <- p

        # --- (Other plot pathways remain unchanged) ---
      } else if (is_prey_age_agg && is_pred_age_agg) {

        message(paste("Generating bar plot for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))

        plot_data_long <- subset_data %>%
          tidyr::pivot_longer(cols = c(Observed, Est), names_to = "Source", values_to = "Proportion") %>%
          dplyr::mutate(Prey_name = species[Prey])

        p_fit <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = factor(Year), y = Proportion, fill = Source)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge") +
          ggplot2::scale_fill_manual(name = "Source", values = c("Observed" = "grey50", "Est" = "red")) +
          ggplot2::labs(x = "Year", y = "Diet Proportion", title = paste("Fit to Aggregated Diet:", species[pred_ind], "on", species[prey_ind])) +
          ggplot2::theme_bw()

        print(p_fit)
        plot_list[[length(plot_list) + 1]] <- p_fit

      } else {

        message(paste("Generating bubble plots for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))
        yrs <- sort(unique(subset_data$Year))

        for(i in 1:length(yrs)) {
          current_yr <- yrs[i]
          comp_tmp_yr <- subset_data %>% dplyr::filter(Year == current_yr)
          if(sum(comp_tmp_yr$Observed, na.rm = TRUE) == 0) next

          title <- paste(pred_legend, "- Prey:", species[prey_ind], "- Year:", current_yr)
          if(current_yr == 0) title <- paste(pred_legend, "- Prey:", species[prey_ind], "(Avg over Years)")

          p_obs <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = Pred_age, y = Prey_age, size = Observed)) +
            ggplot2::geom_point(alpha=0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Observed", size = "Prop.")

          p_est <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = Pred_age, y = Prey_age, size = Est)) +
            ggplot2::geom_point(alpha=0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Estimated", size = "Prop.")

          p_pear <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = Pred_age, y = Prey_age, size = AbsPearson, color = Pearson < 0)) +
            ggplot2::geom_point(alpha=0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Pearson Residuals", size = "Abs(Resid)")

          p1 <- cowplot::plot_grid(p_obs, p_est, p_pear, nrow = 1)
          p1 <- cowplot::ggdraw(p1) + cowplot::draw_label(title, x = 0.5, y = 0.98)

          print(p1)
          plot_list[[length(plot_list) + 1]] <- p1

          if (!is.null(file)) {
            ggplot2::ggsave(paste0(file, "_diet_bubble_Pred", pred_ind, "_Prey", prey_ind, "_Yr", current_yr, ".png"), p1, width = 12, height = 4)
          }
        }
      }
    }
  }
  return(invisible(plot_list))
}

plot_diet_comp2(run_ms_LN_3)

