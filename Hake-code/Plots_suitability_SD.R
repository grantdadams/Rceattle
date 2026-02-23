# 1. SETUP & EXTRACTION
mod      <- run_ms_CSL_Mest
nspp     <- mod$data_list$nspp
nsex_vec <- mod$data_list$nsex          # per-species sex count, length = nspp
sp_names <- mod$data_list$spnames
prey_idx <- 1                           # Hake as prey

dims         <- dim(mod$quantities$suitability)
max_age_pred <- dims[3]
max_age_prey <- dims[4]
nsex_max     <- max(nsex_vec)           # used only for array reconstruction

# 2. EXTRACT FROM SDREPORT
sd_sum       <- summary(mod$sdrep)
suit_indices <- which(rownames(sd_sum) == "suitability_final_yr")
if (length(suit_indices) == 0) stop("'suitability_final_yr' not found in sdreport.")

suit_vals <- sd_sum[suit_indices, 1]
suit_sds  <- sd_sum[suit_indices, 2]

# Sanity check: length should equal nspp*nsex_max * nspp*nsex_max * max_age_pred * max_age_prey
expected_len <- (nspp * nsex_max)^2 * max_age_pred * max_age_prey
if (length(suit_vals) != expected_len) {
  warning(sprintf(
    "Expected %d suitability values but got %d. Check nsex or array dimensions.",
    expected_len, length(suit_vals)
  ))
}

# 3. RESHAPE — column-major order matches TMB/C++ loop order:
#    (pred_sp_sex, prey_sp_sex, pred_age, prey_age)
suit_array    <- array(suit_vals, dim = c(nspp * nsex_max, nspp * nsex_max,
                                          max_age_pred,    max_age_prey))
suit_sd_array <- array(suit_sds,  dim = c(nspp * nsex_max, nspp * nsex_max,
                                          max_age_pred,    max_age_prey))

# 4. BUILD DATA FRAME — recode sex labels so there are only two levels
suit_list <- list()

for (spp in 1:nspp) {

  n_sexes <- nsex_vec[spp]

  for (sex_idx in 1:n_sexes) {

    pred_idx_flat <- ifelse(sex_idx == 1, spp, spp + nspp)

    # Always call sex 1 "Female/Combined" regardless of whether species is sexed
    sex_label <- ifelse(sex_idx == 1, "Female/Combined", "Male")

    sub_ests <- suit_array[pred_idx_flat, prey_idx, , ]
    sub_sds  <- suit_sd_array[pred_idx_flat, prey_idx, , ]

    suit_mean    <- apply(sub_ests, 2, mean, na.rm = TRUE)
    suit_se_mean <- apply(sub_sds,  2, function(x) sqrt(mean(x^2, na.rm = TRUE)))

    suit_list[[paste0(sp_names[spp], "_", sex_label)]] <- data.frame(
      Prey_Age    = seq_len(max_age_prey),
      Predator    = sp_names[spp],
      Sex         = sex_label,
      Suitability = suit_mean,
      SE          = suit_se_mean
    ) %>%
      mutate(
        Lower_95 = pmax(0, Suitability - 1.96 * SE),
        Upper_95 = Suitability + 1.96 * SE
      )
  }
}

suit_df <- bind_rows(suit_list) %>%
  mutate(Sex = factor(Sex, levels = c("Female/Combined", "Male")))

# 5. PLOT
final_colors <- c("Hake" = "#7570B3", "ATF" = "#1B9E77",
                  "Sablefish" = "#E7298A", "CSL" = "#D95F02")

ggplot(suit_df %>% filter(Prey_Age <= 20),
       aes(x = Prey_Age, y = Suitability,
           color = Predator, fill = Predator)) +

  # Ribbon with a faint border so it's visible even when very narrow
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95, group = Sex),
              alpha = 0.25, linewidth = 0.3) +

  geom_line(aes(linetype = Sex), linewidth = 1.2) +

  facet_wrap(~ Predator, scales = "free_y") +

  scale_color_manual(values = final_colors) +
  scale_fill_manual(values  = final_colors) +
  scale_x_continuous(breaks = seq(2, 20, by = 2)) +

  # Only two linetype levels now
  scale_linetype_manual(
    values = c("Female/Combined" = "solid", "Male" = "dashed")
  ) +

  expand_limits(y = 0) +

  labs(
    # title = "Pacific Hake Suitability with 95% Confidence Intervals",
    x        = "Hake Age (Prey)",
    y        = "Suitability Index",
    color    = "Predator",
    fill     = "Predator",
    linetype = "Sex"
  ) +

  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    strip.text      = element_text(face = "bold"),
    panel.spacing   = unit(0.5, "lines")
  )


# 1. Extract log_phi estimates and standard errors
# Filter the summary of the sdreport for the log_phi parameter
phi_summary <- summary(mod$sdrep)
phi_summary <- phi_summary[rownames(phi_summary) == "log_phi", ]
phi_summary
