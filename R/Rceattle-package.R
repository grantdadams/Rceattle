#' @rawNamespace useDynLib(Rceattle, .registration=TRUE); useDynLib(ceattle_v01_11)
#' @keywords internal
#' @importFrom foreach %dopar%
#' @importFrom dplyr filter mutate select arrange group_by summarise summarize rename left_join right_join full_join bind_rows bind_cols distinct pull n lag vars matches
#' @importFrom plyr rbind.fill
#' @importFrom stats sd median nlminb dnorm qlnorm aggregate rgamma rmultinom
#' @importFrom graphics mtext axis text symbols filled.contour persp curve points image
#' @importFrom grDevices col2rgb rgb
#' @importFrom utils write.csv type.convert
#' @importFrom ggplot2 geom_tile scale_y_continuous coord_equal scale_x_continuous element_rect element_text scale_fill_viridis_c geom_contour geom_contour_filled scale_fill_viridis_d geom_line scale_color_viridis_c
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

utils::globalVariables(c(
  # Data column names used in dplyr/ggplot2 non-standard evaluation
  "Year", "Fleet_code", "Catch", "Species", "Length", "Bin", "Sex",
  "Pred", "Pred_age", "Pred_sex", "Sample_size", "Diet_hat",
  "Stomach_proportion_by_weight", "Diet_weights_mcallister",
  "Prey", "Prey_sex", "Prey_age", "stratum_id", "stomach_id",
  "Selectivity_index", "Catchability", "Q_index", "Fleet_type",
  "Selectivity", "Comp_loglike", "CAAL_loglike", "Wt_index", "Tmp_ind",
  "Wt_name", "Age_transition_name", "Age_transition_index",
  "Age", "True_age", "M",
  "Observation", "Log_sd", "Age0_Length1", "Month", "Length_bin",
  "Estimate_survey_sd", "Survey_sd_prior", "Nselages", "Sel_sd_prior",
  "Estimate_q", "Age_first_selected", "Age_max_selected",
  "Age_max_selected_upper",
  "Time_varying_sel", "Time_varying_sel_sd_prior",
  "Sel_curve_pen1", "Sel_curve_pen2",
  # Loop iterators and anonymous function variables
  "i", "x", "spp", "sim", "var",
  # Default argument dataset references
  "ms_run", "ss_run",
  # oce function used without explicit namespace
  "oce.colorsDensity"
))
