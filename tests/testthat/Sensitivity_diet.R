#############################################
### Sensitivity Diet
## Using upper and lower 95 of the bootsrapped average.
#############################################

#diet_upper95 <- Rceattle::read_data(file = "sensitivity/070725_ATF_Hake_model_Neff_upper95.xlsx")
#diet_lower95 <- Rceattle::read_data(file = "sensitivity/070725_ATF_Hake_model_Neff_lower95.xlsx")

diet_upper95 <- Rceattle::read_data(file = "sensitivity/070725_Sablefish_Hake_model_upper95.xlsx")
diet_lower95 <- Rceattle::read_data(file = "sensitivity/070725_Sablefish_Hake_model_lower95.xlsx")

## upper 95 diet estimates
test_data$diet_data <- diet_upper95$diet_data

run_ms_LN_upper <- Rceattle::fit_mod(data_list = test_data,
                                 inits = inits, # Initial parameters from single species ests
                                 map = map,
                                 M1Fun = build_M1(M1_model = 1,
                                                  updateM1 = TRUE,
                                                  M1_use_prior = FALSE,
                                                  M2_use_prior = FALSE),
                                 file = NULL, # Don't save
                                 estimateMode = 0, # estimate
                                 niter = 3, # 3 iterations around population and predation dynamics
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # MSVPA based
                                 loopnum = 5,
                                 phase = TRUE,
                                 suitMode = c(0, 4), # empirical + LN suitability
                                 initMode = 2,
                                 verbose = 1)

run_ms_LN_upper$quantities$jnll_comp
save(run_ms_LN_upper, file = "Sensitivity/models/SBF_run_ms_LN_estM_upper95.Rdata")

## lower 95 diet estimates
test_data$diet_data <- diet_lower95$diet_data

run_ms_LN_lower <- Rceattle::fit_mod(data_list = test_data,
                                     inits = inits, # Initial parameters from single species ests
                                     map = map,
                                     M1Fun = build_M1(M1_model = 1,
                                                      updateM1 = TRUE,
                                                      M1_use_prior = FALSE,
                                                      M2_use_prior = FALSE),
                                     file = NULL, # Don't save
                                     estimateMode = 0, # estimate
                                     niter = 3, # 3 iterations around population and predation dynamics
                                     random_rec = FALSE, # No random recruitment
                                     msmMode = 1, # MSVPA based
                                     loopnum = 5,
                                     phase = TRUE,
                                     suitMode = c(0, 4), # empirical + LN suitability
                                     initMode = 2,
                                     verbose = 1)

run_ms_LN_lower$quantities$jnll_comp
save(run_ms_LN_lower, file = "Sensitivity/models/SBF_run_ms_LN_estM_lower95.Rdata")

models <- list(run_ms_LN_3, run_ms_LN_lower, run_ms_LN_upper)
model_names <- c("Base Model", "Diet_lower95", "Diet_upper95")

models <- list(run_ms_LN_3, run_ms_LN_lower)
model_names <- c("Base Model", "Diet_lower95")

# --------------------------------------------------------------------------
# STEP 1: Master Data Extraction
# This loop builds one clean data frame with all the results we need.
# --------------------------------------------------------------------------
data_frames_list <- list()

for (i in 1:length(models)) {

  model_object <- models[[i]]
  model_name <- model_names[i]

  # Get hindcast years for this model
  yrs <- model_object$data_list$styr:model_object$data_list$endyr

  # --- Loop through both SPECIES (1 and 2) for standard quantities ---
  for (species_index in 1:model_object$data_list$nspp) {

    species_name <- model_object$data_list$spnames[species_index]
    quantities_to_extract <- c("ssb", "biomass", "R")
    variable_names <- c("Spawning Biomass", "Total Biomass", "Recruitment")

    for (j in 1:length(quantities_to_extract)) {
      quantity_name <- quantities_to_extract[j]

      # Extract values for the current species and hindcast years
      value_vec <- model_object$quantities[[quantity_name]][species_index, 1:length(yrs)]

      # Extract the corresponding errors
      error_vector_all <- model_object$sdrep$sd[which(names(model_object$sdrep$value) == quantity_name)]
      num_errors_per_spp <- length(error_vector_all) / model_object$data_list$nspp
      start_index <- (species_index - 1) * num_errors_per_spp + 1
      end_index <- species_index * num_errors_per_spp
      error_vec <- error_vector_all[start_index:end_index][1:length(yrs)]

      # Create the data frame
      temp_df <- data.frame(
        model = model_name,
        species = species_name,
        year = yrs,
        variable = variable_names[j],
        value = value_vec,
        error = error_vec
      )
      data_frames_list <- append(data_frames_list, list(temp_df))
    }
  }

  # --- Add the custom consumption metrics ---
  b_eaten <- model_object$quantities$B_eaten

  # Hake Consumed by ATF
  hake_eaten_by_atf <- sapply(1:length(yrs), function(yr_index) sum(b_eaten[c(2, 4), 1, , , yr_index]))
  temp_df_atf_cons <- data.frame(model = model_name, species = "ATF", year = yrs, variable = "Consumption", value = hake_eaten_by_atf, error = NA)
  data_frames_list <- append(data_frames_list, list(temp_df_atf_cons))

  # Hake Cannibalism (Corrected based on our debugging)
  hake_eaten_by_hake <- sapply(1:length(yrs), function(yr_index) sum(b_eaten[1, 1, , , yr_index]))
  temp_df_hake_cons <- data.frame(model = model_name, species = "Hake", year = yrs, variable = "Consumption", value = hake_eaten_by_hake, error = NA)
  data_frames_list <- append(data_frames_list, list(temp_df_hake_cons))
}

# Combine into the final data frame
sensitivity_results <- do.call("rbind", data_frames_list)


# --------------------------------------------------------------------------
# STEP 2: Plot 1 - Consumption
# --------------------------------------------------------------------------

# Filter the data to ONLY include consumption metrics from 1985 onward
consumption_data <- sensitivity_results %>%
  filter(variable == "Consumption") %>%
  mutate(
    value = value / 1000000,
    panel_title = case_when(
      species == "ATF" ~ "Consumption of Hake by ATF",
      species == "Hake" ~ "Hake Cannibalism"
    )
  )

# Define the desired order for the legend
#model_order <- c("Base Model", "Diet_lower95", "Diet_upper95")
model_order <- c("Base Model", "Diet_lower95")

# Define the specific colors for each model (using a named vector is best practice)
model_colors <- c(
  #"Diet_upper95" = "#440154FF",  # A dark Viridis purple
  "Diet_lower95" = "#21908CFF", # A Viridis teal
  "Base Model" = "#FDE725FF"   # A bright Viridis yellow
)

# Create the consumption plot
consumption_plot <- ggplot(consumption_data, aes(x = year, y = value, color = model)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ panel_title, ncol = 1, scales = "free_y") +
  labs(
    x = "Year",
    y = "Consumption (million tons)",
    title = "Hake Consumption Sensitivity",
    color = "Model Scenario"
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = model_colors)  +
  theme_bw() +
  theme(legend.position = "top")

# Display the plot
print(consumption_plot)

# --------------------------------------------------------------------------
# STEP 3: Plot 2 - Hake Population Dynamics
# --------------------------------------------------------------------------

# Filter data for Hake population dynamics from 1985 onward
hake_popdy_data <- sensitivity_results %>%
  filter(
    species == "Hake",
    variable != "Consumption",
    #year >= 1990
  ) %>%
  mutate(
    value = value / 1000000,
    error = error / 1000000,
    variable_label = case_when(
      variable == "Spawning Biomass" ~ "Spawning Biomass (Mt)",
      variable == "Total Biomass" ~ "Total Biomass (Mt)",
      variable == "Recruitment" ~ "Recruitment (millions)"
    )
  ) %>%
  mutate(
    min = value - (2 * error),
    max = value + (2 * error)
  )


# Create the Hake population plot
hake_plot <- ggplot(hake_popdy_data, aes(x = year, y = value, color = model)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ variable_label, ncol = 1, scales = "free_y") +
  labs(
    x = "Year",
    y = "",
    title = "Hake Population Dynamics Sensitivity (1980-2019)",
    color = "Model Scenario",
    fill = "Model Scenario"
  ) +
  scale_y_continuous(limits = c(0, NA), labels = scales::comma) +
  theme_bw() +
  theme(legend.position = "top")

# Display the plot
print(hake_plot)

plot_biomass(Rceattle = models, model_names = model_names, add_ci = T)
plot_ssb(Rceattle = models, model_names = model_names, add_ci = T)
plot_recruitment(Rceattle = models, model_names = model_names, add_ci = T)

plot_diet_comp(run_ms_LN_upper)
plot_diet_comp(run_ms_LN_lower)

######### check diet comps and Macallister Ianelli
base_model<- run_ms_LN_3
base_model$estimated_params$diet_comp_weights
base_model$data_list$Diet_comp_weights
base_model$data_list$Diet_weights_mcallister
base_model$quantities$jnll
base_model$estimated_params$log_phi
base_model$quantities$vulnerability
base_model$quantities$vulnerability_other
base_model$quantities$jnll_comp

#### plot models ###
models <- list(ms_run, run_ms_LN_upper)
model_names <- c("Cannib", "Diet_upper95")

# Plot biomass trajectory
plot_biomass(Rceattle = models, model_names = model_names) #Now biomass looks alike
plot_biomass(Rceattle = models, model_names = model_names, add_ci = TRUE)
plot_ssb(Rceattle = models, model_names = model_names, add_ci = TRUE)
plot_recruitment(Rceattle = models, model_names = model_names, add_ci = TRUE)

plot_b_eaten_prop(Rceattle = models, model_names = model_names)


##########==================
# --------------------------------------------------------------------------
# STEP: Create Stacked Barplot of Hake Consumption by Predator
# --------------------------------------------------------------------------

# Extract consumption data and reshape for stacked barplot
consumption_stacked_data <- sensitivity_results %>%
  filter(variable == "Consumption" & year >= 1985) %>%
  mutate(
    value_mt = value / 1000000,  # Convert to million tons
    predator = case_when(
      species == "ATF" ~ "ATF",
      species == "Hake" ~ "Hake (Cannibalism)"
    )
  ) %>%
  select(model, year, predator, value_mt)

# Calculate total consumption per year per model
total_consumption <- consumption_stacked_data %>%
  group_by(model, year) %>%
  summarise(total = sum(value_mt, na.rm = TRUE), .groups = "drop")

# Create stacked barplot
consumption_barplot <- ggplot(consumption_stacked_data, aes(x = year, y = value_mt, fill = predator)) +
  geom_col(position = "stack", width = 0.8) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    x = "Year",
    y = "Hake Consumption (Million Tons)",
    title = "Total Hake Consumption by Predator (1985-2019)",
    fill = "Predator"
  ) +
  scale_fill_manual(values = c("Hake (Cannibalism)" = "#E31A1C", "ATF" = "#1F78B4")) +
  scale_y_continuous(limits = c(0, NA), labels = scales::comma) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(consumption_barplot)

# --------------------------------------------------------------------------
# Alternative: Summary barplot by time periods
# --------------------------------------------------------------------------

# Calculate mean consumption by predator for different time periods
consumption_summary <- consumption_stacked_data %>%
  mutate(
    period = case_when(
      year >= 1985 & year <= 1994 ~ "1985-1994",
      year >= 1995 & year <= 2004 ~ "1995-2004",
      year >= 2005 & year <= 2014 ~ "2005-2014",
      year >= 2015 & year <= 2019 ~ "2015-2019",
      TRUE ~ "Other"
    )
  ) %>%
  filter(period != "Other") %>%
  group_by(model, period, predator) %>%
  summarise(mean_consumption = mean(value_mt, na.rm = TRUE), .groups = "drop")

# Create summary barplot by periods
summary_barplot <- ggplot(consumption_summary, aes(x = period, y = mean_consumption, fill = predator)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~ model, ncol = 3) +
  labs(
    x = "Time Period",
    y = "Mean Hake Consumption (Million Tons)",
    title = "Mean Hake Consumption by Predator Across Time Periods",
    fill = "Predator"
  ) +
  scale_fill_manual(values = c("Hake (Cannibalism)" = "#E31A1C", "ATF" = "#1F78B4")) +
  scale_y_continuous(limits = c(0, NA), labels = scales::comma) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the summary plot
print(summary_barplot)

# --------------------------------------------------------------------------
# Data table summary for reference
# --------------------------------------------------------------------------

# Create a summary table of total consumption by model and predator
consumption_table <- consumption_stacked_data %>%
  group_by(model, predator) %>%
  summarise(
    total_consumption_mt = sum(value_mt, na.rm = TRUE),
    mean_annual_consumption_mt = mean(value_mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(model, desc(total_consumption_mt))

print("Summary of Hake Consumption by Predator:")
print(consumption_table)

devtools::document()
