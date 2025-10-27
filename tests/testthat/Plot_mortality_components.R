# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(stringr)
library(RColorBrewer)

plot_mortality_components <- function(Rceattle,
                                      species = NULL,
                                      ages = 1:5,
                                      file = NULL,
                                      minyr = NULL,
                                      incl_proj = FALSE,
                                      width = 10,
                                      height = 8,
                                      ...) {

  # --- 1. Data Extraction & Setup ---

  if(class(Rceattle) != "Rceattle"){
    stop("This function supports a single Rceattle object.")
  }

  data_list <- Rceattle$data_list
  spnames <- data_list$spnames
  nspp <- data_list$nspp

  if(is.null(minyr)){ minyr <- data_list$styr }
  endyr <- if(incl_proj) data_list$projyr else data_list$endyr
  years <- minyr:endyr

  if(is.null(species)){
    species <- 1:nspp
  }

  # --- 2. Extract and Tidy All Mortality Arrays ---

  # Helper function to convert a mortality array to a tidy data frame
  tidy_mortality <- function(mort_array, mort_name) {
    df <- as.data.frame.table(mort_array, responseName = "Mortality", stringsAsFactors = FALSE)
    colnames(df) <- c("Species_Name_Raw", "Sex", "Age_Str", "Year", "Mortality")

    df %>%
      # !! --- THIS IS THE KEY FIX --- !!
      # Create clean Species and Age columns by parsing the raw text
      mutate(
        Year = as.numeric(Year),
        Age = as.numeric(gsub("\\D", "", Age_Str)),
        Species = case_when(
          grepl("Hake", Species_Name_Raw) ~ "Hake",
          grepl("ATF", Species_Name_Raw) ~ "ATF",
          # Add other species from your 'spnames' if needed
          TRUE ~ "Unknown"
        ),
        Type = mort_name
      ) %>%
      # Filter for the primary "female/combined" sex data to match original plot logic
      filter(grepl("Sex combined or females", Sex)) %>%
      # Filter for desired species, ages, and years
      filter(Species %in% spnames[species], Age %in% ages, Year %in% years) %>%
      select(Species, Age, Year, Type, Mortality)
  }

  # Extract M1, M2, and Total M using the corrected helper function
  m1_df <- tidy_mortality(Rceattle$quantities$M1_at_age, "M1 (Natural)")
  m2_df <- tidy_mortality(Rceattle$quantities$M2_at_age, "M2 (Predation)")
  total_m_df <- tidy_mortality(Rceattle$quantities$M_at_age, "Total M")

  combined_mortality_df <- bind_rows(m1_df, m2_df, total_m_df) %>%
    mutate(Type = factor(Type, levels = c("M1 (Natural)", "M2 (Predation)", "Total M")))

  if(nrow(combined_mortality_df) == 0) {
    stop("Data frame is empty after filtering. This is unexpected with the new parsing logic.")
  }

  # --- 3. Plotting with ggplot2 ---

  p <- ggplot(combined_mortality_df, aes(x = Year, y = Mortality, color = Type)) +
    geom_line(linewidth = 1.1) +
    facet_grid(Species ~ Age, scales = "free_y", labeller = label_both) +
    labs(
      title = "Mortality Components by Species and Age",
      x = "Year",
      y = "Mortality Rate",
      color = "Mortality Component"
    ) +
    expand_limits(y = 0) +
    scale_color_manual(values = c("M1 (Natural)" = "cornflowerblue", "M2 (Predation)" = "yellow4", "Total M" = "black")) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )

  if (!is.null(file)) {
    ggsave(
      filename = paste0(file, "_mortality_components.png"),
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 300
    )
  }

  print(p)
  invisible(p)
}

# This will now produce the correct plot
plot_mortality_components(
  Rceattle = run_ms_LN_3,
  species = 1, # Hake
  ages = c(1:5)     # Age 1
)

plot_m_at_age(run_ms_LN_3, age=2, species = 1)

## ACCESS DATA
## We need to extract M2 per predator.
dim(run_ms_LN_3$quantities$M2_at_age) #Pred, sex, age, year
run_ms_LN_3$quantities$M2_at_age[1,1,1,1] #this is predation mortality on hake
run_ms_LN_3$quantities$M2_at_age[2,2,1,1] #this is predations mortality on ATF = 0
#but the predation mortality on hake is combined (hake cannib and ATF)

#accees total mortality
dim(run_ms_LN_3$quantities$M_at_age)
run_ms_LN_3$quantities$M_at_age[1,1,1,1] #this is total mortality on hake (M1 + M2)
run_ms_LN_3$quantities$M_at_age[2,2,1,1] #this is total mortality on ATF (we dont have M2, hence this is = to M1)

#acces residual mortality (M1)
dim(run_ms_LN_3$quantities$M1_at_age)
run_ms_LN_3$quantities$M1_at_age[1,1,1,1] #this is residual mortality on hake (M1)
run_ms_LN_3$quantities$M1_at_age[2,2,1,1] #this is residual mortality on ATF (M1)

#We can get predation by predator in here, NOW this is a 5dim array
dim(run_ms_LN_3$quantities$M2_prop) #pred:sex, prey:sex, pred age, prey age, year
run_ms_LN_3$quantities$M2_prop[1,1,1,1,1]

plot_m2_at_age_prop(run_ms_LN_3)

### Partition M2 by predator
plot_M2_by_predator <- function(Rceattle,
                                species = NULL,
                                ages = 1:5,
                                file = NULL,
                                minyr = NULL,
                                incl_proj = FALSE,
                                width = 10,
                                height = 8,
                                ...) {

  # --- 1. Data Extraction & Setup ---

  if(class(Rceattle) != "Rceattle"){
    stop("This function supports a single Rceattle object.")
  }

  data_list <- Rceattle$data_list
  spnames <- data_list$spnames
  nspp <- data_list$nspp

  if(is.null(minyr)){ minyr <- data_list$styr }
  endyr <- if(incl_proj) data_list$projyr else data_list$endyr
  years <- minyr:endyr

  if(is.null(species)){
    species <- 1:nspp
  }

  target_spnames <- spnames[species]

  # --- 2. Helper Function for Tidying M2_prop ---

  tidy_m2_prop <- function(m2_prop_array, m1_dimnames, all_spnames, target_spnames, target_ages, target_years) {

    species_names <- m1_dimnames[[1]]
    sex_names <- m1_dimnames[[2]]
    age_labels <- m1_dimnames[[3]]
    year_labels <- m1_dimnames[[4]]

    sp_sex_labels <- paste(
      rep(species_names, times = length(sex_names)),
      rep(sex_names, each = length(species_names)),
      sep = "_"
    )

    dn_m2_prop <- list(
      Pred_SpSex = sp_sex_labels,
      Prey_SpSex = sp_sex_labels,
      Pred_Age   = age_labels,
      Prey_Age   = age_labels,
      Year       = year_labels
    )

    tryCatch({
      dimnames(m2_prop_array) <- dn_m2_prop
    }, error = function(e) {
      print("Error applying dimnames. Array dimensions and names do not match.")
      print(paste("M2_prop dim 1 length:", dim(m2_prop_array)[1], " | Generated names length:", length(sp_sex_labels)))
      print(paste("M2_prop dim 3 length:", dim(m2_prop_array)[3], " | Generated names length:", length(age_labels)))
      stop(e)
    })

    df <- as.data.frame.table(m2_prop_array, responseName = "Mortality", stringsAsFactors = FALSE)

    if(nrow(df) == 0) return(data.frame())

    spname_pattern <- paste(all_spnames, collapse = "|")

    df_processed <- df %>%
      mutate(
        Year = as.numeric(Year),
        Prey_Age = as.numeric(gsub("\\D", "", Prey_Age)),
        Pred_Species = stringr::str_extract(Pred_SpSex, spname_pattern),
        Prey_Species = stringr::str_extract(Prey_SpSex, spname_pattern)
      ) %>%
      filter(
        Prey_Species %in% target_spnames,
        Prey_Age %in% target_ages,
        Year %in% target_years,
        grepl("Sex combined or females", Prey_SpSex)
      )

    if(nrow(df_processed) == 0) return(data.frame())

    df_agg <- df_processed %>%
      group_by(Prey_Species, Prey_Age, Year, Pred_Species) %>%
      summarise(Mortality = sum(Mortality, na.rm = TRUE), .groups = 'drop') %>%
      filter(Mortality > 1e-10) %>%
      rename(Species = Prey_Species, Age = Prey_Age) %>%
      mutate(Type = paste0("M2 (from ", Pred_Species, ")")) %>%
      select(Species, Age, Year, Type, Mortality)

    return(df_agg)
  }

  # --- 3. Extract, Tidy, and Combine Data (ONLY M2) ---

  m1_dims <- dimnames(Rceattle$quantities$M1_at_age)

  if(is.null(Rceattle$quantities$M2_prop)) {
    stop("`Rceattle$quantities$M2_prop` not found. Did you re-run the model with the modified TMB code?")
  }

  combined_mortality_df <- tidy_m2_prop(
    m2_prop_array = Rceattle$quantities$M2_prop,
    m1_dimnames = m1_dims,
    all_spnames = spnames,
    target_spnames = target_spnames,
    target_ages = ages,
    target_years = years
  )

  if(nrow(combined_mortality_df) == 0) {
    stop("Data frame is empty after filtering. Check species, ages, and years.")
  }

  type_levels <- sort(unique(combined_mortality_df$Type))

  combined_mortality_df <- combined_mortality_df %>%
    mutate(Type = factor(Type, levels = type_levels))

  # --- 4. Dynamic Color Palette (ONLY M2 colors) ---

  n_m2_types <- length(type_levels)
  all_colors <- RColorBrewer::brewer.pal(max(3, n_m2_types), "Dark2")
  names(all_colors) <- type_levels

  # --- 5. Plotting with ggplot2 ---

  p <- ggplot(combined_mortality_df, aes(x = Year, y = Mortality, color = Type)) +
    geom_line(linewidth = 1.1) +
    facet_grid(Species ~ Age, scales = "free_y", labeller = label_both) +
    labs(
      title = "M2 (Predation Mortality) by Predator and Prey Age",
      x = "Year",
      y = "M2 Rate",
      color = "Predator"
    ) +
    expand_limits(y = 0) +
    scale_color_manual(values = all_colors) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )

  # --- 6. Calculate Mean Tables (NEW SECTION) ---

  # Mean by prey species, prey age, and predator (averaged over years)
  m2_mean_by_age <- combined_mortality_df %>%
    group_by(Species, Age, Type) %>%
    summarise(Mean_Mortality = mean(Mortality, na.rm = TRUE), .groups = 'drop')

  # Mean by prey species and predator (averaged across specified ages & years)
  m2_mean_overall <- combined_mortality_df %>%
    group_by(Species, Type) %>%
    summarise(Mean_Mortality = mean(Mortality, na.rm = TRUE), .groups = 'drop')

  # --- 7. Save Plot, Print Plot, and Return List (UPDATED) ---

  if (!is.null(file)) {
    ggsave(
      filename = paste0(file, "_M2_by_predator.png"),
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 300
    )
  }

  print(p)

  # Return a list containing the plot and the data frames
  invisible(list(
    plot = p,
    m2_timeseries = combined_mortality_df,
    m2_mean_by_age = m2_mean_by_age,
    m2_mean_overall = m2_mean_overall
  ))
}

# Run the function and store the output in a new object 'm2_results'
m2_results <- plot_M2_by_predator(
  Rceattle = run_ms_LN_3,
  species = 1, # Hake
  ages = 1:20)

# The plot will still appear automatically.
# Now you can access the tables from the list:

# 1. View the full time-series data
print(m2_results$m2_timeseries)

# 2. View the mean M2 by predator and prey age
print(m2_results$m2_mean_by_age)
a<- m2_results$m2_mean_by_age

# 3. View the mean M2 by predator, averaged over all specified ages
print(m2_results$m2_mean_overall)
