# --- Forward Mappings (String -> Integer) ---
# - Maintains backwards compatabilitiy
sel_map <- c(
  "Fixed" = 0,
  "Logistic" = 1,
  "NonParametric" = 2,
  "DoubleLogistic" = 3,
  "DescendingLogistic" = 4,
  "Hake" = 5,
  "2DAR1" = 6,
  "3DAR1" = 7
)

q_map <- c(
  "Fixed" = 0,
  "Estimated" = 1,
  "Estimated-with-prior" = 2,
  "Analytical" = 3,
  "PowerEquation" = 4,
  "Environmental" = 5,
  "AR1" = 6
)

comp_loglike_map <- c(
  "MultinomialAFSC" = -1,
  "Multinomial" = 0,
  "DirichletMultinomial" = 1
)

fleet_map <- c(
  "Fishery" = 1,
  "Survey" = 2,
  "Off" = 0
)

# --- Reverse Mappings (Integer -> String) ---
# Automatically invert the maps above.
sel_rev_map <- setNames(names(sel_map), as.character(sel_map))
q_rev_map <- setNames(names(q_map), as.character(q_map))
comp_loglike_map_rev_map <- setNames(names(comp_loglike_map), as.character(comp_loglike_map))
fleet_rev_map <- setNames(names(fleet_map), as.character(fleet_map))
