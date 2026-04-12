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
  "Environmental" = 5,
  "AR1-Deviates" = 6
)

# --- Reverse Mappings (Integer -> String) ---
# Automatically invert the maps above.
# E.g., c("1" = "Logistic", "2" = "Non-parametric", ...)
sel_rev_map <- setNames(names(sel_map), as.character(sel_map))
q_rev_map <- setNames(names(q_map), as.character(q_map))
