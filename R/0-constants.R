# --- Forward Mappings (String -> Integer) ---
# - Maintains backwards compatibility
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

# Initial age-structure mode
# 0 = Free parameters for initial age-structure
# 1 = Equilibrium, no init devs, Finit = 0 (unfished)
# 2 = Equilibrium + init devs, Finit = 0  [default]
# 3 = Non-equilibrium: Finit estimated, init devs included
# 4 = Non-equilibrium: Finit scales R0
initMode_map <- c(
  "FreeParams"           = 0,
  "Equilibrium"          = 1,
  "EquilibriumDev"       = 2,
  "NonEquilibrium"       = 3,
  "NonEquilibriumScaled" = 4
)

# Predator-prey suitability mode (per predator species)
# 1 and 3 are blocked in data_check() until growth-model validation is added
suitMode_map <- c(
  "Empirical"       = 0,
  "GammaLength"     = 1,   # NOT YET AVAILABLE
  "GammaWeight"     = 2,
  "LognormalLength" = 3,   # NOT YET AVAILABLE
  "LognormalWeight" = 4,
  "NormalLength"    = 5,
  "NormalWeight"    = 6
)

# Stock-recruit relationship
srr_map <- c(
  "Mean"            = 0,
  "MeanEnv"         = 1,
  "BevertonHolt"    = 2,
  "BevertonHoltEnv" = 3,
  "Ricker"          = 4,
  "RickerEnv"       = 5
)

# Harvest control rule
hcr_map <- c(
  "NoFishing"    = 0,
  "CMSY"         = 1,
  "ConstantF"    = 2,
  "ConstantFSSB" = 3,
  "ConstantFSPR" = 4,
  "NPFMC"        = 5,
  "PFMC"         = 6,
  "SESSF"        = 7
)

# --- Reverse Mappings (Integer -> String) ---
# Automatically invert the maps above.
sel_rev_map         <- setNames(names(sel_map),         as.character(sel_map))
q_rev_map           <- setNames(names(q_map),           as.character(q_map))
comp_loglike_rev_map <- setNames(names(comp_loglike_map), as.character(comp_loglike_map))
fleet_rev_map       <- setNames(names(fleet_map),       as.character(fleet_map))
initMode_rev_map    <- setNames(names(initMode_map),    as.character(initMode_map))
suitMode_rev_map    <- setNames(names(suitMode_map),    as.character(suitMode_map))
srr_rev_map         <- setNames(names(srr_map),         as.character(srr_map))
hcr_rev_map         <- setNames(names(hcr_map),         as.character(hcr_map))


# # Helper: convert a single string value using a map, pass integers through unchanged
# .conv <- function(x, map) {
#   if (is.character(x) && x %in% names(map)) unname(map[[x]]) else x
# }
#
# # initMode (scalar)
# data_list$initMode <- as.integer(.conv(data_list$initMode, initMode_map))
