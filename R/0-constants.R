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

tv_sel_map <-c(
  "Off" = 0,
  "IID" = 1,
  "AR1" = 2,
  "Block" = 3,
  "RandomWalk" = 4,
  "RandomWalkAscending" = 5
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

tv_q_map <- c(
  "Off" = 0,
  "IID" = 1,
  "AR1" = 2,
  "Block" = 3,
  "RandomWalk" = 4
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
  "NonEquilibrium"       = 2,
  "FishedNonEquilibrium"       = 3,
  "FishedNonEquilibriumScaled" = 4
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
