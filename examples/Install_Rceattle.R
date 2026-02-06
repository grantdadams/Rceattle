
# Install from GitHub ----
# - Recommended
# Install devtools if you don't already have it
install.packages("devtools")

# Install TMB (instructions: https://github.com/kaskr/adcomp/wiki/Download)
install.packages('TMB', type = 'source')
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper") # install tmb helper functions
# Try "TMB::runExample(all = TRUE)" to see if TMB works

# Install dependencies
install.packages("pacman")
pacman::p_load(dplyr,
               ggplot2,
               MASS,
               oce,
               readxl,
               TMB,
               devtools,
               writexl,
               reshape2,
               gplots,
               tidyr,
               testthat,
               foreach,
               R.utils,
               knitr,
               doParallel)

# Install Rceattle
devtools::install_github("grantdadams/Rceattle")


