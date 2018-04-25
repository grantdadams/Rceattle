# This functions runs CEATTLE

Rceattle <- function(){
  #--------------------------------------------------
  # 1. DATA and MODEL PREP
  #--------------------------------------------------
  # Check if require packages are installed and install if not
  if("TMB" %in% rownames(installed.packages()) == FALSE) {install.packages("TMB")}
  # Load required packages
  library(TMB)

  # Compile CEATTLE
  setwd("src")
  compile("CEATTLE_BSAI_v01.cpp")
  dyn.load(dynlib("CEATTLE_BSAI_v01"))
  setwd("../")


}
