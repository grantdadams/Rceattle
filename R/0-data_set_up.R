
#' Function to set up data for package
#' Used internally
set_updata <- function(){
  library(Rceattle)
  BS2017SS <- build_dat(ctlFilename = "asmnt2017_0A_corrected", dat_dir = "~/Documents/GitHub/RceattleRuns/Data/data-raw/BSAI/BS_SS_Files/dat files", nspp = 3)
  BS2017SS$fsh_emp_sel <- BS2017SS$fsh_emp_sel[-1,]
  # "~/Documents/GitHub/RceattleRuns/Data/data-raw/BSAI/BS_MS_5_Loops_Files/dat files"
  # "C:/Users/Grant Adams/Documents/GitHub/RceattleRuns/Data/data-raw/BSAI/BS_MS_5_Loops_Files/dat files"
  BS2017MS <- build_dat(ctlFilename = "asmnt2017_2A_corrected", dat_dir = "~/Documents/GitHub/RceattleRuns/Data/data-raw/BSAI/BS_MS_5_Loops_Files/dat files", nspp = 3) # "C:/Users/Grant Adams/Documents/GitHub/RceattleRuns/Data/data-raw/BSAI/BS_MS_5_Loops_Files/dat files"
  BS2017MS$fsh_emp_sel <- BS2017MS$fsh_emp_sel[-1,]
  usethis::use_data(BS2017MS, BS2017SS, overwrite = TRUE)
}
