library(Rceattle)
###################################################
# ADMB  single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_admb <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3, file_name = "Report/FISH559/Models/ss_admb")
ss_admb_rep <- ss_admb$sdrep

###################################################
# Run in single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run_no_re <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = FALSE, niter = 3, file_name = "Report/FISH559/Models/ss_no_re")
ss_rep_no_re <- ss_run_no_re$rep

###################################################
# Run in single species mode with random recruitment
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run_re <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = TRUE, niter = 3, file_name = "Report/FISH559/Models/ss_re")
ss_rep_re <- ss_run_re$rep


###################################################
# ADMB MS
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_admb <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_3_Loops_Files_Corrected/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3, file_name = "Report/FISH559/Models/ms_admb")


# Loop through number of iterations
load("data/BS_MS_Files/2017_assessment_data_list.RData")
for(i in 3:20){
  ###################################################
  # Run in MS mode using par files
  ms_run_no_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_no_re$final_params, debug = FALSE, random_rec = FALSE, niter = i, file_name = paste0("Report/FISH559/Models/ms_no_re_", i))

  # ###################################################
  # # Run in MS mode using par files and random recruitment
  # ms_run_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ms_run_no_re$final_params, debug = FALSE, random_rec = TRUE, niter = i, file_name = paste0("Report/FISH559/Models/ms_re_", i))
}

