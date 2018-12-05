
###################################################
# ADMB  single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_admb <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3, file_name = "Report/GSS_SS_admb")
ss_admb_rep <- ss_admb$sdrep

###################################################
# Run in single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run_no_re <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = FALSE, niter = 3, file_name = "Report/GSS_SS_no_re")
ss_rep_no_re <- ss_run_no_re$rep

###################################################
# Run in single species mode with random recruitment
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run_re <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = TRUE, niter = 3, file_name = "Report/GSS_SS_re")
ss_rep_re <- ss_run_re$rep


###################################################
# ADMB MS
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_admb <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_3_Loops_Files_Corrected/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3, file_name = "Report/GSS_MS_admb")
ms_admb_rep <- ms_admb$rep
plot(ms_admb$quantities$NByage[,1,1], type = "l", col = 1)
lines(ms_admb$quantities$AvgN[,1,1], lty = 1, col = 2)

###################################################
# Run in MS mode using par files
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run_no_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_no_re$final_params, debug = FALSE, random_rec = FALSE, niter = 3, file_name = "Report/GSS_MS_no_re")
ms_rep_no_re <- ms_run_no_re$rep

###################################################
# Run in MS mode using par files and random recruitment
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_re$final_params, debug = FALSE, random_rec = TRUE, niter = 3, file_name = "Report/GSS_MS_re")
ms_rep_re <- ms_run_re$rep
