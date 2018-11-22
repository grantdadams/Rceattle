load("data/BS_MS_Files/CEATTLE_results.Rdata")
ms_3_loops <- tmp
load("data/BS_MS_5_Loops_Files/CEATTLE_results.Rdata")
ms_5_loops <- tmp
load("data/BS_MS_10_Loops_Files/CEATTLE_results.Rdata")
ms_10_loops <- tmp

ms_3_loops$R_1
ms_5_loops$R_1
ms_10_loops$R_1


load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run_10 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_10_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 10)
ms_rep_10 <- ms_run_10$quantities

ms_run_5 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_5_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 5)
ms_rep_5 <- ms_run_5$quantities

ms_run_3 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3)
ms_rep_3 <- ms_run_3$quantities

pollock_rec <- data.frame(ADMB_3_iter = as.numeric(ms_3_loops$R_1),
                     TMB_3_iter = signif(ms_rep_3$R[1,],6),
                     ADMB_5_iter = as.numeric(ms_5_loops$R_1),
                     TMB_5_iter = signif(ms_rep_5$R[1,],6),
                     ADMB_10_iter = as.numeric(ms_10_loops$R_1),
                     TMB_10_iter = signif(ms_rep_10$R[1,],6))
