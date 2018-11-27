load("data/BS_MS_1_Loops_Files/CEATTLE_results.Rdata")
ms_1_loops <- tmp
load("data/BS_MS_3_Loops_Files/CEATTLE_results.Rdata")
ms_3_loops <- tmp
load("data/BS_MS_5_Loops_Files/CEATTLE_results.Rdata")
ms_5_loops <- tmp
load("data/BS_MS_10_Loops_Files/CEATTLE_results.Rdata")
ms_10_loops <- tmp
load("data/BS_MS_20_Loops_Files/CEATTLE_results.Rdata")
ms_20_loops <- tmp

ms_1_loops$suit_main_1_1
ms_3_loops$suit_main_1_1
ms_5_loops$suit_main_1_1
ms_10_loops$suit_main_1_1


load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run_3 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_3_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3)
ms_rep_3 <- ms_run_3$quantities

ms_run_5 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_5_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 5)
ms_rep_5 <- ms_run_5$quantities

ms_run_10 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_10_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 10)
ms_rep_10 <- ms_run_10$quantities

ms_run_20 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_20_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 10)
ms_rep_20 <- ms_run_10$quantities


year <- 2
pollock_nbyage_yr1 <- data.frame(ADMB_3_iter = as.numeric(ms_3_loops$NByage_1[year,]),
                     TMB_3_iter = signif(ms_rep_3$NByage[year,1:12,1],6),
                     ADMB_5_iter = as.numeric(ms_5_loops$NByage_1[year,]),
                     TMB_5_iter = signif(ms_rep_5$NByage[year,1:12,1],6),
                     ADMB_10_iter = as.numeric(ms_10_loops$NByage_1[year,]),
                     TMB_10_iter = signif(ms_rep_10$NByage[year,1:12,1],6),
                     ADMB_20_iter = as.numeric(ms_20_loops$NByage_1[year,]),
                     TMB_20_iter = signif(ms_rep_20$NByage[year,1:12,1],6))
