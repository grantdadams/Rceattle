library(Rceattle)
data_list_ss <- build_dat(ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "ceattle_v01_02", dat_dir = "data/BSAI/BS_SS_Files/dat files/", nspp = 3)

ss_run <- Rceattle(TMBfilename = "ceattle_v01_02",
                   data_list = data_list_ss,
                   inits = NULL, # Initial parameters = 0
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   random_rec = FALSE, # No random recruitment
                   msmMode = 0, # Single species mode
                   avgnMode = 0,
                   silent = FALSE)



data_list_ms <- build_dat(ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "ceattle_v01_02", dat_dir = "data/BSAI/BS_MS_Files/dat files/", nspp = 3)

ms_run <- Rceattle(TMBfilename = "ceattle_v01_02",
                     data_list = data_list_ms,
                     inits = ss_run$estimated_params, # Initial parameters from ss run
                     file_name = NULL, # Don't save
                     debug = 0, # Estimate
                     random_rec = FALSE, # No random recruitment
                     niter = 10, # Number of iterations around predation/pop dy functions
                     msmMode = 1, # Multi-species holsman mode
                     avgnMode = 0)

plot_biomass(ceattle_list =  list(ms_run, ss_run), model_names = c("MS", "SS"))
plot_recruitment(ceattle_list =  list(ms_run, ss_run), model_names = c("MS", "SS"))

ms_run_3 <- Rceattle(TMBfilename = "ceattle_v01_02",
                   data_list = data_list_ms,
                   inits = ss_run$estimated_params, # Initial parameters from ss run
                   file_name = NULL, # Don't save
                   debug = 1, # Estimate
                   random_rec = FALSE, # No random recruitment
                   niter = 10, # Number of iterations around predation/pop dy functions
                   msmMode = 3, # Multi-species holsman mode
                   avgnMode = 0)
ms_run_3$quantities$NByage[,1:12,1]



ms_run_4 <- Rceattle(TMBfilename = "ceattle_v01_02",
                     data_list = data_list_ms,
                     inits = ss_run$estimated_params, # Initial parameters from ss run
                     file_name = NULL, # Don't save
                     debug = 1, # Estimate
                     random_rec = FALSE, # No random recruitment
                     niter = 10, # Number of iterations around predation/pop dy functions
                     msmMode = 4, # Multi-species holsman mode
                     avgnMode = 0)
ms_run_4$quantities$NByage[,1:12,1]
ms_run_4$quantities$M2[,1:12,1]



ms_run_5 <- Rceattle(TMBfilename = "ceattle_v01_02",
                     data_list = data_list_ms,
                     inits = ss_run$estimated_params, # Initial parameters from ss run
                     file_name = NULL, # Don't save
                     debug = 1, # Estimate
                     random_rec = FALSE, # No random recruitment
                     niter = 10, # Number of iterations around predation/pop dy functions
                     msmMode = 5, # Multi-species holsman mode
                     avgnMode = 0)
ms_run_5$quantities$NByage[,1:12,1]
