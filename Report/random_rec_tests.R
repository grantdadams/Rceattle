library(Rceattle)
data_list_ss <- build_dat(ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "ceattle_v01_02", dat_dir = "data-raw/BSAI/BS_SS_Files/dat files/", nspp = 3)


ss_run <- Rceattle(TMBfilename = "ceattle_v01_02",
                   data_list = data_list_ss,
                   inits = NULL, # Initial parameters = 0
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   random_rec = FALSE, # No random recruitment
                   msmMode = 0, # Single species mode
                   avgnMode = 0,
                   silent = TRUE)

ss_run_re <- Rceattle(TMBfilename = "ceattle_v01_02",
                      data_list = data_list_ss,
                      inits = NULL, # Initial parameters = 0
                      file_name = NULL, # Don't save
                      debug = 0, # Estimate
                      random_rec = TRUE, # No random recruitment
                      msmMode = 0, # Single species mode
                      avgnMode = 0,
                      silent = TRUE)


# Multi-species
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

ms_run_re <- Rceattle(TMBfilename = "ceattle_v01_02",
                      data_list = data_list_ms,
                      inits = ss_run$estimated_params, # Initial parameters from ss run
                      file_name = NULL, # Don't save
                      debug = 0, # Estimate
                      random_rec = TRUE, # No random recruitment
                      niter = 10, # Number of iterations around predation/pop dy functions
                      msmMode = 1, # Multi-species holsman mode
                      avgnMode = 0)




plot_biomass(ceattle_list =  list(ss_run, ss_run_re, ms_run, ms_run_re), model_names = c("SS", "SS_re", "MS", "MS_re"), file_name = "Biomass")
plot_recruitment(ceattle_list =  list(ss_run, ss_run_re, ms_run, ms_run_re), model_names = c("SS", "SS_re", "MS", "MS_re"), ci_col = c(1:4), file_name = "Recruitment")
