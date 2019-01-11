library(Rceattle)

# Fit initial models
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

# Simulation
nsim <- 100
ms_sim <- list()
ms_sim_run <- list()

for(i in 41:nsim){
  ms_sim[[i]] <- sim_mod(ms_run)
  print(paste0("###############"))
  print(paste0("Fitting sim ",i))
  print(paste0("###############"))
  ms_sim_run[[i]] <- Rceattle(TMBfilename = "ceattle_v01_02",
                         data_list = ms_sim[[i]],
                         inits = ms_run$estimated_params, # Initial parameters = 0
                         file_name = NULL, # Don't save
                         debug = 0, # Estimate
                         random_rec = FALSE, # No random recruitment
                         msmMode = 1, # Holsman MS mode
                         avgnMode = 0,
                         silent = TRUE)

}

