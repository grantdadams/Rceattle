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
                   silent = TRUE)

data_list_ms <- build_dat(ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "ceattle_v01_02", dat_dir = "data/BSAI/BS_MS_Files/dat files/", nspp = 3)

ms_run <- Rceattle(TMBfilename = "ceattle_v01_02",
                   data_list = data_list_ms,
                   inits = ss_run$estimated_params, # Initial parameters from ss run
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   random_rec = FALSE, # No random recruitment
                   niter = 10, # Number of iterations around predation/pop dy functions
                   msmMode = 1, # Multi-species holsman mode
                   avgnMode = 0,
                   silent = TRUE)

# # Simulation
# nsim <- 100
# ms_sim <- list()
# ms_sim_run <- list()
#
# for(i in 1:nsim){
#   ms_sim[[i]] <- sim_mod(ms_run)
#   print(paste0("###############"))
#   print(paste0("Fitting sim ",i))
#   print(paste0("###############"))
#   mod <- Rceattle(TMBfilename = "ceattle_v01_02",
#                          data_list = ms_sim[[i]],
#                          inits = ms_run$estimated_params, # Initial parameters = 0
#                          file_name = NULL, # Don't save
#                          debug = 0, # Estimate
#                          random_rec = FALSE, # No random recruitment
#                          msmMode = 1, # Holsman MS mode
#                          avgnMode = 0,
#                          silent = TRUE)
#   mod$opt <- NULL
#   mod$obj <- NULL
#   mod$map <- NULL
#   mod$bounds <- NULL
#   mod$initial_params <- NULL
#
#   ms_sim_run[[i]] <- mod
# }

# save(ms_sim_run, file = "Report/2019_1_Simulation/ms_sim_runs.RData")
# comparison <- compare_sim(operating_mod = ms_run, simulation_mods = ms_sim_run, object = "quantities")


load("Report/2019_1_Simulation/ms_sim_runs.RData")
comparison <- compare_sim(operating_mod = ms_run, simulation_mods = ms_sim_run, object = "quantities")



# Time-series
ms_run_mean <- ms_run
ms_run_mean$quantities <- comparison$Mean

ms_run_median <- ms_run
ms_run_median$quantities <- comparison$Median

plot_biomass( ceattle_list = c(list(ms_run, ms_run_mean, ms_run_median), ms_sim_run, list(ms_run, ms_run_mean, ms_run_median)), line_col = c(1, 2, 3, rep("grey", length(ms_sim_run)), 1, 2, 3), model_names = c("OM", "Mean SM", "Median SM", "SMs"), file_name = "Report/2019_1_Simulation/time")
plot_recruitment( ceattle_list = c(list(ms_run, ms_run_mean, ms_run_median), ms_sim_run, list(ms_run, ms_run_mean, ms_run_median)), line_col = c(1, 2, 3, rep("grey", length(ms_sim_run)), 1, 2, 3), model_names = c("OM", "Mean SM", "Median SM", "SMs"), file_name = "Report/2019_1_Simulation/time")
plot_selectivity( ceattle_list = c(list(ms_run, ms_run_mean, ms_run_median), ms_sim_run, list(ms_run, ms_run_mean, ms_run_median)), line_col = c(1, 2, 3, rep("grey", length(ms_sim_run)), 1, 2, 3), model_names = c("OM", "Mean SM", "Median SM", "SMs"), file_name = "Report/2019_1_Simulation/time")


# MSE
ms_run_mse <- ms_run
ms_run_mse$quantities <- comparison$MSE

plot_biomass( ceattle_list = list(ms_run_mse), line_col = c(1), model_names = c("MSE"), file_name = "Report/2019_1_Simulation/mse")
plot_recruitment( ceattle_list = list(ms_run_mse), line_col = c(1), model_names = c("MSE"), file_name = "Report/2019_1_Simulation/mse")
plot_selectivity( ceattle_list = list(ms_run_mse), line_col = c(1), model_names = c("MSE"), file_name = "Report/2019_1_Simulation/mse")


# MRE
ms_run_mre <- ms_run
ms_run_mre$quantities <- comparison$MRE

plot_biomass( ceattle_list = list(ms_run_mre), line_col = c(1), model_names = c("MRE"), file_name = "Report/2019_1_Simulation/mre")
plot_recruitment( ceattle_list = list(ms_run_mre), line_col = c(1), model_names = c("MRE"), file_name = "Report/2019_1_Simulation/mre")
plot_selectivity( ceattle_list = list(ms_run_mre), line_col = c(1), model_names = c("MRE"), file_name = "Report/2019_1_Simulation/mre")



# CV
ms_run_cv <- ms_run
ms_run_cv$quantities <- comparison$CV

plot_biomass( ceattle_list = list(ms_run_cv), line_col = c(1), model_names = c("CV"), file_name = "Report/2019_1_Simulation/cv")
plot_recruitment( ceattle_list = list(ms_run_cv), line_col = c(1), model_names = c("CV"), file_name = "Report/2019_1_Simulation/cv")
plot_selectivity( ceattle_list = list(ms_run_cv), line_col = c(1), model_names = c("CV"), file_name = "Report/2019_1_Simulation/cv")



