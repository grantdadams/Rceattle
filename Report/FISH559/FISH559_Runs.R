library(Rceattle)
###################################################
# ADMB  single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_admb <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3, file_name = "Report/FISH559/Models/ss_admb")

###################################################
# Run in single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run_no_re <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = FALSE, niter = 3, file_name = "Report/FISH559/Models/ss_no_re")

###################################################
# Run in single species mode with random recruitment
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run_re <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = TRUE, niter = 3, file_name = "Report/FISH559/Models/ss_re")


###################################################
# ADMB MS
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_admb_ssdat <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_10_Loops_Files_Corrected/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 10, file_name = "Report/FISH559/Models/ms_admb_10")


# Loop through number of iterations
load("data/BS_MS_Files/2017_assessment_data_list.RData")
for(i in 1:2){
  ###################################################
  # Run in MS mode using par files
  # ms_run_no_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_no_re$final_params, debug = FALSE, random_rec = FALSE, niter = i, file_name = paste0("Report/FISH559/Models/ms_no_re_5optim_", i))

  ###################################################
  # Run in MS mode using par files and random recruitment
  ms_run_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ms_run_no_re$final_params, debug = FALSE, random_rec = TRUE, niter = i, file_name = paste0("Report/FISH559/Models/ms_re_", i))
}

ms_run_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ms_run_no_re$final_params, debug = FALSE, random_rec = TRUE, niter = 10, file_name = paste0("Report/FISH559/Models/ms_re_", 10))



# Compare other food
# Loop through number of iterations
load("data/BS_MS_Files/2017_assessment_data_list.RData")
other_food <- data_list_ms$other_food
other_food_percent <- seq(0.05, 2, by = 0.05 )

for(i in 20:length(other_food_percent)){
  data_list_ms$other_food <- other_food * other_food_percent[i]

  ###################################################
  # Run in MS mode using par files
  ms_run_no_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_no_re$final_params, debug = FALSE, random_rec = FALSE, niter = 8, file_name = paste0("Report/FISH559/Models/Other Food/ms_no_re_5_optim_other_food_", i))

  # ###################################################
  # # Run in MS mode using par files and random recruitment
  # ms_run_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ms_run_no_re$final_params, debug = FALSE, random_rec = TRUE, niter = i, file_name = paste0("Report/FISH559/Models/ms_re_", i))
}


# Average N
ms_run_no_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_no_re$final_params, debug = FALSE, random_rec = FALSE, niter = 10, file_name = paste0("Report/FISH559/Models/ms_no_re_avgn", 10), AvgN_type = 0)

ms_run_no_re <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = ss_run_no_re$final_params, debug = FALSE, random_rec = FALSE, niter = 10, file_name = paste0("Report/FISH559/Models/ms_no_re_n", 10), AvgN_type = 1)
