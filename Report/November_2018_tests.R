###################################################
# Run in single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = "ceattle.par", debug = FALSE, random_rec = FALSE, niter = 3)
ss_rep <- ss_run$quantities


###################################################
# Run in MS mode using par files
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_3_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 10)
ms_rep_3 <- ms_run$quantities

ms_run$diagnostics


load("data/BS_MS_3_Loops_Files_Corrected/CEATTLE_results.Rdata")
ms_tmp_3 <- tmp
tmb_like <- sum(ms_rep$jnll_comp)

ms_res_3 <- compare_output(ms_rep_3, ms_tmp_3, rel_error = 0.01)
ms_res_3$summary
ms_res_3$likelihood
tmb_nll_comp <- ms_res_3$tmb_like
admb_nll_comp <- ms_res_3$admb_like
tmb_nll_comp
admb_nll_comp
