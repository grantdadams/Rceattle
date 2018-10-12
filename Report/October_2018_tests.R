###################################################
# Run in single species mode
ss_run <- Rceattle(data_list = ss_run$data_list, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", debug = FALSE, inits = FALSE)
ss_rep <- ss_run$rep

# Load previous estimates
load("data/BS_SS_Files/CEATTLE_results.Rdata")
ss_tmp <- tmp

# Compare with current
ss_res <- compare_output(ss_rep, ss_tmp, rel_error = 0.001)
ss_res$summary
ss_res$likelihood
ss_res$tmb_like
ss_res$admb_like

sum(ss_rep$jnll_comp)
ss_tmp$obj_fun


###################################################
# Run in MS mode
ms_run <- Rceattle(data_list = NULL, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", debug = FALSE, inits = FALSE)
ms_rep <- ms_run$rep

# Load previous estimates
load("data/BS_MS_Files/CEATTLE_results.Rdata")
ms_tmp <- tmp

# Compare with current
ms_res <- compare_output(ms_rep, ms_tmp, rel_error = 0.08)
ms_res$summary
ms_res$likelihood
ms_res$tmb_like
ms_res$admb_like

ms_rep$jnll
ms_tmp$obj_fun

ms_rep$M2[1,,1]
ms_tmp$M2_1[1,]

age_comp_like_comparison( data_list = ms_run$data_list, species = 3, ms_rep, ms_tmp, ADMB_TMB = 2)
age_comp_like_comparison( data_list = ms_run$data_list, species = 3, ms_rep, ms_tmp, ADMB_TMB = 1)

