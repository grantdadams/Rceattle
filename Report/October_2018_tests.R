###################################################
# Run in single species mode
version <- "CEATTLE_BSAI_MS_v01"

load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, plot_trajectory = FALSE, random_rec = FALSE, niter = 3)
ss_rep <- ss_run$rep

# Load previous estimates
load("data/BS_SS_Files/CEATTLE_results.Rdata")
ss_tmp <- tmp
tmb_like <- sum(ss_rep$jnll_comp)

# Compare with current
ss_res <- compare_output(ss_rep, ss_tmp, rel_error = 0.0001)
ss_res$summary
ss_res$likelihood
tmb_like_comp <- ss_res$tmb_like
admb_like_cop <- ss_res$admb_like

tmb_like_comp
admb_like_cop
tmb_like
admb_like
ss_tmp$obj_fun


###################################################
# Run in MS mode using par files
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = "ceattle.par", debug = TRUE, plot_trajectory = FALSE, random_rec = FALSE, niter = 3)
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

sum(ms_rep$jnll_comp)
ms_tmp$obj_fun


###################################################
# Run in MS mode using std files
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = "ceattle_est.std", debug = FALSE, plot_trajectory = FALSE, random_rec = FALSE, niter = 3)
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

sum(ms_rep$jnll_comp)
ms_tmp$obj_fun

ms_rep$M2[1,,1]
ms_tmp$M2_1[1,]

###################################################
# Run in MS mode using par files and 10 loops
load("data/BS_MS_10_Loops_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_10_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, plot_trajectory = FALSE, random_rec = FALSE, niter = 10)
ms_rep <- ms_run$rep

# Load previous estimates
load("data/BS_MS_10_Loops_Files/CEATTLE_results.Rdata")
ms_tmp <- tmp

# Compare with current
ms_res <- compare_output(ms_rep, ms_tmp, rel_error = 0.08)
ms_res$summary
ms_res$likelihood
ms_res$tmb_like
ms_res$admb_like

sum(ms_rep$jnll_comp)
ms_tmp$obj_fun

ms_rep$M2[1,,1]
ms_tmp$M2_1[1,]
