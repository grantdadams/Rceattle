###################################################
# Run in single species mode
load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = "ceattle.par", debug = FALSE, random_rec = FALSE, niter = 3)
ss_rep <- ss_run$quantities

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
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 3)
ms_rep_3 <- ms_run$quantities

# Load previous estimates
load("data/BS_MS_Files/CEATTLE_results.Rdata")
ms_tmp_3 <- tmp

# Compare with current
ms_res_3 <- compare_output(ms_rep_3, ms_tmp_3, rel_error = 0.01)
ms_res_3$summary
ms_res_3$likelihood
tmb_nll_comp <- ms_res_3$tmb_like
admb_nll_comp <- ms_res_3$admb_like

tmb_pollock_age1_year1_M2 <- round(ms_rep_3$M2[1,1:12,1], 8)
admb_pollock_age1_year1_M2 <- ms_tmp_3$M2_1[1,]

tmb_jnll <- sum(ms_rep_3$jnll_comp)
admb_jnll <- -1854293.01910412


tmb_nll_comp
admb_nll_comp

tmb_jnll
admb_jnll

tmb_pollock_age1_year1_M2
admb_pollock_age1_year1_M2


###################################################
# Run in MS mode using par files and 5 loops
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_5_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 5)
ms_rep_5 <- ms_run$quantities

# Load previous estimates
load("data/BS_MS_5_Loops_Files/CEATTLE_results.Rdata")
ms_tmp_5 <- tmp

# Compare with current
ms_res_5 <- compare_output(ms_rep_5, ms_tmp_5, rel_error = 0.01)
ms_res_5$summary
ms_res_5$likelihood
tmb_nll_comp <- ms_res_5$tmb_like
admb_nll_comp <- ms_res_5$admb_like

tmb_jnll <- sum(ms_rep_5$jnll_comp)
admb_jnll <- -1854292.98931228


tmb_nll_comp
admb_nll_comp

tmb_jnll
admb_jnll

ms_rep_5$M2[1,1:12,1]
ms_tmp_5$M2_1[1,]


###################################################
# Run in MS mode using par files and 10 loops
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run_10 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_10_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 10)
ms_rep_10 <- ms_run_10$quantities

# Load previous estimates
load("data/BS_MS_10_Loops_Files/CEATTLE_results.Rdata")
ms_tmp_10 <- tmp

# Compare with current
ms_res_10 <- compare_output(ms_rep_10, ms_tmp_10, rel_error = 0.01)
ms_res_10$summary
ms_res_10$likelihood
tmb_nll_comp <- ms_res_10$tmb_like
admb_nll_comp <- ms_res_10$admb_like

tmb_jnll <- sum(ms_rep_10$jnll_comp)
admb_jnll <- -1854292.98946876


tmb_nll_comp
admb_nll_comp

tmb_jnll
admb_jnll

ms_rep_10$M2[1,1:12,1]
ms_tmp_10$M2_1[1,]


###################################################
# Run in MS mode using par files and 20 loops
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run_20 <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_20_Loops_Files/dat files/", inits = "ceattle.par", debug = TRUE, random_rec = FALSE, niter = 20)
ms_rep_20 <- ms_run_20$quantities

# Load previous estimates
load("data/BS_MS_20_Loops_Files/CEATTLE_results.Rdata")
ms_tmp_20 <- tmp

# Compare with current
ms_res_20 <- compare_output(ms_rep_20, ms_tmp_20, rel_error = 0.001)
ms_res_20$summary
ms_res_20$likelihood
tmb_nll_comp <- ms_res_20$tmb_like
admb_nll_comp <- ms_res_20$admb_like

tmb_jnll <- sum(ms_rep_20$jnll_comp)
admb_jnll <- -1854292.98946876


tmb_nll_comp
admb_nll_comp

tmb_jnll
admb_jnll

ms_rep_20$M2[1,1:12,1]
ms_tmp_20$M2_1[1,]
