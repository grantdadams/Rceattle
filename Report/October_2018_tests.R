###################################################
# Run in single species mode
version <- "CEATTLE_BSAI_MS_v01"

load("data/BS_SS_Files/2017_assessment_data_list.RData")
ss_run <- Rceattle(data_list = data_list_ss, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", debug = TRUE, inits = "ceattle.par", plot_trajectory = FALSE, random_rec = FALSE)
ss_rep <- ss_run$rep

# Load previous estimates
load("data/BS_SS_Files/CEATTLE_results.Rdata")
ss_tmp <- tmp

# Compare with current
ss_res <- compare_output(ss_rep, ss_tmp, rel_error = 0.00001)
ss_res$summary
ss_res$likelihood
ss_res$tmb_like
ss_res$admb_like

sum(ss_rep$jnll_comp)
ss_tmp$obj_fun


###################################################
# Run in MS mode using par files
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", debug = TRUE, inits = "ceattle.par")
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

# CHECK 1 -  Check N[1,1,nspp] - CLEAR
exp(ms_run$params$ln_mn_rec + ms_run$params$rec_dev[,1])
ms_rep$NByage[1,1,]
c(ms_tmp$NByage_1[1,1], ms_tmp$NByage_2[1,1], ms_tmp$NByage_3[1,1])

# CHECK 2 - Check N0 - CLEAR
# Species 1
exp(ms_run$params$ln_mn_rec[1] + ms_run$params$rec_dev[1,])
ms_rep$NByage[,1,1]
ms_tmp$NByage_1[,1]

# Species 2
exp(ms_run$params$ln_mn_rec[2] + ms_run$params$rec_dev[2,])
ms_rep$NByage[,1,2]
ms_tmp$NByage_2[,1]

# Species 3
exp(ms_run$params$ln_mn_rec[3] + ms_run$params$rec_dev[3,])
ms_rep$NByage[,1,3]
ms_tmp$NByage_3[,1]

# CHECK 3 - Check F
# Species 1
exp(ms_run$params$ln_mean_F[1] + ms_run$params$F_dev[1,])
ms_rep$F[,1,1]
ms_tmp$F_1[,1]

# Species 2
ms_rep$F[,1,2]
ms_tmp$F_2[,1]

# Species 3
ms_rep$F[,1,3]
ms_tmp$F_3[,1]

# CHECK 4 - Compare M2
# TMB avgN inputs
tmb_est <- compare_pred_function(ms_rep, ms_tmp, TMB = TRUE)
# Species 1
tmb_est$R_M2[,1:12,1]
tmb_est$TMB_M2[,1:12,1]
round(tmb_est$TMB_M2[,1:12,1],6) == round(tmb_est$R_M2[,1:12,1],6)

# Species 2
tmb_est$R_M2[,1:21,2]
tmb_est$TMB_M2[,1:21,2]
round(tmb_est$TMB_M2[,1:21,2],6) == round(tmb_est$R_M2[,1:21,2],6)

# Species 3
tmb_est$R_M2[,1:21,3]
tmb_est$TMB_M2[,1:21,3]
round(tmb_est$TMB_M2[,1:21,3],6) == round(tmb_est$R_M2[,1:21,3],6)


# ADMB avgN inputs
admb_est <- compare_pred_function(ms_rep, ms_tmp, TMB = FALSE)
rounder = 2
# Species 1
admb_est$R_M2[,1:12,1]
ms_tmp$M2_1
round(ms_tmp$M2_1, rounder) == round(admb_est$R_M2[,1:12,1], rounder)

# Species 2
admb_est$R_M2[,1:21,2]
ms_tmp$M2_2
round(ms_tmp$M2_2[,1:12], rounder) == round(admb_est$R_M2[,1:12,2], rounder)

# Species 3
admb_est$R_M2[,1:21,3]
ms_tmp$M2_3
round(ms_tmp$M2_3, rounder) == round(admb_est$R_M2[,1:21,3], rounder)

###################################################
# Run in MS mode using std files
load("data/BS_MS_Files/2017_assessment_data_list.RData")
ms_run <- Rceattle(data_list = data_list_ms, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", debug = TRUE, inits = "ceattle_est.std")
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
