
# Note: TC-Hat and Fsh_age_hat will be off, but like is good


ms_run <- Rceattle(data_list = data_list, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", debug = FALSE)
rep <- ms_run$rep
data_list <-ms_run$data_list

ss_run <- Rceattle(data_list = ss_run$data_list, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", debug = FALSE)
rep <- ss_run$rep

res <- compare_output(rep, tmp)

res[[2]]



opt$objective
rep$jnll
tmp$obj_fun
