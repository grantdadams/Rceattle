# This function extracts data from CEATTLE inputs used by K Holsman's ADMB CEATTLE Model and outputs a list of inputs for the TMB CEATTLE model


to_tmb <- function(){

  # Change wd
  setwd("data")

  # Specify filenames of the documents we are going to use
  datafile_name       = "CEATTLE_results_0.Rdata";

  #---------------------------------------------------------------------
  # 1. Get assessment data
  #---------------------------------------------------------------------
  load(datafile_name)
  assmnt_data <- tmp

  # 1.1 Extract data wanted
  assmntdat_names <- c("nages", "nyrs_tc_biom", "yrs_tc_biom", "tcb_obs", "nyrs_fsh_comp", "yrs_fsh_comp", "fsh_age_type", "fsh_age_bins", "obs_catch", "nyrs_wt_at_age", "yrs_wt_at_age", "wt", "pmature", "M1_base", "nyrs_srv_biom", "yrs_srv_biom", "srv_bio", "nyrs_srv_age", "srv_age_type", "yrs_srv_age", "srv_age_bins", "srv_age_n", "srv_age_obs", "srv_age_sizes", "age_trans_matrix", "n_eit", "yrs_eit", "obs_eit", "obs_eit_age", "nyrs_eit_sel", "yrs_eit_sel", "eit_sel")
  VBGFdat_names <- c("mf_type", "propMorF", "t0", "log_mean_d", "logK", "logH", "Tcoef", "Pcoef")
  stomdat_names <- c

  assmntdat_objno <- grep(paste(assmntdat_names,collapse="|"), names(assmnt_data))

  assmnt_data_extract <- list()
  for(i in 1:length(assmntdat_objno)){
    assmnt_data_extract[[i]] <- assmnt_data[[assmntdat_objno[i]]]
    names(assmnt_data_extract)[i] <- names(assmnt_data)[assmntdat_objno[i]]
  }

  # Combine objects into arrays or lists
  gsub( paste("_", 1:3, sep = ""), "", names(assmnt_data_extract))

  # Change wd back to main
  setwd("../")

  # Return the data to be used in TMB
  return(assmnt_data_extract)
}
