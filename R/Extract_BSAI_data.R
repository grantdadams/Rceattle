# This function extracts data from CEATTLE inputs used by K Holsman's ADMB CEATTLE Model and outputs a list of inputs for the TMB CEATTLE model
#
#
#
#
# Developed by Grant Adams 3-26-18

extract_data <- function(){

  # Load required functions
  source("R/Support Functions/condense_df_fun.R")

  # Change wd
  setwd("data")

  # Specify filenames of the documents we are going to use
  datafile_name       = "assmnt_nov2017_0.dat";
  diet2file_name      = "diet2.dat";
  stomfile_name       = "stomach2017.dat";
  ATFfile_name        = "ATF_Misc.dat";
  Fprofile_datafile   = "F_profile.dat";
  Recfile_name        = "fits_4_CEATTLE/rs_data4CEATTLE_0_0_0.dat";
  retrofile_name      = "retro_data2017_asssmnt.dat";
  futfile_name        = "projection_data2017_assssmnt.dat";
  catch_in_name       = "catch_in.dat";
  setC_name           = "set_catch.dat";
  setF_name           = "setFabcFofl.dat";



  #---------------------------------------------------------------------
  # 1. Get assessment data
  #---------------------------------------------------------------------
  assmnt_data <- list()

  # 1.1 Extract line numbers of inputs
  input_data <- read.delim(datafile_name, header =  T)
  assmnt_dat_names <- c("nages", "nyrs_tc_biom", "yrs_tc_biom", "tcb_obs", "nyrs_fsh_comp", "yrs_fsh_comp", "fsh_age_type", "fsh_age_bins", "obs_catch", "nyrs_wt_at_age", "yrs_wt_at_age", "wt", "pmature", "M1_base", "nyrs_srv_biom", "yrs_srv_biom", "srv_biom", "nyrs_srv_age", "srv_age_type", "yrs_srv_age", "srv_age_bins", "srv_age_n", "srv_age_obs", "srv_age_sizes", "age_trans_matrix", "n_eit", "yrs_eit", "obs_eit", "obs_eit_age", "nyrs_eit_sel", "yrs_eit_sel", "eit_sel")

  VBGF_dat_names <- c("mf_type", "propMorF", "t0", "log_mean_d", "logK", "logH", "Tcoef", "Pcoef")

  assmnt_line_no <- which(substr(input_data[,1], 0 ,2) %in% paste("#", assmnt_dat_names , sep = ""))
  assmnt_line_break <- which(input_data[,1] %in% c("##__________________________________________________________", "# ##__________________________________________________________"))

  # 1.1.2 Get names of the line numbers
  dat_name <- as.matrix(input_data[assmnt_line_break + 1,])[,1]
  dat_name <- gsub("#", "", dat_name)

  # 1.1.3 Remove extra line numbers for obs_catch, wt, srv_age_obs, eit_sel because they are formatted differently.
  # Von bert stuff is also funky

  # 1.2 Extract model inputs
  for(i in 1:length(assmnt_line_no)){
    # 1.2.1 Get row numbers for data
    row_ind_name  <- assmnt_line_no[i] # Get the row number for the object name
    row_ind_first <- assmnt_line_no[i] + 2 # Get the first row number where data are
    row_ind_break <- assmnt_line_break[which(assmnt_line_break > row_ind_first)][1] # Get the row number for the first break after the object
    row_ind_data <- row_ind_first:(row_ind_break-1) # Row numbers where data are

    # 1.2.2 Extract data and remove unwanted columns
    assmnt_data$new_data <- as.matrix(input_data[row_ind_data,])
    colnames(assmnt_data$new_data) = NULL
    assmnt_data$new_data[which(assmnt_data$new_data %in% c("#", ""))] <- NA # Set blank cells to NA
    # Remove columns and rows of only NAs
    assmnt_data$new_data <- assmnt_data$new_data[, colSums(is.na(assmnt_data$new_data)) < nrow(assmnt_data$new_data)]
    if(is.null(nrow(assmnt_data$new_data)) == F){
      assmnt_data$new_data <- assmnt_data$new_data[rowSums(is.na(assmnt_data$new_data)) < ncol(assmnt_data$new_data),]
    }

    # 1.2.3 Convert to numeric if 1 or 2 dimensional
    if(class(assmnt_data$new_data) == "matrix" & is.null(nrow(assmnt_data$new_data)) == F){
      if(nrow(assmnt_data$new_data) == 3){
        assmnt_data$new_data <- condense(assmnt_data$new_data)
      }
    }
    if(class(assmnt_data$new_data) == "character"){
      assmnt_data$new_data <- as.numeric(as.character(assmnt_data$new_data))
    }

    # 1.2.4 Convert to array and make numeric if 3 dimensional
    if(class(assmnt_data$new_data) == "matrix" & is.null(nrow(assmnt_data$new_data)) == F){
      if(nrow(assmnt_data$new_data) > 3){
        new_dat_list <- list()
        dat_breaks <- sort(c(grep("#", assmnt_data$new_data[,1]), nrow(assmnt_data$new_data) + 1)) # Get the row breaks for each species

        # 1.2.4.1 Assign each species' data into a list
        for(j in 1:(length(dat_breaks)-1)){
          dat_lines <- (dat_breaks[j]+1):(dat_breaks[j+1]-1)
          new_dat_list$new_data <- assmnt_data$new_data[dat_lines,]

          # 1.2.4.2 Condense matrix and make numeric if NAs
          new_dat_list$new_data <- condense(new_dat_list$new_data)

          # 1.2.4.3 Name for species
          dat_name <- as.character(paste(assmnt_data$new_data[dat_breaks[j],], collapse = ""))
          dat_name <- gsub("NA", "", dat_name)
          dat_name <- gsub("#", "", dat_name)
          names(new_dat_list)[which(names(new_dat_list) == "new_data" )] <- dat_name
        }
        # 1.2.4.4 Make into array and assign to list
        list_max_dim <- c(max(unlist(lapply(new_dat_list, nrow))), max(unlist(lapply(new_dat_list, ncol))), length(new_dat_list))
        new_dat_array <- array(unlist(new_dat_list), dim = list_max_dim) # Note: ragged arrays will be filled with NA's based on the max size of an object
        dimnames(new_dat_array)[[3]] <- names(new_dat_list)
        assmnt_data$new_data <- new_dat_array
      }
    }

    # 1.2.5 Rename list object to model object
    names(assmnt_data)[which(names(assmnt_data) == "new_data" )] <- dat_name[i]
  }
  # Change wd back to main
  setwd("../")

  # Return the data to be used in TMB
  return(assmnt_data)
}
