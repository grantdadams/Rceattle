# This function extracts data from CEATTLE inputs used by K Holsman's ADMB CEATTLE Model and outputs a list of inputs for the TMB CEATTLE model

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
  assmnt_line_no <- which(substr(input_data[,1], 0 ,2) %in% paste("#", c(LETTERS, letters), sep = ""))
  assmnt_line_break <- which(input_data[,1] %in% c("##__________________________________________________________", "#=============================================================================="))

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
    if(class(assmnt_data$new_data) == "matrix" & nrow(assmnt_data$new_data) == 3){
      assmnt_data$new_data <- condense(assmnt_data$new_data)
    }
    if(class(assmnt_data$new_data) == "character"){
      assmnt_data$new_data <- as.numeric(as.character(assmnt_data$new_data))
    }

    # 1.2.4 Convert to array and make numeric if 3 dimensional
    if(class(assmnt_data$new_data) == "matrix" & nrow(assmnt_data$new_data) > 3){
      new_dat_list <- list()
      dat_breaks <- c(which(is.na(assmnt_data$new_data[,1])), nrow(assmnt_data$new_data) + 1) # Get the row breaks for each species

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
      # Assign to list
      assmnt_data$new_data <- new_dat_list
    }

    # 1.2.5 Rename list object to model object
    dat_name <- as.matrix(input_data[row_ind_name,])[1,1]
    dat_name <- gsub("#", "", dat_name)
    names(assmnt_data)[which(names(assmnt_data) == "new_data" )] <- dat_name
  }











}
