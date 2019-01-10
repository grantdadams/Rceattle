#' Build data list object
#' @description Function to build a \code{data_list} object to be used by Rceattle for GOA groundfish.
#'
#' @param TMBfilename The version of the cpp CEATTLE file found in the src folder
#' @param dat_dir The directory where dat files are stored
#' @param nspp The number of species included in the CEATTLE model
build_dat_goa <- function(TMBfilename = NULL, dat_dir = NULL, nspp = 3) {

  #---------------------------------------------------------------------
  # Step 1 -- Extract data names used in TMB
  #---------------------------------------------------------------------
  cpp_fn <- file(paste("inst/", TMBfilename, ".cpp", sep = ""))
  cpp_file <- readLines(cpp_fn)

  skipp <- grep("MODEL INPUTS", cpp_file) # Line of data files
  nrow <- grep("PARAMETER SECTION", cpp_file) # Last line of data files
  cpp_file <- cpp_file[skipp:nrow]
  data_lines <- grep("DATA_", cpp_file)
  cpp_tmp <- cpp_file[data_lines]
  tt <- strsplit(cpp_tmp, split = c(" ")) # find all the text lines

  dat_names <- c() # Character string of variables used in model
  for (i in 1:length(data_lines)) {
    dat_line <- grep("DATA_", tt[i][[1]])
    dat_call <- paste(tt[i][[1]][ dat_line:(dat_line + 2)], collapse = "")
    dat_names[i] <- sub("\\).*", "", sub(".*\\(", "", dat_call))
  }




  return(dat_list)
}


ss_to_ceattle <- function( dat_file ){

}
