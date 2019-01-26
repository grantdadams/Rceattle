#' Build data list object for TMB from ADMB dat files
#'
#' @description Function to build a \code{data_list} object to be used by Rceattle from ADMB based CEATTLE dat and ctl files for BSAI groundfish. The function first reads the TMB .cpp file speciied by \code{TMBfilename} to identify the names of the data objects used by TMB. The function then searches across the \code{.dat} files in the directory specified by \code{dat_dir} to find the location of each data object. The function then reads in the specific data object from the \code{.dat} file and reformats them to be used by TMB.
#'
#' @param ctlFilename The ADMB control (.ctl) file used for CEATTLE
#' @param TMBfilename The version of the cpp CEATTLE file found in the /src folder
#' @param cpp_directory The directory where the cpp file is found
#' @param dat_dir The directory where ctl and dat files are stored
#' @param nspp The number of species included in the CEATTLE model. Deafualts to 3.
#' @param nselages Number of ages to estimate selectivity. Can either be a single number or vector of length nspp
#'
#' @return A list of data objects used by TMB
#' @export
build_dat <- function(ctlFilename = NULL, TMBfilename = NULL, cpp_directory = NULL, dat_dir = NULL, nspp = 3, nselages = 8) {

  # Get cpp file if not provided
  if(is.null(TMBfilename) | is.null(cpp_directory)){
    cpp_directory <- system.file("executables",package="Rceattle")
    TMBfilename <- "ceattle_v01_02"
  } else{
    cpp_directory <- cpp_directory
    TMBfilename <- TMBfilename
  }

  #---------------------------------------------------------------------
  # Step 1 -- Extract data names used in TMB
  #---------------------------------------------------------------------
  cpp_fn <- file(paste(cpp_directory,"/", TMBfilename, ".cpp", sep = ""))
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

  #---------------------------------------------------------------------
  # Step 2 -- Find location of data in dat files
  #---------------------------------------------------------------------
  ctl_fn <- file(paste(dat_dir, ctlFilename, ".ctl", sep = ""))
  ctl_file <- readLines(ctl_fn)
  skipp <- grep("START filenames", ctl_file) # Line of data files
  nrow <- grep("END filenames", ctl_file) # Last line of data files
  ctl_file <- ctl_file[c(skipp:nrow)]
  dat_files <- ctl_file[-grep("#", ctl_file)]
  dat_files <- dat_files[grep(".dat", dat_files)]
  dat_files <- strsplit(dat_files, "/")
  dat_files <- sapply(dat_files, function(x) tail(x, n = 1))

  dat_loc <- data.frame(dat_name = dat_names, datfile = rep(NA, length(dat_names)))

  for (i in 1:length(dat_files)) {
    fn <- paste(dat_dir, dat_files[i], sep = "")
    if (file.exists(fn)) {
      dat_tmp <- scan(file = fn, what = character(), sep = "\n", quiet = T) # Get values from each line
      dat_tmp <- gsub(" ", "", dat_tmp) # Remove spacing
      dat_tmp <- strsplit(dat_tmp, split = c(";", ":")) # Substring at colon
      dat_tmp <- sapply(dat_tmp, head, n = 1) # Take first element.. usually name

      for (j in 1:length(dat_names)) {
        if (length(grep(paste0("\\b", "#", dat_names[j], "\\b"), dat_tmp)) > 0) {
          dat_loc$datfile[j] <- dat_files[i]
        }
      }
    }
    if (!file.exists(fn)) {
      print(paste("File", dat_files[i], "is not in the directory"))
    }
  }
  dat_loc
  dat_loc <- dat_loc[complete.cases(dat_loc), ]

  #---------------------------------------------------------------------
  # Step 3 -- Extract data from the dat files
  #---------------------------------------------------------------------
  dat_list <- list()
  for (i in 1:nrow(dat_loc)) {
    dat_list[[i]] <- Rceattle:::readdat(fn = paste(dat_dir, dat_loc[i, 2], sep = ""), nm = as.character(dat_loc[i, 1]), nspp = nspp)
    names(dat_list)[i] <- as.character(dat_loc[i, 1])
  }

  #---------------------------------------------------------------------
  # Step 4 -- Clean data remove columns of all NAs
  #---------------------------------------------------------------------
  for (i in 1:length(dat_list)) {
    dat_list[[i]] <- Rceattle:::dim_check(dat_list[[i]])
    dat_list[[i]] <- Rceattle:::remove_na_col(dat_list[[i]])
    dat_list[[i]] <- Rceattle:::list_to_array(dat_list[[i]])
  }

  # Print data included
  not_included <- dat_names[(!(dat_names %in% names(dat_list)))]
  if(length(not_included) > 0){
    print(paste("The following data inputs are not included:", paste(not_included, collapse = ", "), sep = " "))
  } else {
    print("All data inputs items are included.")
  }

  #---------------------------------------------------------------------
  # Steo 5 -- Model configuration
  #---------------------------------------------------------------------
  cpp_file <- readLines(cpp_fn)
  ctl_fn <- file(paste(dat_dir, ctlFilename, ".ctl", sep = ""))
  ctl_file <- readLines(ctl_fn)
  skipp <- grep("MODEL CONFIGURATION", cpp_file) # Line of data files
  nrow <- grep("MODEL INPUTS", cpp_file) # Last line of data files
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
  dat_names <- c(dat_names, "logist_sel_phase")

  for (i in 1:length(dat_names)) {
    skipp <- grep(dat_names[i], ctl_file) # Line of data files
    dat_list$config <- as.numeric(scan(ctl_fn, what = "", flush = F, blank.lines.skip = F, skip = skipp, nlines = 1, quiet = T, sep = ""))
    names(dat_list)[which(names(dat_list) == "config")] <- dat_names[i]
  }

  #---------------------------------------------------------------------
  # Steo 6 -- Convert to TMB configuration
  #---------------------------------------------------------------------
  # Convert selectivity
  for(i in 1:length(dat_list$logist_sel_phase)){
    if(dat_list$logist_sel_phase[i] < 0){
      dat_list$logist_sel_phase[i] <- 1
    }
    if(dat_list$logist_sel_phase[i] > 0){
      dat_list$logist_sel_phase[i] <- 0
    }
  }

  # Rename
  dat_list$srv_sel_type <- dat_list$logist_sel_phase
  dat_list$logist_sel_phase <- NULL


  if(length(nselages) == 1){
    dat_list$nselages <- rep(nselages, dat_list$nspp)
  }
  else if(length(nselages) == dat_list$nspp){
    dat_list$nselages <- nselages
  } else{
    stop("nselages is not of length 1 or nspp")
  }

  return(dat_list)
}
