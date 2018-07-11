build_dat <- function(ctlFilename, TMBfilename, dat_dir) {

  #---------------------------------------------------------------------
  # Step 1 -- Extract data names used in TMB
  #---------------------------------------------------------------------
  cpp_fn <- file(paste("src/", TMBfilename, ".cpp", sep = ""))

  cpp_file <- readLines(cpp_fn)
  skipp <- grep("MODEL INPUTS", cpp_file) # Line of data files
  nrow <- grep("PARAMETER SECTION", cpp_file) # Last line of data files
  cpp_file <- scan(cpp_fn, what = "character", skip = skipp, nlines = (nrow - skipp), flush = T, blank.lines.skip = F, quiet = T)
  nt <- length(cpp_file)
  cpp_tmp <- c()
  data_lines <- grep("DATA_", cpp_file)

  for (i in 1:length(data_lines)) {
    cpp_tmp[i] <- paste(scan(cpp_fn, skip = data_lines[i] + skipp - 1, flush = F, sep = "\t", nlines = 1, quiet = TRUE, what = "character", blank.lines.skip = TRUE), sep = "", collapse = " ")
  }
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
  source("R/Support Functions/readdat_fun.R")
  dat_list <- list()
  for (i in 1:nrow(dat_loc)) {
    dat_list[[i]] <- readdat(paste(dat_dir, dat_loc[i, 2], sep = ""), as.character(dat_loc[i, 1]), nspp = 3)
    names(dat_list)[i] <- as.character(dat_loc[i, 1])
  }

  #---------------------------------------------------------------------
  # Step 4 -- Clean data remove columns of all NAs
  #---------------------------------------------------------------------
  source("R/Support Functions/dim_check.R")
  for (i in 1:length(dat_list)) {
    dat_list[[i]] <- dim_check(dat_list[[i]])
    dat_list[[i]] <- remove_na_col(dat_list[[i]])
    dat_list[[i]] <- list_to_array(dat_list[[i]])
  }


  print(paste("The following items are not included:,", paste(dat_names[(!(dat_names %in% names(dat_list)))], collapse = ", "), sep = " "))


  #---------------------------------------------------------------------
  # Steo 5 -- Model configuration
  #---------------------------------------------------------------------
  cpp_file <- readLines(cpp_fn)
  ctl_fn <- file(paste(dat_dir, ctlFilename, ".ctl", sep = ""))
  ctl_file <- readLines(ctl_fn)
  skipp <- grep("MODEL CONFIGURATION", cpp_file) # Line of data files
  nrow <- grep("MODEL INPUTS", cpp_file) # Last line of data files
  cpp_file <- scan(cpp_fn, what = "character", skip = skipp, nlines = (nrow - skipp), flush = T, blank.lines.skip = F, quiet = T)
  nt <- length(cpp_file)
  cpp_tmp <- c()
  data_lines <- grep("DATA_", cpp_file)

  for (i in 1:length(data_lines)) {
    cpp_tmp[i] <- paste(scan(cpp_fn, skip = data_lines[i] + skipp - 1, flush = F, sep = "\t", nlines = 1, quiet = TRUE, what = "character", blank.lines.skip = TRUE), sep = "", collapse = " ")
  }
  tt <- strsplit(cpp_tmp, split = c(" ")) # find all the text lines

  dat_names <- c() # Character string of variables used in model
  for (i in 1:length(data_lines)) {
    dat_line <- grep("DATA_", tt[i][[1]])
    dat_call <- paste(tt[i][[1]][ dat_line:(dat_line + 2)], collapse = "")
    dat_names[i] <- sub("\\).*", "", sub(".*\\(", "", dat_call))
  }


  for (i in 1:length(dat_names)) {
    skipp <- grep(dat_names[i], ctl_file) # Line of data files
    dat_list$config <- as.numeric(scan(ctl_fn, what = "", flush = F, blank.lines.skip = F, skip = skipp, nlines = 1, quiet = T, sep = ""))
    names(dat_list)[which(names(dat_list) == "config")] <- dat_names[i]
  }


  return(dat_list)
}
