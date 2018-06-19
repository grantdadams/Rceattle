build_params <- function(data_list, nselages = 8) {
  data_list$nselages <- nselages

  #---------------------------------------------------------------------
  # Step 1 -- Extract parameter names and dimensions used in TMB
  #---------------------------------------------------------------------
  cpp_fn <- file(paste("src/", TMBfilename, ".cpp", sep = ""))

  cpp_file <- readLines(cpp_fn)
  nt <- length(cpp_file)
  cpp_tmp <- c()
  param_lines <- grep("PARAMETER_", cpp_file)

  for (i in 1:length(param_lines)) {
    cpp_tmp[i] <- paste(scan(cpp_fn, skip = data_lines[i] - 1, flush = F, sep = "\t", nlines = 1, quiet = TRUE, what = "character", blank.lines.skip = TRUE), sep = "", collapse = " ")
  }
  tt <- strsplit(cpp_tmp, split = c(" ")) # find all the text lines

  param_names <- c() # Character string of variables used in model
  param_dim <- data.frame(matrix(1, nrow = length(param_lines), ncol = 4))

  for (i in 1:length(param_lines)) {
    dat_line <- grep("PARAMETER_", tt[i][[1]])
    dat_call <- paste(tt[[i]][ dat_line:(dat_line + 2)], collapse = "")
    param_names[i] <- sub("\\).*", "", sub(".*\\(", "", dat_call))

    dim_start <- grep("\\[", tt[[i]])
    dim_end <- grep("\\]", tt[[i]])
    dat_call <- paste(tt[[i]][ dim_start:dim_end], collapse = "")
    dat_call <- sub("\\].*", "", sub(".*\\[", "", dat_call))
    dat_call <- unlist(strsplit(dat_call, ","))

    for (j in 1:length(dat_call)) {
      param_dim[i, j] <- dat_call[j]
    }
  }

  #---------------------------------------------------------------------
  # Step 2 -- Build list of parameters size of maximum data dimension
  #---------------------------------------------------------------------
  param_list <- list()

  for(i in 1:nrow(param_dim)){
    for(j in 1:ncol(param_dim)){
      if(is.na(as.numeric(param_dim[i,j]))){
        param_dim[i,j] <- max(data_list[[which(names(data_list) == param_dim[i,j])]], na.rm = T)
      }
    }
    param_dim[i,] <- as.numeric(param_dim[i,])
    if(length(which(param_dim[i,] > 1)) == 1){ # PARAMETER_VECTOR
      param = rep(0, param_dim[i,][which(param_dim[i,] > 1)] )
    }

    if(length(which(param_dim[i,] > 1)) == 2){
      param = matrix(0, nrow = as.numeric(param_dim[i,1]), ncol = as.numeric(param_dim[i,2]) )
    }

    if(length(which(param_dim[i,] > 1)) == 3){
      param = array(0, dim = as.numeric(param_dim[i, 1:3]))
    }
    param_list[[i]] <- param
    names(param_list)[i] <- param_names[i]
  }

  return(param_list)
}
