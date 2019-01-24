#' Build parameter list from cpp file
#'
#' @description Function to read a TMB cpp file and construct parameter list object for Rceattle
#'
#' @param data_list A data_list object created by \code{\link{build_dat}}
#' @param inits Character vector of named initial values from ADMB \code{.std} or \code{.par} files or list of previous parameter estimates from Rceattle model.
#' @param TMBfilename The version of the cpp CEATTLE file.
#' @param cpp_directory The directory where the cpp file is found
#'
#' @return a list of map arguments for each parameter
#' @export
build_params <-
  function(data_list,
           inits = NULL,
           TMBfilename = NULL,
           cpp_directory = "src") {
    closeAllConnections()
    #---------------------------------------------------------------------
    # Step 1 -- Extract parameter names and dimensions used in TMB
    #---------------------------------------------------------------------
    cpp_fn <- file(paste(cpp_directory, "/", TMBfilename, ".cpp", sep = ""))

    cpp_file <- readLines(cpp_fn)
    cpp_file <-
      scan(
        cpp_fn,
        what = "character",
        flush = T,
        blank.lines.skip = F,
        quiet = T
      )
    cpp_tmp <- c()
    search_term <- c("PARAMETER\\(|PARAMETER_")
    param_lines <- grep(search_term, cpp_file, ignore.case = F)

    for (i in 1:length(param_lines)) {
      cpp_tmp[i] <-
        paste(
          scan(
            cpp_fn,
            skip = param_lines[i] - 1,
            flush = F,
            sep = "\t",
            nlines = 1,
            quiet = TRUE,
            what = "character",
            blank.lines.skip = TRUE
          ),
          sep = "",
          collapse = " "
        )
    }
    tt <-
      strsplit(cpp_tmp, split = c(" ")) # find all the text lines

    param_names <- c() # Character string of variables used in model
    param_dim <-
      data.frame(matrix(1, nrow = length(param_lines), ncol = 5))

    for (i in 1:length(param_lines)) {
      dat_line <- grep(search_term, tt[i][[1]], ignore.case = F)
      dat_call <-
        paste(tt[[i]][dat_line:(dat_line + 2)], collapse = "")
      param_names[i] <- sub("\\).*", "", sub(".*\\(", "", dat_call))

      dim_start <- grep("\\[", tt[[i]])
      dim_end <- grep("\\]", tt[[i]])
      dat_call <- paste(tt[[i]][dim_start:dim_end], collapse = "")
      dat_call <- sub("\\].*", "", sub(".*\\[", "", dat_call))
      dat_call <- unlist(strsplit(dat_call, ","))

      for (j in 1:length(dat_call)) {
        param_dim[i, j] <- dat_call[j]
      }
    }

    #---------------------------------------------------------------------
    # Step 2 -- Build list of parameters size of maximum data dimension
    #---------------------------------------------------------------------
    data_list$nspp2 = data_list$nspp + 1
    data_list$nspp_sq = data_list$nspp * data_list$nspp
    data_list$nspp_sq2 = data_list$nspp * (data_list$nspp + 1)
    param_list <- list()

    for (i in 1:nrow(param_dim)) {
      for (j in 1:ncol(param_dim)) {
        if (is.na(as.numeric(param_dim[i, j]))) {
          param_dim[i, j] <-
            max(data_list[[which(names(data_list) == param_dim[i, j])]], na.rm = T)
        }
      }
      param_dim[i,] <- as.numeric(param_dim[i,])

      if (length(which(param_dim[i,] > 1)) == 0) {
        # PARAMETER
        param = 0
      }

      # PARAMETER_VECTOR
      if (length(which(param_dim[i,] > 1)) == 1) {
        param = rep(0, param_dim[i,][which(param_dim[i,] > 1)])
      }

      # PARAMETER_MATRIX
      if (length(which(param_dim[i,] > 1)) == 2) {
        param = matrix(0,
                       nrow = as.numeric(param_dim[i, 1]),
                       ncol = as.numeric(param_dim[i, 2]))
      }

      # 3D PARAMETER_ARRAY
      if (length(which(param_dim[i,] > 1)) == 3) {
        param = array(0, dim = as.numeric(param_dim[i, 1:3]))
      }

      # 4D PARAMETER_ARRAY
      if (length(which(param_dim[i,] > 1)) == 4) {
        param = array(0, dim = as.numeric(param_dim[i, 1:4]))
      }

      param_list[[i]] <- param
      names(param_list)[i] <- param_names[i]
    }


    #---------------------------------------------------------------------
    # Step 3 -- Replace inits with starting values in range
    #---------------------------------------------------------------------
    param_list$log_eit_q <- -6.7025
    param_list$ln_mn_rec <-
      replace(param_list$ln_mn_rec, values = 9)
    param_list$ln_mean_F <-
      replace(param_list$ln_mean_F, values = -.8)
    param_list$log_gam_a <-
      replace(param_list$log_gam_a, values = 0.5)
    param_list$log_gam_b <-
      replace(param_list$log_gam_b, values = -0.5)


    param_list$logH_1 <-
      replace(param_list$logH_1, values = -8.5)
    param_list$logH_1b <-
      replace(param_list$logH_1b, values = 0)
    param_list$logH_1a <-
      replace(param_list$logH_1a, values = -3)

    param_list$logH_2 <-
      replace(param_list$logH_2, values = -9)
    param_list$logH_3 <-
      replace(param_list$logH_3, values = -9)
    param_list$H_4 <-
      replace(param_list$H_4, values = 1)

    # remove last init dev
    param_list$init_dev <- param_list$init_dev[,1:(ncol(param_list$init_dev)-1)]


    #---------------------------------------------------------------------
    # Step 4 -- Replace inits with previous parameters if desired
    #---------------------------------------------------------------------
    if (!is.null(inits)) {

      # If using std file
      if(grepl(".std", inits)){
        # std fild
        std_dat <- read.delim(inits, sep = "")
        std_dat$name <- gsub('[[:digit:]]|\\[|\\]', '', std_dat$name)

        param_names <- names(param_list)

        for (i in 1:length(param_list)) {

          #STD
          if (param_names[i] %in% unique(std_dat$name)) {
            # PARAMETER_VECTOR and PARAMETER
            if (length(which(param_dim[i,] > 1)) <= 1 | length(which(param_dim[i,] == 1)) == ncol(param_dim)) {
              param_list[[param_names[i]]] <-
                replace(param_list[[param_names[i]]],
                        values = std_dat$value[which(std_dat$name == param_names[i])])
            }
            if (length(which(param_dim[i,] > 1)) == 2) {

              # Init devs because odd age distribution
              if (param_names[i] == "init_dev") {
                init_dev <- std_dat$value[which(std_dat$name == param_names[i])]
                init_dev_lines <- c()

                for (j in 1:nrow(param_list[[param_names[i]]])) {
                  init_dev_lines <- c(init_dev_lines, rep(j, data_list$nages[j] - 1))
                  param_list[[param_names[i]]][j, 1:(data_list$nages[j] - 1)] <-
                    replace(param_list[[param_names[i]]][j, 1:(data_list$nages[j] - 1)], values = init_dev[which(init_dev_lines == j)])
                }
              } else{
                param_list[[param_names[i]]] <-
                  matrix(
                    std_dat$value[which(std_dat$name == param_names[i])],
                    byrow = T,
                    ncol = ncol(param_list[[param_names[i]]]),
                    nrow = nrow(param_list[[param_names[i]]])
                  )
              }
            }
          }
        }
      }

      ################################################################################################
      # If using par file
      if(grepl(".par", inits)){
        # std fild
        std_dat <-
          scan(
            inits,
            what = "",
            flush = T,
            blank.lines.skip = F,
            quiet = T
          )

        par_file_names <- c()
        search_term <- c("#")
        par_file_lines <- grep(search_term, std_dat, ignore.case = F)
        par_file_lines <- par_file_lines[-1] # Subtract objective function

        # Get parameter names
        for (i in 1:length(par_file_lines)) {
          par_file_names[i] <-
            paste(
              scan(
                inits,
                skip = par_file_lines[i] - 1,
                flush = F,
                sep = "\t",
                nlines = 1,
                quiet = TRUE,
                what = "character",
                blank.lines.skip = TRUE
              ),
              sep = "",
              collapse = " "
            )
        }

        par_file_names <- gsub("# ", "", par_file_names)
        par_file_names <- gsub(":", "", par_file_names)
        par_file_names_filtered <-  gsub("[0-9]", "", par_file_names)
        par_file_names_filtered <-  gsub("[[:punct:\\_]]", "", par_file_names_filtered)

        # Extract values
        par_list <- list()
        for (i in 1:length(par_file_lines)) {
          par_list[[i]] <-
            scan(inits,what="numeric",flush=F,blank.lines.skip=F,skip=par_file_lines[i],nlines=ifelse(i < length(par_file_lines), par_file_lines[i+1] - par_file_lines[i] - 1, 1), quiet=T,sep="")
          par_list[[i]] <- as.numeric(as.character(par_list[[i]]))
          names(par_list)[i] <- par_file_names[i]
        }

        # param_names <- names(param_list)

        # REPLACE VALUES
        for (i in 1:length(param_list)) {
          if (param_names[i] %in% unique(par_file_names_filtered)) {

            # PARAMETER_VECTOR and PARAMETER
            if (length(which(param_dim[i,] > 1)) <= 1| length(which(param_dim[i,] == 1)) == ncol(param_dim)) {
              param_list[[param_names[i]]] <-
                replace(param_list[[param_names[i]]],
                        values = unlist(par_list[grep( param_names[i], par_file_names)])) #FIXME: this will break if the order of the saved parameters are off
            }

            # PARAMETER_MATRIX
            if (length(which(param_dim[i,] > 1)) == 2) {

              # Init devs because odd age distribution
              if (param_names[i] == "init_dev") {
                init_dev <- unlist(par_list[grep( param_names[i], par_file_names)])
                init_dev_lines <- c()

                for (j in 1:nrow(param_list[[param_names[i]]])) {
                  init_dev_lines <- c(init_dev_lines, rep(j, data_list$nages[j] - 1))
                  param_list[[param_names[i]]][j, 1:(data_list$nages[j] - 1)] <-
                    replace(param_list[[param_names[i]]][j, 1:(data_list$nages[j] - 1)], values = init_dev[which(init_dev_lines == j)])
                }
              } else{
                param_list[[param_names[i]]] <-
                  matrix(
                    unlist(par_list[grep( param_names[i], par_file_names)]),
                    byrow = T,
                    ncol = ncol(param_list[[param_names[i]]]),
                    nrow = nrow(param_list[[param_names[i]]])
                  )
              }
            }
          }
        }
      }
    }
    closeAllConnections()

    return(param_list)
  }
