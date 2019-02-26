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

    data_list$nspp2 = data_list$nspp + 1
    data_list$nspp_sq = data_list$nspp * data_list$nspp
    data_list$nspp_sq2 = data_list$nspp * (data_list$nspp + 1)
    param_list <- list()


    #---------------------------------------------------------------------
    # Step 1 -- Specify parameter names and dimensions used in TMB
    #---------------------------------------------------------------------

    param_list$dummy = 0           # Variable to test derived quantities given input parameters; n = [1]

    # -- 3.1. Recruitment parameters
    param_list$ln_mn_rec = rep(0, data_list$nspp)    # Mean recruitment; n = [1, nspp]
    param_list$ln_rec_sigma = rep(0, data_list$nspp)  # Standard deviation of recruitment deviations; n = [1, nspp]
    param_list$rec_dev = matrix(0, nrow = data_list$nspp, ncol = data_list$nyrs)     # Annual recruitment deviation; n = [nspp, nyrs]


    # -- 3.2. Abundance parameters
    param_list$init_dev = matrix(0, nrow = data_list$nspp, ncol = max(data_list$nages))      # Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this

    # -- 3.3. fishing mortality parameters
    param_list$ln_mean_F = rep(0, data_list$nspp)   # Log mean fishing mortality; n = [1, nspp]
    param_list$F_dev = matrix(0, nrow = data_list$nspp, ncol = data_list$nyrs)     # Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change


    # -- 3.4. Selectivity parameters
    # FIXME - change nspp to n_srv
    # FIXME - change order of selectivity paramters for logistic
    param_list$srv_sel_coff = matrix(0, nrow = sum(data_list$n_srv), ncol = max(data_list$nselages))   # Survey selectivity parameters; n = [nspp, nselages]
    param_list$fsh_sel_coff = matrix(0, nrow = sum(data_list$n_srv), ncol = max(data_list$nselages))  # Fishery age selectivity coef; n = [nspp, nselages]
    param_list$srv_sel_slp = matrix(0, nrow = 2, ncol = sum(data_list$n_srv))  # Survey selectivity paramaters for logistic; n = [2, nspp]
    param_list$srv_sel_inf = matrix(0, nrow = 2, ncol = sum(data_list$n_srv))  # Survey selectivity paramaters for logistic; n = [2, nspp]
    param_list$log_srv_q = data_list$srv_control$log_q_prior   # Survey catchability; n = [sum(n_srv)]


    # -- 3.5. Kinzery predation function parameters
    param_list$logH_1 = matrix(0, nrow = data_list$nspp, ncol = data_list$nspp2)       # Predation functional form; n = [nspp, nspp2]; # FIXME: make matrix; nspp2 = nspp + 1
    param_list$logH_1a = rep(0, data_list$nspp)      # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
    param_list$logH_1b = rep(0, data_list$nspp)      # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix

    param_list$logH_2 = matrix(0, nrow = data_list$nspp, ncol = data_list$nspp)       # Predation functional form; n = [nspp, nspp]
    param_list$logH_3 = matrix(0, nrow = data_list$nspp, ncol = data_list$nspp)       # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
    param_list$H_4 = matrix(0, nrow = data_list$nspp, ncol = data_list$nspp)        # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;

    # 3.6 Gamma selectivity parameters
    param_list$log_gam_a = rep(0, data_list$nspp)    # Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
    param_list$log_gam_b = rep(0, data_list$nspp)   # Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

    #---------------------------------------------------------------------
    # Step 3 -- Replace inits with starting values in range
    #---------------------------------------------------------------------
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
        par_file_names_filtered <- gsub("[0-9]", "", par_file_names)
        par_file_names_filtered <- gsub("[[:punct:\\_]]", "", par_file_names_filtered)

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

    # Start phi at 0.5
    param_list$phi <- replace(param_list$phi, values = rep(log(0.5), length(param_list$phi)))

    closeAllConnections()

    return(param_list)
  }
