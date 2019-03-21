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
#' @param proj_yr The year to project the populations with no fishing. Assumed to be 2100
#' @param stom_tau Stomach content sample size for likelihood. Assumed to be 20.
#'
#' @return A list of data objects used by TMB
#' @export
build_dat <- function(ctlFilename = NULL, TMBfilename = NULL, cpp_directory = NULL, dat_dir = NULL, nspp = 3, nselages = 8, proj_yr = 2100, stom_tau = 20) {

  # Get cpp file if not provided
  if(is.null(TMBfilename) | is.null(cpp_directory)){
    cpp_directory <- system.file("executables",package="Rceattle")
    TMBfilename <- "ceattle_v01_04"
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
  # Step 1 -- Add data names NOT used in TMB
  #---------------------------------------------------------------------
  names_not_in_cpp <- c("nyrs_srv_biom", "yrs_srv_biom", "srv_biom_se",
                        "srv_age_obs", "nyrs_srv_age", "yrs_srv_age", "srv_age_n",
                        "srv_age_type", "srv_age_bins",
                        "n_eit", "yrs_eit", "obs_eit", "eit_sel",
                        "eit_age_n", "obs_eit_age",
                        "fsh_age_bins", "fsh_age_type", "nyrs_fsh_comp", "yrs_fsh_age",
                        "nyrs_fsh_comp" , "yrs_fsh_comp" , "obs_catch", "nyrs_tc_biom" , "tcb_obs", "nyrs_tc_biom", "yrs_tc_biom",
                        "propMorF", "mf_type", "BTempC_retro", "fsh_sel_type")
  names_in_cpp <- dat_names
  dat_names <- c(dat_names, names_not_in_cpp)

  dat_names <- as.character(unique(dat_names))

  #---------------------------------------------------------------------
  # Step 2 -- Find location of data in dat files
  #---------------------------------------------------------------------
  ctl_fn <- file(paste(dat_dir,"/", ctlFilename, ".ctl", sep = ""))
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
    fn <- paste(dat_dir, "/",dat_files[i], sep = "")
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
    dat_list[[i]] <- Rceattle:::readdat(fn = paste(dat_dir, "/", dat_loc[i, 2], sep = ""), nm = as.character(dat_loc[i, 1]), nspp = nspp)
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
    print(paste("The following data inputs are not included in the dat files:", paste(not_included, collapse = ", "), sep = " "))
  } else {
    print("All data inputs items are included.")
  }

  #---------------------------------------------------------------------
  # Steo 5 -- Model configuration
  #---------------------------------------------------------------------
  cpp_file <- readLines(cpp_fn)
  ctl_fn <- file(paste(dat_dir, "/",ctlFilename, ".ctl", sep = ""))
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
  names_in_cpp <- c(names_in_cpp, dat_names)
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
      dat_list$logist_sel_phase[i] <- 2 # Non-parametric
    }
    if(dat_list$logist_sel_phase[i] > 0){
      dat_list$logist_sel_phase[i] <- 1 # Logistic selectivity
    }
  }

  # Rename
  dat_list$srv_sel_type <- dat_list$logist_sel_phase
  dat_list$logist_sel_phase <- NULL

  # Get add in projected year
  # dat_list$proj_yr <- proj_yr
  dat_list$stom_tau <- stom_tau


  if(length(nselages) == 1){
    dat_list$nselages <- rep(nselages, dat_list$nspp)
  } else if(length(nselages) == dat_list$nspp){
    dat_list$nselages <- nselages
  } else{
    stop("nselages is not of length 1 or nspp")
  }

  #---------------------------------------------------------------------
  # Step 7 -- Survey specifications
  #---------------------------------------------------------------------
  # dat_list$n_srv <- c(2,1,1) # Bottom trawl and EIT

  dat_list$srv_control <- data.frame(
    Survey_name = c("BT_Pollock", "BT_Cod", "BT_ATF", "EIT_Pollock"),
    Survey_code = c(1:4),
    Species = c(1:3, 1),
    Selectivity = c(dat_list$srv_sel_type, 0),
    Nselages = rep(NA, 4),
    Comp_type = c(dat_list$srv_age_type, 1),
    Comp_N_bins = c(dat_list$srv_age_bins, ncol(dat_list$obs_eit_age)),
    # Comp_Nyrs = c(dat_list$nyrs_srv_age, dat_list$n_eit),
    Estimate_q = c(0, 0, 0, 1),
    log_q_start = c(0, 0, 0, -6.7025)
  )

  #---------------------------------------------------------------------
  # Step 8 -- Reorganize for survey biomass
  #---------------------------------------------------------------------
  # BT BIOMASS
  dat_list$srv_biom <- data.frame(
    Survey_name = rep(c("BT_Pollock", "BT_Cod", "BT_ATF"), dat_list$nyrs_srv_biom),
    Survey_code = rep(1:3, dat_list$nyrs_srv_biom),
    Species = rep(1:nspp, dat_list$nyrs_srv_biom),
    Year = as.vector(t(dat_list$yrs_srv_biom)),
    Month = rep(rep(6, nspp), dat_list$nyrs_srv_biom),
    Observation = as.vector(t(dat_list$srv_biom)),
    CV = as.vector(t(dat_list$srv_biom_se))
  )

  # SE to CV
  dat_list$srv_biom$CV <- (dat_list$srv_biom$CV/dat_list$srv_biom$Observation)
  dat_list$srv_biom$CV <- sqrt(log((dat_list$srv_biom$CV^2) + 1))

  # Add EIT bit
  eit_biom <- data.frame(
    Survey_name = rep("EIT_Pollock", dat_list$n_eit),
    Survey_code = rep(4, dat_list$n_eit),
    Species = rep(1, dat_list$n_eit),
    Year = dat_list$yrs_eit,
    Month = rep(0, dat_list$n_eit),
    Observation = dat_list$obs_eit,
    CV = rep( 0.2, dat_list$n_eit)
  )
  dat_list$srv_biom <- rbind(dat_list$srv_biom, eit_biom)

  #---------------------------------------------------------------------
  # Step 9 -- Reorganize for survey age comp
  #---------------------------------------------------------------------
  # BT Survey age comp
  dat_list$srv_comp <- data.frame(
    Survey_name = rep(c("BT_Pollock", "BT_Cod", "BT_ATF"), dat_list$nyrs_srv_age),
    Survey_code = rep(1:nspp, dat_list$nyrs_srv_age),
    Species = rep(1:nspp, dat_list$nyrs_srv_age),
    Year = as.vector(t(dat_list$yrs_srv_age)),
    Month = rep(rep(6, nspp), dat_list$nyrs_srv_age),
    Sample_size = as.vector(t(dat_list$srv_age_n))
  )


  Observation <- rbind(dat_list$srv_age_obs[,,1], dat_list$srv_age_obs[,,2], dat_list$srv_age_obs[,,3])
  colnames(Observation) <- paste0("Comp_", 1:ncol(Observation))
  dat_list$srv_comp <- cbind(dat_list$srv_comp, Observation)

  # EIT Survey ag e comp
  eit_comp <- data.frame(
    Survey_name = rep("EIT_Pollock", dat_list$n_eit),
    Survey_code = rep(4, dat_list$n_eit),
    Species = rep(1, dat_list$n_eit),
    Year = dat_list$yrs_eit,
    Month = rep(0, dat_list$n_eit),
    Sample_size = dat_list$eit_age_n
  )

  colnames(dat_list$obs_eit_age) <- paste0("Comp_", 1:ncol(dat_list$obs_eit_age))
  eit_comp <- cbind(eit_comp, dat_list$obs_eit_age)

  # Combine
  dat_list$srv_comp <- plyr::rbind.fill(dat_list$srv_comp, eit_comp)

  #---------------------------------------------------------------------
  # Step 10 -- Reorganize selectivity
  #---------------------------------------------------------------------
  dat_list$emp_srv_sel <- data.frame(
    Survey_name = rep("EIT_Pollock", dat_list$n_eit),
    Survey_code = rep(4, length(dat_list$n_eit)),
    Species = rep(1, length(dat_list$n_eit)),
    Year = dat_list$yrs_eit
  )
  colnames(dat_list$eit_sel) <- paste("Comp_",1:ncol(dat_list$eit_sel), sep = "")
  dat_list$emp_srv_sel <- cbind(dat_list$emp_srv_sel, dat_list$eit_sel[dat_list$yrs_eit - dat_list$styr + 1,])



  #---------------------------------------------------------------------
  # Step 11 -- Fishery specifications
  #---------------------------------------------------------------------
  # dat_list$n_fsh <- c(1,1,1) # Bottom trawl and EIT

  dat_list$fsh_control <- data.frame(
    Fishery_name = c("Pollock", "Cod", "ATF"),
    Fishery_code = c(1:3),
    Species = c(1:3),
    Selectivity = rep(2, 3),
    Nselages = c(dat_list$nselages),
    Comp_type = c(dat_list$fsh_age_type),
    Comp_N_bins = c(dat_list$fsh_age_bins)
  )

  #---------------------------------------------------------------------
  # Step 12 -- Reorganize for fishery biomass
  #---------------------------------------------------------------------
  # BT BIOMASS
  dat_list$catch_biom <- data.frame(
    Fishery_name = rep(c("Pollock", "Cod", "ATF"), dat_list$nyrs_tc_biom),
    Fishery_code = rep(1:3, dat_list$nyrs_tc_biom),
    Species = rep(1:nspp, dat_list$nyrs_tc_biom),
    Year = as.vector(t(dat_list$yrs_tc_biom)),
    Month = rep(rep(0, nspp), dat_list$nyrs_tc_biom),
    Catch_kg = as.vector(t(dat_list$tcb_obs)),
    CV = rep(rep(0.05, nspp), dat_list$nyrs_tc_biom)
    )

  #---------------------------------------------------------------------
  # Step 13 -- Reorganize for fishery age comp
  #---------------------------------------------------------------------
  # Fishery age comp
  dat_list$fsh_comp <- data.frame(
    Fishery_name = rep(c("Pollock", "Cod", "ATF"), dat_list$nyrs_fsh_comp),
    Fishery_code = rep(1:nspp, dat_list$nyrs_fsh_comp),
    Species = rep(1:nspp, dat_list$nyrs_fsh_comp),
    Year = na.exclude(as.vector(t(dat_list$yrs_fsh_comp))),
    Month = rep(rep(0, nspp), dat_list$nyrs_fsh_comp),
    Sample_size = rep(rep(200, nspp), dat_list$nyrs_fsh_comp)
  )


  Observation <- rbind(dat_list$obs_catch[,,1], dat_list$obs_catch[,,2], dat_list$obs_catch[,,3])
  colnames(Observation) <- paste0("Comp_", 1:ncol(Observation))
  Observation <- Observation[rowSums(is.na(Observation)) != ncol(Observation), ]
  dat_list$fsh_comp <- cbind(dat_list$fsh_comp, Observation)

  #---------------------------------------------------------------------
  # Step 14 -- Reorganize selectivity
  #---------------------------------------------------------------------
  dat_list$emp_fsh_sel <- data.frame(
    Fishery_name = NA,
    Fishery_code = NA,
    Species = NA,
    Year = NA,
    Comp_1 = NA
  )


  #---------------------------------------------------------------------
  # Step 15 -- Rename stuff
  #---------------------------------------------------------------------

  dat_list$endyr = 2100
  dat_list$propF = dat_list$propMorF[c(1,2,4),]
  dat_list$BTempC <- dat_list$BTempC_retro

  colnames(dat_list$aLW) <- paste0("Species", 1:nspp)


  dimnames(dat_list$Uobs) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)))

  dimnames(dat_list$UobsWt) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)))

  dimnames(dat_list$UobsAge) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$nages)), paste0('Pred_ln', 1:max( dat_list$nages)))

  dimnames(dat_list$UobsWtAge) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$nages)), paste0('Pred_ln', 1:max( dat_list$nages)))


  dat_list2$avgnMode <- 1
dat_list2$avgnMode <- 0
  dat_list2$random_rec <- FALSE
  dat_list2$niter <- 10
  dat_list2$suitMode <- 0
  dat_list2$est_diet <- 0
  dat_list2$msmMode <- 1
  dat_list2$debug <- TRUE


  #---------------------------------------------------------------------
  # Final Step -- Remove unwanted bits
  #---------------------------------------------------------------------
  names_not_in_cpp <- c(names_not_in_cpp, "BTempC_retro", "propMorF" , "srv_sel_type", "srv_age_type", "fsh_age_bins")

  '%!in%' <- function(x,y)!('%in%'(x,y))
  dat_list2 <- list()

  for(i in 1:length(names_in_cpp)){
    dat_list2[[names_in_cpp[i]]] <-  dat_list[[names_in_cpp[i]]]
  }

  return(dat_list2)
}
