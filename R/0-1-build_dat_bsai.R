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
#' @param proj_yr The year to project the populations with no fishing. Assumed to be 2050
#' @param stom_tau Stomach content sample size for likelihood. Assumed to be 20.
#' @param endyr The end year of the hindcast
#' @param proj_F F vector or single value for projections
#'
#' @return A list of data objects used by TMB
#' @export
build_dat <- function(ctlFilename = NULL, TMBfilename = NULL, cpp_directory = NULL, dat_dir = NULL, nspp = 3, nselages = 8, endyr = 2017, proj_yr = 2050, proj_F = 0, stom_tau = 20) {

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
  names_not_in_cpp <- c("nspp"
    , "nyrs_srv_biom", "yrs_srv_biom", "srv_biom", "srv_biom_se",
                        "srv_age_obs", "nyrs_srv_age", "yrs_srv_age", "srv_age_n",
                        "srv_age_type", "srv_age_bins",
                        "n_eit", "yrs_eit", "obs_eit", "eit_sel",
                        "eit_age_n", "obs_eit_age",
                        "nyrs_fsh_comp" , "yrs_fsh_comp", "fsh_age_type" , "fsh_age_bins", "obs_catch", # FSH Comp
                        "nyrs_tc_biom" , "yrs_tc_biom", "tcb_obs", # Fish biom
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
  dat_list$projyr <- proj_yr
  dat_list$stom_tau <- rep(stom_tau, dat_list$nspp)


  if(length(nselages) == 1){
    dat_list$nselages <- rep(nselages, dat_list$nspp)
  } else if(length(nselages) == dat_list$nspp){
    dat_list$nselages <- nselages
  } else{
    stop("nselages is not of length 1 or nspp")
  }

  dat_list$pop_wt_index <- c(1:3)
  dat_list$ssb_wt_index <- c(1:3)
  dat_list$spawn_month <- rep(0, 3)
  dat_list$pop_alk_index <- c(1:3)
  dat_list$nlengths <- dat_list$fsh_age_bins
  dat_list$endyr <- endyr
  dat_list$nsex <- rep(1, 3)
  dat_list$R_sexr <- rep(NA, 3)

  #---------------------------------------------------------------------
  # Step 7 -- Survey specifications
  #---------------------------------------------------------------------
  # dat_list$n_srv <- c(2,1,1) # Bottom trawl and EIT

  srv_control <- data.frame(
    Fleet_name = c("BT_Pollock", "BT_Cod", "BT_ATF", "EIT_Pollock"),
    Fleet_code = c(4:7),
    Fleet_type = rep(2,4),
    Species = c(1:3, 1),
    Selectivity_index = c(4:7),
    Selectivity = c(dat_list$srv_sel_type, 0),
    Nselages = rep(NA, 4),
    Time_varying_sel = rep(0, 4),
    Sel_sd_prior = rep(0, 4),
    Age_first_selected = rep(1, 4),
    Accumatation_age_lower = rep(1, 4),
    Accumatation_age_upper = c(10, 12, 21, 10),
    Weight1_Numbers2 = rep(1, 4),
    Weight_index = c(1:3,1),
    ALK_index = c(1:3,1),
    Q_index = c(1:4),
    Estimate_q = c(0, 0, 0, 1),
    Log_q_prior = c(0, 0, 0, -6.7025),
    Time_varying_q = rep(0, 4),
    Q_sd_prior = rep(0, 4),
    Estimate_survey_sd = rep(0, 4), # Used to be Estimate_sigma_survey
    Survey_sd_prior = rep(NA, 4), # Used to be Sigma_survey_prior
    Estimate_catch_sd = rep(NA, 4), # Used to be Estimate_sigma_catch
    Catch_sd_prior = rep(NA, 4) # Used to be Sigma_catch_prior
  )


  #---------------------------------------------------------------------
  # Step 11 -- Fishery specifications
  #---------------------------------------------------------------------
  # dat_list$n_fsh <- c(1,1,1) # Bottom trawl and EIT
  dat_list$fleet_control <- data.frame(
    Fleet_name = c("Pollock", "Cod", "ATF"),
    Fleet_code = c(1:3),
    Fleet_type = rep(1,3),
    Species = c(1:3),
    Selectivity_index = 1:3,
    Selectivity = rep(2, 3),
    Nselages = c(dat_list$nselages),
    Time_varying_sel = rep(0, 3),
    Sel_sd_prior = rep(0, 3),
    Age_first_selected = rep(1, 3),
    Accumatation_age_lower = rep(1, 3),
    Accumatation_age_upper = c(10, 12, 21),
    Weight1_Numbers2 = rep(1, 3),
    Weight_index = c(1:3),
    ALK_index = c(1:3),
    Q_index = rep(NA, 3),
    Estimate_q = rep(NA, 3),
    Log_q_prior = rep(NA, 3),
    Time_varying_q = rep(NA, 3),
    Q_sd_prior = rep(NA, 3),
    Estimate_survey_sd = rep(NA, 3),
    Survey_sd_prior = rep(NA, 3),
    Estimate_catch_sd = rep(0, 3), # Used to be Estimate_sigma_catch
    Catch_sd_prior = rep(NA, 3) # Used to be Sigma_catch_prior
  )

  dat_list$fleet_control <- rbind(dat_list$fleet_control, srv_control)

  # Projected fishing mortality
  if(length(proj_F) != nrow(dat_list$fleet_control)){
    proj_F <- rep(proj_F[1], nrow(dat_list$fleet_control))
  }

  dat_list$fleet_control$proj_F <- proj_F

  #---------------------------------------------------------------------
  # Step 8 -- Reorganize for survey biomass
  #---------------------------------------------------------------------
  # BT BIOMASS
  dat_list$srv_biom <- data.frame(
    Fleet_name = rep(c("BT_Pollock", "BT_Cod", "BT_ATF"), dat_list$nyrs_srv_biom),
    Fleet_code = rep(4:6, dat_list$nyrs_srv_biom),
    Species = rep(1:nspp, dat_list$nyrs_srv_biom),
    Year = as.vector(t(dat_list$yrs_srv_biom)),
    Month = rep(rep(6, nspp), dat_list$nyrs_srv_biom),
    Selectivity_block = rep(rep(1, nspp), dat_list$nyrs_srv_biom),
    Q_block = rep(rep(1, nspp), dat_list$nyrs_srv_biom),
    Observation = as.vector(t(dat_list$srv_biom)),
    CV = as.vector(t(dat_list$srv_biom_se))
  )

  # SE to CV
  dat_list$srv_biom$CV <- (dat_list$srv_biom$CV/dat_list$srv_biom$Observation)
  dat_list$srv_biom$CV <- sqrt(log((dat_list$srv_biom$CV^2) + 1))

  # Add EIT bit
  eit_biom <- data.frame(
    Fleet_name = rep("EIT_Pollock", dat_list$n_eit),
    Fleet_code = rep(7, dat_list$n_eit),
    Species = rep(1, dat_list$n_eit),
    Year = dat_list$yrs_eit,
    Month = rep(6, dat_list$n_eit),
    Selectivity_block = rep( 1, dat_list$n_eit),
    Q_block = rep( 1, dat_list$n_eit),
    Observation = dat_list$obs_eit,
    CV = rep( 0.2, dat_list$n_eit)
  )
  dat_list$srv_biom <- rbind(dat_list$srv_biom, eit_biom)

  #---------------------------------------------------------------------
  # Step 9 -- Reorganize for survey age comp
  #---------------------------------------------------------------------
  # BT Survey age comp
  srv_comp <- data.frame(
    Fleet_name = rep(c("BT_Pollock", "BT_Cod", "BT_ATF"), dat_list$nyrs_srv_age),
    Fleet_code = rep(4:6, dat_list$nyrs_srv_age),
    Species = rep(1:nspp, dat_list$nyrs_srv_age),
    Sex = rep(rep(0, nspp), dat_list$nyrs_srv_age),
    Age0_Length1 = rep(c(0,1,1), dat_list$nyrs_srv_age),
    Year = as.vector(t(dat_list$yrs_srv_age)),
    Month = rep(rep(6, nspp), dat_list$nyrs_srv_age),
    Sample_size = as.vector(t(dat_list$srv_age_n))
  )


  Observation <- rbind(dat_list$srv_age_obs[,,1], dat_list$srv_age_obs[,,2], dat_list$srv_age_obs[,,3])
  colnames(Observation) <- paste0("Comp_", 1:ncol(Observation))
  srv_comp <- cbind(srv_comp, Observation)

  # EIT Survey ag e comp
  eit_comp <- data.frame(
    Fleet_name = rep("EIT_Pollock", dat_list$n_eit),
    Fleet_code = rep(7, dat_list$n_eit),
    Species = rep(1, dat_list$n_eit),
    Sex = rep(0, dat_list$n_eit),
    Age0_Length1 = rep(0, dat_list$n_eit),
    Year = dat_list$yrs_eit,
    Month = rep(6, dat_list$n_eit),
    Sample_size = dat_list$eit_age_n
  )

  colnames(dat_list$obs_eit_age) <- paste0("Comp_", 1:ncol(dat_list$obs_eit_age))
  eit_comp <- cbind(eit_comp, dat_list$obs_eit_age)

  # Combine
  srv_comp <- plyr::rbind.fill(srv_comp, eit_comp)

  #---------------------------------------------------------------------
  # Step 10 -- Reorganize selectivity
  #---------------------------------------------------------------------
  dat_list$emp_sel <- data.frame(
    Fleet_name = rep("EIT_Pollock", dat_list$n_eit),
    Fleet_code = rep(7, length(dat_list$n_eit)),
    Species = rep(1, length(dat_list$n_eit)),
    Sex = rep(0, length(dat_list$n_eit)),
    Year = dat_list$yrs_eit
  )
  colnames(dat_list$eit_sel) <- paste("Comp_",1:ncol(dat_list$eit_sel), sep = "")
  dat_list$emp_sel <- cbind(dat_list$emp_sel, dat_list$eit_sel[dat_list$yrs_eit - dat_list$styr + 1,])


  #---------------------------------------------------------------------
  # Step 11 -- Reorganize age transition matrix
  #---------------------------------------------------------------------
  alk_tmp <- dat_list$age_trans_matrix
  age_trans_matrix <- data.frame(ALK_name = c(rep("Pollock", 10), rep("Cod", 12), rep("ATF", 21)),
                                 ALK_index = c(rep(1, 10), rep(2, 12), rep(3, 21)),
                                 Species = c(rep(1, 10), rep(2, 12), rep(3, 21)),
                                 Sex = c(rep(0, 10), rep(0, 12), rep(0, 21)),
                                 Age = c(1:10, 1:12, 1:21))
  age_trans_obs <- rbind(alk_tmp[1:10,,1], alk_tmp[1:12,,2])
  age_trans_obs <-rbind(age_trans_obs, alk_tmp[1:21,,3])
  colnames(age_trans_obs) <- paste0("Length", 1:ncol(age_trans_obs))
  age_trans_matrix <- cbind(age_trans_matrix, age_trans_obs)
  dat_list$age_trans_matrix <- age_trans_matrix

  #---------------------------------------------------------------------
  # Step 12 -- Reorganize for fishery biomass
  #---------------------------------------------------------------------
  # BT BIOMASS
  dat_list$fsh_biom <- data.frame(
    Fleet_name = rep(c("Pollock", "Cod", "ATF"), dat_list$nyrs_tc_biom),
    Fleet_code = rep(1:3, dat_list$nyrs_tc_biom),
    Species = rep(1:nspp, dat_list$nyrs_tc_biom),
    Year = as.vector(t(dat_list$yrs_tc_biom)),
    Month = rep(rep(0, nspp), dat_list$nyrs_tc_biom),
    Selectivity_block = rep(rep(1, nspp), dat_list$nyrs_tc_biom),
    Catch = as.vector(t(dat_list$tcb_obs)),
    CV = rep(rep(0.05, nspp), dat_list$nyrs_tc_biom)
    )

  #---------------------------------------------------------------------
  # Step 13 -- Reorganize for fishery age comp
  #---------------------------------------------------------------------
  # Fishery age comp
  fsh_comp <- data.frame(
    Fleet_name = rep(c("Pollock", "Cod", "ATF"), dat_list$nyrs_fsh_comp),
    Fleet_code = rep(1:nspp, dat_list$nyrs_fsh_comp),
    Species = rep(1:nspp, dat_list$nyrs_fsh_comp),
    Sex = rep(rep(0, nspp), dat_list$nyrs_fsh_comp),
    Age0_Length1 = rep(c(0,1,1), dat_list$nyrs_fsh_comp),
    Year = na.exclude(as.vector(t(dat_list$yrs_fsh_comp))),
    Month = rep(rep(0, nspp), dat_list$nyrs_fsh_comp),
    Sample_size = rep(rep(200, nspp), dat_list$nyrs_fsh_comp)
  )


  Observation <- rbind(dat_list$obs_catch[,,1], dat_list$obs_catch[,,2], dat_list$obs_catch[,,3])
  colnames(Observation) <- paste0("Comp_", 1:ncol(Observation))
  Observation <- Observation[rowSums(is.na(Observation)) != ncol(Observation), ]
  fsh_comp <- cbind(fsh_comp, Observation)

  dat_list$comp_data <- rbind(fsh_comp, srv_comp)




  #---------------------------------------------------------------------
  # Step 15 -- Rename stuff
  #---------------------------------------------------------------------

  dat_list$propF = dat_list$propMorF[c(1,2,4),]
  dat_list$BTempC <- dat_list$BTempC_retro

  colnames(dat_list$aLW) <- paste0("Species", 1:nspp)


  dimnames(dat_list$Uobs) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)))

  dimnames(dat_list$UobsWt) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)), paste0('Pred_ln', 1:max( dat_list$srv_age_bins)))

  dimnames(dat_list$UobsAge) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$nages)), paste0('Pred_ln', 1:max( dat_list$nages)))

  dimnames(dat_list$UobsWtAge) <- list(paste0('Pred', 1:3),paste0('Prey', 1:3), paste0('Pred_ln', 1:max( dat_list$nages)), paste0('Pred_ln', 1:max( dat_list$nages)))


  dat_list$avgnMode <- 1
  dat_list$avgnMode <- 0
  dat_list$random_rec <- FALSE
  dat_list$niter <- 10
  dat_list$suitMode <- 0
  dat_list$est_diet <- 0
  dat_list$msmMode <- 1
  dat_list$debug <- TRUE
  dat_list$minage <- rep(1, dat_list$nspp)
  dat_list$sigma_rec_prior <- rep(sqrt(0.5), dat_list$nspp)
  dat_list$spnames <- c("Pollock", "Cod", "Arrowtooth flounder")

  ###########################
  # Change WT format
  ###########################
  yrs <- dat_list$styr:dat_list$endyr
  wt_new <- data.frame(Wt_name = rep(c("Pollock", "Cod", "ATF"), each = length(yrs)), Wt_index = rep(1:3, each = length(yrs)), Species = rep(1:3, each = length(yrs)), Year = rep(yrs, 3))
  wt_new$Sex = 0

  wt_temp <- rbind(dat_list$wt[,,1], dat_list$wt[,,2])
  wt_temp <- rbind(wt_temp, dat_list$wt[,,3])
  wt_temp <- as.data.frame(wt_temp)
  colnames(wt_temp) <- paste0("Age", 1:ncol(wt_temp))
  wt_new <- cbind(wt_new, wt_temp)
  dat_list$wt <- wt_new

  ###########################
  # Change Mn_LatAge format
  ###########################
  Mn_LatAge_tmp <- data.frame(Species = 1:3, Sex = rep(0, 3))
  colnames(dat_list$Mn_LatAge) <- paste0("Age", 1:ncol(dat_list$Mn_LatAge))
  Mn_LatAge_tmp <- cbind(Mn_LatAge_tmp, dat_list$Mn_LatAge)
  dat_list$Mn_LatAge <- Mn_LatAge_tmp

  ###########################
  # Change M1 format
  ###########################
  M1_tmp <- data.frame(Species = 1:3, Sex = rep(0, 3))
  colnames(dat_list$M1_base) <- paste0("Age", 1:ncol(dat_list$M1_base))
  M1_tmp <- cbind(M1_tmp, dat_list$M1_base)
  dat_list$M1_base <- M1_tmp

  ###########################
  # Adding aging error
  ###########################
  dat_list$age_error <- array(0, dim = c(dat_list$nspp, max(dat_list$nages), max(dat_list$nages)))
  for(sp in 1:dat_list$nspp){
    dat_list$age_error[sp,,] <- diag(1, max(dat_list$nages), max(dat_list$nages))
  }

  dat_list$nlengths <- dat_list$fsh_age_bins



  ###########################
  # Change UobsWtAge and UobsAge format
  ###########################
  UobsAge <- matrix(NA, ncol = 9, nrow = length(dat_list$UobsAge))
  dims <- dim(dat_list$UobsAge)

  ind <- 1
  for (pred in 1:dims[1]) {
    for (prey in 1:dims[2]) {
      for (pred_a in 1:dat_list$nages[pred]) {
        for (prey_a in 1:dat_list$nages[prey]) {
          UobsAge[ind, 1] <- pred
          UobsAge[ind, 2] <- prey
          UobsAge[ind, 3] <- 0
          UobsAge[ind, 4] <- 0
          UobsAge[ind, 5] <- pred_a + dat_list$minage[pred] - 1
          UobsAge[ind, 6] <- prey_a + dat_list$minage[prey] - 1
          UobsAge[ind, 7] <- 0
          UobsAge[ind, 8] <- 20
          UobsAge[ind, 9] <- dat_list$UobsAge[pred, prey, pred_a, prey_a]
          ind = ind + 1
        }
      }
    }
  }
  colnames(UobsAge) <- c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age", "Year", "Sample_size", "Stomach_proportion_by_number")
  dat_list$UobsAge <- UobsAge[which(complete.cases(UobsAge)),]


  UobsWtAge <- matrix(NA, ncol = 9, nrow = length(dat_list$UobsWtAge))
  dims <- dim(dat_list$UobsWtAge)

  ind <- 1
  for (pred in 1:dims[1]) {
    for (prey in 1:dims[2]) {
      for (pred_a in 1:dat_list$nages[pred]) {
        for (prey_a in 1:dat_list$nages[prey]) {
          UobsWtAge[ind, 1] <- pred
          UobsWtAge[ind, 2] <- prey
          UobsWtAge[ind, 3] <- 0
          UobsWtAge[ind, 4] <- 0
          UobsWtAge[ind, 5] <- pred_a + dat_list$minage[pred] - 1
          UobsWtAge[ind, 6] <- prey_a + dat_list$minage[prey] - 1
          UobsWtAge[ind, 7] <- 0
          UobsWtAge[ind, 8] <- 20
          UobsWtAge[ind, 9] <- dat_list$UobsWtAge[pred, prey, pred_a, prey_a]
          ind = ind + 1
        }
      }
    }
  }
  colnames(UobsWtAge) <- c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age",  "Year", "Sample_size", "Stomach_proportion_by_weight")
  dat_list$UobsWtAge <- UobsWtAge[which(complete.cases(UobsWtAge)),]



  #---------------------------------------------------------------------
  # Final Step -- Remove unwanted bits
  #---------------------------------------------------------------------
  names_not_in_cpp <- c(names_not_in_cpp, "BTempC_retro", "propMorF" , "srv_sel_type", "srv_age_type", "fsh_age_bins", "fsh_comp", "srv_comp")

  '%!in%' <- function(x,y)!('%in%'(x,y))
  dat_list2 <- list()

  names_in_cpp <- c(names_in_cpp, "comp_data",
                    "emp_sel", "fleet_control",
                    "fsh_comp", "srv_comp",
                    "fsh_biom", "srv_biom", "proj_F", "minage", "sigma_rec_prior", "spnames", "nsex", "R_sexr", "ssb_wt_index", "spawn_month")

  for(i in 1:length(names_in_cpp)){
    dat_list2[[names_in_cpp[i]]] <-  dat_list[[names_in_cpp[i]]]
  }

  return(dat_list2)
}
