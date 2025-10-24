#' Run a management strategy evaluation
#'
#' @description Runs a forward projecting MSE. Main assumptions are the projected selectivity/catchability, foraging days, and weight-at-age are the same as the terminal year of the hindcast in the operating model. Assumes survey sd is same as average across historic time series, while comp data sample size is same as last year. No implementation error and no observation error for catch!
#'
#' @param om CEATTLE model object exported from \code{\link{Rceattle}}
#' @param em CEATTLE model object exported from \code{\link{Rceattle}}
#' @param nsim Number of simulations to run (default 10)
#' @param start_sim First simulation number to start at. Useful if the code stops at specific seed/sim (default = 1).
#' @param assessment_period Period of years that each assessment is taken
#' @param sampling_period Period of years data sampling is conducted. Single value or vector the same length as the number of fleets.
#' @param simulate_data Include simulated random error proportional to that estimated/provided for the data from the OM.
#' @param regenerate_past Refits the EM to historical/conditioning data prior to the MSE, where the data are generated from the OM with \code{simulate_data = TRUE} or without \code{simulate_data = FALSE} sampling error.
#' @param sample_rec Include resampled recruitment deviates from the"hindcast" in the projection of the OM. Resampled deviates are used rather than sampling from N(0, sigmaR) because initial deviates bias R0 low. If false, uses mean of recruitment deviates.
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier \code{mean rec * (1 + (rec_trend/projection years) * 1:projection years)}. Can be of length 1 or of length nspp. If length 1, all species get the same trend.
#' @param fut_sample future sampling effort relative to last year.  \code{ Log_sd * 1 / fut_sample} for index and \code{ Sample_size * fut_sample} for comps
#' @param cap A cap on the catch in the projection. Can be a single number applied to all species (proportional to recommended catch) or vector of length \code{nspp} applied to each species. Default = NULL
#' @param catch_mult A multiplier for the catch in the projection. Can be a single number or vector of length nspp. Default = NULL
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param file (Optional) Filename where each OM simulation with EMs will be saved. If NULL, no files are saved.
#' @param dir (Optional) Directory where each OM simulation is saved
#' @param seed seed for the simulation
#' @param regenerate_seed seed for regenerating data
#' @param timeout length of time (minutes) estimation will run before stopping a sim (default 999 minutes)
#'
#' @return A list of operating models (differ by simulated recruitment determined by \code{nsim}) and estimation models fit to each operating model (differ by terminal year).
#' @export
#'
#'
run_mse <- function(om = ms_run, em = ss_run, nsim = 10, start_sim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, regenerate_past = FALSE, sample_rec = TRUE, rec_trend = 0, fut_sample = 1, cap = NULL, catch_mult = NULL, seed = 666, regenerate_seed = seed, loopnum = 1, file = NULL, dir = NULL, timeout = 999, endyr = NA){

  # om = om; em = em; nsim = 1; start_sim = 1; assessment_period = 1; sampling_period = 1; simulate_data = TRUE; regenerate_past = FALSE; sample_rec = FALSE; rec_trend = 0; fut_sample = 1; cap = NULL; catch_mult = NULL; seed = 666; regenerate_seed = seed; loopnum = 1; file = NULL; dir = NULL; endyr = NA; timeout = 999

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # MSE SETUP ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  '%!in%' <- function(x,y)!('%in%'(x,y))
  library(dplyr)
  set.seed(regenerate_seed)

  Rceattle_OM_list <- list()
  Rceattle_EM_list <- list()

  # - Set om to project from R0
  om$data_list$proj_mean_rec = 0 # - Sample rec devs assuming this down the line
  #FIXME: SB0 for equilibrium HCRs will have to be adjusted

  # - Adjust cap
  if(!is.null(cap)){
    if(!length(cap) %in% c(1, om$data_list$nspp)){
      stop("cap is not length 1 or length nspp")
    }
  }

  if(!is.null(catch_mult)){
    if(length(catch_mult) == 1){
      catch_mult = rep(catch_mult, om$data_list$nspp)
    }

    if(length(catch_mult) != om$data_list$nspp){
      stop("catch_mult is not length 1 or length nspp")
    }
  }

  if(sum(om$data_list$fleet_control$proj_F_prop) == 0){
    stop("F prop per fllet 'proj_F_prop' is zero")
  }

  # - Years for simulations
  hind_yrs <- (em$data_list$styr) : em$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  om_proj_yrs <- (om$data_list$endyr + 1) : om$data_list$projyr
  om_proj_nyrs <- length(om_proj_yrs)

  em_proj_yrs <- (em$data_list$endyr + 1) : em$data_list$projyr
  em_proj_nyrs <- length(em_proj_yrs)
  nflts = nrow(om$data_list$fleet_control)

  # - N sel ages for sel coff dev
  if(all(is.na(om$data_list$fleet_control$Nselages))){
    nselages_om = dim(om$estimated_params$sel_coff_dev)[3]
  } else {
    nselages_om <- max(om$data_list$fleet_control$Nselages, na.rm = TRUE)
  }

  if(all(is.na(em$data_list$fleet_control$Nselages))){
    nselages_em = dim(em$estimated_params$sel_coff_dev)[3]
  } else {
    nselages_em <- max(em$data_list$fleet_control$Nselages, na.rm = TRUE)
  }

  # - Assessment period
  assess_yrs <- seq(from = om$data_list$endyr + assessment_period, to =  min(c(om$data_list$projyr, em$data_list$projyr, endyr), na.rm = TRUE),  by = assessment_period)

  # - Data sampling period
  if(length(sampling_period)==1){
    sampling_period = rep(sampling_period, nflts)

  }

  if(nflts != nrow(em$data_list$fleet_control)){
    stop("OM and EM fleets do not match or sampling period length is mispecified")
  }
  if(nflts != length(sampling_period)){
    stop("Sampling period length is mispecified, does not match number of fleets")
  }

  # - Set up years of data we are sampling for each fleet
  sample_yrs <- lapply(sampling_period, function(x) seq(from = em$data_list$endyr + x, to = em$data_list$projyr,  by = x))
  fleet_id <- sample_yrs
  for(i in 1:length(sample_yrs)){
    fleet_id[[i]] <- replace(fleet_id[[i]], values = i)
  }
  sample_yrs = data.frame(Fleet_code = unlist(fleet_id), Year = unlist(sample_yrs))


  # * Filter arbitrary "future" data ----
  # -- index_data
  om$data_list$index_data <- om$data_list$index_data %>%
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$index_data <- em$data_list$index_data %>%
    dplyr::filter(abs(Year) <= em$data_list$endyr)

  # -- comp_data
  om$data_list$comp_data <- om$data_list$comp_data %>%
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$comp_data <- em$data_list$comp_data %>%
    dplyr::filter(abs(Year) <= em$data_list$endyr)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Regenerate past data from OM and refit EM ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if(regenerate_past){

    # - Simulate index and comp data and updatae EM
    sim_dat <- sim_mod(om, simulate = FALSE)

    em$data_list$index_data <- sim_dat$index_data
    em$data_list$comp_data <- sim_dat$comp_data

    # Restimate
    em <- fit_mod(
      data_list = em$data_list,
      inits = em$estimated_params,
      map =  NULL,
      bounds = NULL,
      file = NULL,
      estimateMode = ifelse(em$data_list$estimateMode < 3, 0, em$data_list$estimateMode), # Run hindcast and projection, otherwise debug
      HCR = build_hcr(HCR = em$data_list$HCR, # Tier3 HCR
                      DynamicHCR = em$data_list$DynamicHCR,
                      Ftarget = em$data_list$Ftarget,
                      Flimit = em$data_list$Flimit,
                      Ptarget = em$data_list$Ptarget,
                      Plimit = em$data_list$Plimit,
                      Alpha = em$data_list$Alpha,
                      Pstar = em$data_list$Pstar,
                      Sigma = em$data_list$Sigma,
                      Fmult = em$data_list$Fmult,
                      HCRorder = em$data_list$HCRorder
      ),
      recFun = build_srr(srr_fun = em$data_list$srr_fun,
                         srr_pred_fun  = em$data_list$srr_pred_fun ,
                         proj_mean_rec  = em$data_list$proj_mean_rec ,
                         srr_meanyr = em$data_list$srr_meanyr,
                         srr_hat_styr = em$data_list$srr_hat_styr,
                         srr_hat_endyr = em$data_list$srr_hat_endyr,
                         srr_est_mode  = em$data_list$srr_est_mode ,
                         srr_prior  = em$data_list$srr_prior,
                         srr_prior_sd   = em$data_list$srr_prior_sd,
                         Bmsy_lim = em$data_list$Bmsy_lim,
                         srr_indices = em$data_list$srr_indices),
      M1Fun =     build_M1(M1_model = em$data_list$M1_model,
                           M1_re = em$data_list$M1_re,
                           updateM1 = FALSE,
                           M1_use_prior = em$data_list$M1_use_prior,
                           M2_use_prior = em$data_list$M2_use_prior,
                           M_prior = em$data_list$M_prior,
                           M_prior_sd = em$data_list$M_prior_sd,
                           M1_indices = em$data_list$M1_indices),
      random_rec = em$data_list$random_rec,
      niter = em$data_list$niter,
      msmMode = em$data_list$msmMode,
      avgnMode = em$data_list$avgnMode,
      suitMode = em$data_list$suitMode,
      suit_styr = em$data_list$suit_styr,
      suit_endyr = em$data_list$suit_endyr,
      initMode = em$data_list$initMode,
      phase = NULL,
      loopnum = loopnum,
      getsd = FALSE,
      verbose = 0)

    # Update avg F given model fit to regenerated data
    if(em$data_list$HCR == 2){

      # - Get avg F
      avg_F <- exp(em$estimated_params$ln_F) # Average F from last 5 years
      avg_F <- rowMeans(avg_F[,(ncol(avg_F)-4) : ncol(avg_F)])
      avg_F <- data.frame(avg_F = avg_F, spp = em$data_list$fleet_control$Species)
      avg_F <- avg_F %>%
        group_by(spp) %>%
        summarise(avg_F = sum(avg_F)) %>%
        arrange(spp)

      # - Update model
      em <- Rceattle::fit_mod(data_list = em$data_list,
                              inits = em$estimated_params,
                              estimateMode = 2, # Don't estimate
                              HCR = build_hcr(HCR = 2, # Input F
                                              Ftarget = avg_F$avg_F,
                                              Ptarget = em$data_list$Ptarget,
                                              Plimit = em$data_list$Plimit
                              ),
                              recFun = build_srr(srr_fun = em$data_list$srr_fun,
                                                 srr_pred_fun  = em$data_list$srr_pred_fun ,
                                                 proj_mean_rec  = em$data_list$proj_mean_rec ,
                                                 srr_meanyr = em$data_list$srr_meanyr,
                                                 srr_hat_styr = em$data_list$srr_hat_styr,
                                                 srr_hat_endyr = em$data_list$srr_hat_endyr,
                                                 srr_est_mode  = em$data_list$srr_est_mode ,
                                                 srr_prior  = em$data_list$srr_prior,
                                                 srr_prior_sd   = em$data_list$srr_prior_sd,
                                                 Bmsy_lim = em$data_list$Bmsy_lim,
                                                 srr_indices = em$data_list$srr_indices),
                              M1Fun =     build_M1(M1_model = em$data_list$M1_model,
                                                   M1_re = em$data_list$M1_re,
                                                   updateM1 = FALSE,
                                                   M1_use_prior = em$data_list$M1_use_prior,
                                                   M2_use_prior = em$data_list$M2_use_prior,
                                                   M_prior = em$data_list$M_prior,
                                                   M_prior_sd = em$data_list$M_prior_sd,
                                                   M1_indices = em$data_list$M1_indices),
                              random_rec = em$data_list$random_rec,
                              niter = em$data_list$niter,
                              msmMode = em$data_list$msmMode,
                              avgnMode = em$data_list$avgnMode,
                              suitMode = em$data_list$suitMode,
                              suit_styr = em$data_list$suit_styr,
                              suit_endyr = em$data_list$suit_endyr,
                              initMode = em$data_list$initMode,
                              loopnum = loopnum,
                              getsd = FALSE,
                              verbose = 0)
    }
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Expand OM data-dim ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # -- index_data
  proj_srv <- om$data_list$index_data %>%
    dplyr::group_by(Fleet_code) %>%
    dplyr::slice(rep(n(),  om_proj_nyrs)) %>%
    dplyr::mutate(Year = -om_proj_yrs)
  proj_srv$Log_sd <- proj_srv$Log_sd * 1/fut_sample
  proj_srv$Observation <- NA
  om$data_list$index_data  <- rbind(om$data_list$index_data, proj_srv)
  om$data_list$index_data <- dplyr::arrange(om$data_list$index_data, Fleet_code, abs(Year))

  # -- Nbyage
  if(nrow(om$data_list$NByageFixed) > 0){
    proj_nbyage <- om$data_list$NByageFixed %>%
      dplyr::group_by(Species, Sex) %>%
      dplyr::slice(rep(n(),  om_proj_nyrs)) %>%
      dplyr::mutate(Year = om_proj_yrs)
    proj_nbyage <- proj_nbyage[which(om_proj_yrs %!in% om$data_list$NByageFixed$Year),] # Subset rows already forcasted
    om$data_list$NByageFixed  <- rbind(om$data_list$NByageFixed, proj_nbyage)
    om$data_list$NByageFixed <- dplyr::arrange(om$data_list$NByageFixed, Species, Year)
  }

  # -- comp_data
  proj_comp <- om$data_list$comp_data %>%
    dplyr::group_by(Fleet_code, Sex) %>%
    dplyr::slice(rep(n(),  om_proj_nyrs)) %>%
    dplyr::mutate(Year = -om_proj_yrs)
  proj_comp$Sample_size <- proj_comp$Sample_size * fut_sample # Adjust future sampling effort
  proj_comp <- proj_comp %>%
    dplyr::mutate_at(vars(matches("Comp_")), ~ 1)
  om$data_list$comp_data  <- rbind(om$data_list$comp_data, proj_comp)
  om$data_list$comp_data <- dplyr::arrange(om$data_list$comp_data, Fleet_code, abs(Year))

  # -- emp_sel - Use terminal year
  if(nrow(om$data_list$emp_sel) > 0){
    proj_emp_sel <- om$data_list$emp_sel %>%
      dplyr::group_by(Fleet_code, Sex) %>%
      dplyr::slice(rep(n(),  om_proj_nyrs)) %>%
      dplyr::mutate(Year = om_proj_yrs)
    om$data_list$emp_sel  <- rbind(om$data_list$emp_sel, proj_emp_sel)
    om$data_list$emp_sel <- dplyr::arrange(om$data_list$emp_sel, Fleet_code, Year)
  }

  # -- weight
  #FIXME ignores forecasted growth
  proj_wt <- om$data_list$weight %>%
    dplyr::group_by(Wt_index , Sex) %>%
    dplyr::slice(rep(n(),  om_proj_nyrs)) %>%
    dplyr::mutate(Year = om_proj_yrs)
  om$data_list$weight  <- rbind(om$data_list$weight, proj_wt)
  om$data_list$weight <- dplyr::arrange(om$data_list$weight, Wt_index, Year)

  # -- Pyrs
  if(nrow(om$data_list$Pyrs) > 0){
    proj_Pyrs <- om$data_list$Pyrs %>%
      dplyr::group_by(Species, Sex) %>%
      dplyr::slice(rep(n(),  om_proj_nyrs)) %>%
      dplyr::mutate(Year = om_proj_yrs)
    om$data_list$Pyrs  <- rbind(om$data_list$Pyrs, proj_Pyrs)
    om$data_list$Pyrs <- dplyr::arrange(om$data_list$Pyrs, Species, Year)
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Expand EM data-dim ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  #FIXME - assuming same as terminal year of hindcast
  # -- EM emp_sel - Use terminal year
  proj_emp_sel <- em$data_list$emp_sel %>%
    dplyr::group_by(Fleet_code, Sex) %>%
    dplyr::slice(rep(n(),  em_proj_nyrs)) %>%
    dplyr::mutate(Year = em_proj_yrs)
  em$data_list$emp_sel  <- rbind(em$data_list$emp_sel, proj_emp_sel)
  em$data_list$emp_sel <- dplyr::arrange(em$data_list$emp_sel, Fleet_code, Year)

  # -- EM weight
  proj_wt <- em$data_list$weight %>%
    dplyr::group_by(Wt_index , Sex) %>%
    dplyr::slice(rep(n(),  em_proj_nyrs)) %>%
    dplyr::mutate(Year = em_proj_yrs)
  em$data_list$weight  <- rbind(em$data_list$weight, proj_wt)
  em$data_list$weight <- dplyr::arrange(em$data_list$weight, Wt_index, Year)

  # -- EM Pyrs
  if(nrow(em$data_list$Pyrs) > 0){
    proj_Pyrs <- em$data_list$Pyrs %>%
      dplyr::group_by(Species, Sex) %>%
      dplyr::slice(rep(n(),  em_proj_nyrs)) %>%
      dplyr::mutate(Year = em_proj_yrs)
    em$data_list$Pyrs  <- rbind(em$data_list$Pyrs, proj_Pyrs)
    em$data_list$Pyrs <- dplyr::arrange(em$data_list$Pyrs, Species, Year)
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # DO THE MSE ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  ### Set up parallel processing
  library(foreach)
  library(doParallel)

  cores = detectCores() - 6
  registerDoParallel(cores)

  sim_list <- foreach(sim = start_sim:nsim) %dopar% {
    library(Rceattle)
    library(dplyr)

    set.seed(seed = seed + sim) # setting unique seed for each simulation
    kill_sim <- list(kill_sim = FALSE, failure = NA)

    # Set models objects
    sim_list <- list(EM = list())# , OM = list())
    sim_list$EM[[1]] <- em
    # sim_list$OM[[1]] <- om

    em_use <- em
    om_use <- om

    # Sample recruitment
    om_use <- Rceattle::sample_rec(om_use, sample_rec = sample_rec, update_model = FALSE, rec_trend = rec_trend)

    # Run through assessment years
    for(k in 1:length(assess_yrs)){

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 1. Get recommended catch from the EM-HCR ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      new_years <- om_proj_yrs[which(om_proj_yrs <= assess_yrs[k] & om_proj_yrs > om_use$data_list$endyr)]

      # - Get projected catch data from EM
      new_catch_data <- em_use$data_list$catch_data
      dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
      new_catch_data$Catch[dat_fill_ind] <- em_use$quantities$catch_hat[dat_fill_ind]

      # * Catch multiplier ----
      if(!is.null(catch_mult)){
        new_catch_data$Catch[dat_fill_ind] <- new_catch_data$Catch[dat_fill_ind] * catch_mult[new_catch_data$Species[dat_fill_ind]]
      }

      # * Apply cap ----
      if(!is.null(cap)){
        # Applied across species
        if(length(cap) == 1){
          new_catch_data$Catch[dat_fill_ind] <- ifelse(sum(new_catch_data$Catch[dat_fill_ind]) > cap,
                                                       cap * new_catch_data$Catch[dat_fill_ind]/sum(new_catch_data$Catch[dat_fill_ind]),
                                                       new_catch_data$Catch[dat_fill_ind]) # FIXME: does not work for assessments that don't occur annually
        } else { # Species-specific
          new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap[new_catch_data$Species[dat_fill_ind]], cap[new_catch_data$Species[dat_fill_ind]], new_catch_data$Catch[dat_fill_ind])
        }
      }

      # * Exploitable biomass limit ----
      # - If projected catch > exploitable biomass in OM, reduce to exploitable biomass
      exploitable_biomass_data <- om_use$data_list$catch_data
      exploitable_biomass_data$Catch[dat_fill_ind] <- om_use$quantities$max_catch_hat[dat_fill_ind]

      new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > exploitable_biomass_data$Catch[dat_fill_ind],
                                                   exploitable_biomass_data$Catch[dat_fill_ind],
                                                   new_catch_data$Catch[dat_fill_ind])

      new_catch_switch <- sum(new_catch_data$Catch[dat_fill_ind]) #Switch to turn off re-running OM if new catch = 0

      # - Update catch data in OM and EM
      om_use$data_list$catch_data <- new_catch_data
      em_use$data_list$catch_data <- new_catch_data


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 2. Update the OM ----
      # - Estimate Fdev and update dynamics
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # - Update endyr of OM
      nyrs_hind <- om_use$data_list$endyr - om_use$data_list$styr + 1
      om_use$data_list$endyr <- assess_yrs[k]

      # * Update parameters ----
      # -- ln_F
      om_use$estimated_params$ln_F <- cbind(om_use$estimated_params$ln_F, matrix(0, nrow= nrow(om_use$estimated_params$ln_F), ncol = length(new_years)))

      # -- M1_dev
      #FIXME - simulate
      # om_use$estimated_params$ln_M1_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- om_use$estimated_params$ln_M1_dev[,,,nyrs_hind]

      # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
      om_use$estimated_params$index_q_dev <- cbind(om_use$estimated_params$index_q_dev, matrix(om_use$estimated_params$index_q_dev[,ncol(om_use$estimated_params$index_q_dev)], nrow= nrow(om_use$estimated_params$index_q_dev), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      ln_sel_slp_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
      sel_coff_dev = array(0, dim = c(nflts, 2, nselages_om, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for non-parameteric

      ln_sel_slp_dev[,,,1:nyrs_hind] <- om_use$estimated_params$ln_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_inf_dev
      sel_coff_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_coff_dev

      ln_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- ln_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

      om_use$estimated_params$ln_sel_slp_dev <- ln_sel_slp_dev
      om_use$estimated_params$sel_inf_dev <- sel_inf_dev
      om_use$estimated_params$sel_coff_dev <- sel_coff_dev


      # * Update map ----
      # -(Only new parameter we are estimating in OM is the ln_F of the new years)
      om_use$map <- build_map(
        data_list = om_use$data_list,
        params = om_use$estimated_params,
        debug = TRUE,
        random_rec = om_use$data_list$random_rec)
      om_use$map$mapFactor$dummy <- as.factor(NA); om_use$map$mapList$dummy <- NA


      # -- Estimate terminal F for catch
      new_f_yrs <- (ncol(om_use$map$mapList$ln_F) - length(new_years) + 1) : ncol(om_use$map$mapList$ln_F) # - Years of new F
      f_fleets <- om_use$data_list$fleet_control$Fleet_code[which(om_use$data_list$fleet_control$Fleet_type == 1)] # Fleet rows for F
      om_use$map$mapList$ln_F[f_fleets,new_f_yrs] <- replace(om_use$map$mapList$ln_F[f_fleets,new_f_yrs], values = 1:length(om_use$map$mapList$ln_F[f_fleets,new_f_yrs]))

      # -- Map out Fdev for years with 0 catch to very low number
      zero_catch <- om_use$data_list$catch_data %>%
        dplyr::filter(Year <= om_use$data_list$endyr &
                        Catch == 0) %>%
        dplyr::mutate(Year = Year - om_use$data_list$styr + 1) %>%
        dplyr::select(Fleet_code, Year) %>%
        as.matrix()
      om_use$estimated_params$ln_F[zero_catch] <- -999
      om_use$map$mapList$ln_F[zero_catch] <- NA
      om_use$map$mapFactor$ln_F <- factor(om_use$map$mapList$ln_F)
      rm(zero_catch)

      # -- Set estimate mode
      estimate_mode_base <- om_use$data_list$estimateMode
      estimate_mode_use <- ifelse(
        new_catch_switch == 0, 3, # Run in debug mode if catch is 0 for all species
        ifelse(
          estimate_mode_base < 3, 1, # Estimate hindcast only if estimating
          estimate_mode_base)
      )

      if(new_catch_switch == 0){
        om_use$map = NULL
      }

      # * Fit OM with new catch data ----
      kill_sim <- tryCatch({
        R.utils::withTimeout({
          suppressWarnings(
            om_use <- fit_mod(
              data_list = om_use$data_list,
              inits = om_use$estimated_params,
              map = om_use$map,
              bounds = NULL,
              file = NULL,
              estimateMode = estimate_mode_use,
              random_rec = om_use$data_list$random_rec,
              niter = om_use$data_list$niter,
              msmMode = om_use$data_list$msmMode,
              avgnMode = om_use$data_list$avgnMode,
              suitMode = om_use$data_list$suitMode,
              initMode = om_use$data_list$initMode,
              suit_styr = om$data_list$suit_styr,     # This stays the same as original OM to maintain constant suitability
              suit_endyr = om$data_list$suit_endyr,   # This stays the same as original OM to maintain constant suitability
              HCR = build_hcr(HCR = om_use$data_list$HCR,
                              DynamicHCR = om_use$data_list$DynamicHCR,
                              Ftarget = om_use$data_list$Ftarget,
                              Flimit = om_use$data_list$Flimit,
                              Ptarget = om_use$data_list$Ptarget,
                              Plimit = om_use$data_list$Plimit,
                              Alpha = om_use$data_list$Alpha,
                              Pstar = om_use$data_list$Pstar,
                              Sigma = om_use$data_list$Sigma,
                              Fmult = om_use$data_list$Fmult,
                              HCRorder = om_use$data_list$HCRorder
              ),
              recFun = build_srr(srr_fun = om_use$data_list$srr_fun,
                                 srr_pred_fun = om_use$data_list$srr_pred_fun ,
                                 proj_mean_rec = TRUE, # Use mean R for RPs
                                 srr_meanyr = om$data_list$srr_meanyr, # This stays the same as original OM
                                 srr_hat_styr = om$data_list$srr_hat_styr,
                                 srr_hat_endyr = om$data_list$srr_hat_endyr,
                                 srr_est_mode  = om_use$data_list$srr_est_mode ,
                                 srr_prior = om_use$data_list$srr_prior,
                                 srr_prior_sd = om_use$data_list$srr_prior_sd,
                                 Bmsy_lim = om_use$data_list$Bmsy_lim,
                                 srr_indices = om_use$data_list$srr_indices),
              M1Fun = build_M1(M1_model = om_use$data_list$M1_model,
                               M1_re = om_use$data_list$M1_re,
                               updateM1 = FALSE,  # Dont update M1 from data, fix at previous parameters
                               M1_use_prior = om_use$data_list$M1_use_prior,
                               M2_use_prior = om_use$data_list$M2_use_prior,
                               M_prior = om_use$data_list$M_prior,
                               M_prior_sd = om_use$data_list$M_prior_sd,
                               M1_indices = om_use$data_list$M1_indices),
              loopnum = loopnum,
              phase = FALSE,
              getsd = FALSE,
              verbose = 0)
          )
          return(list(kill_sim = FALSE, failure = NA))
        },
        timeout = 60*timeout)
      },
      error = function(e){
        return(list(kill_sim = TRUE, failure = "OM"))
      },
      TimeoutException = function(e){
        return(list(kill_sim = TRUE, failure = "OM"))
      })

      if(kill_sim$kill_sim){
        break()
      }

      # -- Set estimate mode back to original
      om_use$data_list$estimateMode <- estimate_mode_base


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 3. Get actual catch from OM ----
      # - Maybe the OM can't support the TAC
      # (should be minimal with the < exploitable biomass check)
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      # - Get realized catch data from OM
      new_catch_data <- om_use$data_list$catch_data
      dat_fill_ind <- which(new_catch_data$Year %in% new_years)
      new_catch_data$Catch[dat_fill_ind] <- om_use$quantities$catch_hat[dat_fill_ind] # Catch from OM

      # - Update catch data in OM and EM
      om_use$data_list$catch_data <- new_catch_data
      em_use$data_list$catch_data <- new_catch_data


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 4. Simulate data from OM ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # - Simulate new survey and comp data
      sim_dat <- Rceattle::sim_mod(om_use, simulate = simulate_data)

      years_include <- sample_yrs[which(sample_yrs$Year > em_use$data_list$endyr & sample_yrs$Year <= assess_yrs[k]),]

      # -- Add newly simulated survey data to EM and OM
      # - Get simulated survey data
      new_index_data <- sim_dat$index_data %>%
        dplyr::filter(abs(Year) %in% years_include$Year &
                        Fleet_code %in% years_include$Fleet_code) %>%
        dplyr::mutate(Year = -Year)

      # - Add to EM and OM
      om_use$data_list$index_data <- om_use$data_list$index_data %>%
        dplyr::filter(!(abs(Year) %in% years_include$Year &
                          Fleet_code %in% years_include$Fleet_code)) %>%
        rbind(new_index_data %>%
                dplyr::mutate(Year = -abs(Year))) %>%
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$index_data <- em_use$data_list$index_data %>%
        rbind(new_index_data) %>%
        dplyr::arrange(Fleet_code, abs(Year))


      # -- Add newly simulated comp data to EM & OM
      # - Simulated comp data
      new_comp_data <- sim_dat$comp_data %>%
        dplyr::filter(abs(Year) %in% years_include$Year &
                        Fleet_code %in% years_include$Fleet_code) %>%
        dplyr::mutate(Year = -Year)

      new_comp_data$Sample_size <- new_comp_data$Sample_size * as.numeric(rowSums(dplyr::select(new_comp_data, dplyr::contains("Comp_"))) > 0) # Set sample size to 0 if catch is 0
      new_comp_data <- new_comp_data %>%
        dplyr::mutate_at(dplyr::vars(dplyr::contains("Comp_")), ~ .x + 1 * (Sample_size == 0)) # Set all values to 1 if catch is 0

      # - Add to EM and OM
      om_use$data_list$comp_data <- om_use$data_list$comp_data %>%
        dplyr::filter(!(abs(Year) %in% years_include$Year &
                          Fleet_code %in% years_include$Fleet_code)) %>%
        rbind(new_comp_data %>%
                dplyr::mutate(Year = -abs(Year))) %>%
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$comp_data <- em_use$data_list$comp_data %>%
        rbind(new_comp_data) %>%
        dplyr::arrange(Fleet_code, abs(Year))

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 5. Update EM and HCR ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # Update end year and re-estimate
      em_use$data_list$endyr <- assess_yrs[k]

      # Update parameter size and use previous estimates
      # -- ln_F
      em_use$estimated_params$ln_F <- cbind(em_use$estimated_params$ln_F, matrix(0, nrow= nrow(em_use$estimated_params$ln_F), ncol = length(new_years)))

      # # -- ln_M1_dev
      # em_use$estimated_params$ln_M1_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- em_use$estimated_params$ln_M1_dev[,,,nyrs_hind]

      # -- Time-varying survey catachbilitiy - Assume last year - filled by columns
      em_use$estimated_params$index_q_dev <- cbind(em_use$estimated_params$index_q_dev, matrix(em_use$estimated_params$index_q_dev[,ncol(em_use$estimated_params$index_q_dev)], nrow= nrow(em_use$estimated_params$index_q_dev), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      ln_sel_slp_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
      sel_coff_dev = array(0, dim = c(nflts, 2, nselages_em, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for non-parameteric

      ln_sel_slp_dev[,,,1:nyrs_hind] <- em_use$estimated_params$ln_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_inf_dev
      sel_coff_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_coff_dev

      # - Initialize new years with last year
      ln_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- ln_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

      em_use$estimated_params$ln_sel_slp_dev <- ln_sel_slp_dev
      em_use$estimated_params$sel_inf_dev <- sel_inf_dev
      em_use$estimated_params$sel_coff_dev <- sel_coff_dev


      # Restimate
      kill_sim <- tryCatch({
        R.utils::withTimeout({
          suppressWarnings(
            em_use <- fit_mod(
              data_list = em_use$data_list,
              inits = em_use$estimated_params,
              map =  NULL,
              bounds = NULL,
              file = NULL,
              estimateMode = ifelse(em_use$data_list$estimateMode < 3, 0, em_use$data_list$estimateMode), # Run hindcast and projection, otherwise debug
              HCR = build_hcr(HCR = em_use$data_list$HCR, # Tier3 HCR
                              DynamicHCR = em_use$data_list$DynamicHCR,
                              Ftarget = em_use$data_list$Ftarget,
                              Flimit = em_use$data_list$Flimit,
                              Ptarget = em_use$data_list$Ptarget,
                              Plimit = em_use$data_list$Plimit,
                              Alpha = em_use$data_list$Alpha,
                              Pstar = em_use$data_list$Pstar,
                              Sigma = em_use$data_list$Sigma,
                              Fmult = em_use$data_list$Fmult,
                              HCRorder = em$data_list$HCRorder
              ),
              recFun = build_srr(srr_fun = em_use$data_list$srr_fun,
                                 srr_pred_fun = em_use$data_list$srr_pred_fun,
                                 proj_mean_rec = em_use$data_list$proj_mean_rec,
                                 srr_meanyr = em_use$data_list$endyr, # Update end year
                                 srr_hat_styr = em_use$data_list$srr_hat_styr,
                                 srr_hat_endyr = em_use$data_list$srr_hat_endyr,
                                 srr_est_mode  = em_use$data_list$srr_est_mode ,
                                 srr_prior = em_use$data_list$srr_prior,
                                 srr_prior_sd = em_use$data_list$srr_prior_sd,
                                 Bmsy_lim = em_use$data_list$Bmsy_lim,
                                 srr_indices = em_use$data_list$srr_indices),
              M1Fun =     build_M1(M1_model = em_use$data_list$M1_model,
                                   M1_re = em_use$data_list$M1_re,
                                   updateM1 = FALSE,
                                   M1_use_prior = em_use$data_list$M1_use_prior,
                                   M2_use_prior = em_use$data_list$M2_use_prior,
                                   M_prior = em_use$data_list$M_prior,
                                   M_prior_sd = em_use$data_list$M_prior_sd,
                                   M1_indices = em_use$data_list$M1_indices),
              random_rec = em_use$data_list$random_rec,
              niter = em_use$data_list$niter,
              msmMode = em_use$data_list$msmMode,
              avgnMode = em_use$data_list$avgnMode,
              suitMode = em_use$data_list$suitMode,
              suit_styr = em_use$data_list$suit_styr,
              suit_endyr = em_use$data_list$suit_endyr,
              initMode = em_use$data_list$initMode,
              phase = FALSE,
              loopnum = loopnum,
              getsd = FALSE,
              verbose = 0)
          )
          return(list(kill_sim = FALSE, failure = NA))
        },
        timeout = 60*timeout)
      },
      error = function(ex) {
        return(list(kill_sim = TRUE, failure = "EM"))
      },
      TimeoutException = function(ex) {
        return(list(kill_sim = TRUE, failure = "EM"))
      })

      if(is.null(em_use)){
        return(list(kill_sim = TRUE, failure = "EM"))
        break()
      }
      # plot_biomass(list(em_use, om_use), model_names = c("EM", "OM"))
      # End year of assessment

      # - Remove unneeded bits for memory
      em_use$initial_params <- NULL
      em_use$bounds <- NULL
      em_use$map <- NULL
      em_use$phase_params <- NULL
      em_use$obj <- NULL
      em_use$opt <- NULL
      em_use$sdrep <- NULL
      em_use$quantities[names(em_use$quantities) %!in% c("catch_hat",
                                                         "ln_catch_sd",
                                                         "index_hat",
                                                         "ln_index_sd",
                                                         "ssb_depletion",
                                                         "biomass_depletion",
                                                         "biomass",
                                                         "ssb",
                                                         "ssb_depletion",
                                                         "BO",
                                                         "SB0",
                                                         "SBF",
                                                         "F_spp",
                                                         "R",
                                                         "M1_at_age",
                                                         "M_at_age",
                                                         "avg_rec",
                                                         "DynamicB0",
                                                         "DynamicSB0",
                                                         "DynamicSBF",
                                                         "SPR0",
                                                         "SPRlimit",
                                                         "SPRtarget",
                                                         "Ftarget",
                                                         "B_eaten",
                                                         "B_eaten_as_prey",
                                                         "Flimit")] <- NULL

      sim_list$EM[[k+1]] <- em_use
      message(paste0("Sim ",sim, " - EM Year ", assess_yrs[k], " COMPLETE"))
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 6. End year loop ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    }


    # - Rename models
    sim_list$use_sim <- !kill_sim$kill_sim
    sim_list$failure = kill_sim$failure
    sim_list$OM <- om_use # OM
    sim_list$OM_no_F <- remove_F(om_use) # OM with no Fishing
    if(!kill_sim$kill_sim){
      names(sim_list$EM) <- c("EM", paste0("OM_Sim_",sim,". EM_yr_", assess_yrs))
    }

    # - Save
    if(!is.null(dir)){
      dir.create(file.path(getwd(), dir), showWarnings = FALSE, recursive = TRUE)
      saveRDS(sim_list, file = paste0(dir, "/", file, "EMs_from_OM_Sim_",sim, ".rds"))
      sim_list <- NULL
    } else{
      sim_list # Return simlist
    }

  } # End sim loop

  # When you're done, clean up the cluster
  stopImplicitCluster()
  names(sim_list) <- paste0("Sim_", start_sim:nsim)

  if(is.null(dir)){
    return(sim_list)
  }
}
