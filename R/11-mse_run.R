#' Run a management strategy evaluation
#'
#' @description Runs a forward projecting MSE. Main assumptions are the projected selectivity/catchability, foraging days, and weight-at-age are the same as the terminal year of the hindcast in the operating model. Assumes survey sd is same as average across historic time series, while comp data sample size is same as last year. No implementation error!
#'
#' @param om CEATTLE model object exported from \code{\link{Rceattle}}
#' @param em CEATTLE model object exported from \code{\link{Rceattle}}
#' @param nsim Number of simulations to run (default 10)
#' @param assessment_period Period of years that each assessment is taken
#' @param sampling_period Period of years data sampling is conducted. Single value or vector the same length as the number of fleets.
#' @param simulate Include simulated random error proportional to that estimated/provided.
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier mean rec * (1 + (rec_trend/projection years) * 1:projection years)
#' @param cap A cap on the catch in the projection. Can be a single number or vector. Default = NULL
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param file (Optional) Filename where each OM simulation with EMs will be saved. If NULL, no files are saved.
#'
#' @return A list of operating models (differ by simulated recruitment determined by \code{nsim}) and estimation models fit to each operating model (differ by terminal year).
#' @export
#'
#' @examples
mse_run <- function(om = ms_run, em = ss_run, nsim = 10, assessment_period = 1, sampling_period = 1, simulate = TRUE, rec_trend = 0, cap = NULL, seed = 666, loopnum = 1, file = NULL){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  library(dplyr)
  set.seed(seed)

  Rceattle_OM_list <- list()
  Rceattle_EM_list <- list()

  # - Adjust cap
  if(!is.null(cap)){
    if(length(cap) == 1){
      cap = rep(cap, om$data_list$nspp)
    }

    if(length(cap) != om$data_list$nspp){
      stop("cap is not length 1 or length nspp")
    }
  }

  # - Years for simulations
  proj_yrs <- (em$data_list$endyr + 1) : em$data_list$projyr
  proj_nyrs <- length(proj_yrs)
  nflts = nrow(om$data_list$fleet_control)
  nselages_om <- max(om$data_list$fleet_control$Nselages, na.rm = TRUE)
  nselages_em <- max(em$data_list$fleet_control$Nselages, na.rm = TRUE)

  # - Assessment period
  assess_yrs <- seq(from = em$data_list$endyr + assessment_period, to = em$data_list$projyr,  by = assessment_period)

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

  # Update data-files in OM so we can fill in updated years
  # -- srv_biom
  proj_srv <- om$data_list$srv_biom %>%
    group_by(Fleet_code) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = -proj_yrs)
  om$data_list$srv_biom  <- rbind(om$data_list$srv_biom, proj_srv)
  om$data_list$srv_biom <- dplyr::arrange(om$data_list$srv_biom, Fleet_code, abs(Year))

  # -- Nbyage
  if(nrow(om$data_list$NByageFixed) > 0){
    proj_nbyage <- om$data_list$NByageFixed %>%
      group_by(Species, Sex) %>%
      slice(rep(n(),  proj_nyrs)) %>%
      mutate(Year = proj_yrs)
    proj_nbyage <- proj_nbyage[which(proj_yrs %!in% om$data_list$NByageFixed$Year),] # Subset rows already forcasted
    om$data_list$NByageFixed  <- rbind(om$data_list$NByageFixed, proj_nbyage)
    om$data_list$NByageFixed <- dplyr::arrange(om$data_list$NByageFixed, Species, Year)
  }

  # -- comp_data
  proj_comp <- om$data_list$comp_data %>%
    group_by(Fleet_code, Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = -proj_yrs)
  om$data_list$comp_data  <- rbind(om$data_list$comp_data, proj_comp)
  om$data_list$comp_data <- dplyr::arrange(om$data_list$comp_data, Fleet_code, abs(Year))

  # -- emp_sel - Use terminal year
  if(nrow(om$data_list$emp_sel) > 0){
    proj_emp_sel <- om$data_list$emp_sel %>%
      group_by(Fleet_code, Sex) %>%
      slice(rep(n(),  proj_nyrs)) %>%
      mutate(Year = proj_yrs)
    om$data_list$emp_sel  <- rbind(om$data_list$emp_sel, proj_emp_sel)
    om$data_list$emp_sel <- dplyr::arrange(om$data_list$emp_sel, Fleet_code, Year)
  }

  # -- wt
  #FIXME ignrores forecasted growth
  proj_wt <- om$data_list$wt %>%
    group_by(Wt_index , Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  om$data_list$wt  <- rbind(om$data_list$wt, proj_wt)
  om$data_list$wt <- dplyr::arrange(om$data_list$wt, Wt_index, Year)

  # -- Pyrs
  proj_Pyrs <- om$data_list$Pyrs %>%
    group_by(Species, Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  om$data_list$Pyrs  <- rbind(om$data_list$Pyrs, proj_Pyrs)
  om$data_list$Pyrs <- dplyr::arrange(om$data_list$Pyrs, Species, Year)


  #--------------------------------------------------
  # Update data in EM
  #FIXME - assuming same as terminal year of hindcast
  # -- EM emp_sel - Use terminal year
  proj_emp_sel <- em$data_list$emp_sel %>%
    group_by(Fleet_code, Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  em$data_list$emp_sel  <- rbind(em$data_list$emp_sel, proj_emp_sel)
  em$data_list$emp_sel <- dplyr::arrange(em$data_list$emp_sel, Fleet_code, Year)

  # -- EM wt
  proj_wt <- em$data_list$wt %>%
    group_by(Wt_index , Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  em$data_list$wt  <- rbind(em$data_list$wt, proj_wt)
  em$data_list$wt <- dplyr::arrange(em$data_list$wt, Wt_index, Year)

  # -- EM Pyrs
  proj_Pyrs <- em$data_list$Pyrs %>%
    group_by(Species, Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  em$data_list$Pyrs  <- rbind(em$data_list$Pyrs, proj_Pyrs)
  em$data_list$Pyrs <- dplyr::arrange(em$data_list$Pyrs, Species, Year)

  #--------------------------------------------------
  # Do the MSE
  for(sim in 1:nsim){

    # Set models
    Rceattle_EM_list[[sim]] <- list()
    Rceattle_EM_list[[sim]][[1]] <- em
    em_use <- em
    om_use <- om

    # Replace simulate future rec devs
    if(simulate){
      for(sp in 1:om_use$data_list$nspp){
        rec_dev <- rnorm( length(om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1]),
                          mean = log((1+(rec_trend/proj_nyrs) * 1:proj_nyrs)),
                          sd = exp(om_use$estimated_params$ln_rec_sigma[sp]))

        om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1] <- replace(
          om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1],
          values = rnorm( length(om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1]),
                          mean = 0,
                          sd = exp(om_use$estimated_params$ln_rec_sigma[sp])) # Assumed value from penalized likelihood
        )
      }
    }



    # Run through model
    for(k in 1:(length(assess_yrs))){

      # ------------------------------------------------------------
      # 1. OBSERVATION MODEL
      # ------------------------------------------------------------
      new_years <- proj_yrs[which(proj_yrs <= assess_yrs[k] & proj_yrs > om_use$data_list$endyr)]

      # - Get projected catch data from EM
      new_catch_data <- em_use$data_list$fsh_biom
      dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
      new_catch_data$Catch[dat_fill_ind] <- em_use$quantities$fsh_bio_hat[dat_fill_ind]
      if(!is.null(cap)){
        new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap[new_catch_data$Species[dat_fill_ind]], cap[new_catch_data$Species[dat_fill_ind]], new_catch_data$Catch[dat_fill_ind])
      }

      # - Update catch data in OM and EM
      om_use$data_list$fsh_biom <- new_catch_data
      em_use$data_list$fsh_biom <- new_catch_data

      # - Update endyr of OM
      nyrs_hind <- om_use$data_list$endyr - om_use$data_list$styr + 1
      om_use$data_list$endyr <- assess_yrs[k]

      # - Update parameters
      # -- F_dev
      om_use$estimated_params$F_dev <- cbind(om_use$estimated_params$F_dev, matrix(0, nrow= nrow(om_use$estimated_params$F_dev), ncol = length(new_years)))

      # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
      om_use$estimated_params$ln_srv_q_dev <- cbind(om_use$estimated_params$ln_srv_q_dev, matrix(om_use$estimated_params$ln_srv_q_dev[,ncol(om_use$estimated_params$ln_srv_q_dev)], nrow= nrow(om_use$estimated_params$ln_srv_q_dev), ncol = length(new_years)))

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


      # - Update map (Only new parameter we are estimating in OM is the F_dev of the new years)
      om_use$map <- build_map(
        data_list = om_use$data_list,
        params = om_use$estimated_params,
        debug = TRUE,
        random_rec = om_use$data_list$random_rec)
      om_use$map$mapFactor$dummy <- as.factor(NA); om_use$map$mapList$dummy <- NA


      # -- Estimate terminal F for catch
      new_f_yrs <- (ncol(om_use$map$mapList$F_dev) - length(new_years) + 1) : ncol(om_use$map$mapList$F_dev) # - Years of new F
      f_fleets <- om_use$data_list$fleet_control$Fleet_code[which(om_use$data_list$fleet_control$Fleet_type == 1)] # Fleet rows for F
      om_use$map$mapList$F_dev[f_fleets,new_f_yrs] <- replace(om_use$map$mapList$F_dev[f_fleets,new_f_yrs], values = 1:length(om_use$map$mapList$F_dev[f_fleets,new_f_yrs]))

      # -- Map out Fdev for years with 0 catch to very low number
      fsh_biom <- om_use$data_list$fsh_biom
      fsh_ind <- fsh_biom$Fleet_code[which(fsh_biom$Catch == 0)]
      yr_ind <- fsh_biom$Year[which(fsh_biom$Catch == 0)] - om_use$data_list$styr + 1
      om_use$map$mapList$F_dev[fsh_ind, yr_ind] <- NA
      om_use$map$mapFactor$F_dev <- factor( om_use$map$mapList$F_dev)

      # - Fit OM with new catch data
      om_use <- fit_mod(
        data_list = om_use$data_list,
        inits = om_use$estimated_params,
        map =  NULL,
        bounds = NULL,
        file = NULL,
        estimateMode = ifelse(om_use$data_list$estimateMode < 3, 1, om_use$data_list$estimateMode), # Estimate hindcast only if estimating
        random_rec = om_use$data_list$random_rec,
        niter = om_use$data_list$niter,
        msmMode = om_use$data_list$msmMode,
        avgnMode = om_use$data_list$avgnMode,
        minNByage = om_use$data_list$minNByage,
        suitMode = om_use$data_list$suitMode,
        suityr = om$data_list$endyr,
        loopnum = loopnum,
        phase = NULL,
        getsd = FALSE,
        verbose = 0)


      # ------------------------------------------------------------
      # 2. ESTIMATION MODEL
      # ------------------------------------------------------------
      # - Simulate new survey and comp data
      sim_dat <- sim_mod(om_use, simulate = simulate)

      years_include <- sample_yrs[which(sample_yrs$Year > em_use$data_list$endyr & sample_yrs$Year <= assess_yrs[k]),]

      # -- Add newly simulated survey data to EM
      new_srv_biom <- sim_dat$srv_biom[which(abs(sim_dat$srv_biom$Year) %in% years_include$Year & sim_dat$srv_biom$Fleet_code %in% years_include$Fleet_code),]
      new_srv_biom$Year <- -new_srv_biom$Year
      em_use$data_list$srv_biom <- rbind(em_use$data_list$srv_biom, new_srv_biom)
      em_use$data_list$srv_biom <- em_use$data_list$srv_biom[
        with(em_use$data_list$srv_biom, order(Fleet_code, abs(Year))),]

      # -- Add newly simulated comp data to EM
      new_comp_data <- sim_dat$comp_data[which(abs(sim_dat$comp_data$Year) %in% years_include$Year & sim_dat$comp_data$Fleet_code %in% years_include$Fleet_code),]
      new_comp_data$Year <- -new_comp_data$Year
      em_use$data_list$comp_data <- rbind(em_use$data_list$comp_data, new_comp_data)
      em_use$data_list$comp_data <- em_use$data_list$comp_data[
        with(em_use$data_list$comp_data, order(Fleet_code, abs(Year))),]

      # Update end year and re-estimate
      em_use$data_list$endyr <- assess_yrs[k]

      # Update parameters
      # -- F_dev
      em_use$estimated_params$F_dev <- cbind(em_use$estimated_params$F_dev, matrix(0, nrow= nrow(em_use$estimated_params$F_dev), ncol = length(new_years)))

      # -- Time-varying survey catachbilitiy - Assume last year - filled by columns
      em_use$estimated_params$ln_srv_q_dev <- cbind(em_use$estimated_params$ln_srv_q_dev, matrix(em_use$estimated_params$ln_srv_q_dev[,ncol(em_use$estimated_params$ln_srv_q_dev)], nrow= nrow(em_use$estimated_params$ln_srv_q_dev), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      ln_sel_slp_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
      sel_coff_dev = array(0, dim = c(nflts, 2, nselages_om, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for non-parameteric

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
      em_use <- fit_mod(
        data_list = em_use$data_list,
        inits = em_use$estimated_params,
        map =  NULL,
        bounds = NULL,
        file = NULL,
        estimateMode = ifelse(em_use$data_list$estimateMode < 3, 0, em_use$data_list$estimateMode), # Run hindcast and projection, otherwise debug
        HCR = build_hcr(HCR = em_use$data_list$HCR, # Tier3 HCR
                        DynamicHCR = em_use$data_list$DynamicHCR,
                        FsprTarget = em_use$data_list$FsprTarget,
                        FsprLimit = em_use$data_list$FsprLimit,
                        Ptarget = em_use$data_list$Ptarget,
                        Plimit = em_use$data_list$Plimit,
                        Alpha = em_use$data_list$Alpha,
                        Pstar = em_use$data_list$Pstar,
                        Sigma = em_use$data_list$Sigma
        ),
        random_rec = em_use$data_list$random_rec,
        niter = em_use$data_list$niter,
        msmMode = em_use$data_list$msmMode,
        avgnMode = em_use$data_list$avgnMode,
        minNByage = em_use$data_list$minNByage,
        suitMode = em_use$data_list$suitMode,
        phase = NULL,
        loopnum = loopnum,
        getsd = FALSE,
        verbose = 0)
      # plot_biomass(list(em_use, om_use), model_names = c("EM", "OM"))
      # End year of assessment

      # - Remove unneeded bits for memory reasons
      em_use$initial_params <- NULL
      em_use$bounds <- NULL
      em_use$map <- NULL
      em_use$obj <- NULL
      em_use$opt <- NULL
      em_use$sdrep <- NULL
      em_use$quantities[names(em_use$quantities) %!in% c("fsh_bio_hat",
                                                         "biomass",
                                                         "F",
                                                         "F_tot",
                                                         "mn_rec"  ,
                                                         "biomassSSB" ,
                                                         "R",
                                                         "srv_log_sd_hat",
                                                         "SB0",
                                                         "DynamicSB0",
                                                         "SPR0",
                                                         "SPRlimit",
                                                         "SPRtarget",
                                                         "DynamicNbyageSPR",
                                                         "DynamicSPR0",
                                                         "DynamicSPRlimit",
                                                         "DynamicSPRtarget",
                                                         "proj_F",
                                                         "Ftarget",
                                                         "Flimit",
                                                         "FlimitSPR",
                                                         "FtargetSPR",
                                                         "DynamicFlimitSPR",
                                                         "DynamicFtargetSPR")] <- NULL

      Rceattle_EM_list[[sim]][[k+1]] <- em_use
      message(paste0("Sim ",sim, " - EM Year ", assess_yrs[k], " COMPLETE"))
    }

    # Save models
    Rceattle_OM_list[[sim]] <- om_use
    names(Rceattle_EM_list[[sim]]) <- c("EM", paste0("OM_Sim_",sim,". EM_projyr_", assess_yrs))

    if(!is.null(file)){
      saveRDS(list(OM = om_use, EM = Rceattle_EM_list[[sim]]), file = paste0(file, "EMs_from_OM_Sim_",sim, ".rds"))
      Rceattle_EM_list[[sim]] <- NULL
    }
  }

  # - Name them
  names(Rceattle_OM_list) <- paste0("OM_Sim_",1:nsim)
  names(Rceattle_EM_list) <- paste0("OM_Sim_",1:nsim)

  return(list(OM_list = Rceattle_OM_list, EM_list = Rceattle_EM_list, OM = om, EM = em))
}
