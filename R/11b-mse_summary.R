#' Function to load .RDs files from MSE runs
#'
#' @param dir Directory used to save files from \code{\link{mse_run}}
#' @param file file name used to save files from \code{\link{mse_run}}
#'
#' @return list of MSE simulations/run
#' @export
#'
load_mse <- function(dir = NULL, file = NULL){

  # - Get file names
  mse_files <- list.files(path = dir, pattern = paste0(file, "EMs_from_OM_Sim_"))
  mse_order <- as.numeric(gsub(".rds", "", sapply(strsplit(mse_files, "EMs_from_OM_Sim_"), "[[", 2)))
  mse_files <- mse_files[order(mse_order)]

  ### Set up parallel processing
  library(foreach)
  library(doParallel)

  cores = (detectCores()/2)-1
  registerDoParallel(cores)

  mse <- foreach(i = 1:length(mse_files),
                 .combine = "c") %dopar% {
                   mse_tmp <- list(readRDS(file = paste0(dir,"/", mse_files[i])))

                   for(em in 2:length(mse_tmp[[1]]$EM)){
                     mse_tmp[[1]]$EM[[em]]$data_list$wt <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$emp_sel <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$age_trans_matrix <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$age_error <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$NByageFixed <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$aLW <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$UobsWtAge <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$Pyrs <- NULL
                     mse_tmp[[1]]$EM[[em]]$data_list$aLW <- NULL
                     mse_tmp[[1]]$EM[[em]]$estimated_params <- NULL
                   }

                   # Only use these bits for OM
                   mse_tmp[[1]]$OM$initial_params <- NULL
                   mse_tmp[[1]]$OM$bounds <- NULL
                   mse_tmp[[1]]$OM$map <- NULL
                   mse_tmp[[1]]$OM$obj <- NULL
                   mse_tmp[[1]]$OM$opt <- NULL
                   mse_tmp[[1]]$OM$sdrep <- NULL
                   mse_tmp[[1]]$OM$quantities[!names(mse_tmp[[1]]$OM$quantities) %in% c("fsh_bio_hat",
                                                                                        "fsh_log_sd_hat",
                                                                                        "srv_bio_hat",
                                                                                        "srv_log_sd_hat",
                                                                                        "depletion",
                                                                                        "depletionSSB",
                                                                                        "biomass",
                                                                                        "biomassSSB",
                                                                                        "BO",
                                                                                        "SB0",
                                                                                        "SBF",
                                                                                        "F_spp",
                                                                                        "R",
                                                                                        "M1",
                                                                                        "M",
                                                                                        "mean_rec",
                                                                                        "DynamicB0",
                                                                                        "DynamicSB0",
                                                                                        "DynamicSBF",
                                                                                        "SPR0",
                                                                                        "SPRlimit",
                                                                                        "SPRtarget",
                                                                                        "Ftarget",
                                                                                        "Flimit")] <- NULL


                   # Only use these bits for OM no F
                   mse_tmp[[1]]$OM_no_F$initial_params <- NULL
                   mse_tmp[[1]]$OM_no_F$bounds <- NULL
                   mse_tmp[[1]]$OM_no_F$map <- NULL
                   mse_tmp[[1]]$OM_no_F$obj <- NULL
                   mse_tmp[[1]]$OM_no_F$opt <- NULL
                   mse_tmp[[1]]$OM_no_F$sdrep <- NULL
                   mse_tmp[[1]]$OM_no_F$quantities[!names(mse_tmp[[1]]$OM_no_F$quantities) %in% c("fsh_bio_hat",
                                                                                                  "fsh_log_sd_hat",
                                                                                                  "srv_bio_hat",
                                                                                                  "srv_log_sd_hat",
                                                                                                  "depletion",
                                                                                                  "depletionSSB",
                                                                                                  "biomass",
                                                                                                  "biomassSSB",
                                                                                                  "BO",
                                                                                                  "SB0",
                                                                                                  "SBF",
                                                                                                  "F_spp",
                                                                                                  "R",
                                                                                                  "M1",
                                                                                                  "M",
                                                                                                  "mean_rec",
                                                                                                  "DynamicB0",
                                                                                                  "DynamicSB0",
                                                                                                  "DynamicSBF",
                                                                                                  "SPR0",
                                                                                                  "SPRlimit",
                                                                                                  "SPRtarget",
                                                                                                  "Ftarget",
                                                                                                  "Flimit")] <- NULL

                   # - Return
                   mse_tmp[[1]]$name <- mse_files[i]
                   mse_tmp
                 }
  # mse <- lapply(mse_files, function(x) readRDS(file = paste0(dir,"/", x)))
  closeAllConnections()

  names(mse) <- paste0("Sim_", 1:length(mse))
  return(mse)
}


# Average catch
# Interannual catch variation
# SSB/SSB40% from single species model
# Probability of being lower that SSB20% from single species model


#' Management strategy evaluation performance metric summary
#'
#' @param mse MSE runs from \code{\link{mse_run}} or \code{\link{load_mse}}
#'
#' @return Alist of two data.frames with MSE summary statistics of performance metrics including:
#' data.frame 1
#' 1.	Average annual catch across projection years and simulations per fleet and across fleets
#' 2.	Average interannual variation in catch (IAV) across projection years (n) per fleet and across fleets
#' 3.	% of years in which the fishery is closed across simulations (s)
#' 4.	Average relative mean squared error in estimate of spawning biomass in the terminal year across simulations
#' 5. % of years in which the population is perceived as undergoing overfishing as determined from F_Limit across simulations via \code{\link{build_hcr}} in the EM
#' 6.	% of years in which the population is perceived to be overfished  as determined from B_Limit across simulations via \code{\link{build_hcr}} in the EM
#' 7. % of years in which the population is undergoing overfishing as determined from the “true” F_Limit across simulations via \code{\link{build_hcr}} in the OM
#' 8. % of years in which the population is overfished as determined from the “true” B_Limit across simulations via \code{\link{build_hcr}} in the OM
#' 9.	Average ratio of spawning biomass over B_target in the terminal year across simulations in the OM
#' 10-14. Terminal biomass, SSB, SSB depletion (relative to equilibrium), SSB depletion (relative to dynamic SB0)
#'
#' @export
#'
mse_summary <- function(mse){

  ############################################
  ## Set up
  ############################################
  library(dplyr)

  # HCR Switches (make length of nspp if not)
  extend_length <- function(x, nspp){
    if(length(x) == nspp){ return(x)}
    else {return(rep(x, nspp))}
  }

  ## OM dimensions ----
  # - determined from OM sim 1
  # - should be the same as for the EM
  msmMode <- mse$Sim_1$OM$data_list$msmMode
  nspp <- mse$Sim_1$OM$data_list$nspp
  nsex <- mse$Sim_1$OM$data_list$nsex
  flt_type <- mse$Sim_1$OM$data_list$fleet_control$Fleet_type
  flts <- mse$Sim_1$OM$data_list$fleet_control$Fleet_code[which(flt_type == 1)]
  nflts = length(flts)
  flt_spp <- mse$Sim_1$OM$data_list$fleet_control$Species
  styr <- mse$Sim_1$OM$data_list$styr
  endyr <- mse$Sim_1$EM$EM$data_list$endyr
  projyr <- mse$Sim_1$EM$EM$data_list$projyr
  projyrs <- (endyr+1):projyr

  ## HCR ----
  HCR <- mse$Sim_1$EM[[1]]$data_list$HCR
  DynamicHCR <- mse$Sim_1$EM[[1]]$data_list$DynamicHCR
  Ptarget <- extend_length(mse$Sim_1$EM[[1]]$data_list$Ptarget, nspp)
  Plimit <- extend_length(mse$Sim_1$EM[[1]]$data_list$Plimit, nspp)
  Alpha <- extend_length(mse$Sim_1$EM[[1]]$data_list$Alpha, nspp)
  # -- HCR = 0: No catch - Params off
  # -- HCR = 1: Constant catch - Params off
  # -- HCR = 2: Constant input F - Params off
  # -- HCR = 3: F that acheives X% of SSB0 in the end of the projection - Ftarget on
  # -- HCR = 4: Constant target Fspr - Ftarget on
  # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 6: PFMC Cat 1 - Flimit on
  # -- HCR = 7: SESSF Tier 1 - Flimit and Ftarget on

  ## MSE specifications
  nsim <- length(mse)

  ## MSE Output ----
  # - Catch is by fleet
  mse_summary <- data.frame(matrix(NA, nrow = nflts+nspp+1, ncol = 21))
  colnames(mse_summary) <- c("Species",
                             "Fleet_name",
                             "Fleet_code",
                             "Average Catch",
                             "Catch IAV",
                             "P(Closed)",
                             "Avg SSB Relative MSE",
                             "Avg terminal SSB Relative MSE",
                             "EM: P(Fy > Flimit)",
                             "EM: P(SSB < SSBlimit)",
                             "OM: P(Fy > Flimit)",
                             "OM: P(SSB < SSBlimit)",
                             "EM: P(Fy > Flimit) but OM: P(Fy < Flimit)",
                             "EM: P(Fy < Flimit) but OM: P(Fy > Flimit)",
                             "EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)",
                             "EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)",
                             # "OM: Recovery Time",
                             "OM: Terminal B",
                             "OM: Terminal SSB",
                             "OM: Terminal SSB Depletion",
                             "OM: Terminal SSB Depletion (Dynamic)",
                             "OM: Average SSB Depletion")
  mse_summary$Fleet_name <- c(rep(NA, nspp), mse$Sim_1$OM$data_list$fleet_control$Fleet_name[flts], "All")
  mse_summary$Fleet_code <- c(rep(NA, nspp), mse$Sim_1$OM$data_list$fleet_control$Fleet_code[flts], "All")
  mse_summary$Species <- c(mse$Sim_1$OM$data_list$spnames, mse$Sim_1$OM$data_list$fleet_control$Species[flts], "All")


  ## Catch performance metrics by fleet ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  for(i in 1:nflts){
    flt = flts[i]

    # * Mean catch ----
    mse_summary$`Average Catch`[i+nspp] <- mean(
      sapply(mse, function(x)
        x$OM$data_list$fsh_biom %>%
          filter(Fleet_code == flt & Year %in% projyrs) %>%
          pull(Catch)
      ), na.rm = TRUE)

    # * Catch IAV ----
    catch_list_tmp <- lapply(mse, function(x)
      x$OM$data_list$fsh_biom %>%
        filter(Fleet_code == flt & Year %in% projyrs) %>%
        pull(Catch)
    )

    # -- Average across simulations by fleet
    mse_summary$`Catch IAV`[i+nspp] = 0
    for(sim in 1:nsim){
      iav_tmp <- sum((lag(catch_list_tmp[[sim]], 1) - catch_list_tmp[[sim]])^2, na.rm = TRUE)/(length(projyrs) - 1) # Squared difference
      iav_tmp <- sqrt(iav_tmp) / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
      mse_summary$`Catch IAV`[i+nspp] <- mse_summary$`Catch IAV`[i+nspp] + iav_tmp
    }

    # * P(Closed) ----
    mse_summary$`P(Closed)`[i+nspp] <- sum(
      sapply(catch_list_tmp, function(x)
        length(which(x < 1)) # Using less than 1 here just in case super small catches and fishery is effectively close
        /length(x)))/nsim

  }


  ## Catch performance metrics by species ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  for(sp in 1:nspp){

    # * Mean catch ----
    mse_summary$`Average Catch`[sp] <- mean(
      sapply(mse, function(x)
        x$OM$data_list$fsh_biom %>%
          filter(Species == sp & Year %in% projyrs) %>%
          pull(Catch)
      ), na.rm = TRUE)

    # - Catch IAV ----
    catch_list_tmp <- suppressMessages(lapply(mse, function(x)
      x$OM$data_list$fsh_biom %>%
        filter(Species == sp & Year %in% projyrs) %>%
        group_by(Year) %>%
        summarise(Catch = sum(Catch)) %>%
        pull(Catch)
    )) # Sum catch across species

    # - Average across simulations
    mse_summary$`Catch IAV`[sp] <- 0 # Initialize
    for(sim in 1:nsim){
      iav_tmp <- sum((lag(catch_list_tmp[[sim]], 1) - catch_list_tmp[[sim]])^2, na.rm = TRUE)/(length(projyrs) - 1) # Squared difference
      iav_tmp <- sqrt(iav_tmp) / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
      mse_summary$`Catch IAV`[sp] <- mse_summary$`Catch IAV`[sp] + iav_tmp
    }

    # * P(Closed) ----
    mse_summary$`P(Closed)`[sp] <- sum(
      sapply(catch_list_tmp, function(x)
        length(which(x < 1)) # Using less than 1 here just in case super small catches and fishery is effectively close
        /length(x)))/nsim
  }


  ## Catch performance metrics across species ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  catch_list_tmp <- suppressMessages(lapply(mse, function(x)
    x$OM$data_list$fsh_biom %>%
      filter(Year %in% projyrs) %>%
      group_by(Year) %>%
      summarise(Catch = sum(Catch)) %>%
      pull(Catch)
  )
  ) # Sum catch across species

  # * Mean catch ----
  mse_summary$`Average Catch`[nspp + nflts + 1] <- mean(unlist(catch_list_tmp), na.rm = TRUE)

  # * Catch IAV ----
  # -- Average across simulations
  mse_summary$`Catch IAV`[nspp + nflts + 1] <- 0 # Initialize
  for(sim in 1:nsim){
    iav_tmp <- sum((lag(catch_list_tmp[[sim]], 1) - catch_list_tmp[[sim]])^2, na.rm = TRUE)/(length(projyrs) - 1) # Squared difference
    iav_tmp <- sqrt(iav_tmp) / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
    mse_summary$`Catch IAV`[nspp + nflts + 1] <- mse_summary$`Catch IAV`[nspp + nflts + 1] + iav_tmp
  }


  ## Conservation performance metrics ----
  # - Avg terminal SSB MSE
  # - EM: P(Fy > Flimit)
  # - EM: P(SSB < SSBlimit)
  # - OM: P(Fy > Flimit)
  # - OM: P(SSB < SSBlimit)
  # - EM: P(Fy > Flimit) but OM: P(Fy < Flimit)
  # - EM: P(Fy < Flimit) but OM: P(Fy > Flimit)
  # - EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)
  # - EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)
  # - OM: Terminal Depletion Relative to equilibrium SB0
  # - OM: Terminal Depletion Relative to dynamic SB0
  # - EM: Average age-1 M
  # - EM: Variance of age-1 M

  # -- Tier 3 for single-species models
  # - Produces vectors of Flimits given depletion and input Flimit (Fspr)
  # - Note, it doesnt have Plimit because thats for cod
  flimit_tier3_fun <- function(depletionSSB, biomassSSB, SBF, plimit, alpha, Flimit){
    tier3_flimit <- c()
    for(i in 1:length(biomassSSB)){

      # Tier-3 HCR
      if(biomassSSB[i] >= SBF[i]){
        tier3_flimit[i] = Flimit
      }else if(biomassSSB[i] < SBF[i] & biomassSSB[i] > alpha * SBF[i]){
        tier3_flimit[i] = Flimit * (biomassSSB[i]/SBF[i] - alpha)/(1-alpha)
      }else{
        tier3_flimit[i] = 0
      }

      # If below 20%
      if(depletionSSB[i] < plimit){
        tier3_flimit[i] = 0
      }
    }
    return(Flimit)
  }


  for(sp in 1:nspp){

    ## Perceived status
    spp_rows <- which(flt_spp == sp) # FIXME: May want to select only fisheries (will bug if not survey)

    em_f_flimit <- c()
    em_sb_sblimit <- c()

    for(sim in 1:length(mse)){
      for(em in 2:length(mse[[sim]]$EM)){ # First EM is "conditioned" model

        # Terminal year of intermediate assessment
        end_yr_col <- mse[[sim]]$EM[[em]]$data_list$endyr - styr+1

        # * EM: P(F > Flimit) ----
        em_f_flimit <- c(em_f_flimit,
                         mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col]
                         > mse[[sim]]$EM[[em]]$quantities$Flimit[sp]
        )

        # - Tier 3
        if(HCR == 5 & DynamicHCR == FALSE){ # Adjust Tier 3
          em_f_flimit <- c(em_f_flimit,
                           mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col] >
                             flimit_tier3_fun(
                               depletionSSB = mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp,end_yr_col],
                               biomassSSB = mse[[sim]]$EM[[em]]$quantities$biomassSSB[sp,end_yr_col],
                               SBF = mse[[sim]]$EM[[em]]$quantities$SBF[sp,length(projyrs)],
                               Plimit[sp], Alpha[sp], mse[[sim]]$EM[[em]]$quantities$Flimit[sp]
                             )
          )
        }

        # Dynamic Tier 3
        if(HCR == 5 & DynamicHCR == TRUE){ # Adjust Tier 3
          em_f_flimit <- c(em_f_flimit,
                           mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col] >
                             flimit_tier3_fun(
                               depletionSSB = mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp,end_yr_col],
                               biomassSSB = mse[[sim]]$EM[[em]]$quantities$biomassSSB[sp,end_yr_col],
                               SBF = mse[[sim]]$EM[[em]]$quantities$DynamicSBF[sp,length(projyrs)],
                               Plimit[sp], Alpha[sp], mse[[sim]]$EM[[em]]$quantities$Flimit[sp]
                             )
          )
        }


        # * EM: P(SSB < SSBlimit) ----
        # Update over fished definition for the following because Plimit is used for something else
        # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
        # -- HCR = 6: PFMC Cat 1 - Flimit on
        #FIXME: convert to SBF when using ricker in EM
        if(HCR == 2){ # Avg F SPR based
          em_sb_sblimit <- c(em_sb_sblimit,
                             mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.5 * 0.35) # 0.5 * SB35%
        }else if(HCR == 4){ # New England SPR based
          em_sb_sblimit <- c(em_sb_sblimit,
                             mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.5 * 0.4) # 0.5 * SB40%
        }else if(HCR == 5){ # Tier 3 SPR Based
          em_sb_sblimit <- c(em_sb_sblimit,
                             mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.5 * 0.35) # 0.5 * SB35%
        }else if(HCR == 6){ # Cat 1 Depletion based
          if(Ptarget[sp] == 0.25){
            em_sb_sblimit <- c(em_sb_sblimit,
                               mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.125) # 0.125 * SB0
          }
          if(Ptarget[sp] == 0.4){
            em_sb_sblimit <- c(em_sb_sblimit,
                               mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.25) # 0.25 * SB0
          }
        }else if(HCR == 7){ # Tier 1 Depletion based
          em_sb_sblimit <- c(em_sb_sblimit,
                             mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < Plimit[sp])
        } else { # Otherwise depletion based
          em_sb_sblimit <- c(em_sb_sblimit,
                             mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < Plimit[sp])
        }
      }
    }

    ## Perceived status
    # - EM: P(F > Flimit)
    mse_summary$`EM: P(Fy > Flimit)`[sp] <- sum(em_f_flimit)/length(em_f_flimit)

    # - EM: P(SSB < SSBlimit)
    mse_summary$`EM: P(SSB < SSBlimit)`[sp] <- sum(em_sb_sblimit)/length(em_sb_sblimit)



    ## Actual status
    # * OM: P(F > Flimit) ----
    om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] > (x$OM$quantities$Flimit[sp]))

    # - Tier 3
    if(mse$Sim_1$OM$data_list$msmMode == 0 & HCR == 5 & DynamicHCR == FALSE){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] >
                              flimit_tier3_fun(
                                depletionSSB = x$OM$quantities$depletionSSB[sp,(projyrs - styr + 1)],
                                biomassSSB = x$OM$quantities$biomassSSB[sp,(projyrs - styr + 1)],
                                SBF = x$OM$quantities$SBF[sp,(projyrs - styr + 1)],
                                Plimit[sp], Alpha[sp], Flimit = x$OM$quantities$Flimit[sp]
                              )
      )
    }

    # - Dynamic Tier 3
    if(mse$Sim_1$OM$data_list$msmMode == 0 & HCR == 5 & DynamicHCR == TRUE){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] >
                              flimit_tier3_fun(
                                depletionSSB = x$OM$quantities$depletionSSB[sp,(projyrs - styr + 1)],
                                biomassSSB = x$OM$quantities$biomassSSB[sp,(projyrs - styr + 1)],
                                SBF = x$OM$quantities$DynamicSBF[sp,(projyrs - styr + 1)],
                                Plimit[sp], Alpha[sp], Flimit = x$OM$quantities$Flimit[sp]
                              )
      )
    }

    om_f_flimit <- unlist(om_f_flimit)
    mse_summary$`OM: P(Fy > Flimit)`[sp] <- sum(om_f_flimit)/length(om_f_flimit)


    # * OM: P(SSB < SSBlimit) ----
    # - Multi-species
    om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < x$OM$data_list$Plimit[sp])

    # Update over fished definition for the following because Plimit is used for something else
    # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
    # -- HCR = 6: PFMC Cat 1 - Flimit on
    if(mse$Sim_1$OM$data_list$msmMode == 0){

      # Default is depletion based
      om_sb_sblimit <-lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < x$OM$data_list$Plimit[sp])

      # - Avg F SPR based
      if(HCR == 2){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)] < 0.5 * x$OM$quantities$SBF[sp, length(projyrs)]) # 0.5 * SB35%
      }

      # - New England SPR based
      if(HCR == 4 & !DynamicHCR){
        om_sb_sblimit <-lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)] < 0.5 * x$OM$quantities$SBF[sp, length(projyrs)]) # 0.5 * SB40%
      }else if(HCR == 4 & DynamicHCR){ # Dynamic New England SPR based
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)] < 0.5 * x$OM$quantities$DynamicSBF[sp, length(projyrs)]) # 0.5 * SB40%
      }

      # - Tier 3 SPR Based
      if(HCR == 5 & !DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)] < 0.5 * x$OM$quantities$SBF[sp, length(projyrs)]) # 0.5 * SB35%
      }else if(HCR == 5 & DynamicHCR){ # Dynamic Tier 3 SPR Based
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)] < 0.5 * x$OM$quantities$DynamicSBF[sp, length(projyrs)]) # 0.5 * SB35%
      }

      # - Cat 1 Depletion based
      if(HCR == 6){
        if(Ptarget[sp] == 0.25){
          om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < 0.125) # 0.125 * SB0
        }
        if(Ptarget[sp] == 0.4){
          om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < 0.25) # 0.25 * SB0
        }
      }
    }

    om_sb_sblimit <- unlist(om_sb_sblimit)
    mse_summary$`OM: P(SSB < SSBlimit)`[sp] <- sum(om_sb_sblimit)/length(om_sb_sblimit)


    ## * Perceived status relative to actual status ----
    # - EM: P(Fy > Flimit) but OM: P(Fy < Flimit)
    mse_summary$`EM: P(Fy > Flimit) but OM: P(Fy < Flimit)`[sp] <- sum(em_f_flimit == 1 & om_f_flimit == 0)/length(om_f_flimit)

    # - EM: P(Fy < Flimit) but OM: P(Fy > Flimit)
    mse_summary$`EM: P(Fy < Flimit) but OM: P(Fy > Flimit)`[sp] <- sum(em_f_flimit == 0 & om_f_flimit == 1)/length(om_f_flimit)

    # - EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)
    mse_summary$`EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)`[sp] <- sum(em_sb_sblimit == 1 & om_sb_sblimit == 0)/length(om_sb_sblimit)

    # - EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)
    mse_summary$`EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)`[sp] <- sum(em_sb_sblimit == 0 & om_sb_sblimit == 1)/length(om_sb_sblimit)


    # * Bias in terminal SSB ----
    # - last projection year
    terminal_ssb_om <- sapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyr - styr + 1)])
    terminal_ssb_em <- lapply(mse, function(x) x$EM[[length(x$EM)]]$quantities$biomassSSB[sp, (projyr - styr + 1)])
    mse_summary$`Avg terminal SSB Relative MSE`[sp] = mean((unlist(terminal_ssb_em) -  unlist(terminal_ssb_om))^2 / unlist(terminal_ssb_om)^2, na.rm = TRUE)

    # * Bias in terminal SSB ----
    # - across all projection years
    ssb_om <- lapply(mse, function(x) x$OM$quantities$biomass[sp, (projyrs - styr + 1)])
    terminal_ssb_em_all <- list()
    for(i in 1:length(mse)){
      terminal_ssb_em_all[[i]] <- sapply(mse[[i]]$EM[-1], function(x) x$quantities$biomassSSB[sp, (x$data_list$endyr - styr + 1)])
    }
    mse_summary$`Avg SSB Relative MSE`[sp] = mean((unlist(terminal_ssb_em_all) -  unlist(ssb_om))^2 / unlist(ssb_om)^2, na.rm = TRUE)

    # * OM: Terminal B, SSB, depletion ----
    terminal_b_om <- sapply(mse, function(x) x$OM$quantities$biomass[sp, (projyr - styr + 1)])
    terminal_ssb_om <- sapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyr - styr + 1)])

    if(mse$Sim_1$OM$data_list$msmMode == 0){ # Take dynamic SB0 for multi-species model from OM projected with no F
      terminal_sb0_om <- sapply(mse, function(x) x$OM$quantities$SB0[sp, (projyr - styr + 1)])
      terminal_dynamic_sb0_om <- sapply(mse, function(x) x$OM$quantities$DynamicSB0[sp, (projyr - styr + 1)])
    }

    if(mse$Sim_1$OM$data_list$msmMode > 0){ # Take dynamic SB0 for multi-species model from OM projected with no F
      terminal_sb0_om <- sapply(mse, function(x) x$OM$quantities$SB0[sp]) # FIXME: SBO is adjusted in wrapper function
      terminal_dynamic_sb0_om <- sapply(mse, function(x) x$OM_no_F$quantities$biomassSSB[sp, (projyr - styr + 1)])
    }

    mse_summary$`OM: Terminal B`[sp] <- mean(terminal_b_om)
    mse_summary$`OM: Terminal SSB`[sp] <- mean(terminal_ssb_om)
    mse_summary$`OM: Terminal SSB Depletion`[sp] <- mean(terminal_ssb_om/terminal_sb0_om)
    mse_summary$`OM: Terminal SSB Depletion (Dynamic)`[sp] <- mean(terminal_ssb_om/terminal_dynamic_sb0_om)


    # - OM: Average SSB depletion
    sb_depletion <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)])
    sb_depletion <- unlist(sb_depletion)
    mse_summary$`OM: Average SSB Depletion`[sp] <- mean(sb_depletion)
  }

  return(mse_summary = mse_summary)
}
