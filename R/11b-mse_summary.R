#' Function to load .RDs files from MSE runs
#'
#' @param dir Directory used to save files from \code{\link{mse_run}}
#' @param file file name used to save files from \code{\link{mse_run}}
#'
#' @return list of MSE simulations/run
#' @export
#'
load_mse <- function(dir = NULL, file = NULL){
  mse_files <- list.files(path = dir, pattern = paste0(file, "EMs_from_OM_Sim_"))
  mse_order <- as.numeric(gsub(".rds", "", sapply(strsplit(mse_files, "EMs_from_OM_Sim_"), "[[", 2)))
  mse_files <- mse_files[order(mse_order)]
  mse = list()
  for(i in 1:length(mse_files)){
    mse[[i]] <- readRDS(file = paste0(dir,"/", mse_files[i]))
    mse[[i]]$OM$bounds <- NULL
    for(em in 2:length(mse[[i]]$EM)){
      mse[[i]]$EM[[em]]$data_list$wt <- NULL
      mse[[i]]$EM[[em]]$data_list$emp_sel <- NULL
      mse[[i]]$EM[[em]]$data_list$age_trans_matrix <- NULL
      mse[[i]]$EM[[em]]$data_list$age_error <- NULL
      mse[[i]]$EM[[em]]$data_list$NByageFixed <- NULL
      mse[[i]]$EM[[em]]$data_list$aLW <- NULL
      mse[[i]]$EM[[em]]$data_list$UobsWtAge <- NULL
      mse[[i]]$EM[[em]]$data_list$Pyrs <- NULL
      mse[[i]]$EM[[em]]$data_list$aLW <- NULL
      mse[[i]]$EM[[em]]$estimated_params <- NULL
    }
  }
  # mse <- lapply(mse_files, function(x) readRDS(file = paste0(dir,"/", x)))
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
#'
#' data.frame 2
#' 4.	Average relative mean squared error in estimate of spawning biomass in the terminal year across simulations
#' 5. % of years in which the population is perceived as undergoing overfishing as determined from F_Limit across simulations via \code{\link{build_hcr}} in the EM
#' 6.	% of years in which the population is perceived to be overfished  as determined from B_Limit across simulations via \code{\link{build_hcr}} in the EM
#' 7. % of years in which the population is undergoing overfishing as determined from the “true” F_Limit across simulations via \code{\link{build_hcr}} in the OM
#' 8. % of years in which the population is overfished as determined from the “true” B_Limit across simulations via \code{\link{build_hcr}} in the OM
#' 9.	Average ratio of spawning biomass over B_target in the terminal year across simulations in the OM
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
  projyrs_ind <- projyrs - endyr

  ## HCR ----
  HCR <- mse$Sim_1$EM[[1]]$data_list$HCR
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
  mse_summary <- data.frame(matrix(NA, nrow = nflts+nspp+1, ncol = 16))
  colnames(mse_summary) <- c("Species",
                             "Fleet_name",
                             "Fleet_code",
                             "Average Catch",
                             "Catch IAV",
                             "P(Closed)",
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
                             "OM: Terminal SSB Depletion")
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
      iav_tmp <- iav_tmp / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
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
      iav_tmp <- iav_tmp / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
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
    iav_tmp <- iav_tmp / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
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
  # - OM: Terminal SSB/SSBtarget
  # - EM: Average age-1 M
  # - EM: Variance of age-1 M

  # -- Tier 3 for single-species models
  # - Produces vectors of Flimits given depletion and input Flimit (Fspr)
  # - Note, it doesnt have Plimit because thats for cod
  flimit_vec_tier3_fun <- function(depletion, ptarget, alpha, Flimit){
    dynamic_tier3_flimit <- c()
    for(i in 1:length(depletion)){
      if(depletion[i] >= ptarget){
        dynamic_tier3_flimit[i] = Flimit
      }else if(depletion[i] < ptarget & depletion[i] > alpha * ptarget){
        dynamic_tier3_flimit[i] = Flimit * (depletion[i]/0.4 - alpha)/(1-alpha)
      }else{
        dynamic_tier3_flimit[i] = 0
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
        # - Tier 3
        if(HCR == 5){ # Adjust Tier 3
          #FIXME: using Ptarget rather than SPR
          if(mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] >= Ptarget[sp]){ # Above target
            em_f_flimit <- c(em_f_flimit,
                             mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col]
                             > mse[[sim]]$EM[[em]]$quantities$Flimit[sp]
            )
          }else if(mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < Ptarget[sp] & mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] > Alpha[sp] * Ptarget[sp]){ # Below target, but above limit

            em_f_flimit <- c(em_f_flimit,

                             mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col]
                             > mse[[sim]]$EM[[em]]$quantities$Flimit[sp] * (mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col]/Ptarget[sp] - Alpha[sp])/(1-Alpha[sp])
            )
          }else{ # Below limit
            em_f_flimit <- c(em_f_flimit,
                             mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col]
                             > 0)
          }
        } else{
          em_f_flimit <- c(em_f_flimit,
                           mse[[sim]]$EM[[em]]$quantities$F_spp[sp,end_yr_col]
                           > mse[[sim]]$EM[[em]]$quantities$Flimit[sp]
          )
        }


        # * EM: P(SSB < SSBlimit) ----
        # Update over fished definition for the following because Plimit is used for something else
        # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
        # -- HCR = 6: PFMC Cat 1 - Flimit on

        if(HCR %in% c(2,5)){
          em_sb_sblimit <- c(em_sb_sblimit,
                             mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.5 * 0.35) # 0.5 * SB35%
        }else if(HCR == 6){
          if(Ptarget[sp] == 0.25){
            em_sb_sblimit <- c(em_sb_sblimit,
                               mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.125) # 0.125 * SB0
          }

          if(Ptarget[sp] == 0.4){
            em_sb_sblimit <- c(em_sb_sblimit,
                               mse[[sim]]$EM[[em]]$quantities$depletionSSB[sp, end_yr_col] < 0.25) # 0.125 * SB0
          }
        }else{
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

    if(mse$Sim_1$OM$data_list$msmMode == 0 & HCR == 5){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] > flimit_vec_tier3_fun(
        depletion = x$OM$quantities$depletionSSB[sp,(projyrs - styr + 1)],
        ptarget = x$OM$data_list$Ptarget[sp],
        alpha  = x$OM$data_list$Alpha[sp],
        Flimit = x$OM$quantities$Flimit[sp]
      )
      )
    }

    om_f_flimit <- unlist(om_f_flimit)
    mse_summary$`OM: P(Fy > Flimit)`[sp] <- sum(om_f_flimit)/length(om_f_flimit)


    # * OM: P(SSB < SSBlimit) ----
    om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < x$OM$data_list$Plimit[sp])

    # Update over fished definition for the following because Plimit is used for something else
    # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
    # -- HCR = 6: PFMC Cat 1 - Flimit on
    if(mse$Sim_1$OM$data_list$msmMode == 0){
      if(HCR %in% c(2,5)){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < 0.5*0.35) # 0.5 * SB35%
      }

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


    ## Perceived status relative to actual status
    # - EM: P(Fy > Flimit) but OM: P(Fy < Flimit)
    mse_summary$`EM: P(Fy > Flimit) but OM: P(Fy < Flimit)`[sp] <- sum(em_f_flimit == 1 & om_f_flimit == 0)/length(om_f_flimit)

    # - EM: P(Fy < Flimit) but OM: P(Fy > Flimit)
    mse_summary$`EM: P(Fy < Flimit) but OM: P(Fy > Flimit)`[sp] <- sum(em_f_flimit == 0 & om_f_flimit == 1)/length(om_f_flimit)

    # - EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)
    mse_summary$`EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)`[sp] <- sum(em_sb_sblimit == 1 & om_sb_sblimit == 0)/length(om_sb_sblimit)

    # - EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)
    mse_summary$`EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)`[sp] <- sum(em_sb_sblimit == 0 & om_sb_sblimit == 1)/length(om_sb_sblimit)


    ## Average recovery time
    # om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)] < x$OM$data_list$Plimit)=
    # mse_summary$`OM: Recovery Time`[sp] <-

    ## Bias in terminal SSB
    sb_em <- lapply(mse, function(x) x$EM[[length(x$EM)]]$quantities$biomassSSB[sp, (projyrs - styr + 1)])
    sb_om <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)])

    mse_summary$`Avg terminal SSB Relative MSE`[sp] = mean((unlist(sb_em) -  unlist(sb_om))^2 / unlist(sb_om)^2, na.rm = TRUE)

    # - OM: Terminal SSB depletion
    sb_depletion <- lapply(mse, function(x) x$OM$quantities$depletionSSB[sp, (projyrs - styr + 1)])
    sb_depletion <- unlist(sb_depletion)
    mse_summary$`OM: Terminal SSB Depletion`[sp] <- mean(sb_depletion)
  }

  return(mse_summary = mse_summary)
}
