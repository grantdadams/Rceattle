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
  mse <- lapply(mse_files, function(x) readRDS(file = paste0(dir,"/", x)))
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
#' 4.	Average mean squared error in estimate of spawning biomass in the terminal year across simulations
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

  ## Operating model dimensions
  # - determined from OM sim 1
  # - should be the same as for the EM
  nspp <- mse$Sim_1$OM$data_list$nspp
  flt_type <- mse$Sim_1$OM$data_list$fleet_control$Fleet_type
  flts <- mse$Sim_1$OM$data_list$fleet_control$Fleet_code[which(flt_type == 1)]
  nflts = length(flts)
  flt_spp <- mse$Sim_1$OM$data_list$fleet_control$Species
  styr <- mse$Sim_1$OM$data_list$styr
  endyr <- mse$Sim_1$OM$data_list$suityr
  projyr <- mse$Sim_1$OM$data_list$projyr
  projyrs <- (endyr+1):projyr
  projyrs_ind <- projyrs - endyr

  ## MSE specifications
  nsim <- length(mse)

  ## MSE Output
  # - Catch is by fleet
  catch_summary_stats <- data.frame(matrix(0, nrow = nflts+nspp+1, ncol = 6))
  colnames(catch_summary_stats) <- c("Species", "Fleet_name", "Fleet_code","Average Catch", "Catch IAV", "% Years closed")
  catch_summary_stats$Fleet_name <- c(mse$Sim_1$OM$data_list$fleet_control$Fleet_name[flts], mse$Sim_1$OM$data_list$spnames, "All")
  catch_summary_stats$Fleet_code <- c(mse$Sim_1$OM$data_list$fleet_control$Fleet_code[flts], mse$Sim_1$OM$data_list$spnames, "All")
  catch_summary_stats$Species <- c(mse$Sim_1$OM$data_list$fleet_control$Species[flts], mse$Sim_1$OM$data_list$spnames, "All")

  # - Biomass is by species
  biomass_summary_stats <-  data.frame(matrix(0, nrow = nspp, ncol = 9))
  colnames(biomass_summary_stats) <-   c("Species", "Avg terminal SSB MSE", "EM: P(Fy > Flimit)", "EM: P(SSB < SSBlimit)", "EM: P(SSB < Dynamic SSBlimit)", "OM: P(Fy > Flimit)", "OM: P(SSB < SSBlimit)", "OM: P(SSB < Dynamic SSBlimit)", "OM: Terminal SSB/SSBtarget")
  biomass_summary_stats$Species <- c(mse$Sim_1$OM$data_list$spnames)


  ############################################
  ## Catch performance metrics by fleet
  ############################################
  # - Average catch
  # - Catch IAV
  # - % Years closed
  for(i in 1:nflts){
    flt = flts[i]

    # - Mean catch by fleet
    catch_summary_stats$`Average Catch`[flt] <- mean(sapply(mse, function(x)
      x[[1]]$data_list$fsh_biom$Catch[which(x[[1]]$data_list$fsh_biom$Fleet_code == flt &
                                              x[[1]]$data_list$fsh_biom$Year %in% projyrs)]), na.rm = TRUE)

    # - Catch IAV by fleet
    catch_list_tmp <- lapply(mse, function(x)
      x[[1]]$data_list$fsh_biom$Catch[which(x[[1]]$data_list$fsh_biom$Fleet_code == flt &
                                              x[[1]]$data_list$fsh_biom$Year %in% projyrs)])

    # -- Average across simulations by fleet
    for(sim in 1:nsim){
      catch_summary_stats$`Catch IAV`[flt] <- catch_summary_stats$`Catch IAV`[flt] + (sum((catch_list_tmp[[sim]][projyrs_ind[-1]] - catch_list_tmp[[sim]][projyrs_ind[-length(projyrs_ind)]])^2, na.rm = TRUE)/(length(projyrs_ind) - 1) / (sum(catch_list_tmp[[sim]][projyrs_ind], na.rm = TRUE)/ length(projyrs_ind)))/nsim
    }

    # - % Years closed by fleet
    catch_summary_stats$`% Years closed`[flt] <- mean(sapply(mse, function(x)
      length(x[[1]]$data_list$fsh_biom$Catch[which(x[[1]]$data_list$fsh_biom$Fleet_code == flt &
                                                     x[[1]]$data_list$fsh_biom$Year %in% projyrs &
                                                     x[[1]]$data_list$fsh_biom$Catch < 1)]) # Using less than 1 here just in case super small catches and fishery is effectively close
      /length(projyrs) * 100))
  }


  ############################################
  ## Catch performance metrics by species
  # - Average catch
  # - Catch IAV
  # - % Years closed
  for(sp in 1:nspp){

    # - Mean catch by species
    catch_summary_stats$`Average Catch`[sp + nflts] <- mean(sapply(mse, function(x)
      x[[1]]$data_list$fsh_biom$Catch[which(x[[1]]$data_list$fsh_biom$Species == sp &
                                              x[[1]]$data_list$fsh_biom$Year %in% projyrs)]), na.rm = TRUE)

    # - Catch IAV by species
    catch_list_tmp <- lapply(mse, function(x)
      x[[1]]$data_list$fsh_biom %>%
        filter(Species == sp & Year > endyr) %>%
        group_by(Year) %>%
        summarise(Catch = sum(Catch))) # Sum catch across species

    # -- Average across simulations
    for(sim in 1:nsim){
      catch_summary_stats$`Catch IAV`[sp+nflts] <- catch_summary_stats$`Catch IAV`[sp+nflts] + (sum((catch_list_tmp[[sim]]$Catch[projyrs_ind[-1]] - catch_list_tmp[[sim]]$Catch[projyrs_ind[-length(projyrs_ind)]])^2, na.rm = TRUE)/(length(projyrs_ind) - 1) / (sum(catch_list_tmp[[sim]]$Catch[projyrs_ind], na.rm = TRUE)/ length(projyrs_ind)))/nsim
    }

    # - % Years closed by species
    catch_summary_stats$`% Years closed`[sp+nflts] <- mean(sapply(catch_list_tmp, function(x)
      length(which(x$Catch < 1)) # Using less than 1 here just in case super small catches and fishery is effectively close
      /length(x$Catch) * 100))
  }


  ############################################
  ## Catch performance metrics across species
  # - Average catch
  # - Catch IAV
  # - % Years closed
  catch_list_tmp <- lapply(mse, function(x)
    x[[1]]$data_list$fsh_biom %>%
      filter(Year > endyr) %>%
      group_by(Year) %>%
      summarise(Catch = sum(Catch))) # Sum catch across species

  # - Mean catch across species
  catch_summary_stats$`Average Catch`[nspp + nflts + 1] <- mean(sapply(catch_list_tmp, function(x) mean(x$Catch)), na.rm = TRUE)

  # - Catch IAV across species
  # -- Average across simulations
  for(sim in 1:nsim){
    catch_summary_stats$`Catch IAV`[nspp + nflts + 1] <- catch_summary_stats$`Catch IAV`[nspp + nflts + 1] + (sum((catch_list_tmp[[sim]]$Catch[projyrs_ind[-1]] - catch_list_tmp[[sim]]$Catch[projyrs_ind[-length(projyrs_ind)]])^2, na.rm = TRUE)/(length(projyrs_ind) - 1) / (sum(catch_list_tmp[[sim]]$Catch[projyrs_ind], na.rm = TRUE)/ length(projyrs_ind)))/nsim
  }

  # - % Years closed across species
  catch_summary_stats$`% Years closed`[nspp + nflts + 1] <- NA


  ############################################
  ## Conservation performance metrics
  ############################################
  #FIXME currently based on EM rather than OM
  # - Avg terminal SSB MSE
  # - EM: P(Fy > Flimit)
  # - EM: P(SSB < SSBlimit)
  # - EM: P(SSB < Dynamic SSBlimit)
  # - OM: P(Fy > Flimit)
  # - OM: P(SSB < SSBlimit)
  # - OM: P(SSB < Dynamic SSBlimit)
  # - OM: Terminal SSB/SSBtarget
  for(sp in 1:nspp){

    ## Perceived status
    spp_rows <- which(flt_spp == sp) # FIXME: May want to select only fisheries (will bug if not survey)

    flimit_ratio_tmp <- c()
    sb_sblimit_tmp <- c()
    sb_sblimit_dynamic_tmp <- c()

    for(sim in 1:length(mse)){
      for(em in 2:length(mse[[sim]]$EM)){ # First EM is conditioned model
        # - EM: P(F > Flimit)
        end_yr_col <- mse[[sim]]$EM[[em]]$data_list$endyr - styr+1
        # -- Get terminal F by species
        flimit_ratio_tmp <- c(flimit_ratio_tmp,
                              sum(exp(mse[[sim]]$EM[[em]]$estimated_params$ln_mean_F + mse[[sim]]$EM[[em]]$estimated_params$F_dev)[spp_rows,end_yr_col])
                              / mse[[sim]]$EM[[em]]$quantities$Flimit[sp,end_yr_col]
        )

        # - EM: P(SSB < SSBlimit)
        sb_sblimit_tmp <- c(sb_sblimit_tmp,
                            mse[[sim]]$EM[[em]]$quantities$biomassSSB[sp, (projyrs - styr + 1)]/ (mse[[sim]]$EM[[em]]$quantities$SB0[sp] * mse[[sim]]$EM[[em]]$data_list$Plimit)
        )

        # - EM: P(SSB < Dynamic SSBlimit)
        sb_sblimit_dynamic_tmp <- c(sb_sblimit_dynamic_tmp,
                                    mse[[sim]]$EM[[em]]$quantities$biomassSSB[sp, (projyrs - styr + 1)]/ (mse[[sim]]$EM[[em]]$quantities$DynamicSB0[sp] * mse[[sim]]$EM[[em]]$data_list$Plimit)
        )
      }
    }

    ## Summarize
    # - EM: P(F > Flimit)
    biomass_summary_stats$`EM: P(Fy > Flimit)`[sp] <- length(which(flimit_ratio_tmp > 1))/length(flimit_ratio_tmp)

    # - EM: P(SSB < SSBlimit)
    biomass_summary_stats$`EM: P(SSB < SSBlimit)`[sp] <- length(which(sb_sblimit_tmp < 1))/length(sb_sblimit_tmp)

    # - EM: P(SSB < Dynamic SSBlimit)
    biomass_summary_stats$`EM: P(SSB < Dynamic SSBlimit)`[sp] <- length(which(sb_sblimit_tmp < 1))/length(sb_sblimit_tmp)

    ## Actual status
    # - OM: P(F > Flimit)
    flimit_ratio_tmp <- lapply(mse, function(x) colSums(exp(x$OM$estimated_params$ln_mean_F
                                                            + x$OM$estimated_params$F_dev)[spp_rows,])
                               / (x$OM$quantities$Flimit[sp,]))
    flimit_ratio_tmp <- unlist(flimit_ratio_tmp)
    biomass_summary_stats$`OM: P(Fy > Flimit)`[sp] <- length(which(flimit_ratio_tmp > 1))/length(flimit_ratio_tmp)


    # - OM: P(SSB < SSBlimit)
    sb_sblimit_tmp <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)]/ (x$OM$quantities$SB0[sp] * x$OM$data_list$Plimit))
    sb_sblimit_tmp <- unlist(sb_sblimit_tmp)
    biomass_summary_stats$`OM: P(SSB < SSBlimit)`[sp] <- length(which(sb_sblimit_tmp < 1))/length(sb_sblimit_tmp)


    # - OM: P(SSB < Dynamic SSBlimit)
    sb_sblimit_tmp <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)]/ (x$OM$quantities$DynamicSB0[sp,projyrs - styr + 1] * x$OM$data_list$Plimit))
    sb_sblimit_tmp <- unlist(sb_sblimit_tmp)
    biomass_summary_stats$`OM: P(SSB < Dynamic SSBlimit)`[sp] <- length(which(sb_sblimit_tmp < 1))/length(sb_sblimit_tmp)

    # - Bias in terminal SSB
    sb_em_tmp <- lapply(mse, function(x) x$EM[[length(x$EM)]]$quantities$biomassSSB[sp, (projyrs - styr + 1)])
    sb_om_tmp <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)])

    biomass_summary_stats$`Avg terminal SSB MSE`[sp] = mean((unlist(sb_em_tmp) -  unlist(sb_om_tmp))^2, na.rm = TRUE)

    # - OM: Terminal SSB/SSBtarget status
    sb_sbtarget_tmp <- lapply(mse, function(x) x$OM$quantities$biomassSSB[sp, (projyrs - styr + 1)]/ (x$OM$quantities$SB0[sp] * x$OM$data_list$Ptarget))
    sb_sbtarget_tmp <- unlist(sb_sbtarget_tmp)
    biomass_summary_stats$`OM: Terminal SSB/SSBtarget`[sp] <- mean(sb_sbtarget_tmp)
  }

  return(list(biomass_summary_stats = biomass_summary_stats, catch_summary_stats = catch_summary_stats))
}
