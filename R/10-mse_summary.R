

#' Management strategy evaluation performance metric summary
#'
#' @param mse MSE runs from \code{\link{run_mse}} or \code{\link{load_mse}}
#' @param om_only If TRUE, report only operating-model (true) status and skip the
#'   estimation-model (EM) perception metrics.
#'
#' @return A named list, one element per entity dimension so each metric lives
#'   only where it applies (no NA padding):
#'   * `species` -- a data.frame with one row per species (keyed by `Species`)
#'     of the conservation / status metrics: per-species average catch, catch
#'     inter-annual variability (IAV), and P(Closed) = the probability the
#'     fishery is closed (catch ~ 0); the relative mean-squared error of
#'     terminal and average SSB; the estimation-model- and operating-model-
#'     perceived overfishing / overfished probabilities (via
#'     \code{\link{build_hcr}}) and the probability each status is misclassified
#'     (estimation model disagrees with the operating-model truth); and terminal
#'     biomass, SSB, dynamic SB0, SSB depletion (equilibrium and dynamic),
#'     average SSB depletion, and SSB-collapse counts.
#'   * `fleet` -- a data.frame with one row per fishery fleet (keyed by
#'     `Fleet_code` / `Fleet_name`) of average catch, catch IAV, and P(Closed).
#'   * `total` -- a named numeric of the across-fleet total average catch and
#'     catch IAV.
#'   * `meta` -- run provenance (`nsim`, `nspp`, `nflts`, `HCR`, projection-year
#'     range).
#'
#'   All metrics are averaged across projection years and simulations.
#'
#' @export
#'
mse_summary <- function(mse, om_only = FALSE){

  ############################################
  ## Set up
  ############################################
  # HCR Switches (make length of nspp if not)
  extend_length <- function(x, nspp){
    if(length(x) == nspp){ return(x)}
    else {return(rep(x, nspp))}
  }

  ## Drop simulations that did not run to completion ----
  # run_mse() returns only a failure marker for those (no OM, no EMs), so there
  # is nothing to summarise and averaging over them would mix a partial
  # trajectory into every performance metric.
  failed <- vapply(mse, function(x) isFALSE(x$use_sim), logical(1))
  if (any(failed)) {
    warning(sum(failed), " of ", length(mse), " simulations did not complete ",
            "and are excluded: ",
            paste(names(mse)[failed], collapse = ", "), call. = FALSE)
    mse <- mse[!failed]
  }
  if (!length(mse)) {
    stop("No simulations completed; there is nothing to summarise.", call. = FALSE)
  }

  ## Simulations whose unfished reference run failed ----
  # The unfished run is a separate fit and can fail on its own, in which case
  # run_mse() keeps the simulation but leaves OM_no_F NULL. Such a simulation is
  # valid for everything except the no-F comparisons, so it is excluded from
  # those specifically -- numerator and denominator both. Left in,
  # `sum(NULL < 1000) > 0` answers FALSE, counting a missing reference run as
  # "did not collapse" and biasing the collapse metrics low.
  has_no_f <- vapply(mse, function(x) !is.null(x$OM_no_F), logical(1))
  mse_no_f <- mse[has_no_f]
  if (any(!has_no_f)) {
    warning(sum(!has_no_f), " of ", length(mse), " simulations have no unfished ",
            "reference run; the 'OM no F' and 'SSB Collapse from F' metrics are ",
            "over the remaining ", length(mse_no_f), ".", call. = FALSE)
  }

  ## OM dimensions ----
  # - determined from OM sim 1
  # - should be the same as for the EM
  msmMode <- mse[[1]]$OM$data_list$msmMode
  nspp <- mse[[1]]$OM$data_list$nspp
  nsex <- mse[[1]]$OM$data_list$nsex
  flt_type <- mse[[1]]$OM$data_list$fleet_control$Fleet_type
  flts <- mse[[1]]$OM$data_list$fleet_control$Fleet_code[which(flt_type == "Fishery")]
  nflts = length(flts)
  flt_spp <- mse[[1]]$OM$data_list$fleet_control$Species
  styr <- mse[[1]]$OM$data_list$styr
  endyr <- mse[[1]]$EM[[1]]$data_list$endyr
  projyr <- min(c(mse[[1]]$EM[[1]]$data_list$projyr, mse[[1]]$OM$data_list$projyr))
  projyrs <- (endyr+1):projyr

  ## HCR ----
  # - Normalize HCR to integer code. build_hcr() accepts either integer or
  #   string alias (e.g. "NPFMC"), and fit_mod()/run_mse() may store either
  #   form depending on processing path. All downstream comparisons in this
  #   function use integer codes (HCR == 5, etc.), so coerce once here.
  .normalize_hcr <- function(x) {
    if (is.numeric(x) || (is.character(x) && !is.na(suppressWarnings(as.integer(x))))) {
      return(as.integer(x))
    }
    idx <- match(as.character(x), names(hcr_map))
    if (is.na(idx)) {
      stop("Unrecognized HCR value: '", x, "'. Expected one of ",
           paste(names(hcr_map), collapse = ", "),
           " or integer code 0:", max(hcr_map), ".")
    }
    as.integer(unname(hcr_map[idx]))
  }
  HCR <- .normalize_hcr(mse[[1]]$EM[[1]]$data_list$HCR)
  DynamicHCR <- as.logical(mse[[1]]$EM[[1]]$data_list$DynamicHCR)
  Ptarget <- extend_length(mse[[1]]$EM[[1]]$data_list$Ptarget, nspp)
  Plimit <- extend_length(mse[[1]]$EM[[1]]$data_list$Plimit, nspp)
  Alpha <- extend_length(mse[[1]]$EM[[1]]$data_list$Alpha, nspp)
  # -- HCR = 0: No catch - Params off
  # -- HCR = 1: Constant catch - Params off
  # -- HCR = 2: Constant input F - Params off
  # -- HCR = 3: F that achieves X% of SSB0 in the end of the projection - Ftarget on
  # -- HCR = 4: Constant target Fspr - Ftarget on
  # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 6: PFMC Cat 1 - Flimit on
  # -- HCR = 7: SESSF Tier 1 - Flimit and Ftarget on

  ## MSE specifications
  nsim <- length(mse)

  ## MSE Output ----
  # - Catch is by fleet
  mse_summary <- data.frame(matrix(NA, nrow = nflts+nspp+1, ncol = 25))
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
                             "OM: Terminal Dynamic SB0",
                             "OM: Terminal SSB Depletion",
                             "OM: Terminal SSB Depletion (Dynamic)",
                             "OM: Average SSB Depletion",
                             "OM no F: SSB Collapse",
                             "OM: SSB Collapse",
                             "OM: SSB Collapse from F")
  mse_summary$Fleet_name <- c(rep("Conservation metric", nspp), mse[[1]]$OM$data_list$fleet_control$Fleet_name[flts], "All")
  mse_summary$Fleet_code <- c(rep(NA, nspp), mse[[1]]$OM$data_list$fleet_control$Fleet_code[flts], "All")
  mse_summary$Species <- c(mse[[1]]$OM$data_list$spnames,
                           mse[[1]]$OM$data_list$spnames[mse[[1]]$OM$data_list$fleet_control$Species[flts]],
                           "All")


  ## Catch performance metrics by fleet ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  for(i in 1:nflts){
    flt = flts[i]

    # * Mean catch ----
    mse_summary$`Average Catch`[i+nspp] <- mean(
      sapply(mse, function(x)
        x$OM$data_list$catch_data |>
          filter(Fleet_code == flt & Year %in% projyrs) |>
          pull(Catch)
      ), na.rm = TRUE)

    # * Catch IAV ----
    catch_list_tmp <- lapply(mse, function(x)
      x$OM$data_list$catch_data |>
        dplyr::filter(Fleet_code == flt & Year %in% projyrs) |>
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
        x$OM$data_list$catch_data |>
          filter(Species == sp & Year %in% projyrs) |>
          pull(Catch)
      ), na.rm = TRUE)

    # - Catch IAV ----
    catch_list_tmp <- suppressMessages(lapply(mse, function(x)
      x$OM$data_list$catch_data |>
        filter(Species == sp & Year %in% projyrs) |>
        group_by(Year) |>
        summarise(Catch = sum(Catch)) |>
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
    x$OM$data_list$catch_data |>
      filter(Year %in% projyrs) |>
      group_by(Year) |>
      summarise(Catch = sum(Catch)) |>
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
  # - EM: Average age-1 M_at_age
  # - EM: Variance of age-1 M_at_age

  # -- Tier 3 for single-species models
  # - Produces vectors of Flimits given depletion and input Flimit (Fspr)
  # - Note, it doesnt have Plimit because thats for cod
  flimit_tier3_fun <- function(ssb_depletion, ssb, SBF, plimit, alpha, Flimit){
    tier3_flimit <- c()
    for(i in 1:length(ssb)){

      # Tier-3 HCR
      if(ssb[i] >= SBF[i]){
        tier3_flimit[i] = Flimit
      }else if(ssb[i] < SBF[i] & ssb[i] > alpha * SBF[i]){
        tier3_flimit[i] = Flimit * (ssb[i]/SBF[i] - alpha)/(1-alpha)
      }else{
        tier3_flimit[i] = 0
      }

      # If below 20%
      if(ssb_depletion[i] < plimit){
        tier3_flimit[i] = 0
      }
    }
    return(tier3_flimit)
  }


  # Helper: SSB-limit threshold (depletion or absolute) given HCR
  # - returns list(value, is_depletion) so callers can compare against
  #   either ssb_depletion or absolute ssb depending on HCR
  ssb_limit_thresh <- function(HCR, DynamicHCR, Ptarget_sp, Plimit_sp, SBF_val, DynamicSBF_val){
    if(HCR == 2)                                  list(val = 0.5 * SBF_val,        is_dep = FALSE)
    else if(HCR == 4 & !DynamicHCR)               list(val = 0.5 * SBF_val,        is_dep = FALSE)
    else if(HCR == 4 &  DynamicHCR)               list(val = 0.5 * DynamicSBF_val, is_dep = FALSE)
    else if(HCR == 5 & !DynamicHCR)               list(val = 0.5 * SBF_val,        is_dep = FALSE)
    else if(HCR == 5 &  DynamicHCR)               list(val = 0.5 * DynamicSBF_val, is_dep = FALSE)
    else if(HCR == 6 & Ptarget_sp == 0.4)         list(val = 0.25,                 is_dep = TRUE)
    else if(HCR == 6)                             list(val = 0.5 * Plimit_sp,      is_dep = TRUE)
    else                                          list(val = Plimit_sp,            is_dep = TRUE)
  }

  for(sp in 1:nspp){

    ## Perceived status (EM) and OM at matching assess years ----
    if(!om_only){
      spp_rows <- which(flt_spp == sp) # FIXME: May want to select only fisheries (will bug if not survey)

      em_f_flimit          <- c()
      em_sb_sblimit        <- c()
      om_f_flimit_assess   <- c()
      om_sb_sblimit_assess <- c()
      em_ssb_assess        <- c()
      om_ssb_assess        <- c()

      for(sim in 1:length(mse)){
        em_list <- mse[[sim]]$EM
        om_q    <- mse[[sim]]$OM$quantities

        for(em in 2:length(em_list)){ # First EM is "conditioned" model

          # Terminal year of intermediate assessment
          em_endyr_sim <- em_list[[em]]$data_list$endyr
          end_yr_col   <- em_endyr_sim - styr + 1
          em_q         <- em_list[[em]]$quantities

          # * EM: P(F > Flimit) ----
          # - Tier-3 aware: replace base Flimit with depletion-adjusted Flimit when HCR == 5
          em_Flimit <- em_q$Flimit[sp]
          if(HCR == 5){
            em_SBF_yr <- if(DynamicHCR) em_q$DynamicSBF[sp, end_yr_col] else em_q$SBF[sp, end_yr_col]
            em_Flimit <- flimit_tier3_fun(
              ssb_depletion = em_q$ssb_depletion[sp, end_yr_col],
              ssb           = em_q$ssb[sp, end_yr_col],
              SBF           = em_SBF_yr,
              plimit = Plimit[sp], alpha = Alpha[sp],
              Flimit = em_q$Flimit[sp]
            )
          }
          em_f_flimit <- c(em_f_flimit, em_q$F_spp[sp, end_yr_col] > em_Flimit)

          # * EM: P(SSB < SSBlimit) ----
          # FIXME: EM uses fixed-depletion proxy for HCR 2 (0.5*0.35);
          #        OM uses 0.5*SBF. Acceptable when EM is mean-recruit, inconsistent under Ricker EM.
          em_thresh <- if(HCR == 2)                              0.5 * 0.35
                       else if(HCR == 4)                         0.5 * Ptarget[sp]
                       else if(HCR == 5)                         0.5 * Plimit[sp]
                       else if(HCR == 6 & Ptarget[sp] == 0.40)   0.25
                       else if(HCR == 6)                         0.5 * Ptarget[sp]
                       else                                      Plimit[sp]
          em_sb_sblimit <- c(em_sb_sblimit, em_q$ssb_depletion[sp, end_yr_col] < em_thresh)

          # * OM at the SAME assess year (for cross-tab and SSB MSE) ----
          # - OM F > Flimit (Tier-3 aware for single-species mode)
          om_Flimit <- om_q$Flimit[sp]
          if(msmMode == 0 & HCR == 5){
            om_SBF_yr <- if(DynamicHCR) om_q$DynamicSBF[sp, end_yr_col] else om_q$SBF[sp, end_yr_col]
            om_Flimit <- flimit_tier3_fun(
              ssb_depletion = om_q$ssb_depletion[sp, end_yr_col],
              ssb           = om_q$ssb[sp, end_yr_col],
              SBF           = om_SBF_yr,
              plimit = Plimit[sp], alpha = Alpha[sp],
              Flimit = om_q$Flimit[sp]
            )
          }
          om_f_flimit_assess <- c(om_f_flimit_assess, om_q$F_spp[sp, end_yr_col] > om_Flimit)

          # - OM SSB < SSBlimit (uses absolute SBF for SPR-based HCRs, depletion otherwise)
          if(msmMode == 0){
            om_thr <- ssb_limit_thresh(HCR, DynamicHCR, Ptarget[sp], mse[[sim]]$OM$data_list$Plimit[sp],
                                       SBF_val        = om_q$SBF[sp, end_yr_col],
                                       DynamicSBF_val = om_q$DynamicSBF[sp, end_yr_col])
            om_below <- if(om_thr$is_dep) om_q$ssb_depletion[sp, end_yr_col] < om_thr$val
                        else              om_q$ssb[sp,            end_yr_col] < om_thr$val
          } else {
            om_below <- om_q$ssb_depletion[sp, end_yr_col] < mse[[sim]]$OM$data_list$Plimit[sp]
          }
          om_sb_sblimit_assess <- c(om_sb_sblimit_assess, om_below)

          # - Paired SSB values for relative MSE across assessments
          em_ssb_assess <- c(em_ssb_assess, em_q$ssb[sp, end_yr_col])
          om_ssb_assess <- c(om_ssb_assess, om_q$ssb[sp, end_yr_col])
        }
      }

      ## Perceived status (averages over all assess-year observations)
      mse_summary$`EM: P(Fy > Flimit)`[sp]    <- mean(em_f_flimit,   na.rm = TRUE)
      mse_summary$`EM: P(SSB < SSBlimit)`[sp] <- mean(em_sb_sblimit, na.rm = TRUE)
    }


    ## Actual status across all projection years ----
    # * OM: P(F > Flimit) ----
    om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] > (x$OM$quantities$Flimit[sp]))

    # - Tier 3
    if(msmMode == 0 & HCR == 5 & DynamicHCR == FALSE){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] >
                              flimit_tier3_fun(
                                ssb_depletion = x$OM$quantities$ssb_depletion[sp,(projyrs - styr + 1)],
                                ssb = x$OM$quantities$ssb[sp,(projyrs - styr + 1)],
                                SBF = x$OM$quantities$SBF[sp,(projyrs - styr + 1)],
                                Plimit[sp], Alpha[sp], Flimit = x$OM$quantities$Flimit[sp]
                              )
      )
    }

    # - Dynamic Tier 3
    if(msmMode == 0 & HCR == 5 & DynamicHCR == TRUE){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] >
                              flimit_tier3_fun(
                                ssb_depletion = x$OM$quantities$ssb_depletion[sp,(projyrs - styr + 1)],
                                ssb = x$OM$quantities$ssb[sp,(projyrs - styr + 1)],
                                SBF = x$OM$quantities$DynamicSBF[sp,(projyrs - styr + 1)],
                                Plimit[sp], Alpha[sp], Flimit = x$OM$quantities$Flimit[sp]
                              )
      )
    }

    om_f_flimit <- unlist(om_f_flimit)
    mse_summary$`OM: P(Fy > Flimit)`[sp] <- mean(om_f_flimit, na.rm = TRUE)


    # * OM: P(SSB < SSBlimit) ----
    # - default (multi-species or HCR with no SPR ref): depletion-based
    om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, (projyrs - styr + 1)] < x$OM$data_list$Plimit[sp])

    if(msmMode == 0){
      proj_idx <- projyrs - styr + 1

      # - Avg F SPR based (HCR 2)
      if(HCR == 2){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$SBF[sp, proj_idx])
      }

      # - New England SPR based (HCR 4)
      if(HCR == 4 & !DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$SBF[sp, proj_idx])
      } else if(HCR == 4 & DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$DynamicSBF[sp, proj_idx])
      }

      # - Tier 3 SPR based (HCR 5)
      if(HCR == 5 & !DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$SBF[sp, proj_idx])
      } else if(HCR == 5 & DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$DynamicSBF[sp, proj_idx])
      }

      # - Cat 1 depletion based (HCR 6)
      if(HCR == 6){
          om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, proj_idx] < Ptarget[sp] * 0.5)
        if(Ptarget[sp] == 0.4){
          om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, proj_idx] < 0.25)
        }
      }
    }

    om_sb_sblimit <- unlist(om_sb_sblimit)
    mse_summary$`OM: P(SSB < SSBlimit)`[sp] <- mean(om_sb_sblimit, na.rm = TRUE)


    ## * Perceived status relative to actual status (paired at assess years) ----
    if(!om_only){
      mse_summary$`EM: P(Fy > Flimit) but OM: P(Fy < Flimit)`[sp]      <- mean(em_f_flimit   == 1 & om_f_flimit_assess   == 0, na.rm = TRUE)
      mse_summary$`EM: P(Fy < Flimit) but OM: P(Fy > Flimit)`[sp]      <- mean(em_f_flimit   == 0 & om_f_flimit_assess   == 1, na.rm = TRUE)
      mse_summary$`EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)`[sp] <- mean(em_sb_sblimit == 1 & om_sb_sblimit_assess == 0, na.rm = TRUE)
      mse_summary$`EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)`[sp] <- mean(em_sb_sblimit == 0 & om_sb_sblimit_assess == 1, na.rm = TRUE)
    }


    # * Bias in SSB ----
    if(!om_only){
      # - terminal projection year (last EM vs OM)
      terminal_ssb_om <- sapply(mse, function(x) x$OM$quantities$ssb[sp, (projyr - styr + 1)])
      terminal_ssb_em <- sapply(mse, function(x) x$EM[[length(x$EM)]]$quantities$ssb[sp, (projyr - styr + 1)])
      mse_summary$`Avg terminal SSB Relative MSE`[sp] <- mean((terminal_ssb_em - terminal_ssb_om)^2 / terminal_ssb_om^2, na.rm = TRUE)

      # - across all intermediate assessments (EM terminal SSB vs OM SSB at the same year)
      mse_summary$`Avg SSB Relative MSE`[sp] <- mean((em_ssb_assess - om_ssb_assess)^2 / om_ssb_assess^2, na.rm = TRUE)
    }

    # * OM: Terminal B, SSB, depletion ----
    terminal_b_om <- sapply(mse, function(x) x$OM$quantities$biomass[sp, (projyr - styr + 1)])
    terminal_ssb_om <- sapply(mse, function(x) x$OM$quantities$ssb[sp, (projyr - styr + 1)])

    if(mse[[1]]$OM$data_list$msmMode == 0){ # Take dynamic SB0 for multi-species model from OM projected with no F
      terminal_sb0_om <- sapply(mse, function(x) x$OM$quantities$SB0[sp, (projyr - styr + 1)])
      terminal_dynamic_sb0_om <- sapply(mse, function(x) x$OM$quantities$DynamicSB0[sp, (projyr - styr + 1)])
    }

    if(mse[[1]]$OM$data_list$msmMode > 0){ # Take dynamic SB0 for multi-species model from OM projected with no F
      terminal_sb0_om <- sapply(mse, function(x) x$OM$quantities$SB0[sp]) # FIXME: SBO is adjusted in wrapper function
      terminal_dynamic_sb0_om <- if (length(mse_no_f)) {
        sapply(mse_no_f, function(x) x$OM_no_F$quantities$ssb[sp, (projyr - styr + 1)])
      } else NA_real_
    }

    mse_summary$`OM: Terminal B`[sp] <- mean(terminal_b_om)
    mse_summary$`OM: Terminal SSB`[sp] <- mean(terminal_ssb_om)
    mse_summary$`OM: Terminal Dynamic SB0`[sp] <- mean(terminal_dynamic_sb0_om)
    mse_summary$`OM: Terminal SSB Depletion`[sp] <- mean(terminal_ssb_om/terminal_sb0_om)
    mse_summary$`OM: Terminal SSB Depletion (Dynamic)`[sp] <- mean(terminal_ssb_om/terminal_dynamic_sb0_om)


    # - OM: Average SSB depletion
    sb_depletion <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, (projyrs - styr + 1)])
    sb_depletion <- unlist(sb_depletion)
    mse_summary$`OM: Average SSB Depletion`[sp] <- mean(sb_depletion)

    # * OM: Collapse ----
    mse_summary$`OM no F: SSB Collapse`[sp] <- if (length(mse_no_f)) {
      sum(sapply(mse_no_f, function(x) sum(x$OM_no_F$quantities$ssb[sp, (projyrs - styr + 1)] < 1000) > 0))
    } else NA_integer_
    mse_summary$`OM: SSB Collapse`[sp] <- sum(sapply(mse, function(x) sum(x$OM$quantities$ssb[sp, (projyrs - styr + 1)] < 1000) > 0))

    # -- OM without F is above cutoff, but OM with F is below
    mse_summary$`OM: SSB Collapse from F`[sp] <- if (length(mse_no_f)) {
      sum(sapply(mse_no_f,
                 function(x) sum((x$OM$quantities$ssb[sp, (projyrs - styr + 1)] < 1000) *
                                   (x$OM_no_F$quantities$ssb[sp, (projyrs - styr + 1)] > 1000)) > 0))
    } else NA_integer_
  }

  ############################################
  ## Reshape into per-entity tidy frames
  ############################################
  # The metrics above were accumulated into a single stacked frame (species
  # rows, then fishery-fleet rows, then an "All" row), so every column was
  # NA-padded outside the entity it applies to and the same column name (e.g.
  # "Average Catch") meant a species value, a fleet value, or a total depending
  # on the row. Return them split by entity instead: one frame per dimension,
  # each with a labelled key column and no padding.
  metric_cols   <- setdiff(colnames(mse_summary),
                           c("Species", "Fleet_name", "Fleet_code"))
  fleet_metrics <- c("Average Catch", "Catch IAV", "P(Closed)")
  species_rows  <- seq_len(nspp)
  fleet_rows    <- nspp + seq_len(nflts)
  total_row     <- nspp + nflts + 1L

  species <- data.frame(
    Species = mse_summary$Species[species_rows],
    mse_summary[species_rows, metric_cols, drop = FALSE],
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)

  fleet <- data.frame(
    # Fleet_code was character in the stacked frame only because of the "All"
    # total row (now in `total`); the per-fleet frame carries the integer codes.
    Fleet_code = as.integer(mse_summary$Fleet_code[fleet_rows]),
    Fleet_name = mse_summary$Fleet_name[fleet_rows],
    mse_summary[fleet_rows, fleet_metrics, drop = FALSE],
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)

  total <- unlist(mse_summary[total_row, c("Average Catch", "Catch IAV")])

  list(
    species = species,
    fleet   = fleet,
    total   = total,
    meta    = list(nsim = nsim, nspp = nspp, nflts = nflts, HCR = HCR,
                   proj_years = c(first = min(projyrs), last = max(projyrs)))
  )
}


#' Check which saved MSE simulation files can be loaded
#'
#' @param dir Directory used to save files from \code{\link{run_mse}}
#' @param file file name used to save files from \code{\link{run_mse}}
#'
#' @return list of MSE simulations/run
#' @export
#'
check_mse <- function(dir = NULL, file = NULL){

  # - Get file names
  mse_files <- list.files(path = dir, pattern = paste0(file, "EMs_from_OM_Sim_"))
  mse_order <- as.numeric(gsub(".rds", "", sapply(strsplit(mse_files, "EMs_from_OM_Sim_"), "[[", 2)))
  mse_files <- mse_files[order(mse_order)]

  ### Set up parallel processing (PSOCK cluster -- cross-platform,
  ### avoids the foreach::%dopar% rlang deparser interaction).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  use_parallel <- length(mse_files) > 1L && cores > 1L

  read_one <- function(i) {
    mse_tmp <- readRDS(file = paste0(dir, "/", mse_files[i]))
    check_tmp <- data.frame(Dir = dir, File = mse_files[i], Use = NA, Ind = i)
    if (!is.null(mse_tmp$use_sim)) {
      check_tmp$Use <- mse_tmp$use_sim
    }
    check_tmp
  }

  if (use_parallel) {
    cl <- parallel::makeCluster(min(cores, length(mse_files)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("dir", "mse_files"), envir = environment())
    rows <- parallel::parLapply(cl, seq_along(mse_files), read_one)
  } else {
    rows <- lapply(seq_along(mse_files), read_one)
  }
  check_df <- do.call(rbind, rows)

  return(check_df)
}




#' Load saved MSE simulation runs
#'
#' @param dir Directory used to save files from \code{\link{run_mse}}
#' @param file file name used to save files from \code{\link{run_mse}}
#' @param exclude index of MSE simulations not to load
#' @param include_em whether the EMs should be loaded or not (default = TRUE)
#'
#' @return list of MSE simulations/run
#' @export
#'
load_mse <- function(dir = NULL, file = NULL, exclude = NULL, include_em = TRUE){

  # - Get file names
  mse_files <- list.files(path = dir, pattern = paste0(file, "EMs_from_OM_Sim_"))
  mse_order <- as.numeric(gsub(".rds", "", sapply(strsplit(mse_files, "EMs_from_OM_Sim_"), "[[", 2)))
  mse_files <- mse_files[order(mse_order)]


  if(!is.null(exclude)){
    mse_files <- mse_files[-exclude]
  }

  ### Set up parallel processing
  # cores = (detectCores()/2)-1
  # doParallel::registerDoParallel(cores)
  mse_tmp <- list()
  for(i in 1:length(mse_files)){
    # mse <- foreach(i = 1:length(mse_files),
    #                .combine = "c") %dopar% {
    # mse_tmp <- list(readRDS(file = paste0(dir,"/", mse_files[i])))
    mse_tmp[[i]] <- (readRDS(file = paste0(dir,"/", mse_files[i])))

    # - Remove EM
    if(!include_em){
      mse_tmp[[i]]$EM <- NULL
    }

    # - Clean EM data
    if(length(mse_tmp[[i]]$EM) > 1){
      for(em in 2:length(mse_tmp[[i]]$EM)){
        mse_tmp[[i]]$EM[[em]]$data_list$weight <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$emp_sel <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$age_trans_matrix <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$age_error <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$NByageFixed <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$aLW <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$diet_data <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$ration_data <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$aLW <- NULL
        mse_tmp[[i]]$EM[[em]]$estimated_params <- NULL
      }
    }


    # Only use these bits for OM
    mse_tmp[[i]]$OM$initial_params <- NULL
    mse_tmp[[i]]$OM$bounds <- NULL
    mse_tmp[[i]]$OM$map <- NULL
    mse_tmp[[i]]$OM$obj <- NULL
    mse_tmp[[i]]$OM$opt <- NULL
    mse_tmp[[i]]$OM$sdrep <- NULL
    mse_tmp[[i]]$OM$quantities[!names(mse_tmp[[i]]$OM$quantities) %in% .mse_keep_quantities] <- NULL


    # Only use these bits for OM no F
    mse_tmp[[i]]$OM_no_F$initial_params <- NULL
    mse_tmp[[i]]$OM_no_F$bounds <- NULL
    mse_tmp[[i]]$OM_no_F$map <- NULL
    mse_tmp[[i]]$OM_no_F$obj <- NULL
    mse_tmp[[i]]$OM_no_F$opt <- NULL
    mse_tmp[[i]]$OM_no_F$sdrep <- NULL
    mse_tmp[[i]]$OM_no_F$quantities[!names(mse_tmp[[i]]$OM_no_F$quantities) %in% .mse_keep_quantities] <- NULL

    # - Return
    mse_tmp[[i]]$name <- mse_files[i]
    #  mse_tmp
  }
  # mse <- lapply(mse_files, function(x) readRDS(file = paste0(dir,"/", x)))
  closeAllConnections()

  names(mse_tmp) <- paste0("Sim_", 1:length(mse_tmp))
  return(mse_tmp)
}

