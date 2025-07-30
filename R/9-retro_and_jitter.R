#' Retrospective peels
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param peels the number of retrospective peels to use in the calculation of rho and for model estimation
#' @param rescale TRUE/FALSE wether to subset and rescale environmental predictors for the range of peel years.
#' @param nyrs_forecast Number of forecast years to calculate Mohn's Rho in addition to terminal year
#'
#' @return a list of 1. list of Rceattle models and 2. vector of Mohn's rho for each species
#'
#' @examples
#' data(BS2017SS) # ?BS2017SS for more information on the data
#' data('BS2017MS') # Note: the only difference is the residual mortality is lower
#'
#' ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
#'                             inits = NULL, # Initial parameters = 0
#'                             file = NULL, # Don't save
#'                             debug = 0, # Estimate
#'                             random_rec = FALSE, # No random recruitment
#'                             msmMode = 0, # Single species mode
#'                             silent = TRUE)
#'
#' retro <- retrospective(ss_run, peels = 10)
#' @export
retrospective <- function(Rceattle = NULL, peels = NULL, rescale = FALSE, nyrs_forecast = 3) {
  if (class(Rceattle) != "Rceattle") {
    stop("Object is not of class 'Rceattle'")
  }

  library(dplyr)
  # Get objects
  mod_list <- list(Rceattle)
  endyr <- Rceattle$data_list$endyr
  styr <- Rceattle$data_list$styr
  nyrs <- length(styr:endyr)
  projyr <- Rceattle$data_list$projyr
  nyrs_proj <- projyr - styr + 1

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Run retrospective peels ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  ind <- 2
  for (i in 1:peels) {

    # * Get end year of peel ----
    data_list <- Rceattle$data_list
    endyr_peel <- endyr - i
    data_list$endyr_peel <- endyr_peel
    nyrs_peel <- endyr_peel - styr + 1
    peel_prj_yrs <- (endyr_peel+1):endyr
    nyrs_proj_peel <- length(peel_prj_yrs)


    # * Turn off data after endyr_peel ----
    # - Catch is retained to update model dynamics
    data_list$index_data <- data_list$index_data %>%
      dplyr::filter(Year <= endyr_peel)
    #dplyr::mutate(Year = ifelse(Year > endyr_peel, -Year, Year))

    data_list$comp_data <- data_list$comp_data %>%
      dplyr::filter(Year <= endyr_peel)
    #dplyr::mutate(Year = ifelse(Year > endyr_peel, - Year, Year))

    data_list$diet_data <- data_list$diet_data %>%
      dplyr::filter(Year <= endyr_peel)
    #dplyr::mutate(Year = ifelse(Year > endyr_peel, - Year, Year))


    # * Assume weight/pyrs/emp_sel is same as last year of peel ----
    # -- weight
    #FIXME ignores forecasted growth
    data_list$weight <- data_list$weight %>%
      dplyr::filter(Year <= endyr_peel)

    proj_wt <- data_list$weight %>%
      group_by(Wt_index , Sex) %>%
      slice(rep(n(),  nyrs_proj_peel)) %>%
      mutate(Year = peel_prj_yrs)
    data_list$weight  <- rbind(data_list$weight, proj_wt)
    data_list$weight <- dplyr::arrange(data_list$weight, Wt_index, Year)

    # -- emp_sel
    data_list$emp_sel <- data_list$emp_sel %>%
      dplyr::filter(Year <= endyr_peel)
    if(nrow(data_list$emp_sel) > 0){
      proj_emp_sel <- data_list$emp_sel %>%
        group_by(Fleet_code, Sex) %>%
        slice(rep(n(),  nyrs_proj_peel)) %>%
        mutate(Year = peel_prj_yrs)
      data_list$emp_sel  <- rbind(data_list$emp_sel, proj_emp_sel)
      data_list$emp_sel <- dplyr::arrange(data_list$emp_sel, Fleet_code, Year)
    }

    # -- Pyrs
    data_list$Pyrs <- data_list$Pyrs %>%
      dplyr::filter(Year <= endyr_peel)

    if(nrow(data_list$Pyrs)){
      proj_Pyrs <- data_list$Pyrs %>%
        group_by(Species, Sex) %>%
        slice(rep(n(),  nyrs_proj_peel)) %>%
        mutate(Year = peel_prj_yrs)
      data_list$Pyrs  <- rbind(data_list$Pyrs, proj_Pyrs)
      data_list$Pyrs <- dplyr::arrange(data_list$Pyrs, Species, Year)
    }


    # * Rescale environmental predictors ----
    if(rescale){
      data_list$env_data <- data_list$env_data %>%
        dplyr::filter(Year <= endyr_peel)
      data_list$env_data[,2:ncol(data_list$env_data)]<-scale(data_list$env_data[,2:ncol(data_list$env_data)])
    }


    # * Adjust parameters ----
    # - Assume selectivity is same as last year of peel
    # - F remains on to fit to catch
    inits <- Rceattle$estimated_params
    inits$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- 0
    inits$ln_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- inits$ln_M1_dev[,,,nyrs_peel]
    inits$index_q_dev[,(nyrs_peel+1):nyrs] <- inits$index_q_dev[,nyrs_peel]
    inits$ln_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- inits$ln_sel_slp_dev[,,,nyrs_peel]
    inits$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_inf_dev[,,,nyrs_peel]
    inits$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_coff_dev[,,,nyrs_peel]

    # * Adjust map size ----
    map <- Rceattle$map
    map$mapList$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- NA
    map$mapFactor$rec_dev <- factor(map$mapList$rec_dev)

    map$mapList$ln_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- NA
    map$mapFactor$ln_M1_dev <- factor(map$mapList$ln_M1_dev)

    map$mapList$index_q_dev[,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$index_q_dev <- factor(map$mapList$index_q_dev)

    map$mapList$ln_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$ln_sel_slp_dev <- factor(map$mapList$ln_sel_slp_dev)

    map$mapList$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_inf_dev <- factor(map$mapList$sel_inf_dev)

    map$mapList$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)

    # -- Map out Fdev for years with 0 catch to very low number
    zero_catch <- data_list$catch_data %>%
      dplyr::filter(Year <= endyr &
                      Catch == 0) %>%
      dplyr::mutate(Year = Year - styr + 1) %>%
      select(Fleet_code, Year) %>%
      as.matrix()
    inits$ln_F[zero_catch] <- -999
    map$mapList$ln_F[zero_catch] <- NA
    map$mapFactor$ln_F <- factor(map$mapList$ln_F)
    rm(zero_catch)


    # * Refit ----
    newmod <- suppressWarnings(
      Rceattle::fit_mod(
        data_list = data_list,
        inits = inits,
        map =  map,
        bounds = NULL,
        file = NULL,
        estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
        HCR = build_hcr(HCR = data_list$HCR,
                        DynamicHCR = data_list$DynamicHCR,
                        Ftarget = data_list$Ftarget,
                        Flimit = data_list$Flimit,
                        Ptarget = data_list$Ptarget,
                        Plimit = data_list$Plimit,
                        Alpha = data_list$Alpha,
                        Pstar = data_list$Pstar,
                        Sigma = data_list$Sigma,
                        Fmult = data_list$Fmult,
                        HCRorder = data_list$HCRorder
        ),
        recFun = build_srr(srr_fun = data_list$srr_fun,
                           srr_pred_fun  = data_list$srr_pred_fun ,
                           proj_mean_rec  = data_list$proj_mean_rec ,
                           srr_meanyr = min(data_list$srr_meanyr, endyr_peel), # Update end year if less than srr_meanyr
                           srr_hat_styr = data_list$srr_hat_styr,
                           srr_hat_endyr = min(data_list$srr_hat_endyr, endyr_peel),
                           srr_est_mode  = data_list$srr_est_mode ,
                           srr_prior  = data_list$srr_prior,
                           srr_prior_sd   = data_list$srr_prior_sd,
                           Bmsy_lim = data_list$Bmsy_lim,
                           srr_indices = data_list$srr_indices),
        M1Fun =     build_M1(M1_model = data_list$M1_model,
                             M1_re = data_list$M1_re,
                             updateM1 = FALSE,  # Dont update M1 from data, fix at previous parameters
                             M1_use_prior = data_list$M1_use_prior,
                             M2_use_prior = data_list$M2_use_prior,
                             M_prior = data_list$M_prior,
                             M_prior_sd = data_list$M_prior_sd,
                             M1_indices = data_list$M1_indices),
        random_rec = data_list$random_rec,
        niter = data_list$niter,
        msmMode = data_list$msmMode,
        avgnMode = data_list$avgnMode,
        suitMode = data_list$suitMode,
        suit_styr = data_list$suit_styr,
        suit_endyr = min(data_list$suit_endyr, endyr_peel),   # Update to end year if less than suit_endyr
        initMode = data_list$initMode,
        phase = TRUE, # Phasing or else the parameters dont wanna move
        loopnum = data_list$loopnum,
        getsd = TRUE,
        verbose = 0)
    )

    # * Forecast ----
    # Adjust forecased rec_dev in new mod for bias and refit
    inits <- newmod$estimated_params
    for(sp in 1:newmod$data_list$nspp){

      # -- where SR curve is estimated directly
      if(newmod$data_list$srr_fun == newmod$data_list$srr_pred_fun){
        rec_dev <- log(mean(newmod$quantities$R[sp,1:nyrs_peel]))  - log(newmod$quantities$R0[sp])
      }

      # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
      if(newmod$data_list$srr_fun != newmod$data_list$srr_pred_fun){
        rec_dev <- log(mean((log(newmod$quantities$R) - log(newmod$quantities$R_hat))[sp, 1:nyrs_peel])) # - Scale mean rec for rec trend

      }

      # - Update OM with devs
      inits$rec_dev[sp, (peel_prj_yrs - styr + 1)] <- replace(
        inits$rec_dev[sp, (peel_prj_yrs - styr + 1)],
        values =  rec_dev)
    }

    newmod <- suppressWarnings(
      Rceattle::fit_mod(
        data_list = data_list,
        inits = inits,
        map =  map,
        bounds = NULL,
        file = NULL,
        estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
        HCR = build_hcr(HCR = data_list$HCR,
                        DynamicHCR = data_list$DynamicHCR,
                        Ftarget = data_list$Ftarget,
                        Flimit = data_list$Flimit,
                        Ptarget = data_list$Ptarget,
                        Plimit = data_list$Plimit,
                        Alpha = data_list$Alpha,
                        Pstar = data_list$Pstar,
                        Sigma = data_list$Sigma,
                        Fmult = data_list$Fmult,
                        HCRorder = data_list$HCRorder
        ),
        recFun = build_srr(srr_fun = data_list$srr_fun,
                           srr_pred_fun  = data_list$srr_pred_fun ,
                           proj_mean_rec  = data_list$proj_mean_rec ,
                           srr_meanyr = min(data_list$srr_meanyr, endyr_peel), # Update end year if less than srr_meanyr
                           srr_hat_styr = data_list$srr_hat_styr,
                           srr_hat_endyr = min(data_list$srr_hat_endyr, endyr_peel),
                           srr_est_mode  = data_list$srr_est_mode ,
                           srr_prior  = data_list$srr_prior,
                           srr_prior_sd   = data_list$srr_prior_sd,
                           Bmsy_lim = data_list$Bmsy_lim,
                           srr_indices = data_list$srr_indices),
        M1Fun =     build_M1(M1_model = data_list$M1_model,
                             M1_re = data_list$M1_re,
                             updateM1 = FALSE,  # Dont update M1 from data, fix at previous parameters
                             M1_use_prior = data_list$M1_use_prior,
                             M2_use_prior = data_list$M2_use_prior,
                             M_prior = data_list$M_prior,
                             M_prior_sd = data_list$M_prior_sd,
                             M1_indices = data_list$M1_indices),
        random_rec = data_list$random_rec,
        niter = data_list$niter,
        msmMode = data_list$msmMode,
        avgnMode = data_list$avgnMode,
        suitMode = data_list$suitMode,
        suit_styr = data_list$suit_styr,
        suit_endyr = min(data_list$suit_endyr, endyr_peel),   # Update to end year if less than suit_endyr
        initMode = data_list$initMode,
        phase = TRUE, # Phasing or else the parameters dont wanna move
        loopnum = data_list$loopnum,
        getsd = TRUE,
        verbose = 0)
    )

    # gc()
    #
    # map$mapFactor <- map$mapFactor[names(newmod$map$mapFactor)]
    # check <- c()
    # check_na <- c()
    # for(j in 1:length(map$mapList)){
    #   check[j] <- sum(map$mapFactor[[j]] != newmod$map$mapFactor[[j]], na.rm = TRUE)
    #   check_na[j] <- sum(is.na(map$mapFactor[[j]]) != is.na(newmod$map$mapFactor[[j]]), na.rm = TRUE)
    # }

    # Refit model If converged
    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        mod_list[[ind]] <- newmod
        ind <- ind + 1
      }
    }
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Calculate Mohs rho ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # * Data frame to save ----
  objects <- c("biomass", "ssb", "R", "F_spp")
  mohns <- data.frame(matrix(0, nrow = length(objects) * (nyrs_forecast+1), ncol = 3 + data_list$nspp))
  colnames(mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  # * Loop through peels ----
  for (i in 1:(length(mod_list) - 1)) {
    endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
    ind <- 1

    # * Loop output ----
    for (j in 1:length(objects)) {
      term_quantities <- mod_list[[1]]$quantities[[objects[j]]]
      retro_quantities <- mod_list[[i + 1]]$quantities[[objects[j]]]

      # * Loop forecast years ----
      for(yr in 0:nyrs_forecast){

        # If data exist for forecast (save)
        if(endyr_peel + yr <= endyr){

          # * Get full and peeled models ----
          base <- term_quantities[, nyrs_peel + yr]
          peel <- retro_quantities[, nyrs_peel + yr]
          rel_error <- ((peel - base)/base)

          # * Save and sum relative error ----
          mohns[ind, 1] <- objects[j]         # Object
          mohns[ind, 2] <- yr                 # Year
          mohns[ind, 3] <- mohns[ind, 3] + 1  # N
          mohns[ind, 4:(data_list$nspp + 3) ] <- mohns[j, 4:(data_list$nspp + 3)] + rel_error # Relative error
        }
        ind = ind+1
      }
    }
  }

  # * Divide N ----
  mohns[, 4:(data_list$nspp + 3) ] <- mohns[, 4:(data_list$nspp + 3)]/mohns[, 3]



  # * Beta coefficients ----
  objects <- colnames(Rceattle$estimated_params$beta_rec_pars)
  beta_mohns <- data.frame(matrix(0, nrow = length(objects), ncol = 3 + data_list$nspp))
  colnames(beta_mohns) <- c("Object", "Forecast year", "N", data_list$spnames)


  # * Loop through peels ----
  for (i in 1:(length(mod_list) - 1)) {
    endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
    ind <- 1

    # * Loop output ----
    for (j in 1:length(objects)) {
      base <- mod_list[[1]]$estimated_params$beta_rec_pars[,j]
      peel <- mod_list[[i + 1]]$estimated_params$beta_rec_pars[,j]
      rel_error <- ((peel - base)/base)

      # * Save and sum relative error ----
      beta_mohns[j, 1] <- objects[j]        # Object
      beta_mohns[j, 2] <- 0                 # Year
      beta_mohns[j, 3] <- beta_mohns[j, 3] + 1   # N
      beta_mohns[j, 4:(data_list$nspp + 3) ] <- beta_mohns[j, 4:(data_list$nspp + 3)] + rel_error # Relative error
    }
  }

  # * Divide N ----
  beta_mohns[, 4:(data_list$nspp + 3) ] <- beta_mohns[, 4:(data_list$nspp + 3)]/beta_mohns[, 3]


  return(list(Rceattle_list = rev(mod_list), mohns = rbind(mohns, beta_mohns)))
}




#' Jitter analysis
#'
#' @description Run's the Rceattle model at initial values that are +- N(0, 1) from the initial parameters.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param njitter the number of jitters to run
#' @param phase as in \code{\link{fit_mod}} default = FALSE
#' @param seed random number seed
#'
#' @return a list of Rceattle models
#'
#' @examples
#' data(BS2017SS) # ?BS2017SS for more information on the data
#' data('BS2017MS') # Note: the only difference is the residual mortality is lower
#'
#' ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
#'                             inits = NULL, # Initial parameters = 0
#'                             file = NULL, # Don't save
#'                             debug = 0, # Estimate
#'                             random_rec = FALSE, # No random recruitment
#'                             msmMode = 0, # Single species mode
#'                             silent = TRUE)
#'
#' jitters <- jitter(ss_run, njitter = 10)
#' @export
jitter <- function(Rceattle = NULL, njitter = 50, phase = FALSE, seed = 123) {
  if (class(Rceattle) != "Rceattle") {
    stop("Object is not of class 'Rceattle'")
  }

  set.seed(seed)

  # Run jitters ----
  mod_list <- list()
  ind = 1
  for (i in 1:njitter) {


    # * Adjust initial values ----
    inits <- Rceattle$initial_params
    mapList <- Rceattle$map$mapList
    data_list <- Rceattle$data_list

    for(j in 1:length(inits)){
      par <- names(inits)[j]
      inits[[j]] <- replace(inits[[j]],
                            values = ifelse(is.na(as.numeric(mapList[[par]])),
                                            as.numeric(inits[[j]]),
                                            as.numeric(inits[[j]]) + rnorm(length(as.numeric(inits[[j]])), 0, 1))
      )
    }


    # * Refit ----
    newmod <-
      suppressMessages(
        suppressWarnings(
          Rceattle::fit_mod(
            data_list = data_list,
            inits = inits,
            map =  NULL,
            bounds = NULL,
            file = NULL,
            estimateMode = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode), # Run hindcast and projection, otherwise debug
            HCR = build_hcr(HCR = data_list$HCR,
                            DynamicHCR = data_list$DynamicHCR,
                            Ftarget = data_list$Ftarget,
                            Flimit = data_list$Flimit,
                            Ptarget = data_list$Ptarget,
                            Plimit = data_list$Plimit,
                            Alpha = data_list$Alpha,
                            Pstar = data_list$Pstar,
                            Sigma = data_list$Sigma,
                            Fmult = data_list$Fmult,
                            HCRorder = data_list$HCRorder
            ),
            recFun = build_srr(srr_fun = data_list$srr_fun,
                               srr_pred_fun  = data_list$srr_pred_fun ,
                               proj_mean_rec  = data_list$proj_mean_rec ,
                               srr_meanyr = min(data_list$srr_meanyr, data_list$endyr), # Update end year if less than srr_meanyr
                               srr_hat_styr = data_list$srr_hat_styr,
                               srr_hat_endyr = data_list$srr_hat_endyr,
                               srr_est_mode  = data_list$srr_est_mode ,
                               srr_prior  = data_list$srr_prior,
                               srr_prior_sd   = data_list$srr_prior_sd,
                               Bmsy_lim = data_list$Bmsy_lim,
                               srr_indices = data_list$srr_indices),
            M1Fun =     build_M1(M1_model = data_list$M1_model,
                                 M1_re = data_list$M1_re,
                                 updateM1 = FALSE,  # Dont update M1 from data, fix at previous parameters
                                 M1_use_prior = data_list$M1_use_prior,
                                 M2_use_prior = data_list$M2_use_prior,
                                 M_prior = data_list$M_prior,
                                 M_prior_sd = data_list$M_prior_sd,
                                 M1_indices = data_list$M1_indices),
            random_rec = data_list$random_rec,
            niter = data_list$niter,
            msmMode = data_list$msmMode,
            avgnMode = data_list$avgnMode,
            suitMode = data_list$suitMode,
            suit_styr = data_list$suit_styr,
            suit_endyr = min(data_list$suit_endyr, data_list$endyr),   # Update to end year if less than suit_endyr
            initMode = data_list$initMode,
            phase = phase,
            loopnum = data_list$loopnum,
            getsd = TRUE,
            verbose = 0)
        )
      )

    # Refit model If converged
    if (!is.null(newmod$opt$Convergence_check)) {
      if (newmod$opt$Convergence_check != "The model is definitely not converged") {
        mod_list[[ind]] <- newmod
        ind <- ind + 1
      }
    }
  }


  # Plot ----
  jnll <- sapply(mod_list, function(x) x$quantities$jnll)
  # plot(x = 1:length(jnll), y = jnll)


  # Return ----
  return(list(Rceattle_list = mod_list, nll = jnll))
}

