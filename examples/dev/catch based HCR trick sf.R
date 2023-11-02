################################################
# Data
################################################
library(Rceattle)
library(dplyr)
# data("BS2017SS") # Single-species data. ?BS2017SS for more information on the data
# data("BS2017MS") # Multi-species data. Note: the only difference is the residual mortality (M1_base) is lower

data <- read_data("examples/dev/hake_intrasp_230912.xlsx")



################################################
# No fishing models - for reference
################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
# data$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ss_run <- Rceattle::fit_mod(data_list = data,
                            inits = NULL, # Initial parameters = 0
                            M1Fun = Rceattle::build_M1(M1_model = 1,
                                                       M1_use_prior = TRUE,
                                                       M1_prior_mean = 0.2,
                                                       M1_prior_sd = .1),
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)

# - Multi-species
# For the a multispecies model we from the single species parameters.
# BS2017MS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ms_run <- Rceattle::fit_mod(data_list = data,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            M1Fun = Rceattle::build_M1(M1_model = 1,
                                                       M1_use_prior = TRUE,
                                                       M1_prior_mean = 0.2,
                                                       M1_prior_sd = .1),
                            file = NULL, # Don't save
                            phase = "default",
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

##############################################################################
# Fish at F-SPR40$ or F to get B40% (equivalent in single-species model)
##############################################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
# BS2017SS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ss_run_F <- Rceattle::fit_mod(data_list = data,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              M1Fun = Rceattle::build_M1(M1_model = 1,
                                                         M1_use_prior = TRUE,
                                                         M1_prior_mean = 0.2,
                                                         M1_prior_sd = .1),
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              HCR = build_hcr(HCR = 4, # Constant Fspr 40%
                                              FsprTarget = 0.4, # 0.75 * F40%
                              ),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 1)

# - Multi-species
# For the a multispecies model we from the single species parameters.
# BS2017MS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ms_run_F <- Rceattle::fit_mod(data_list = data,
                              inits = ms_run$estimated_params, # Initial parameters from single species ests
                              M1Fun = Rceattle::build_M1(M1_model = 1,
                                                         M1_use_prior = TRUE,
                                                         M1_prior_mean = 0.2,
                                                         M1_prior_sd = .1),
                              HCR = build_hcr(HCR = 3, # F that acheives X% of SSB0
                                              FsprTarget = 0.4
                              ),
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              niter = 3, # 3 iterations around population and predation dynamics
                              random_rec = FALSE, # No random recruitment
                              msmMode = 1, # MSVPA based
                              suitMode = 0, # empirical suitability
                              verbose = 1)

##############################################################################
# Project model at 40-10 F based HCR
##############################################################################
# - Single-species
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
# BS2017SS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ss_run_4010_F <- Rceattle::fit_mod(data_list = data,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              M1Fun = Rceattle::build_M1(M1_model = 1,
                                                         M1_use_prior = TRUE,
                                                         M1_prior_mean = 0.2,
                                                         M1_prior_sd = .1),
                              estimateMode = 0, # Estimate
                              HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                                        FsprLimit = 0.4, # F40%
                                                        Ptarget = 0.4, # Target is 40% B0
                                                        Plimit = 0.1, # No fishing when SB<SB10
                                                        Pstar = 0.5,
                                                        Sigma = 0.5),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = NULL,
                              verbose = 1)

# - Multi-species
# For the a multispecies model we from the single species parameters.
# BS2017MS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
ms_run_4010_F <- Rceattle::fit_mod(data_list = data,
                              inits = ms_run$estimated_params, # Initial parameters from single species ests
                              M1Fun = Rceattle::build_M1(M1_model = 1,
                                                         M1_use_prior = TRUE,
                                                         M1_prior_mean = 0.2,
                                                         M1_prior_sd = .1),
                              HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                                        FsprLimit = 0.4, # F40%
                                                        Ptarget = 0.4, # Target is 40% B0
                                                        Plimit = 0.1, # No fishing when SB<SB10
                                                        Pstar = 0.5,
                                                        Sigma = 0.5),
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              niter = 3, # 3 iterations around population and predation dynamics
                              random_rec = FALSE, # No random recruitment
                              msmMode = 1, # MSVPA based
                              suitMode = 0, # empirical suitability
                              verbose = 1)


##############################################################################
# Project model at 40-10 catch based HCR (see function below)
##############################################################################
ss_run_4010_catch <- catch_hcr(model = ss_run_F)

ms_run_F$data_list$MSSB0 <- ms_run$quantities$SB0[, ncol(ms_run$quantities$biomassSSB)]
ms_run_4010_catch <- catch_hcr(model = ms_run_F)

##############################################################################
# Plot and compare
##############################################################################
# - Single species
plot_ssb(list(ss_run, ss_run_F, ss_run_4010_F, ss_run_4010_catch), model_names = c("No F", "F SPR40%", "40-10 F", "40-10 Catch"), incl_proj = TRUE)

plot_catch(list(ss_run, ss_run_F, ss_run_4010_F, ss_run_4010_catch), model_names = c("No F", "F SPR40%", "40-10 F", "40-10 Catch"), incl_proj = TRUE)

plot_catch(ss_run_F, incl_proj = TRUE) # does catch in the projection vary?

# Multi-species
plot_ssb(list(ms_run, ms_run_F, ms_run_4010_F, ms_run_4010_catch), incl_proj = TRUE)

plot_catch(list(ms_run, ms_run_F, ms_run_4010_F, ms_run_4010_catch), incl_proj = TRUE)

plot_catch(ms_run_4010_catch, incl_proj = TRUE) # does catch in the projection vary?



#' Project sloping catch-based HCR
#'
#' @description Projects a CEATTLE model using a sloping catch-based HCR. CEATTLE does not do catch-based HCRs internally so this must be done for projections. Projects the model until \code{projyr} in the data. The projection uses mean(R) for the projection.
#'
#' @param model CEATTLE model object exported from \code{\link{Rceattle}}
#' @param ptarget Target depletion (start of HCR) defaults to 0.4
#' @param plimit Limit depletion that defaults to 0.1
#' @param assessment_period Period of years that each assessment is taken
#' @param cap A cap on the catch in the projection. Can be a single number or vector. Default = NULL
#'
#' @return The model with an updated catch series
#' @export
#'
#'
catch_hcr <- function(model = model, ptarget = 0.4, plimit = 0.1, assessment_period = 1, cap = NULL){
  # model = ss_run
  # ptarget = 0.4
  # plimit = 0.1
  # assessment_period = 1
  # cap = NULL

  # - Adjust cap
  if(!is.null(cap)){
    if(length(cap) == 1){
      cap = rep(cap, model$data_list$nspp)
    }

    if(length(cap) != model$data_list$nspp){
      stop("cap is not length 1 or length nspp")
    }
  }

  # - Years for simulations
  endyr_base <- model$data_list$endyr
  hind_yrs <- (model$data_list$styr) : model$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  proj_yrs <- (model$data_list$endyr + 1) : model$data_list$projyr
  proj_nyrs <- length(proj_yrs)
  nflts = nrow(model$data_list$fleet_control)
  nselages <- max(model$data_list$fleet_control$Nselages, na.rm = TRUE)

  # - Assessment period
  assess_yrs <- seq(from = model$data_list$endyr + assessment_period, to = model$data_list$projyr,  by = assessment_period)

  # Assign multi-species SB0 and B0 if used
  model$data_list$MSSB0 <- model$quantities$SB0[,ncol(model$quantities$SB0)]
  model$data_list$MSB0 <- model$quantities$B0[,ncol(model$quantities$B0)]

  #--------------------------------------------------
  # Update data-files for updated years ----
  #--------------------------------------------------

  # -- Nbyage
  if(nrow(model$data_list$NByageFixed) > 0){
    proj_nbyage <- model$data_list$NByageFixed %>%
      group_by(Species, Sex) %>%
      slice(rep(n(),  proj_nyrs)) %>%
      mutate(Year = proj_yrs)
    proj_nbyage <- proj_nbyage[which(proj_yrs %!in% model$data_list$NByageFixed$Year),] # Subset rows already forcasted
    model$data_list$NByageFixed  <- rbind(model$data_list$NByageFixed, proj_nbyage)
    model$data_list$NByageFixed <- dplyr::arrange(model$data_list$NByageFixed, Species, Year)
  }

  # -- emp_sel - Use terminal year
  if(nrow(model$data_list$emp_sel) > 0){
    proj_emp_sel <- model$data_list$emp_sel %>%
      group_by(Fleet_code, Sex) %>%
      slice(rep(n(),  proj_nyrs)) %>%
      mutate(Year = proj_yrs)
    model$data_list$emp_sel  <- rbind(model$data_list$emp_sel, proj_emp_sel)
    model$data_list$emp_sel <- dplyr::arrange(model$data_list$emp_sel, Fleet_code, Year)
  }

  # -- wt
  #FIXME ignrores forecasted growth
  proj_wt <- model$data_list$wt %>%
    group_by(Wt_index , Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  model$data_list$wt  <- rbind(model$data_list$wt, proj_wt)
  model$data_list$wt <- dplyr::arrange(model$data_list$wt, Wt_index, Year)

  # -- Pyrs
  proj_Pyrs <- model$data_list$Pyrs %>%
    group_by(Species, Sex) %>%
    slice(rep(n(),  proj_nyrs)) %>%
    mutate(Year = proj_yrs)
  model$data_list$Pyrs  <- rbind(model$data_list$Pyrs, proj_Pyrs)
  model$data_list$Pyrs <- dplyr::arrange(model$data_list$Pyrs, Species, Year)


  #--------------------------------------------------
  # Do the projection ----
  #--------------------------------------------------

  # Replace future rec devs with mean-rec
  for(sp in 1:model$data_list$nspp){
    rec_dev <- log(mean(model$quantities$R[sp,1:hind_nyrs]))  - log(model$quantities$R0[sp]) # - Scale mean rec for rec trend

    # - Update model with devs
    model$estimated_params$rec_dev[sp,proj_yrs - model$data_list$styr + 1] <- replace(
      model$estimated_params$rec_dev[sp,proj_yrs - model$data_list$styr + 1],
      values =  rec_dev)
  }

  model_base <- model

  # Run through assessment years
  for(k in 1:length(assess_yrs)){

    # ------------------------------------------------------------
    # 1. GET RECOMMENDED TAC FROM EM-HCR ----
    # ------------------------------------------------------------
    new_years <- proj_yrs[which(proj_yrs <= assess_yrs[k] & proj_yrs > model$data_list$endyr)]
    nyrs_hind <- model$data_list$endyr - model$data_list$styr + 1

    # - Get projected catch data
    new_catch_data <- model$data_list$fsh_biom
    dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
    fspp <- new_catch_data$Species[dat_fill_ind]
    catch_proj <- model_base$quantities$fsh_bio_hat[dat_fill_ind]

    # Terminal SSB and SB0
    SB0 <- model$quantities$SB0[,ncol(model$quantities$SB0)]
    SBterm <- model$quantities$biomassSSB[,hind_nyrs]

    # Apply 40-10 HCR ----
    for(i in 1:length(fspp)){
      # SB < SB0 * Ptarget
      if(SBterm[fspp[i]] < ptarget * SB0[fspp[i]]){
        catch_proj[i] <- catch_proj[i] * (SB0[fspp[i]] * ptarget * (SBterm[fspp[i]] - SB0[fspp[i]] * plimit)) / (SBterm[fspp[i]] * SB0[fspp[i]] * (ptarget - plimit))

        print(paste0("Scaling catch for ", model$data_list$spnames[fspp[i]]))
      }
      # SB < SB0 * Plimit
      if(SBterm[fspp[i]] < plimit * SB0[fspp[i]]){
        catch_proj[i] <- 0
      }
    }

    # - Fill in to OM
    new_catch_data$Catch[dat_fill_ind] <- catch_proj

    # - Apply cap
    if(!is.null(cap)){
      new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap[new_catch_data$Species[dat_fill_ind]], cap[new_catch_data$Species[dat_fill_ind]], new_catch_data$Catch[dat_fill_ind])
    }

    # - Update catch data in OM and EM
    model$data_list$fsh_biom <- new_catch_data

    # ------------------------------------------------------------
    # 2. UPDATE MODEL ----
    # ------------------------------------------------------------
    # - Update endyr
    model$data_list$endyr <- assess_yrs[k]

    # - Update parameters
    # -- F_dev
    model$estimated_params$F_dev <- cbind(model$estimated_params$F_dev, matrix(0, nrow= nrow(model$estimated_params$F_dev), ncol = length(new_years)))

    # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
    model$estimated_params$ln_srv_q_dev <- cbind(model$estimated_params$ln_srv_q_dev, matrix(model$estimated_params$ln_srv_q_dev[,ncol(model$estimated_params$ln_srv_q_dev)], nrow= nrow(model$estimated_params$ln_srv_q_dev), ncol = length(new_years)))

    # -- Time-varing selectivity - Assume last year - filled by columns
    ln_sel_slp_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
    sel_inf_dev = array(0, dim = c(2, nflts, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic
    sel_coff_dev = array(0, dim = c(nflts, 2, nselages, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for non-parameteric

    ln_sel_slp_dev[,,,1:nyrs_hind] <- model$estimated_params$ln_sel_slp_dev
    sel_inf_dev[,,,1:nyrs_hind] <- model$estimated_params$sel_inf_dev
    sel_coff_dev[,,,1:nyrs_hind] <- model$estimated_params$sel_coff_dev

    ln_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- ln_sel_slp_dev[,,,nyrs_hind]
    sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
    sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

    model$estimated_params$ln_sel_slp_dev <- ln_sel_slp_dev
    model$estimated_params$sel_inf_dev <- sel_inf_dev
    model$estimated_params$sel_coff_dev <- sel_coff_dev


    # - Update map (Only new parameter we are estimating in model is the F_dev of the new years)
    model$map <- build_map(
      data_list = model$data_list,
      params = model$estimated_params,
      debug = TRUE,
      random_rec = model$data_list$random_rec)
    model$map$mapFactor$dummy <- as.factor(NA); model$map$mapList$dummy <- NA


    # -- Estimate terminal F for updated catch data from HCR
    new_f_yrs <- (ncol(model$map$mapList$F_dev) - length(new_years) + 1) : ncol(model$map$mapList$F_dev) # - Years of new F
    f_fleets <- model$data_list$fleet_control$Fleet_code[which(model$data_list$fleet_control$Fleet_type == 1)] # Fleet rows for F
    model$map$mapList$F_dev[f_fleets,new_f_yrs] <- replace(model$map$mapList$F_dev[f_fleets,new_f_yrs], values = 1:length(model$map$mapList$F_dev[f_fleets,new_f_yrs]))

    # -- Map out Fdev for years with 0 catch to very low number
    fsh_biom <- model$data_list$fsh_biom
    fsh_ind <- fsh_biom$Fleet_code[which(fsh_biom$Catch == 0)]
    yr_ind <- fsh_biom$Year[which(fsh_biom$Catch == 0)] - model$data_list$styr + 1
    model$map$mapList$F_dev[fsh_ind, yr_ind] <- NA
    model$map$mapFactor$F_dev <- factor( model$map$mapList$F_dev)

    # - Fit OM with new catch data
    model <- fit_mod(
      data_list = model$data_list,
      inits = model$estimated_params,
      map =  model$map,
      bounds = NULL,
      file = NULL,
      estimateMode = ifelse(model$data_list$estimateMode < 3, 1, model$data_list$estimateMode), # Estimate hindcast only if estimating
      random_rec = model$data_list$random_rec,
      niter = model$data_list$niter,
      msmMode = model$data_list$msmMode,
      avgnMode = model$data_list$avgnMode,
      minNByage = model$data_list$minNByage,
      suitMode = model$data_list$suitMode,
      initMode = model$data_list$initMode,
      meanyr = model$data_list$meanyr,
      HCR = build_hcr(HCR = model$data_list$HCR, # Tier3 HCR
                      DynamicHCR = model$data_list$DynamicHCR,
                      FsprTarget = model$data_list$FsprTarget,
                      FsprLimit = model$data_list$FsprLimit,
                      Ptarget = model$data_list$Ptarget,
                      Plimit = model$data_list$Plimit,
                      Alpha = model$data_list$Alpha,
                      Pstar = model$data_list$Pstar,
                      Sigma = model$data_list$Sigma,
                      Fmult = model$data_list$Fmult
      ),
      recFun = build_srr(srr_fun = model$data_list$srr_fun,
                         srr_pred_fun = model$data_list$srr_pred_fun ,
                         proj_mean_rec = model$data_list$proj_mean_rec,
                         srr_est_mode  = model$data_list$srr_est_mode ,
                         srr_prior_mean = model$data_list$srr_prior_mean,
                         srr_prior_sd = model$data_list$srr_prior_sd),
      M1Fun = build_M1(M1_model= model$data_list$M1_model,
                       updateM1 = FALSE,
                       M1_use_prior = TRUE,
                       M2_use_prior = TRUE,
                       M1_prior_mean = 0.2,
                       M1_prior_sd = 0.1), # Dont update M1 from data, fix at previous parameters
      loopnum = 3,
      phase = NULL,
      getsd = FALSE,
      verbose = 0)

  }

  # Retun model and set endyr to value from start
  model$data_list$endyr <- endyr_base
  return(model)
}

