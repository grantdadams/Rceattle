#' Run TMB using phases
#'
#' This function runs TMB with ADMB-like phasing of parameter estimation.
#' Function with normal inputs, passed via “...”, plus two additional arguments, “phase”
#' Optimizer by default is nlminb
#' phase is a tagged list where missing elements are populated with a vector of 1s, and non-missing elements are integers, and where the optimizer loops through values of phase while progressively changing map to turn on parameters
#'
#' @param  data A list to be passed to TMB
#' @param  parameters A list of parameters of the model
#' @param  map a list of map object from the model
#' @param  random A character vector of names of parameters that are random effects
#' @param  phases A list of the phases for the parameters of the model (same structure as your parameter list)
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param  model_name A string describing the model name. Must be the name of your .cpp file
#' @return A list of parameter estimates and their standard errors
#' @author Gavin Fay https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R
#' @export
#'
TMBphase <- function(data, parameters, map, random, phases, model_name,
                     silent, use_gradient = TRUE,
                     control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)) {

  # function to fill list component with a factor
  fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}

  #loop over phases
  for (phase_cur in 1:max(unlist(phases))) {
    #phase_cur <- 1

    # work out the map for this phase
    # if phases for parameters is less than the current phase
    # then map will contain a factor filled with NAs
    map_use <- map
    j <- 0
    for (i in 1:length(parameters)) {
      if (phases[[i]]>phase_cur) {
        map_val <- which(names(map_use) %in% names(phases)[i])
        map_use[[map_val]] <- fill_vals(map[[map_val]],NA)
      }
    }

    # remove the random effects if they are not estimated
    random_use <- random

    # initialize the parameters at values in previous phase
    params_use <- parameters
    if (phase_cur>1) params_use <- obj$env$parList(opt$par)

    # Fit the model
    obj <- TMB::MakeADFun(data,parameters =  params_use,random=random_use,DLL=model_name,map=map_use, silent = silent)

    if(use_gradient){
      opt <- nlminb(obj$par,obj$fn,obj$gr, control = control)
    }else{
      opt <- nlminb(obj$par,obj$fn)
    }
    last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))

    # write.csv(phase_cur, file = paste0("Phase",phase_cur))
    #close phase loop
  }

  return(last_par)
}

#' Function to set phasing order
#'
#' @returns list of parameter names with associated phase
#' @export
#'
set_phases <- function(){
  phaseList = list(
    dummy = 1,
    ln_pop_scalar = 4, # Scalar for input numbers-at-age
    rec_pars = 1, # Stock-recruit parameters or log(mean rec) if no stock-recruit relationship
    beta_rec_pars = 3,
    R_ln_sd = 2, # Variance for annual recruitment deviats
    rec_dev = 2, # Annual recruitment deviats
    init_dev = 2, # Age specific initial age-structure deviates or parameters
    # sex_ratio_ln_sd = 3, # Variance of sex ratio (usually fixed)
    ln_M1 = 4, #  Estimated natural or residual mortality
    ln_Flimit = 3, # Estimated F limit
    ln_Ftarget = 3, # Estimated F target
    ln_Finit = 3, # Estimated fishing mortality for non-equilibrium initial age-structure
    proj_F_prop = 1, # Fixed fleet-specific proportion of Flimit and Ftarget apportioned within each species
    ln_F = 1, # Annual fleet specific fishing mortality
    index_ln_q = 3, # Survey catchability
    index_q_dev = 5, # Annual survey catchability deviates (if time-varying)
    index_q_ln_sd = 4, # Prior SD for survey catchability deviates
    index_q_beta = 4, # Regression coefficients for environmental linkage
    index_q_rho = 4, # AR1 correlation parameter
    index_q_dev_ln_sd = 4, # SD for annual survey catchability deviates (if time-varying)
    sel_coff = 3, # Non-parametric selectivity coefficients
    sel_coff_dev = 4, # Annual deviates for non-parametric selectivity coefficients
    ln_sel_slp = 3, # Slope parameters for logistic forms of selectivity
    sel_inf = 3, # Asymptote parameters for logistic forms of selectivity
    ln_sel_slp_dev = 5, # Annual deviates for slope parameters for logistic forms of selectivity (if time-varying)
    sel_inf_dev = 5, # Annual deviates for asymptote parameters for logistic forms of selectivity (if time-varying)
    sel_dev_ln_sd = 4, # SD for annual selectivity deviates (if time-varying)
    sel_curve_pen = 4, # Penalty for non-parametric selectivity
    index_ln_sd = 2, # Log SD for survey lognormal index likelihood (usually input)
    catch_ln_sd = 2, # Log SD for lognormal catch likelihood (usually input)
    comp_weights = 5, # Weights for multinomial comp likelihood
    # ,logH_1 = 6,  # Functional form parameter (not used in MSVPA functional form)
    # logH_1a = 6, # Functional form parameter (not used in MSVPA functional form)
    # logH_1b = 6, # Functional form parameter (not used in MSVPA functional form)
    # logH_2 = 6, # Functional form parameter (not used in MSVPA functional form)
    # logH_3 = 6, # Functional form parameter (not used in MSVPA functional form)
    # H_4 = 6, # Functional form parameter (not used in MSVPA functional form)
    log_gam_a = 5, # Suitability parameter
    log_gam_b = 5, # Suitability parameter
    log_phi = 5 # Suitability parameter
  )


  # debugphase = list(
  #   dummy = 1,
  #   ln_pop_scalar = 5, # Scalar for input numbers-at-age
  #   rec_pars = 1, # Stock-recruit parameters or log(mean rec) if no stock-recruit relationship
  #   R_ln_sd = 4, # Variance for annual recruitment deviats
  #   rec_dev = 2, # Annual recruitment deviats
  #   init_dev = 3, # Age specific initial age-structure deviates or parameters
  #   sex_ratio_ln_sd = 3, # Variance of sex ratio (usually fixed)
  #   ln_M1 = 4, #  Estimated natural or residual mortality
  #   ln_mean_F = 6, # Mean fleet-specific fishing mortality
  #   ln_Flimit = 15, # Estimated F limit
  #   ln_Ftarget = 15, # Estimated F target
  #   ln_Finit = 7, # Estimated fishing mortality for non-equilibrium initial age-structure
  #   proj_F_prop = 14, # Fixed fleet-specific proportion of Flimit and Ftarget apportioned within each species
  #   F_dev = 7, # Annual fleet specific fishing mortality deviates
  #   index_ln_q = 10, # Survey catchability
  #   index_q_dev = 11, # Annual survey catchability deviates (if time-varying)
  #   index_q_ln_sd = 15, # Prior SD for survey catchability deviates
  #   index_q_dev_ln_sd = 15, # SD for annual survey catchability deviates (if time-varying)
  #   sel_coff = 8, # Non-parametric selectivity coefficients
  #   sel_coff_dev = 11, # Annual deviates for non-parametric selectivity coefficients
  #   ln_sel_slp = 9, # Slope parameters for logistic forms of selectivity
  #   sel_inf = 9, # Asymptote parameters for logistic forms of selectivity
  #   ln_sel_slp_dev = 11, # Annual deviates for slope parameters for logistic forms of selectivity (if time-varying)
  #   sel_inf_dev = 11, # Annual deviates for asymptote parameters for logistic forms of selectivity (if time-varying)
  #   sel_dev_ln_sd = 12, # SD for annual selectivity deviates (if time-varying)
  #   sel_curve_pen = 13, # Penalty for non-parametric selectivity
  #   index_ln_sd = 14, # Log SD for survey lognormal index likelihood (usually input)
  #   catch_ln_sd = 14, # Log SD for lognormal catch likelihood (usually input)
  #   comp_weights = 15, # Weights for multinomial comp likelihood
  #   logH_1 = 15,  # Functional form parameter (not used in MSVPA functional form)
  #   logH_1a = 15, # Functional form parameter (not used in MSVPA functional form)
  #   logH_1b = 15, # Functional form parameter (not used in MSVPA functional form)
  #   logH_2 = 15, # Functional form parameter (not used in MSVPA functional form)
  #   logH_3 = 15, # Functional form parameter (not used in MSVPA functional form)
  #   H_4 = 15, # Functional form parameter (not used in MSVPA functional form)
  #   log_gam_a = 15, # Suitability parameter (not used in MSVPA style)
  #   log_gam_b = 15, # Suitability parameter (not used in MSVPA style)
  #   log_phi = 15 # Suitability parameter (not used in MSVPA style)
  # )

  return(phaseList)
}
