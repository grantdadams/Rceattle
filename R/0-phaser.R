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
#' @examples
#'  setwd("~/Dropbox/ADMB/TMBphase/R")
#'  Y<-scan('thetalog.dat', skip=3, quiet=TRUE)
#'  data <- list(Y=Y)
#'  parameters <- list(
#'    X=data$Y*0,
#'    logr0=0,
#'    logtheta=0,
#'    logK=6,
#'    logQ=0,
#'    logR=0
#'  )
#' parameters$logQ <- -3
#'  random <- "X"
#'  model_name <- "thetalog"
#'  phases <- list(
#'    X=2,
#'    logr0=1,
#'    logtheta=1,
#'    logK=1,
#'    logQ=2,
#'    logR=1
#'  )
#'  TMBphase(data, parameters, random, model_name, optimizer = "nlminb")

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

    #remove the random effects if they are not estimated
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
