#' Run TMB using phases
#'
#' This function runs TMB with ADMB-like phasing of parameter estimation.
#' Phasing within an optimizer
#' Function with normal inputs, passed via “...”, plus two additional arguments, “phase” and “optimizer”
#' Optimizer by default is nlminb
#' phase is a tagged list where missing elements are populated with a vector of 1s, and non-missing elements are integers, and where the optimizer loops through values of phase while progressively changing map to turn on parameters
#'
#' @param  data A list to be passed to TMB
#' @param  parameters A list of parameters of the model
#' @param  map a list of map object from the model
#' @param  random A character vector of names of parameters that are random effects
#' @param  phases A list of the phases for the parameters of the model (same structure as your parameter list)
#' @param cpp_directory a string describing the model cpp directory.
#' @param  model_name A string describing the model name. Must be the name of your .cpp file
#' @param  optimizer The optimizer to use. Default is nlminb (This is not currently active)
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

TMBphase <- function(data, parameters, map, random, phases, cpp_directory, model_name,
                     optimizer = "nlminb", silent) {

  # function to fill list component with a factor
  fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}

  # # compile the model - NOTE: redundant
  # old_wd <- getwd()
  # setwd(cpp_directory)
  # TMB::compile(paste0(model_name, ".cpp"))
  # dyn.load(TMB::dynlib(paste0(model_name)), silent = TRUE)
  # setwd(old_wd)

  #loop over phases
  for (phase_cur in 1:max(unlist(phases))) {
    #phase_cur <- 1

    #work out the map for this phase
    # if phases for parameters is less than the current phase
    # then map will contain a factor filled with NAs
    map_use <- map
    j <- 0
    for (i in 1:length(parameters)) {
      if (phases[[i]]>phase_cur) {
        map_use[[i]] <- fill_vals(map[[i]],NA)
      }
    }
    #map_use

    #remove the random effects if they are not estimated
    random_use <- random

    # initialize the parameters at values in previous phase
    params_use <- parameters
    if (phase_cur>1) params_use <- obj$env$parList(opt$par)

    # Fit the model
    obj <- TMB::MakeADFun(data,parameters =  params_use,random=random_use,DLL=model_name,map=map_use, silent = silent)
    #obj <- TMB::MakeADFun(data,params_use,DLL=DLL_use,map=map_use)
    TMB::newtonOption(obj,smartsearch=FALSE)
    #obj$fn()
    #obj$gr()
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    last_par = suppressWarnings(obj$env$parList(obj$env$last.par.best))

    #close phase loop
  }

  return(last_par)
}
