#' Run TMB using phases
#'
#' @description This function runs TMB with ADMB-like phasing of parameter estimation.
#' Phasing within an optimizer
#' Function with normal inputs, passed via “...”, plus two additional arguments, “phase” and “optimizer”
#' Optimizer by default is nlminb
#' phase is a tagged list where missing elements are populated with a vector of 1s, and non-missing elements are integers, and where the optimizer loops through values of phase while progressively changing map to turn on parameters
#' Adapted from code developed by Gavin Fay - https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R
#'
#' @param  data A list to be passed to TMB
#' @param  parameters A list of parameters of the model
#' @param  random A character vector of names of parameters that are random effects
#' @param  phases A list of the phases for the parameters of the model (same structure as your parameter list)
#' @param  model_name A string describing the model name. Must be the name of your .cpp file
#' @param  optimizer The optimizer to use. Default is nlminb (This is not currently active)
#' @return A list of parameter estimates and their standard errors
#' @author Gavin Fay
#' @export
#'


TMBphase <- function(data, parameters, random, phases, model_name,
                     optimizer = "nlminb") {

  # function to fill list component with a factor
  fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}

  # compile the model
  TMB::compile(paste0(model_name,".cpp"))
  dyn.load(TMB::dynlib(model_name))
  DLL_use <- model_name

  #loop over phases
  for (phase_cur in 1:max(unlist(phases))) {
    #phase_cur <- 1

    #work out the map for this phase
    # if phases for parameters is less than the current phase
    # then map will contain a factor filled with NAs
    map_use <- list()
    j <- 0
    for (i in 1:length(parameters)) {
      if (phases[[i]]>phase_cur) {
        j <- j+1
        map_use[[j]] <- fill_vals(parameters[[i]],NA)
        names(map_use)[j] <- names(parameters)[i]
      }
    }
    #map_use

    #remove the random effects if they are not estimated
    random_use <- random[!random%in%names(map_use)]

    # initialize the parameters at values in previous phase
    params_use <- parameters
    if (phase_cur>1) params_use <- obj$env$parList(opt$par)

    # Fit the model
    obj <- TMB::MakeADFun(data,params_use,random=random_use,DLL=DLL_use,map=map_use)
    #obj <- TMB::MakeADFun(data,params_use,DLL=DLL_use,map=map_use)
    TMB::newtonOption(obj,smartsearch=FALSE)
    #obj$fn()
    #obj$gr()
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    rep <- TMB::sdreport(obj)
    rep

    #close phase loop
  }

  return(rep)

  #close function TMBphase
}
