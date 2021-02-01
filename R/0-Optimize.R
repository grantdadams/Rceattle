
#' Optimize a TMB model
#'
#' \code{Optimize} runs a TMB model and generates standard diagnostics
#' https://raw.githubusercontent.com/kaskr/TMB_contrib_R/master/TMBhelper/R/Optimize.R
#' @param obj The compiled TMB object
#' @param startpar Starting values for fixed effects
#' @param lower lower bounds on fixed effects
#' @param upper upper bounds on fixed effects
#' @param getsd Boolean whether to run standard error calculation
#' @param control list of options to pass to \code{nlminb}
#' @param bias.correct Boolean whether to do epsilon bias-correction
#' @param bias.correct.control tagged list of options for epsilon bias-correction, where \code{vars_to_correct} is a character-vector of ADREPORT variables that should be bias-corrected
#' @param savedir directory to save results (if \code{savedir=NULL}, then results aren't saved)
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param n sample sizes (if \code{n!=Inf} then \code{n} is used to calculate BIC and AICc)
#' @param ... list of settings to pass to \code{sdreport}
#'
#' @return the standard output from \code{nlminb}, except with additional diagnostics and timing info, and a new slot containing the output from \code{sdreport}

#' @examples
#' TMBhelper::Optimize( Obj ) # where Obj is a compiled TMB object

#' @export
Optimize = function(obj, fn = obj$fn, gr = obj$gr, startpar = obj$par, lower = rep(-Inf, length(startpar)), upper = rep(Inf,
                                                                                                                        length(startpar)), getsd = TRUE, control = list(eval.max = 1e+06, iter.max = 1e+06, trace = 0), bias.correct = FALSE, bias.correct.control = list(sd = FALSE,
                                                                                                                                                                                                                                                                          split = NULL, nsplit = NULL, vars_to_correct = NULL), savedir = NULL, loopnum = 5, newtonsteps = 0, n = Inf, ...) {

  # Specify default values for bias.correct.control, which is internally called `BS.control`
  BS.control = list(sd = FALSE, split = NULL, nsplit = NULL, vars_to_correct = NULL)

  # Replace defaults for `BS.control` with provided values (if any)
  for (i in 1:length(bias.correct.control)) {
    if (tolower(names(bias.correct.control)[i]) %in% tolower(names(BS.control))) {
      BS.control[[match(tolower(names(bias.correct.control)[i]), tolower(names(BS.control)))]] = bias.correct.control[[i]]
    }
  }

  # Run first time
  start_time = Sys.time()
  parameter_estimates = nlminb(start = startpar, objective = fn, gradient = gr, control = control, lower = lower, upper = upper)

  # Re-run to further decrease final gradient
  for (i in seq(2, loopnum, length = max(0, loopnum - 1))) {
    Temp = parameter_estimates[c("iterations", "evaluations")]
    parameter_estimates = nlminb(start = parameter_estimates$par, objective = fn, gradient = gr, control = control, lower = lower,
                                 upper = upper)
    parameter_estimates[["iterations"]] = parameter_estimates[["iterations"]] + Temp[["iterations"]]
    parameter_estimates[["evaluations"]] = parameter_estimates[["evaluations"]] + Temp[["evaluations"]]
  }

  ## Run some Newton steps
  for (i in seq_len(newtonsteps)) {
    g = as.numeric(gr(parameter_estimates$par))
    h = optimHess(parameter_estimates$par, fn = fn, gr = gr)
    parameter_estimates$par = parameter_estimates$par - solve(h, g)
    parameter_estimates$objective = fn(parameter_estimates$par)
  }

  # Exclude difficult-to-interpret messages
  parameter_estimates = parameter_estimates[c("par", "objective", "iterations", "evaluations")]

  # Add diagnostics
  parameter_estimates[["time_for_MLE"]] = Sys.time() - start_time
  parameter_estimates[["max_gradient"]] = max(abs(gr(parameter_estimates$par)))
  parameter_estimates[["Convergence_check"]] = ifelse(parameter_estimates[["max_gradient"]] < 1e-04, "There is no evidence that the model is not converged",
                                                      paste0("The maximum gradient is ", round(parameter_estimates[["max_gradient"]], 6)))
  parameter_estimates[["number_of_coefficients"]] = c(Total = length(unlist(obj$env$parameters)), Fixed = length(startpar),
                                                      Random = length(unlist(obj$env$parameters)) - length(startpar))
  parameter_estimates[["AIC"]] = TMBhelper::TMBAIC(opt = parameter_estimates)
  if (n != Inf) {
    parameter_estimates[["AICc"]] = TMBhelper::TMBAIC(opt = parameter_estimates, n = n)
    parameter_estimates[["BIC"]] = TMBhelper::TMBAIC(opt = parameter_estimates, p = log(n))
  }
  parameter_estimates[["diagnostics"]] = data.frame(Param = names(startpar), starting_value = startpar, MLE = parameter_estimates$par, final_gradient = as.vector(gr(parameter_estimates$par)))

  # Get standard deviations
  if (getsd == TRUE) {
    # Compute hessian
    sd_time = Sys.time()
    h = optimHess(parameter_estimates$par, fn = fn, gr = gr)
    # Check for problems
    if (is.character(try(chol(h), silent = TRUE))) {
      warning("Hessian is not positive definite, so standard errors are not available")
      if (!is.null(savedir)) {
        capture.output(parameter_estimates, file = file.path(savedir, "parameter_estimates.txt"))
      }
      return(list(opt = parameter_estimates, h = h))
    }
    # Compute standard errors
    if (bias.correct == FALSE | is.null(BS.control[["vars_to_correct"]])) {
      if (!is.null(BS.control[["nsplit"]])) {
        if (BS.control[["nsplit"]] == 1)
          BS.control[["nsplit"]] = NULL
      }
      parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, par.fixed = parameter_estimates$par, hessian.fixed = h, bias.correct = bias.correct,
                                                  bias.correct.control = BS.control[c("sd", "split", "nsplit")], ...)
    } else {
      if ("ADreportIndex" %in% names(obj$env)) {
        Which = as.vector(unlist(obj$env$ADreportIndex()[BS.control[["vars_to_correct"]]]))
      } else {
        # Run first time to get indices
        parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, par.fixed = parameter_estimates$par, hessian.fixed = h,
                                                    bias.correct = FALSE, ...)
        # Determine indices
        Which = which(rownames(summary(parameter_estimates[["SD"]], "report")) %in% BS.control[["vars_to_correct"]])
      }
      # Split up indices
      if (!is.null(BS.control[["nsplit"]]) && BS.control[["nsplit"]] > 1) {
        Which = split(Which, cut(seq_along(Which), BS.control[["nsplit"]]))
      }
      Which = Which[sapply(Which, FUN = length) > 0]
      if (length(Which) == 0)
        Which = NULL
      # Repeat SD with indexing
      message(paste0("Bias correcting ", length(Which), " derived quantities"))
      parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, par.fixed = parameter_estimates$par, hessian.fixed = h, bias.correct = TRUE,
                                                  bias.correct.control = list(sd = BS.control[["sd"]], split = Which, nsplit = NULL), ...)
    }
    # Update
    parameter_estimates[["Convergence_check"]] = ifelse(parameter_estimates$SD$pdHess == TRUE, parameter_estimates[["Convergence_check"]],
                                                        "The model is definitely not converged")
    parameter_estimates[["time_for_sdreport"]] = Sys.time() - sd_time
  }
  parameter_estimates[["time_for_run"]] = Sys.time() - start_time

  # Save results
  if (!is.null(savedir)) {
    save(parameter_estimates, file = file.path(savedir, "parameter_estimates.RData"))
    capture.output(parameter_estimates, file = file.path(savedir, "parameter_estimates.txt"))
  }

  # Print warning to screen
  if (parameter_estimates[["Convergence_check"]] != "There is no evidence that the model is not converged") {
    message("#########################")
    message(parameter_estimates[["Convergence_check"]])
    message("#########################")
  }

  # Return stuff
  return(parameter_estimates)
}
