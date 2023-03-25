##' Specify the stock-recruit relationship (SRR) for Rceattle
##'
##' @param srr_fun Stock recruit function to be used for Rceattle (see below). Default = 0
##' @param proj_mean_rec Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = SRR(omega, rec_devs)
##' @param srr_use_prior TRUE or FALSE to use normally distributed prior on stock recruit parameter.
##' @param srr_prior_mean mean for normally distributed prior for stock-recruit parameter
##' @param srr_prior_sd Prior standard deviation for stock-recruit parameter
##' @details
##'
##' **Stock recruitment relationships currently implemented in Rceattle:**
##'
##' \code{srr_fun = 0} No stock recruit relationship. Recruitment is a function of \deqn{R_y = exp(R0 + R_{dev,y})}. Sensu Alaska, steepness = 0.99.
##'
##' \code{srr_fun = 1} Beverton-holt stock-recruitment relationship
##'   \deqn{R_y = \frac{\alpha * SB_{y-minage}}{1+\beta * SB_{y-minage}}}. Prior is on alpha.
##'
##' \code{srr_fun = 2} Beverton-holt stock-recruitment relationship with environmental covariates impacting larval survival rate
##'   \deqn{R_y = \frac{\alpha * e^{X * \Beta} * SB_{y-minage}}{1+\beta * SB_{y-minage}}}. Prior is on alpha.
##'
##' \code{srr_fun = 3} Ricker stock-recruitment relationship.
##'   \deqn{R_y = \alpha * SB_{y-minage}} * exp(-beta * SB_{y-minage})}. Prior is on alpha.
##'
##'
##' @return A \code{list} containing the stock recruitment relationship settings
##' @export
##'
build_srr <- function(srr_fun = 0,
                      proj_mean_rec = TRUE,
                      srr_use_prior = FALSE,
                      srr_prior_mean = 0.40,
                      srr_prior_sd = 0.35){
  list(srr_fun = srr_fun,
       proj_mean_rec = proj_mean_rec,
       srr_use_prior = srr_use_prior,
       srr_prior_mean = srr_prior_mean,
       srr_prior_sd = srr_prior_sd
       )
}
