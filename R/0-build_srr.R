##' Specify the stock-recruit relationship (SRR) for Rceattle
##'
##' @param srr_fun Stock recruit function to be used for hindcast estimation of Rceattle (see below). Default = 0
##' @param srr_pred_fun stock recruit function for projection, reference points, and penalties to be used for Rceattle (see below) when \code{srr_fun == 0}. It treats treat the stock-recruit curve as an additional penalty onto the annualy estimated recruitment from the hindcast (sensu AMAK and Jim Ianelli's pollock model). If \code{srr_fun > 0} then \code{srr_pred_fun <- srr_fun} and no additional penalty is included.
##' @param proj_mean_rec Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = SRR(omega, rec_devs)
##' @param srr_meanyr Integer. The last year used to calculate mean recruitment, starting at \code{styr}. Defaults to $endyr$ in $data_list$. Used for MSE runs where mean recruitment is held at the value estimated from the years used to condition the OM, but F is estimated for years beyond those used to condition the OM to account for projected catch.
##' @param R_hat_yr Integer. The last year used to calculate predicted recruitment for estimating recruitment functions as additional penalties sensu AMAK and Jim Ianelli's pollock model, starting at \code{styr}. Defaults to $endyr$ in $data_list$. Useful if environmental data used to condition stock-recruit relationships is not available until end-year, but projections are desired.
##' @param srr_est_mode Switch to determine estimation mode. 0 = fix alpha to prior mean, 1 = freely estimate alpha and beta, 2 = use lognormally distributed prior for alpha.
##' @param srr_prior_mean mean for normally distributed prior for stock-recruit parameter
##' @param srr_prior_sd Prior standard deviation for stock-recruit parameter
##' @param srr_env_indices vector or single index indicating the columns of \code{env_data} to use in a environmentally driven stock recruit curve.
##' @param Bmsy_lim Upper limit for ricker based SSB-MSY (e.g 1/Beta). Will add a likelihood penalty if beta is estimated above this limit.
##' @description
##'
##' **Stock recruitment relationships currently implemented in Rceattle:**
##'
##' \code{srr_fun = 0} No stock recruit relationship. Recruitment is a function of R0 and annual deviates (i.e. steepness = 0.99).
##'  \deqn{R_y = exp(R0 + R_{dev,y})}
##'
##' \code{srr_fun = 1} Environmentally driven recruitment without stock recruit relationship
##'  \deqn{R_y = exp(R0 + R_{dev,y} + X * \beta_X)}
##'
##' \code{srr_fun = 2} Beverton-holt stock-recruitment relationship
##'   \deqn{R_y = \frac{\alpha_{srr} * SB_{y-minage}}{1+\beta_{srr} * SB_{y-minage}}}
##'
##' \code{srr_fun = 3} Beverton-holt stock-recruitment relationship with environmental covariates impacting larval survival rate and prior is on alpha.
##'   \deqn{R_y = \frac{\alpha_{srr} * e^{X * \beta_X} * SB_{y-minage}}{1+\beta_{srr} * SB_{y-minage}}}
##'
##' \code{srr_fun = 4} Ricker stock-recruitment relationship
##'   \deqn{R_y = \alpha_{srr} * SB_{y-minage} * exp(-\beta_{srr} * SB_{y-minage})}
##'
##' \code{srr_fun = 5} Ricker stock-recruitment relationship with environmental covariates impacting larval survival rate and prior is on alpha.
##'   \deqn{R_y = \alpha_{srr} e^{X * \beta_X} * SB_{y-minage} * exp(-\beta_{srr} * SB_{y-minage})}
##'
##'
##' @return A \code{list} containing the stock recruitment relationship settings
##' @export
##'
build_srr <- function(srr_fun = 0,
                      srr_pred_fun = srr_fun,
                      proj_mean_rec = TRUE,
                      srr_meanyr = NULL,
                      R_hat_yr = NULL,
                      srr_est_mode = 1,
                      srr_prior_mean = 4,
                      srr_prior_sd = 1,
                      srr_env_indices = 1,
                      Bmsy_lim = -999){

  # Set pred/RP/penalty to same as SR curve if SR fun > 0
  if(srr_fun > 0){
    srr_pred_fun = srr_fun
  }

  if(!srr_pred_fun %in% c(3,4)){
    Bmsy_lim = -999
  }

  list(srr_fun = srr_fun,
       srr_pred_fun = srr_pred_fun,
       proj_mean_rec = proj_mean_rec,
       srr_meanyr = srr_meanyr,
       R_hat_yr = R_hat_yr,
       srr_est_mode = srr_est_mode,
       srr_prior_mean = srr_prior_mean,
       srr_prior_sd = srr_prior_sd,
       srr_env_indices = srr_env_indices,
       Bmsy_lim = Bmsy_lim
  )
}



#' Define M1 specifications
#'
#' @param M1_model M1 set-up. 0 = use fixed natural mortality from M1_base in data, 1 = estimate sex- and age-invariant M1, 2 = sex-specific (two-sex model), age-invariant M1, 3 =   estimate sex- and age-specific M1.
#' @param updateM1 If using initial parameters, use M1 from data instead.
#' @param M1_use_prior Have M1 come from a lognormal prior
#' @param M2_use_prior Include M2 in prior for multi-species models
#' @param M1_prior_mean Mean of M prior on natural scale
#' @param M1_prior_sd SD of lognormal M prior.
#'
#' @return A list of switches for defining the M1 model
#' @export
#'
build_M1 <- function(M1_model = 0, #FIXME est_M1 from data
                     updateM1 = FALSE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE,
                     M1_prior_mean = 0.40,
                     M1_prior_sd = 0.35){
  list(
    M1_model= M1_model,
    updateM1 = updateM1,
    M1_use_prior = M1_use_prior,
    M2_use_prior = M2_use_prior,
    M1_prior_mean = M1_prior_mean,
    M1_prior_sd = M1_prior_sd
  )
}
