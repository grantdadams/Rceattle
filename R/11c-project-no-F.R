#' Function to update hindcast and set F to 0. Useful for determining dynamic reference points for multi-species models under climate-change
#'
#' @param Rceattle
#'
#' @return
#' @export
#'
#' @examples
remove_F <- function(Rceattle){

  # * Years for F = 0 ----
  # - don't want hindcast or it will bias suitability in Multi-species models
  proj_years <- (Rceattle$data_list$suit_endyr+1):Rceattle$data_list$projyr - Rceattle$data_list$styr + 1
  fdevs_cols <- 1:ncol(Rceattle$estimated_params$ln_F)
  fdevs_change <- which(fdevs_cols %in% proj_years)

  # * Set F to 0 ----
  Rceattle$estimated_params$ln_F[,fdevs_change] <- replace(Rceattle$estimated_params$ln_F[,fdevs_change], values = -999)

  # * Update fit ----
  estMode <- Rceattle$data_list$estimateMode
  Rceattle <- fit_mod(
    data_list = Rceattle$data_list,
    inits = Rceattle$estimated_params,
    map =  NULL,
    bounds = NULL,
    file = NULL,
    estimateMode = 3,
    HCR = build_hcr(HCR = Rceattle$data_list$HCR, # Tier3 HCR
                    DynamicHCR = Rceattle$data_list$DynamicHCR,
                    Ftarget = Rceattle$data_list$Ftarget,
                    Flimit = Rceattle$data_list$Flimit,
                    Ptarget = Rceattle$data_list$Ptarget,
                    Plimit = Rceattle$data_list$Plimit,
                    Alpha = Rceattle$data_list$Alpha,
                    Pstar = Rceattle$data_list$Pstar,
                    Sigma = Rceattle$data_list$Sigma,
                    Fmult = Rceattle$data_list$Fmult,
                    HCRorder = Rceattle$data_list$HCRorder
    ),
    recFun = build_srr(srr_fun = Rceattle$data_list$srr_fun,
                       srr_pred_fun  = Rceattle$data_list$srr_pred_fun ,
                       proj_mean_rec  = Rceattle$data_list$proj_mean_rec ,
                       srr_meanyr = min(Rceattle$data_list$srr_meanyr, Rceattle$data_list$endyr), # Update end year if less than srr_meanyr
                       srr_hat_styr = Rceattle$data_list$srr_hat_styr,
                       srr_hat_endyr = Rceattle$data_list$srr_hat_endyr,
                       srr_est_mode  = Rceattle$data_list$srr_est_mode ,
                       srr_prior  = Rceattle$data_list$srr_prior,
                       srr_prior_sd   = Rceattle$data_list$srr_prior_sd,
                       Bmsy_lim = Rceattle$data_list$Bmsy_lim,
                       srr_env_indices = Rceattle$data_list$srr_env_indices),
    M1Fun =     build_M1(M1_model= Rceattle$data_list$M1_model,
                         updateM1 = FALSE,
                         M1_use_prior = Rceattle$data_list$M1_use_prior,
                         M2_use_prior = Rceattle$data_list$M2_use_prior,
                         M_prior = Rceattle$data_list$M_prior,
                         M_prior_sd = Rceattle$data_list$M_prior_sd),
    random_rec = Rceattle$data_list$random_rec,
    niter = Rceattle$data_list$niter,
    msmMode = Rceattle$data_list$msmMode,
    avgnMode = Rceattle$data_list$avgnMode,
    suitMode = Rceattle$data_list$suitMode,
    suit_styr = Rceattle$data_list$suit_styr,
    suit_endyr = min(Rceattle$data_list$suit_endyr, Rceattle$data_list$endyr),   # Update to end year if less than suit_endyr
    initMode = Rceattle$data_list$initMode,
    phase = FALSE,
    loopnum = Rceattle$data_list$loopnum,
    getsd = TRUE,
    verbose = 0)

  Rceattle$data_list$estimateMode <- estMode
  return(Rceattle)

}

