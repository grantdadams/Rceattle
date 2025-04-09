#' Updates a model projection with a recruitment trend
#'
#' @param Rceattle Rceattle model
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier \code{mean rec * (1 + (rec_trend/projection years) * 1:projection years)}. Can be of length 1 or of length nspp. If length 1, all species get the same trend.
#'
#' @return An updated Rceattle model where the projection has a linear recruitment trend
#' @export
#'
#' @examples
project_trend <- function(Rceattle, rec_trend = 0, sample_rec = FALSE){

  # - Years for projection
  hind_yrs <- (Rceattle$data_list$styr) : Rceattle$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  proj_yrs <- (Rceattle$data_list$endyr + 1) : Rceattle$data_list$projyr
  proj_nyrs <- length(proj_yrs)

  # - Adjust rec trend
  if(length(rec_trend)==1){
    rec_trend = rep(rec_trend, Rceattle$data_list$nspp)
  }

  # - Replace future rec devs
  for(sp in 1:Rceattle$data_list$nspp){

    # -- Projection uses mean R of hindcast
    if(Rceattle$data_list$proj_mean_rec == 1){
      rec_dev <- log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
    }

    # -- Projection uses exp(ln_R0)
    if(Rceattle$data_list$proj_mean_rec == 0){
      rec_dev <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs]) * (1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs))  - Rceattle$estimated_params$ln_mean_rec[sp] # - Scale mean rec for rec trend
    }

    # - Update OM with devs
    Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1] <- replace(
      Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1],
      values =  rec_dev)
  }

  # - Update
  estMode <- Rceattle$data_list$estimateMode
  Rceattle <- fit_mod(
    data_list = Rceattle$data_list,
    inits = Rceattle$estimated_params,
    map =  Rceattle$map,
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
                    Sigma = Rceattle$data_list$Sigma
    ),
    recFun = build_srr(srr_fun = Rceattle$data_list$srr_fun,
                       srr_pred_fun  = Rceattle$data_list$srr_pred_fun ,
                       proj_mean_rec  = Rceattle$data_list$proj_mean_rec ,
                       srr_est_mode  = Rceattle$data_list$srr_est_mode ,
                       srr_prior_mean  = Rceattle$data_list$srr_prior_mean,
                       srr_prior_sd   = Rceattle$data_list$srr_prior_sd ),
    M1Fun =     build_M1(M1_model= Rceattle$data_list$M1_model,
                         updateM1 = FALSE,
                         M1_use_prior = Rceattle$data_list$M1_use_prior,
                         M2_use_prior = Rceattle$data_list$M2_use_prior,
                         M1_prior_mean = Rceattle$data_list$M1_prior_mean,
                         M1_prior_sd = Rceattle$data_list$M1_prior_sd),
    random_rec = Rceattle$data_list$random_rec,
    niter = Rceattle$data_list$niter,
    msmMode = Rceattle$data_list$msmMode,
    avgnMode = Rceattle$data_list$avgnMode,
    minNByage = Rceattle$data_list$minNByage,
    suitMode = Rceattle$data_list$suitMode,
    initMode = Rceattle$data_list$initMode,
    phase = NULL,
    loopnum = 3,
    getsd = FALSE,
    verbose = 0)

  Rceattle$data_list$estimateMode <- estMode
  return(Rceattle)
}


#' Function to update hindcast removing, setting F to 0 and rec to mean recruitment
#'
#' @param Rceattle
#'
#' @return
#' @export
#'
#' @examples
remove_rec_dev_and_F <- function(Rceattle, rec_trend = 0){

  # - Years for projection
  hind_yrs <- (Rceattle$data_list$styr) : Rceattle$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  proj_yrs <- (Rceattle$data_list$endyr + 1) : Rceattle$data_list$projyr
  proj_nyrs <- length(proj_yrs)

  Rceattle$data_list$proj_mean_rec <- FALSE

  # - Adjust rec trend
  if(length(rec_trend)==1){
    rec_trend = rep(rec_trend, Rceattle$data_list$nspp)
  }

  # - Replace hindcast rec devs
  for(sp in 1:Rceattle$data_list$nspp){
    ## HINDCAST
    # Get dev from mean rec to ln_R0 for hindcast
    rec_dev <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs])) - Rceattle$estimated_params$ln_mean_rec[sp] # - Scale mean rec for rec trend

    # - Update OM with devs
    Rceattle$estimated_params$rec_dev[sp,1:hind_nyrs] <- replace(
      Rceattle$estimated_params$rec_dev[sp,1:hind_nyrs],
      values =  rec_dev)

    ## PROJECTION
    # -- Projection uses exp(ln_R0)
    rec_dev <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs]) * (1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs))  - Rceattle$estimated_params$ln_mean_rec[sp] # - Scale mean rec for rec trend


    # - Update OM with devs
    Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1] <- replace(
      Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1],
      values =  rec_dev)
  }

  # Set fishing to 0
  Rceattle$estimated_params$ln_F <- replace(Rceattle$estimated_params$ln_F, values = -999)

  # - Update
  estMode <- Rceattle$data_list$estimateMode
  Rceattle <- fit_mod(
    data_list = Rceattle$data_list,
    inits = Rceattle$estimated_params,
    map =  Rceattle$map,
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
                    Sigma = Rceattle$data_list$Sigma
    ),
    recFun = build_srr(srr_fun = Rceattle$data_list$srr_fun,
                       srr_pred_fun  = Rceattle$data_list$srr_pred_fun ,
                       proj_mean_rec  = Rceattle$data_list$proj_mean_rec ,
                       srr_est_mode  = Rceattle$data_list$srr_est_mode ,
                       srr_prior_mean  = Rceattle$data_list$srr_prior_mean,
                       srr_prior_sd   = Rceattle$data_list$srr_prior_sd ),
    M1Fun =     build_M1(M1_model= Rceattle$data_list$M1_model,
                         updateM1 = FALSE,
                         M1_use_prior = Rceattle$data_list$M1_use_prior,
                         M2_use_prior = Rceattle$data_list$M2_use_prior,
                         M1_prior_mean = Rceattle$data_list$M1_prior_mean,
                         M1_prior_sd = Rceattle$data_list$M1_prior_sd),
    random_rec = Rceattle$data_list$random_rec,
    niter = Rceattle$data_list$niter,
    msmMode = Rceattle$data_list$msmMode,
    avgnMode = Rceattle$data_list$avgnMode,
    minNByage = Rceattle$data_list$minNByage,
    suitMode = Rceattle$data_list$suitMode,
    initMode = Rceattle$data_list$initMode,
    phase = NULL,
    loopnum = 3,
    getsd = FALSE,
    verbose = 0)

  Rceattle$data_list$estimateMode <- estMode
  return(Rceattle)

}



#' Function to project the model until proj yr under no rec, then again under n proj years given rec trend
#'
#' @param Rceattle
#'
#' @return
#' @export
#'
#' @examples
equilibrate_and_project <- function(Rceattle, rec_trend = 0){

  # - Years for projection
  hind_yrs <- (Rceattle$data_list$styr) : Rceattle$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  proj_yrs <- (Rceattle$data_list$endyr + 1) : Rceattle$data_list$projyr
  proj_nyrs <- length(proj_yrs)

  new_projyr <- Rceattle$data_list$projyr + proj_nyrs
  new_projyrs <- (Rceattle$data_list$projyr + 1):new_projyr

  Rceattle$data_list$projyr <- new_projyr

  Rceattle$data_list$proj_mean_rec <- FALSE

  # - Adjust rec trend
  if(length(rec_trend)==1){
    rec_trend = rep(rec_trend, Rceattle$data_list$nspp)
  }

  # Lengthen rec devs and F paramers
  Rceattle$estimated_params$rec_dev <- cbind(Rceattle$estimated_params$rec_dev, matrix(NA, nrow = Rceattle$data_list$nspp, ncol = proj_nyrs))
  Rceattle$estimated_params$ln_Flimit <- cbind(Rceattle$estimated_params$ln_Flimit, matrix(Rceattle$estimated_params$ln_Flimit[,ncol(Rceattle$estimated_params$ln_Flimit)], nrow = Rceattle$data_list$nspp, ncol = proj_nyrs))
  Rceattle$estimated_params$ln_Ftarget <- cbind(Rceattle$estimated_params$rec_dev, matrix(Rceattle$estimated_params$ln_Ftarget[,ncol(Rceattle$estimated_params$ln_Ftarget)], nrow = Rceattle$data_list$nspp, ncol = proj_nyrs))

  # - Replace hindcast rec devs
  for(sp in 1:Rceattle$data_list$nspp){

    ## PROJECTION
    # -- Projection uses exp(ln_R0)
    rec_dev_old_proj <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs]) * (1+(0/proj_nyrs) * 1:proj_nyrs))  - Rceattle$estimated_params$ln_mean_rec[sp] # - Scale mean rec for rec trend

    rec_dev_new_proj <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs]) * (1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs))  - Rceattle$estimated_params$ln_mean_rec[sp] # - Scale mean rec for rec trend


    # - Update OM with devs
    Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1] <- replace(
      Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1],
      values =  rec_dev_old_proj)

    Rceattle$estimated_params$rec_dev[sp,new_projyrs - Rceattle$data_list$styr + 1] <- replace(
      Rceattle$estimated_params$rec_dev[sp,new_projyrs - Rceattle$data_list$styr + 1],
      values =  rec_dev_new_proj)
  }

  # - Update
  estMode <- Rceattle$data_list$estimateMode
  Rceattle <- fit_mod(
    data_list = Rceattle$data_list,
    inits = Rceattle$estimated_params,
    map =  Rceattle$map,
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
                    Sigma = Rceattle$data_list$Sigma
    ),
    recFun = build_srr(srr_fun = Rceattle$data_list$srr_fun,
                       srr_pred_fun  = Rceattle$data_list$srr_pred_fun ,
                       proj_mean_rec  = Rceattle$data_list$proj_mean_rec ,
                       srr_est_mode  = Rceattle$data_list$srr_est_mode ,
                       srr_prior_mean  = Rceattle$data_list$srr_prior_mean,
                       srr_prior_sd   = Rceattle$data_list$srr_prior_sd ),
    M1Fun =     build_M1(M1_model= Rceattle$data_list$M1_model,
                         updateM1 = FALSE,
                         M1_use_prior = Rceattle$data_list$M1_use_prior,
                         M2_use_prior = Rceattle$data_list$M2_use_prior,
                         M1_prior_mean = Rceattle$data_list$M1_prior_mean,
                         M1_prior_sd = Rceattle$data_list$M1_prior_sd),
    random_rec = Rceattle$data_list$random_rec,
    niter = Rceattle$data_list$niter,
    msmMode = Rceattle$data_list$msmMode,
    avgnMode = Rceattle$data_list$avgnMode,
    minNByage = Rceattle$data_list$minNByage,
    suitMode = Rceattle$data_list$suitMode,
    initMode = Rceattle$data_list$initMode,
    phase = NULL,
    loopnum = 3,
    getsd = FALSE,
    verbose = 0)

  Rceattle$data_list$estimateMode <- estMode
  return(Rceattle)

}
