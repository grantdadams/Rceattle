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
    estimateMode = 3, # Update quantities, not parameters
    random_rec = Rceattle$data_list$random_rec,
    niter = Rceattle$data_list$niter,
    msmMode = Rceattle$data_list$msmMode,
    avgnMode = Rceattle$data_list$avgnMode,
    minNByage = Rceattle$data_list$minNByage,
    suitMode = Rceattle$data_list$suitMode,
    meanyr = Rceattle$data_list$endyr,
    updateM1 = FALSE, # Dont update M1 from data, fix at previous parameters
    loopnum = 2,
    phase = NULL,
    getsd = FALSE,
    verbose = 0)

  Rceattle$data_list$estimateMode <- estMode
  return(Rceattle)
}
