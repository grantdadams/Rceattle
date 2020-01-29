#' Run a management strategy evaluation
#'
#' @param operating_model CEATTLE model object exported from \code{\link{Rceattle}}
#' @param estimation_model CEATTLE model object exported from \code{\link{Rceattle}}
#' @param nsim Number of simulations to run (default 10)
#' @param assessment_period Period of years that each assessment is taken
#' @param sampling_period Period of years data sampling is conducted
#' @param simulate Include simulated random error proportional to that estimated/provided.
#' @param cap A cap on the catch in the projection. Can be a single number or vector. Default = NULL
#'
#' @return A list of
#' @export
#'
#' @examples
mse_run <- function(operating_model = ss_run, estimation_model = ss_run, nsim = 10, assessment_period = 1, sampling_period = 1, simulate = TRUE, cap = NULL, seed = 666){
  '%!in%' <- function(x,y)!('%in%'(x,y))

  set.seed(seed)

  Rceattle_OM_list <- list()
  Rceattle_EM_list <- list()

  # Update data-files in OM so we can fill in updated years
  # -- srv_biom
  for(flt in (unique(operating_model$data_list$srv_biom$Fleet_code))){
    srv_biom_sub <- operating_model$data_list$srv_biom[which(operating_model$data_list$srv_biom$Fleet_code == flt),]
    yrs_proj <- (operating_model$data_list$endyr + 1):operating_model$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_srv_biom <- data.frame(Fleet_name = rep(srv_biom_sub$Fleet_name[1], nyrs_proj),
                                Fleet_code = rep(flt, nyrs_proj),
                                Species = rep(srv_biom_sub$Species[1], nyrs_proj),
                                Year = -yrs_proj,
                                Month = rep(srv_biom_sub$Month[length(srv_biom_sub$Month)], nyrs_proj),
                                Selectivity_block = rep(srv_biom_sub$Selectivity_block[length(srv_biom_sub$Selectivity_block)], nyrs_proj),
                                Q_block = rep(srv_biom_sub$Selectivity_block[length(srv_biom_sub$Selectivity_block)], nyrs_proj),
                                Observation = rep(NA, nyrs_proj),
                                CV = rep(mean(srv_biom_sub$CV, na.rm = TRUE), nyrs_proj)) # FIXME: not sure what the best is to put here, using mean
    operating_model$data_list$srv_biom <- rbind(operating_model$data_list$srv_biom, proj_srv_biom)
  }
  operating_model$data_list$srv_biom <- operating_model$data_list$srv_biom[
    with(operating_model$data_list$srv_biom, order(Fleet_code, Year)),]


  # -- comp_data
  for(flt in (unique(operating_model$data_list$comp_data$Fleet_code))){
    comp_data_sub <- operating_model$data_list$comp_data[which(operating_model$data_list$comp_data$Fleet_code == flt),]
    yrs_proj <- (operating_model$data_list$endyr + 1):operating_model$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_comp_data <- data.frame(Fleet_name = rep(comp_data_sub$Fleet_name[1], nyrs_proj),
                                 Fleet_code = rep(flt, nyrs_proj),
                                 Species = rep(comp_data_sub$Species[1], nyrs_proj),
                                 Sex = rep(comp_data_sub$Sex[length(comp_data_sub$Sex)], nyrs_proj),
                                 Age0_Length1 = rep(comp_data_sub$Age0_Length1[length(comp_data_sub$Age0_Length1)], nyrs_proj),
                                 Year = -yrs_proj,
                                 Month = rep(comp_data_sub$Month[length(comp_data_sub$Month)], nyrs_proj),
                                 Sample_size = rep(comp_data_sub$Sample_size[length(comp_data_sub$Sample_size)], nyrs_proj))
    proj_comps <- data.frame(matrix(0, ncol = ncol(comp_data_sub) - ncol(proj_comp_data), nrow = nyrs_proj))
    colnames(proj_comps) <-paste0("Comp_", 1:ncol(proj_comps))
    proj_comp_data <-cbind(proj_comp_data, proj_comps)

    operating_model$data_list$comp_data <- rbind(operating_model$data_list$comp_data, proj_comp_data)
  }
  operating_model$data_list$comp_data <- operating_model$data_list$comp_data[
    with(operating_model$data_list$comp_data, order(Fleet_code, Year)),]

  # -- emp_sel - Use terminal year
  for(flt in (unique(operating_model$data_list$emp_sel$Fleet_code))){
    emp_sel_sub <- operating_model$data_list$emp_sel[which(operating_model$data_list$emp_sel$Fleet_code == flt),]
    yrs_proj <- (operating_model$data_list$endyr + 1):operating_model$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_emp_sel <- data.frame(Fleet_name = rep(emp_sel_sub$Fleet_name[1], nyrs_proj),
                               Fleet_code = rep(flt, nyrs_proj),
                               Species = rep(emp_sel_sub$Species[1], nyrs_proj),
                               Sex = rep(emp_sel_sub$Sex[length(emp_sel_sub$Sex)], nyrs_proj),
                               Year = yrs_proj)
    proj_comps <- data.frame(matrix(emp_sel_sub[nrow(emp_sel_sub), (ncol(proj_emp_sel) + 1) : ncol(emp_sel_sub)], ncol = ncol(emp_sel_sub) - ncol(proj_emp_sel), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Comp_", 1:ncol(proj_comps))
    proj_emp_sel <-cbind(proj_emp_sel, proj_comps)

    operating_model$data_list$emp_sel <- rbind(operating_model$data_list$emp_sel, proj_emp_sel)
  }
  operating_model$data_list$emp_sel <- operating_model$data_list$emp_sel[
    with(operating_model$data_list$emp_sel, order(Fleet_code, Year)),]

  # -- wt
  for(flt in (unique(operating_model$data_list$wt$Wt_index))){
    wt_sub <- operating_model$data_list$wt[which(operating_model$data_list$wt$Wt_index == flt),]
    yrs_proj <- (operating_model$data_list$endyr + 1):operating_model$data_list$projyr
    yrs_proj <- yrs_proj[which(yrs_proj %!in% wt_sub$Year)]
    nyrs_proj <- length(yrs_proj)
    proj_wt <- data.frame(Wt_name = rep(wt_sub$Wt_name[1], nyrs_proj),
                          Wt_index = rep(flt, nyrs_proj),
                          Species = rep(wt_sub$Species[1], nyrs_proj),
                          Year = yrs_proj,
                          Sex = rep(wt_sub$Sex[length(wt_sub$Sex)], nyrs_proj)
    )
    proj_comps <- data.frame(matrix(wt_sub[nrow(wt_sub), (ncol(proj_wt)+1) : ncol(wt_sub)], ncol = ncol(wt_sub) - ncol(proj_wt), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_wt <-cbind(proj_wt, proj_comps)
    colnames(proj_wt) <- colnames(operating_model$data_list$wt)
    operating_model$data_list$wt <- rbind(operating_model$data_list$wt, proj_wt)
  }
  operating_model$data_list$wt <- operating_model$data_list$wt[
    with(operating_model$data_list$wt, order(Wt_index, Year)),]


  # -- Pyrs
  for(flt in (unique(operating_model$data_list$wt$Species))){
    Pyrs_sub <- operating_model$data_list$Pyrs[which(operating_model$data_list$Pyrs$Species == flt),]
    yrs_proj <- (operating_model$data_list$endyr + 1):operating_model$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_Pyrs <- data.frame(Species = rep(Pyrs_sub$Species[1], nyrs_proj),
                            Sex = rep(Pyrs_sub$Sex[length(Pyrs_sub$Sex)], nyrs_proj),
                            Year = yrs_proj)

    proj_comps <- data.frame(matrix(Pyrs_sub[nrow(Pyrs_sub), (ncol(proj_Pyrs)+1) : ncol(Pyrs_sub)], ncol = ncol(Pyrs_sub) - ncol(proj_Pyrs), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_Pyrs <-cbind(proj_Pyrs, proj_comps)
    colnames(proj_Pyrs) <- colnames(operating_model$data_list$Pyrs)
    operating_model$data_list$Pyrs <- rbind(operating_model$data_list$Pyrs, proj_Pyrs)
  }
  operating_model$data_list$Pyrs <- operating_model$data_list$Pyrs[
    with(operating_model$data_list$Pyrs, order(Species, Year)),]


  # Update emp sel in EM
  # -- emp_sel - Use terminal year
  for(flt in (unique(estimation_model$data_list$emp_sel$Fleet_code))){
    emp_sel_sub <- estimation_model$data_list$emp_sel[which(estimation_model$data_list$emp_sel$Fleet_code == flt),]
    yrs_proj <- (estimation_model$data_list$endyr + 1):estimation_model$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_emp_sel <- data.frame(Fleet_name = rep(emp_sel_sub$Fleet_name[1], nyrs_proj),
                               Fleet_code = rep(flt, nyrs_proj),
                               Species = rep(emp_sel_sub$Species[1], nyrs_proj),
                               Sex = rep(emp_sel_sub$Sex[length(emp_sel_sub$Sex)], nyrs_proj),
                               Year = yrs_proj)
    proj_comps <- data.frame(matrix(emp_sel_sub[nrow(emp_sel_sub), (ncol(proj_emp_sel) + 1) : ncol(emp_sel_sub)], ncol = ncol(emp_sel_sub) - ncol(proj_emp_sel), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Comp_", 1:ncol(proj_comps))
    proj_emp_sel <-cbind(proj_emp_sel, proj_comps)

    estimation_model$data_list$emp_sel <- rbind(estimation_model$data_list$emp_sel, proj_emp_sel)
  }
  estimation_model$data_list$emp_sel <- estimation_model$data_list$emp_sel[
    with(estimation_model$data_list$emp_sel, order(Fleet_code, Year)),]


  # -- wt
  for(flt in (unique(estimation_model$data_list$wt$Wt_index))){
    wt_sub <- estimation_model$data_list$wt[which(estimation_model$data_list$wt$Wt_index == flt),]
    yrs_proj <- (estimation_model$data_list$endyr + 1):estimation_model$data_list$projyr
    yrs_proj <- yrs_proj[which(yrs_proj %!in% wt_sub$Year)]
    nyrs_proj <- length(yrs_proj)
    proj_wt <- data.frame(Wt_name = rep(wt_sub$Wt_name[1], nyrs_proj),
                          Wt_index = rep(flt, nyrs_proj),
                          Species = rep(wt_sub$Species[1], nyrs_proj),
                          Year = yrs_proj,
                          Sex = rep(wt_sub$Sex[length(wt_sub$Sex)], nyrs_proj)
    )
    proj_comps <- data.frame(matrix(wt_sub[nrow(wt_sub), (ncol(proj_wt)+1) : ncol(wt_sub)], ncol = ncol(wt_sub) - ncol(proj_wt), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_wt <-cbind(proj_wt, proj_comps)
    colnames(proj_wt) <- colnames(operating_model$data_list$wt)
    estimation_model$data_list$wt <- rbind(estimation_model$data_list$wt, proj_wt)
  }
  estimation_model$data_list$wt <- estimation_model$data_list$wt[
    with(estimation_model$data_list$wt, order(Wt_index, Year)),]


  # -- Pyrs
  for(flt in (unique(estimation_model$data_list$wt$Species))){
    Pyrs_sub <- estimation_model$data_list$Pyrs[which(estimation_model$data_list$Pyrs$Species == flt),]
    yrs_proj <- (estimation_model$data_list$endyr + 1):estimation_model$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_Pyrs <- data.frame(Species = rep(Pyrs_sub$Species[1], nyrs_proj),
                            Sex = rep(Pyrs_sub$Sex[length(Pyrs_sub$Sex)], nyrs_proj),
                            Year = yrs_proj)

    proj_comps <- data.frame(matrix(Pyrs_sub[nrow(Pyrs_sub), (ncol(proj_Pyrs)+1) : ncol(Pyrs_sub)], ncol = ncol(Pyrs_sub) - ncol(proj_Pyrs), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_Pyrs <-cbind(proj_Pyrs, proj_comps)
    colnames(proj_Pyrs) <- colnames(operating_model$data_list$Pyrs)
    estimation_model$data_list$Pyrs <- rbind(estimation_model$data_list$Pyrs, proj_Pyrs)
  }
  estimation_model$data_list$Pyrs <- estimation_model$data_list$Pyrs[
    with(estimation_model$data_list$Pyrs, order(Species, Year)),]



  # Years for simulations
  proj_yrs <- (estimation_model$data_list$endyr + 1) : estimation_model$data_list$projyr
  proj_nyrs <- length(proj_yrs)
  assess_yrs <- seq(from = estimation_model$data_list$endyr + assessment_period, to = estimation_model$data_list$projyr,  by = assessment_period)
  sample_yrs <- seq(from = estimation_model$data_list$endyr + sampling_period, to = estimation_model$data_list$projyr,  by = sampling_period)


  # Do the MSE
  for(sim in 1:nsim){

    # Set models
    Rceattle_EM_list[[sim]] <- list()
    Rceattle_EM_list[[sim]][[1]] <- estimation_model
    estimation_model_use <- estimation_model
    operating_model_use <- operating_model



    # Replace future rec_devs with numbers
    if(simulate){
      for(sp in 1:operating_model_use$data_list$nspp){
        operating_model_use$estimated_params$rec_dev[sp,proj_yrs - operating_model_use$data_list$styr + 1] <- replace(
          operating_model_use$estimated_params$rec_dev[sp,proj_yrs - operating_model_use$data_list$styr + 1],
          values = rnorm( length(operating_model_use$estimated_params$rec_dev[sp,proj_yrs - operating_model_use$data_list$styr + 1]),
                          mean = 0,
                          sd = exp(operating_model_use$estimated_params$ln_rec_sigma[sp])) # Assumed value from penalized likelihood
        )
      }
    }



    # Run through model
    for(k in 1:length(assess_yrs)){

      # Get projected catch data from EM
      new_catch_data <- estimation_model_use$data_list$fsh_biom
      new_years <- proj_yrs[which(proj_yrs <= assess_yrs[k] & proj_yrs > operating_model_use$data_list$endyr)]
      dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
      new_catch_data$Catch[dat_fill_ind] <- estimation_model_use$quantities$fsh_bio_hat[dat_fill_ind]
      if(!is.null(cap)){
        new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap, cap, new_catch_data$Catch[dat_fill_ind])
      }

      # Update catch data in OM and EM
      operating_model_use$data_list$fsh_biom <- new_catch_data
      estimation_model_use$data_list$fsh_biom <- new_catch_data

      # Update endyr of OM
      nyrs_hind <- operating_model_use$data_list$endyr - operating_model_use$data_list$styr + 1
      operating_model_use$data_list$endyr <- assess_yrs[k]

      # Update parameters
      # -- F_dev
      operating_model_use$estimated_params$F_dev <- cbind(operating_model_use$estimated_params$F_dev, matrix(0, nrow= nrow(operating_model_use$estimated_params$F_dev), ncol = length(new_years)))

      # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
      operating_model_use$estimated_params$ln_srv_q_dev <- cbind(operating_model_use$estimated_params$ln_srv_q_dev, matrix(operating_model_use$estimated_params$ln_srv_q_dev[,ncol(operating_model_use$estimated_params$ln_srv_q_dev)], nrow= nrow(operating_model_use$estimated_params$ln_srv_q_dev), ncol = length(new_years)))
      operating_model_use$estimated_params$ln_srv_q_dev_re <- cbind(operating_model_use$estimated_params$ln_srv_q_dev_re, matrix(operating_model_use$estimated_params$ln_srv_q_dev_re[,ncol(operating_model_use$estimated_params$ln_srv_q_dev_re)], nrow= nrow(operating_model_use$estimated_params$ln_srv_q_dev_re), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      n_selectivities <- nrow(operating_model_use$data_list$fleet_control)

      sel_slp_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]
      sel_inf_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]

      sel_slp_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effect deviations paramaters for logistic; n = [2, nspp]
      sel_inf_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effectdeviations paramaters for logistic; n = [2, nspp]

      sel_slp_dev[,,,1:nyrs_hind] <- operating_model_use$estimated_params$sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- operating_model_use$estimated_params$sel_inf_dev
      sel_slp_dev_re[,,,1:nyrs_hind] <- operating_model_use$estimated_params$sel_slp_dev_re
      sel_inf_dev_re[,,,1:nyrs_hind] <- operating_model_use$estimated_params$sel_inf_dev_re

      sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_slp_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_slp_dev_re[,,,nyrs_hind]
      sel_inf_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev_re[,,,nyrs_hind]

      operating_model_use$estimated_params$sel_slp_dev <- sel_slp_dev
      operating_model_use$estimated_params$sel_inf_dev <- sel_inf_dev
      operating_model_use$estimated_params$sel_slp_dev_re <- sel_slp_dev_re
      operating_model_use$estimated_params$sel_inf_dev_re <- sel_inf_dev_re


      # Update map (Only new parameter we are estimating in OM is the F_dev of the new years)
      operating_model_use$map <- Rceattle::build_map(operating_model_use$data_list, params = operating_model_use$estimated_params, debug = operating_model_use$data_list$debug, random_rec = FALSE)

      for (i in 1:length(operating_model_use$map[[2]])) {
        operating_model_use$map[[2]][[i]] <- replace(operating_model_use$map[[2]][[i]], values = rep(NA, length(operating_model_use$map[[2]][[i]])))
      }

      new_f_ind <- (ncol(operating_model_use$map[[2]]$F_dev) - length(new_years) + 1) : ncol(operating_model_use$map[[2]]$F_dev)
      operating_model_use$map[[2]]$F_dev[,new_f_ind] <- replace(operating_model_use$map[[2]]$F_dev[,new_f_ind], values = 1:length(operating_model_use$map[[2]]$F_dev[,new_f_ind]))


      # -- Map out surveys
      for (i in 1:nrow(operating_model_use$data_list$fleet_control)) {
        # Turn of F and F dev if not estimating of it is a Survey
        if (operating_model_use$data_list$fleet_control$Fleet_type[i] %in% c(0, 2)) {
          operating_model_use$map[[2]]$F_dev[i, ] <- NA
        }
      }

      # -- Map out Fdev for years with 0 catch to very low number
      fsh_biom <- operating_model_use$data_list$fsh_biom
      fsh_ind <- fsh_biom$Fleet_code[which(fsh_biom$Catch == 0)]
      yr_ind <- fsh_biom$Year[which(fsh_biom$Catch == 0)] - operating_model_use$data_list$styr + 1
      operating_model_use$map[[2]]$F_dev[fsh_ind, yr_ind] <- NA

      for (i in 1:length( operating_model_use$map[[2]])) {
        operating_model_use$map[[1]][[i]] <- factor( operating_model_use$map[[2]][[i]])
      }

      operating_model_use$estimated_params$ln_FSPR <- replace(operating_model_use$estimated_params$ln_FSPR, values = rep(-10, length(operating_model_use$estimated_params$ln_FSPR)))

      # Fit OM with new catch data
      operating_model_use <- fit_mod(TMBfilename = operating_model_use$TMBfilename,
                                         cpp_directory = operating_model_use$cpp_directory,
                                         data_list = operating_model_use$data_list,
                                         inits = operating_model_use$estimated_params,
                                         map =  operating_model_use$map,
                                         bounds = NULL,
                                         file = NULL,
                                         debug = operating_model_use$data_list$debug,
                                         random_rec = operating_model_use$data_list$random_rec,
                                         niter = operating_model_use$data_list$niter,
                                         msmMode = operating_model_use$data_list$msmMode,
                                         avgnMode = operating_model_use$data_list$avgnMode,
                                         minNByage = operating_model_use$data_list$debug,
                                         suitMode = operating_model_use$data_list$suitMode,
                                         phase = NULL,
                                         silent = TRUE,
                                         getsd = FALSE,
                                         recompile = FALSE)

      # Update standard error from survey for simulation
      operating_model_use$data_list$srv_biom$CV <- operating_model_use$quantities$srv_cv_hat

      # Simulate new survey and comp data
      sim_dat <- sim_mod(operating_model_use, simulate = simulate)

      years_include <- sample_yrs[which(sample_yrs > estimation_model_use$data_list$endyr & sample_yrs <= assess_yrs[k])]

      # -- Add survey data to EM
      new_srv_biom <- sim_dat$srv_biom[which(abs(sim_dat$srv_biom$Year) %in% years_include),]
      new_srv_biom$Year <- -new_srv_biom$Year
      estimation_model_use$data_list$srv_biom <- rbind(estimation_model_use$data_list$srv_biom, new_srv_biom)
      estimation_model_use$data_list$srv_biom <- estimation_model_use$data_list$srv_biom[
        with(estimation_model_use$data_list$srv_biom, order(Fleet_code, Year)),]

      # -- Add comp data to EM
      new_comp_data <- sim_dat$comp_data[which(abs(sim_dat$comp_data$Year) %in% years_include),]
      new_comp_data$Year <- -new_comp_data$Year
      estimation_model_use$data_list$comp_data <- rbind(estimation_model_use$data_list$comp_data, new_comp_data)

      estimation_model_use$data_list$comp_data <- estimation_model_use$data_list$comp_data[
        with(estimation_model_use$data_list$comp_data, order(Fleet_code, Year)),]

      # Update end year and re-estimate
      estimation_model_use$data_list$endyr <- assess_yrs[k]


      # Update parameters
      # -- F_dev
      estimation_model_use$estimated_params$F_dev <- cbind(estimation_model_use$estimated_params$F_dev, matrix(0, nrow= nrow(estimation_model_use$estimated_params$F_dev), ncol = length(new_years)))

      # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
      estimation_model_use$estimated_params$ln_srv_q_dev <- cbind(estimation_model_use$estimated_params$ln_srv_q_dev, matrix(estimation_model_use$estimated_params$ln_srv_q_dev[,ncol(estimation_model_use$estimated_params$ln_srv_q_dev)], nrow= nrow(estimation_model_use$estimated_params$ln_srv_q_dev), ncol = length(new_years)))
      estimation_model_use$estimated_params$ln_srv_q_dev_re <- cbind(estimation_model_use$estimated_params$ln_srv_q_dev_re, matrix(estimation_model_use$estimated_params$ln_srv_q_dev_re[,ncol(estimation_model_use$estimated_params$ln_srv_q_dev_re)], nrow= nrow(estimation_model_use$estimated_params$ln_srv_q_dev_re), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      n_selectivities <- nrow(estimation_model_use$data_list$fleet_control)

      sel_slp_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]
      sel_inf_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]

      sel_slp_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effect deviations paramaters for logistic; n = [2, nspp]
      sel_inf_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effectdeviations paramaters for logistic; n = [2, nspp]

      sel_slp_dev[,,,1:nyrs_hind] <- estimation_model_use$estimated_params$sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- estimation_model_use$estimated_params$sel_inf_dev
      sel_slp_dev_re[,,,1:nyrs_hind] <- estimation_model_use$estimated_params$sel_slp_dev_re
      sel_inf_dev_re[,,,1:nyrs_hind] <- estimation_model_use$estimated_params$sel_inf_dev_re

      sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_slp_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_slp_dev_re[,,,nyrs_hind]
      sel_inf_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev_re[,,,nyrs_hind]

      estimation_model_use$estimated_params$sel_slp_dev <- sel_slp_dev
      estimation_model_use$estimated_params$sel_inf_dev <- sel_inf_dev
      estimation_model_use$estimated_params$sel_slp_dev_re <- sel_slp_dev_re
      estimation_model_use$estimated_params$sel_inf_dev_re <- sel_inf_dev_re


      # Restimate
      estimation_model_use <- fit_mod(TMBfilename = estimation_model_use$TMBfilename,
                                      cpp_directory = estimation_model_use$cpp_directory,
                                      data_list = estimation_model_use$data_list,
                                      inits = estimation_model_use$estimated_params,
                                      map =  NULL,
                                      bounds = NULL,
                                      file = NULL,
                                      debug = estimation_model_use$data_list$debug,
                                      random_rec = estimation_model_use$data_list$random_rec,
                                      niter = estimation_model_use$data_list$niter,
                                      msmMode = estimation_model_use$data_list$msmMode,
                                      avgnMode = estimation_model_use$data_list$avgnMode,
                                      minNByage = estimation_model_use$data_list$debug,
                                      suitMode = estimation_model_use$data_list$suitMode,
                                      phase = "default",
                                      silent = TRUE,
                                      getsd = FALSE,
                                      recompile = FALSE)
      #plot_biomass(list(estimation_model_use, operating_model_use), model_names = c("EM", "OM"))
      # End year of assessment

      estimation_model_use$initial_params <- NULL
      estimation_model_use$bounds <- NULL
      estimation_model_use$map <- NULL
      estimation_model_use$obj <- NULL
      estimation_model_use$opt <- NULL
      estimation_model_use$sdrep <- NULL
      estimation_model_use$quantities[[names(estimation_model_use$quantities) %!in% c("catch_hat", "biomass", "mn_rec", "SB0", "biomassSSB" , "R", "srv_cv_hat")]] <- NULL

      Rceattle_EM_list[[sim]][[k+1]] <- estimation_model_use
    }

    plot_biomass(list(estimation_model_use, operating_model_use), model_names = c("EM", "OM"))

    # Save models
    Rceattle_OM_list[[sim]] <- operating_model_use
  }

  return(list(Rceattle_OM_list = Rceattle_OM_list, Rceattle_EM_list = Rceattle_EM_list))
}
