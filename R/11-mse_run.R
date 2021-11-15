#' Run a management strategy evaluation
#'
#' @description Runs a forward projecting MSE. Main assumptions are the projected selectivity/catchability, foraging days, and weight-at-age are the same as the terminal year of the hindcast in the operating model. Assumes survey sd is same as average across historic time series, while comp data sample size is same as last year. No implementation error!
#'
#' @param om CEATTLE model object exported from \code{\link{Rceattle}}
#' @param em CEATTLE model object exported from \code{\link{Rceattle}}
#' @param nsim Number of simulations to run (default 10)
#' @param assessment_period Period of years that each assessment is taken
#' @param sampling_period Period of years data sampling is conducted
#' @param simulate Include simulated random error proportional to that estimated/provided.
#' @param cap A cap on the catch in the projection. Can be a single number or vector. Default = NULL
#'
#' @return A list of operating models (differ by simulated recruitment determined by \code{nsim}) and estimation models fit to each operating model (differ by terminal year).
#' @export
#'
#' @examples
mse_run <- function(om = ms_run, em = ss_run, nsim = 10, assessment_period = 1, sampling_period = 1, simulate = TRUE, cap = NULL, seed = 666){
  '%!in%' <- function(x,y)!('%in%'(x,y))

  set.seed(seed)

  Rceattle_OM_list <- list()
  Rceattle_EM_list <- list()

  # - Adjust cap
  if(!is.null(cap)){
    if(length(cap) == 1){
      cap = rep(cap, om$data_list$nspp)
    }

    if(length(cap) != om$data_list$nspp){
      stop("cap is not length 1 or length nspp")
    }
  }

  # Update data-files in OM so we can fill in updated years
  # -- srv_biom
  for(flt in (unique(om$data_list$srv_biom$Fleet_code))){
    sub_rows <- which(om$data_list$srv_biom$Fleet_code == flt)
    srv_biom_sub <- om$data_list$srv_biom[sub_rows,]
    yrs_proj <- (om$data_list$endyr + 1):om$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_srv_biom <- data.frame(Fleet_name = rep(srv_biom_sub$Fleet_name[1], nyrs_proj),
                                Fleet_code = rep(flt, nyrs_proj),
                                Species = rep(srv_biom_sub$Species[1], nyrs_proj),
                                Year = -yrs_proj,
                                Month = rep(srv_biom_sub$Month[length(srv_biom_sub$Month)], nyrs_proj),
                                Selectivity_block = rep(srv_biom_sub$Selectivity_block[length(srv_biom_sub$Selectivity_block)], nyrs_proj),
                                Q_block = rep(srv_biom_sub$Selectivity_block[length(srv_biom_sub$Selectivity_block)], nyrs_proj),
                                Observation = rep(NA, nyrs_proj),
                                Log_sd = rep(mean(om$quantities$srv_log_sd_hat[sub_rows], na.rm = TRUE), nyrs_proj)) # FIXME: not sure what the best is to put here, using mean
    om$data_list$srv_biom <- rbind(om$data_list$srv_biom, proj_srv_biom)
  }
  om$data_list$srv_biom <- om$data_list$srv_biom[
    with(om$data_list$srv_biom, order(Fleet_code, abs(Year))),]
  om$quantities$srv_log_sd_hat <- om$data_list$srv_biom$Log_sd
  om$quantities$srv_bio_hat <- om$data_list$srv_biom$Observation


  # -- Nbyage
  for(spp in 1:om$data_list$nspp){
    if(om$data_list$estDynamics[spp] > 0){
      for(sex in 1:om$data_list$nsex){
        sub_rows <- which(om$data_list$NByageFixed$Species == spp &
                            om$data_list$NByageFixed$Sex == sex)
        NByageFixed_sub <- om$data_list$NByageFixed[sub_rows,]
        yrs_proj <- (om$data_list$endyr + 1):om$data_list$projyr
        yrs_proj <- yrs_proj[which(yrs_proj %!in% NByageFixed_sub$Year)] # What projection years are not present
        nyrs_proj <- length(yrs_proj)

        proj_NByageFixed <- data.frame(Species_name = rep(comp_data_sub$Fleet_name[1], nyrs_proj),
                                       Species  = rep(flt, nyrs_proj),
                                     Sex = rep(comp_data_sub$Sex[length(comp_data_sub$Sex)], nyrs_proj),
                                     Year = -yrs_proj) # Negative year for predicting
        proj_nbyage <- data.frame(matrix(0, ncol = ncol(NByageFixed_sub) - ncol(proj_NByageFixed), nrow = nyrs_proj))
        colnames(proj_nbyage) <-paste0("Age", 1:ncol(proj_comps))
        proj_NByageFixed <-cbind(proj_NByageFixed, proj_nbyage)
        om$data_list$NByageFixed <- rbind(om$data_list$NByageFixed, proj_NByageFixed)
      }
    }
  }
  om$data_list$NByageFixed <- om$data_list$NByageFixed[
    with(om$data_list$NByageFixed, order(Species, abs(Year))),]


  # -- comp_data
  #FIXME may not work if male/females composition is separate
  for(flt in (unique(om$data_list$comp_data$Fleet_code))){
    comp_data_sub <- om$data_list$comp_data[which(om$data_list$comp_data$Fleet_code == flt),]
    yrs_proj <- (om$data_list$endyr + 1):om$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_comp_data <- data.frame(Fleet_name = rep(comp_data_sub$Fleet_name[1], nyrs_proj),
                                 Fleet_code = rep(flt, nyrs_proj),
                                 Species = rep(comp_data_sub$Species[1], nyrs_proj),
                                 Sex = rep(comp_data_sub$Sex[length(comp_data_sub$Sex)], nyrs_proj),
                                 Age0_Length1 = rep(comp_data_sub$Age0_Length1[length(comp_data_sub$Age0_Length1)], nyrs_proj),
                                 Year = -yrs_proj, # Negative year for predicting
                                 Month = rep(comp_data_sub$Month[length(comp_data_sub$Month)], nyrs_proj),
                                 Sample_size = rep(comp_data_sub$Sample_size[length(comp_data_sub$Sample_size)], nyrs_proj))
    proj_comps <- data.frame(matrix(0, ncol = ncol(comp_data_sub) - ncol(proj_comp_data), nrow = nyrs_proj))
    colnames(proj_comps) <-paste0("Comp_", 1:ncol(proj_comps))
    proj_comp_data <-cbind(proj_comp_data, proj_comps)

    om$data_list$comp_data <- rbind(om$data_list$comp_data, proj_comp_data)
  }
  om$data_list$comp_data <- om$data_list$comp_data[
    with(om$data_list$comp_data, order(Fleet_code, abs(Year))),]

  # -- emp_sel - Use terminal year
  for(flt in (unique(om$data_list$emp_sel$Fleet_code))){
    emp_sel_sub <- om$data_list$emp_sel[which(om$data_list$emp_sel$Fleet_code == flt),]
    yrs_proj <- (om$data_list$endyr + 1):om$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_emp_sel <- data.frame(Fleet_name = rep(emp_sel_sub$Fleet_name[1], nyrs_proj),
                               Fleet_code = rep(flt, nyrs_proj),
                               Species = rep(emp_sel_sub$Species[1], nyrs_proj),
                               Sex = rep(emp_sel_sub$Sex[length(emp_sel_sub$Sex)], nyrs_proj),
                               Year = yrs_proj)
    proj_comps <- data.frame(matrix(emp_sel_sub[nrow(emp_sel_sub), (ncol(proj_emp_sel) + 1) : ncol(emp_sel_sub)], ncol = ncol(emp_sel_sub) - ncol(proj_emp_sel), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Comp_", 1:ncol(proj_comps))
    proj_emp_sel <-cbind(proj_emp_sel, proj_comps)

    om$data_list$emp_sel <- rbind(om$data_list$emp_sel, proj_emp_sel)
  }
  om$data_list$emp_sel <- om$data_list$emp_sel[
    with(om$data_list$emp_sel, order(Fleet_code, Year)),]

  # -- wt
  for(flt in (unique(om$data_list$wt$Wt_index))){
    wt_sub <- om$data_list$wt[which(om$data_list$wt$Wt_index == flt),]
    yrs_proj <- (om$data_list$endyr + 1):om$data_list$projyr
    yrs_proj <- yrs_proj[which(yrs_proj %!in% wt_sub$Year)]
    nyrs_proj <- length(yrs_proj)
    proj_wt <- data.frame(Wt_name = rep(wt_sub$Wt_name[1], nyrs_proj),
                          Wt_index = rep(flt, nyrs_proj),
                          Species = rep(wt_sub$Species[1], nyrs_proj),
                          Sex = rep(wt_sub$Sex[length(wt_sub$Sex)], nyrs_proj),
                          Year = yrs_proj
    )
    proj_comps <- data.frame(matrix(wt_sub[nrow(wt_sub), (ncol(proj_wt)+1) : ncol(wt_sub)], ncol = ncol(wt_sub) - ncol(proj_wt), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_wt <-cbind(proj_wt, proj_comps)
    colnames(proj_wt) <- colnames(om$data_list$wt)
    om$data_list$wt <- rbind(om$data_list$wt, proj_wt)
  }
  om$data_list$wt <- om$data_list$wt[
    with(om$data_list$wt, order(Wt_index, Year)),]


  # -- Pyrs
  for(flt in (unique(om$data_list$wt$Species))){
    Pyrs_sub <- om$data_list$Pyrs[which(om$data_list$Pyrs$Species == flt),]
    yrs_proj <- (om$data_list$endyr + 1):om$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_Pyrs <- data.frame(Species = rep(Pyrs_sub$Species[1], nyrs_proj),
                            Sex = rep(Pyrs_sub$Sex[length(Pyrs_sub$Sex)], nyrs_proj),
                            Year = yrs_proj)

    proj_comps <- data.frame(matrix(Pyrs_sub[nrow(Pyrs_sub), (ncol(proj_Pyrs)+1) : ncol(Pyrs_sub)], ncol = ncol(Pyrs_sub) - ncol(proj_Pyrs), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_Pyrs <-cbind(proj_Pyrs, proj_comps)
    colnames(proj_Pyrs) <- colnames(om$data_list$Pyrs)
    om$data_list$Pyrs <- rbind(om$data_list$Pyrs, proj_Pyrs)
  }
  om$data_list$Pyrs <- om$data_list$Pyrs[
    with(om$data_list$Pyrs, order(Species, Year)),]


  # Update data in EM
  #FIXME - assuming same as terminal year of hindcast
  # -- emp_sel - Use terminal year
  for(flt in (unique(em$data_list$emp_sel$Fleet_code))){
    emp_sel_sub <- em$data_list$emp_sel[which(em$data_list$emp_sel$Fleet_code == flt),]
    yrs_proj <- (em$data_list$endyr + 1):em$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_emp_sel <- data.frame(Fleet_name = rep(emp_sel_sub$Fleet_name[1], nyrs_proj),
                               Fleet_code = rep(flt, nyrs_proj),
                               Species = rep(emp_sel_sub$Species[1], nyrs_proj),
                               Sex = rep(emp_sel_sub$Sex[length(emp_sel_sub$Sex)], nyrs_proj),
                               Year = yrs_proj)
    proj_comps <- data.frame(matrix(emp_sel_sub[nrow(emp_sel_sub), (ncol(proj_emp_sel) + 1) : ncol(emp_sel_sub)], ncol = ncol(emp_sel_sub) - ncol(proj_emp_sel), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Comp_", 1:ncol(proj_comps))
    proj_emp_sel <-cbind(proj_emp_sel, proj_comps)

    em$data_list$emp_sel <- rbind(em$data_list$emp_sel, proj_emp_sel)
  }
  em$data_list$emp_sel <- em$data_list$emp_sel[
    with(em$data_list$emp_sel, order(Fleet_code, Year)),]


  # -- wt
  for(flt in (unique(em$data_list$wt$Wt_index))){
    wt_sub <- em$data_list$wt[which(em$data_list$wt$Wt_index == flt),]
    yrs_proj <- (em$data_list$endyr + 1):em$data_list$projyr
    yrs_proj <- yrs_proj[which(yrs_proj %!in% wt_sub$Year)]
    nyrs_proj <- length(yrs_proj)
    proj_wt <- data.frame(Wt_name = rep(wt_sub$Wt_name[1], nyrs_proj),
                          Wt_index = rep(flt, nyrs_proj),
                          Species = rep(wt_sub$Species[1], nyrs_proj),
                          Sex = rep(wt_sub$Sex[length(wt_sub$Sex)], nyrs_proj),
                          Year = yrs_proj
    )
    proj_comps <- data.frame(matrix(wt_sub[nrow(wt_sub), (ncol(proj_wt)+1) : ncol(wt_sub)], ncol = ncol(wt_sub) - ncol(proj_wt), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_wt <-cbind(proj_wt, proj_comps)
    colnames(proj_wt) <- colnames(om$data_list$wt)
    em$data_list$wt <- rbind(em$data_list$wt, proj_wt)
  }
  em$data_list$wt <- em$data_list$wt[
    with(em$data_list$wt, order(Wt_index, Year)),]


  # -- Pyrs
  for(flt in (unique(em$data_list$wt$Species))){
    Pyrs_sub <- em$data_list$Pyrs[which(em$data_list$Pyrs$Species == flt),]
    yrs_proj <- (em$data_list$endyr + 1):em$data_list$projyr
    nyrs_proj <- length(yrs_proj)
    proj_Pyrs <- data.frame(Species = rep(Pyrs_sub$Species[1], nyrs_proj),
                            Sex = rep(Pyrs_sub$Sex[length(Pyrs_sub$Sex)], nyrs_proj),
                            Year = yrs_proj)

    proj_comps <- data.frame(matrix(Pyrs_sub[nrow(Pyrs_sub), (ncol(proj_Pyrs)+1) : ncol(Pyrs_sub)], ncol = ncol(Pyrs_sub) - ncol(proj_Pyrs), nrow = nyrs_proj, byrow = TRUE))
    colnames(proj_comps) <-paste0("Age", 1:ncol(proj_comps))
    proj_Pyrs <-cbind(proj_Pyrs, proj_comps)
    colnames(proj_Pyrs) <- colnames(om$data_list$Pyrs)
    em$data_list$Pyrs <- rbind(em$data_list$Pyrs, proj_Pyrs)
  }
  em$data_list$Pyrs <- em$data_list$Pyrs[
    with(em$data_list$Pyrs, order(Species, Year)),]



  # Years for simulations
  proj_yrs <- (em$data_list$endyr + 1) : em$data_list$projyr
  proj_nyrs <- length(proj_yrs)
  assess_yrs <- seq(from = em$data_list$endyr + assessment_period, to = em$data_list$projyr,  by = assessment_period)
  sample_yrs <- seq(from = em$data_list$endyr + sampling_period, to = em$data_list$projyr,  by = sampling_period)


  # Do the MSE
  for(sim in 1:nsim){

    # Set models
    Rceattle_EM_list[[sim]] <- list()
    Rceattle_EM_list[[sim]][[1]] <- em
    em_use <- em
    om_use <- om



    # Replace future rec_devs with numbers
    if(simulate){
      for(sp in 1:om_use$data_list$nspp){
        om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1] <- replace(
          om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1],
          values = rnorm( length(om_use$estimated_params$rec_dev[sp,proj_yrs - om_use$data_list$styr + 1]),
                          mean = 0,
                          sd = exp(om_use$estimated_params$ln_rec_sigma[sp])) # Assumed value from penalized likelihood
        )
      }
    }



    # Run through model
    for(k in 1:(length(assess_yrs))){

      # ------------------------------------------------------------
      # 1. OBSERVATION MODEL
      # ------------------------------------------------------------
      new_years <- proj_yrs[which(proj_yrs <= assess_yrs[k] & proj_yrs > om_use$data_list$endyr)]

      # - Get projected catch data from EM
      new_catch_data <- em_use$data_list$fsh_biom
      dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
      new_catch_data$Catch[dat_fill_ind] <- em_use$quantities$fsh_bio_hat[dat_fill_ind]
      if(!is.null(cap)){
        new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap[new_catch_data$Species[dat_fill_ind]], cap[new_catch_data$Species[dat_fill_ind]], new_catch_data$Catch[dat_fill_ind])
      }

      # - Update catch data in OM and EM
      om_use$data_list$fsh_biom <- new_catch_data
      em_use$data_list$fsh_biom <- new_catch_data

      # - Update endyr of OM
      nyrs_hind <- om_use$data_list$endyr - om_use$data_list$styr + 1
      om_use$data_list$endyr <- assess_yrs[k]

      # - Update parameters
      # -- F_dev
      om_use$estimated_params$F_dev <- cbind(om_use$estimated_params$F_dev, matrix(0, nrow= nrow(om_use$estimated_params$F_dev), ncol = length(new_years)))

      # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
      om_use$estimated_params$ln_srv_q_dev <- cbind(om_use$estimated_params$ln_srv_q_dev, matrix(om_use$estimated_params$ln_srv_q_dev[,ncol(om_use$estimated_params$ln_srv_q_dev)], nrow= nrow(om_use$estimated_params$ln_srv_q_dev), ncol = length(new_years)))

      #FIXME: update random effects q if used again
      # om_use$estimated_params$ln_srv_q_dev_re <- cbind(om_use$estimated_params$ln_srv_q_dev_re, matrix(om_use$estimated_params$ln_srv_q_dev_re[,ncol(om_use$estimated_params$ln_srv_q_dev_re)], nrow= nrow(om_use$estimated_params$ln_srv_q_dev_re), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      n_selectivities <- nrow(om_use$data_list$fleet_control)

      ln_sel_slp_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]
      sel_inf_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]

      #FIXME: update random effects sel if used again
      # sel_slp_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effect deviations paramaters for logistic; n = [2, nspp]
      # sel_inf_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effectdeviations paramaters for logistic; n = [2, nspp]

      ln_sel_slp_dev[,,,1:nyrs_hind] <- om_use$estimated_params$ln_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_inf_dev

      #FIXME: update random effects sel if used again
      # sel_slp_dev_re[,,,1:nyrs_hind] <- om_use$estimated_params$sel_slp_dev_re
      # sel_inf_dev_re[,,,1:nyrs_hind] <- om_use$estimated_params$sel_inf_dev_re

      ln_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- ln_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      #FIXME: update random effects sel if used again
      # sel_slp_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_slp_dev_re[,,,nyrs_hind]
      # sel_inf_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev_re[,,,nyrs_hind]

      om_use$estimated_params$ln_sel_slp_dev <- ln_sel_slp_dev
      om_use$estimated_params$sel_inf_dev <- sel_inf_dev
      #FIXME: update random effects sel if used again
      # om_use$estimated_params$sel_slp_dev_re <- sel_slp_dev_re
      # om_use$estimated_params$sel_inf_dev_re <- sel_inf_dev_re


      # - Update map (Only new parameter we are estimating in OM is the F_dev of the new years)
      om_use$map <- build_map(
        data_list = om_use$data_list,
        params = om_use$estimated_params,
        debug = om_use$data_list$debug,
        random_rec = om_use$data_list$random_rec)

      # -- Fill in with NA's
      for (i in 1:length(om_use$map[[2]])) {
        om_use$map[[2]][[i]] <- replace(om_use$map[[2]][[i]], values = rep(NA, length(om_use$map[[2]][[i]])))
      }

      # -- Estimate terminal F for catch
      new_f_ind <- (ncol(om_use$map[[2]]$F_dev) - length(new_years) + 1) : ncol(om_use$map[[2]]$F_dev)
      om_use$map[[2]]$F_dev[,new_f_ind] <- replace(om_use$map[[2]]$F_dev[,new_f_ind], values = 1:length(om_use$map[[2]]$F_dev[,new_f_ind]))


      # --- Turn off F for surveys
      for (i in 1:nrow(om_use$data_list$fleet_control)) {
        # Turn of F and F dev if not estimating of it is a Survey
        if (om_use$data_list$fleet_control$Fleet_type[i] %in% c(0, 2)) {
          om_use$map[[2]]$F_dev[i, ] <- NA
        }
      }

      # -- Map out Fdev for years with 0 catch to very low number
      fsh_biom <- om_use$data_list$fsh_biom
      fsh_ind <- fsh_biom$Fleet_code[which(fsh_biom$Catch == 0)]
      yr_ind <- fsh_biom$Year[which(fsh_biom$Catch == 0)] - om_use$data_list$styr + 1
      om_use$map[[2]]$F_dev[fsh_ind, yr_ind] <- NA

      for (i in 1:length( om_use$map[[2]])) {
        om_use$map[[1]][[i]] <- factor( om_use$map[[2]][[i]])
      }

      # om_use$estimated_params$ln_FSPR <- replace(om_use$estimated_params$ln_FSPR, values = rep(-10, length(om_use$estimated_params$ln_FSPR)))

      # - Fit OM with new catch data
      om_use <- fit_mod(
        data_list = om_use$data_list,
        inits = om_use$estimated_params,
        map =  om_use$map,
        bounds = NULL,
        file = NULL,
        debug = om_use$data_list$debug,
        random_rec = om_use$data_list$random_rec,
        niter = om_use$data_list$niter,
        msmMode = om_use$data_list$msmMode,
        avgnMode = om_use$data_list$avgnMode,
        minNByage = om_use$data_list$minNByage,
        suitMode = om_use$data_list$suitMode,
        suityr = om$data_list$endyr,
        phase = NULL,
        getsd = FALSE,
        verbose = 0)

      # ------------------------------------------------------------
      # 2. ESTIMATION MODEL
      # ------------------------------------------------------------
      # - Simulate new survey and comp data
      sim_dat <- sim_mod(om_use, simulate = simulate)

      years_include <- sample_yrs[which(sample_yrs > em_use$data_list$endyr & sample_yrs <= assess_yrs[k])]

      # -- Add newly simulated survey data to EM
      new_srv_biom <- sim_dat$srv_biom[which(abs(sim_dat$srv_biom$Year) %in% years_include),]
      new_srv_biom$Year <- -new_srv_biom$Year
      em_use$data_list$srv_biom <- rbind(em_use$data_list$srv_biom, new_srv_biom)
      em_use$data_list$srv_biom <- em_use$data_list$srv_biom[
        with(em_use$data_list$srv_biom, order(Fleet_code, abs(Year))),]

      # -- Add newly simulated comp data to EM
      new_comp_data <- sim_dat$comp_data[which(abs(sim_dat$comp_data$Year) %in% years_include),]
      new_comp_data$Year <- -new_comp_data$Year
      em_use$data_list$comp_data <- rbind(em_use$data_list$comp_data, new_comp_data)
      em_use$data_list$comp_data <- em_use$data_list$comp_data[
        with(em_use$data_list$comp_data, order(Fleet_code, abs(Year))),]

      # Update end year and re-estimate
      em_use$data_list$endyr <- assess_yrs[k]

      # Update parameters
      # -- F_dev
      em_use$estimated_params$F_dev <- cbind(em_use$estimated_params$F_dev, matrix(0, nrow= nrow(em_use$estimated_params$F_dev), ncol = length(new_years)))

      # -- Time-varying survey catachbilitiy - Assume last year - filled by columns
      em_use$estimated_params$ln_srv_q_dev <- cbind(em_use$estimated_params$ln_srv_q_dev, matrix(em_use$estimated_params$ln_srv_q_dev[,ncol(em_use$estimated_params$ln_srv_q_dev)], nrow= nrow(em_use$estimated_params$ln_srv_q_dev), ncol = length(new_years)))
      # em_use$estimated_params$ln_srv_q_dev_re <- cbind(em_use$estimated_params$ln_srv_q_dev_re, matrix(em_use$estimated_params$ln_srv_q_dev_re[,ncol(em_use$estimated_params$ln_srv_q_dev_re)], nrow= nrow(em_use$estimated_params$ln_srv_q_dev_re), ncol = length(new_years)))

      # -- Time-varing selectivity - Assume last year - filled by columns
      n_selectivities <- nrow(em_use$data_list$fleet_control)

      ln_sel_slp_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]
      sel_inf_dev = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity deviations paramaters for logistic; n = [2, nspp]

      # sel_slp_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effect deviations paramaters for logistic; n = [2, nspp]
      # sel_inf_dev_re = array(0, dim = c(2, n_selectivities, 2, nyrs_hind + length(new_years)))  # selectivity random effectdeviations paramaters for logistic; n = [2, nspp]

      ln_sel_slp_dev[,,,1:nyrs_hind] <- em_use$estimated_params$ln_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_inf_dev
      # sel_slp_dev_re[,,,1:nyrs_hind] <- em_use$estimated_params$sel_slp_dev_re
      # sel_inf_dev_re[,,,1:nyrs_hind] <- em_use$estimated_params$sel_inf_dev_re

      # - Initialize next year with terminal year
      ln_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- ln_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      # sel_slp_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_slp_dev_re[,,,nyrs_hind]
      # sel_inf_dev_re[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev_re[,,,nyrs_hind]

      em_use$estimated_params$ln_sel_slp_dev <- ln_sel_slp_dev
      em_use$estimated_params$sel_inf_dev <- sel_inf_dev
      # em_use$estimated_params$sel_slp_dev_re <- sel_slp_dev_re
      # em_use$estimated_params$sel_inf_dev_re <- sel_inf_dev_re


      # Restimate
      em_use <- fit_mod(
        data_list = em_use$data_list,
        inits = em_use$estimated_params,
        map =  NULL,
        bounds = NULL,
        file = NULL,
        debug = em_use$data_list$debug,
        random_rec = em_use$data_list$random_rec,
        niter = em_use$data_list$niter,
        msmMode = em_use$data_list$msmMode,
        avgnMode = em_use$data_list$avgnMode,
        minNByage = em_use$data_list$minNByage,
        suitMode = em_use$data_list$suitMode,
        phase = NULL,
        getsd = FALSE,
        verbose = 0)
      # plot_biomass(list(em_use, om_use), model_names = c("EM", "OM"))
      # End year of assessment

      # - Remove unneeded bits for memory reasons
      em_use$initial_params <- NULL
      em_use$bounds <- NULL
      em_use$map <- NULL
      em_use$obj <- NULL
      em_use$opt <- NULL
      em_use$sdrep <- NULL
      em_use$quantities[names(em_use$quantities) %!in% c("fsh_bio_hat", "biomass", "F", "F_tot", "mn_rec", "SB0", "SB40", "F40_tot" , "F35_tot" , "biomassSSB" , "R", "srv_log_sd_hat", "FSPR")] <- NULL

      Rceattle_EM_list[[sim]][[k+1]] <- em_use
      message(paste0("Sim ",sim, " - EM Year ", assess_yrs[k], " COMPLETE"))
    }

    # Save models
    Rceattle_OM_list[[sim]] <- om_use
    names(Rceattle_EM_list[[sim]]) <- c("EM", paste0("OM_Sim_",sim,". EM_projyr_", assess_yrs))
  }

  # - Name them
  names(Rceattle_OM_list) <- paste0("OM_Sim_",1:nsim)
  names(Rceattle_EM_list) <- paste0("OM_Sim_",1:nsim)

  return(list(OM_list = Rceattle_OM_list, EM_list = Rceattle_EM_list, OM = om, EM = em))
}
