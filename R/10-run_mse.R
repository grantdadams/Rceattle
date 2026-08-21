# Derived quantities kept when an MSE operating/estimation model object is
# slimmed for storage; every other quantity is dropped to keep saved runs small.
# Shared by run_mse() and mse_summary() so the retained set stays in one place.
.mse_keep_quantities <- c(
  "catch_hat", "catch_sd", "index_hat", "index_sd",
  "log_catch_sd", "log_index_sd",          # deprecated spellings, dropped next release
  "ssb_depletion", "biomass_depletion", "biomass", "ssb",
  "BO", "SB0", "SBF", "F_spp", "R",
  "M1_at_age", "M_at_age", "avg_rec",
  "DynamicB0", "DynamicSB0", "DynamicSBF",
  "SPR0", "SPRlimit", "SPRtarget", "Ftarget",
  "B_eaten", "B_eaten_as_prey", "Flimit"
)

# Shortening the operating model's projection horizon -------------------------
#
# Between assessments the operating model only has to reach one assessment step
# past its terminal year: run_mse() reads `max_catch_hat` in the upcoming
# assessment year to cap the TAC at exploitable biomass, and nothing looks
# further ahead than that. `projyr` is what sizes the AD tape, so refitting on
# the shorter horizon keeps each refit proportional to the years realized so
# far rather than to the whole projection.
#
# clean_data() filters every data frame to styr:projyr, so the shorter horizon
# also drops the projection-year placeholder rows run_mse() appends once up
# front. Those rows are only ever rewritten for years at or before the current
# assessment year -- inside the shortened horizon -- so the dropped ones still
# hold the values they were created with and are restored verbatim.

# Parameter blocks build_params() dimensions by styr:projyr, with the position
# of their year dimension. Every other block is hindcast-length.
.mse_proj_param_yrdim <- c(rec_dev = 2L, log_M1_dev = 4L)

# Data frames carried out to projyr -- exactly the ones clean_data() filters to
# styr:projyr, so exactly the ones a shortened horizon truncates.
.mse_proj_tables <- c("index_data", "comp_data", "caal_data", "catch_data",
                      "diet_data", "NByageFixed", "emp_sel", "weight",
                      "ration_data")

# Between assessments the loop reads two reported quantities off the operating
# model, both indexed by catch_data's rows: max_catch_hat, to cap the next TAC
# at exploitable biomass, and catch_hat, for the catch the operating model
# actually took. Those are lifted back onto the full-length table; the rest of
# the fit's quantities are left as the shortened refit produced them. Nothing
# else reads them -- sim_mod() works off the un-restored fit, a completed
# simulation returns the terminal full-horizon fit, and a failed one returns no
# model at all.
.mse_proj_catch_quantities <- c("catch_hat", "max_catch_hat")

# Replace the year dimension of an array with `keep`, whatever its rank.
.mse_slice_year_dim <- function(x, yr_dim, keep) {
  idx <- lapply(dim(x), seq_len)
  idx[[yr_dim]] <- keep
  do.call(`[`, c(list(x), idx, list(drop = FALSE)))
}

.mse_trim_proj_params <- function(params, nyrs_keep) {
  for (nm in names(.mse_proj_param_yrdim)) {
    if (is.null(params[[nm]])) next
    params[[nm]] <- .mse_slice_year_dim(params[[nm]],
                                        .mse_proj_param_yrdim[[nm]],
                                        seq_len(nyrs_keep))
  }
  params
}

# Write the shortened blocks back over the head of the full-horizon ones, so
# the projection years the refit did not cover keep their existing values --
# for rec_dev those are the recruitment deviations sample_rec() drew for this
# simulation.
.mse_restore_proj_params <- function(params, params_full) {
  for (nm in names(.mse_proj_param_yrdim)) {
    short <- params[[nm]]
    full  <- params_full[[nm]]
    if (is.null(short) || is.null(full)) next
    yr_dim <- .mse_proj_param_yrdim[[nm]]
    idx <- lapply(dim(full), seq_len)
    idx[[yr_dim]] <- seq_len(dim(short)[yr_dim])
    params[[nm]] <- do.call(`[<-`, c(list(full), idx, list(value = short)))
  }
  params
}

# Put a shortened-horizon fit back on the full projection horizon: data frames
# regain their future rows, projection-length parameter blocks regain their
# future slices, and the two catch quantities the loop reads are re-expanded to
# match (NA in the years the refit did not cover).
.mse_restore_om_horizon <- function(fit, data_list_full, params_full) {
  short_projyr <- fit$data_list$projyr

  for (nm in .mse_proj_tables) {
    full  <- data_list_full[[nm]]
    short <- fit$data_list[[nm]]
    if (is.null(full) || is.null(short) || !nrow(full)) next

    keep <- abs(full$Year) <= short_projyr
    # The refit sees the full-horizon table and only filters it, so its rows are
    # the kept rows in the same order and its columns are unchanged. Check
    # rather than assume: either mismatch would silently misalign the rows or
    # columns written back below, and with them every row-indexed quantity.
    # Compare years by value -- clean_data() may return them as integer where
    # the source table held doubles.
    yrs_full  <- as.numeric(abs(full$Year))[keep]
    yrs_short <- as.numeric(abs(short$Year))
    if (length(yrs_full) != length(yrs_short) || any(yrs_full != yrs_short) ||
        !identical(names(full), names(short))) {
      stop("run_mse(): '", nm, "' did not survive the shortened operating-model ",
           "horizon unchanged (", sum(keep), " rows expected, ", nrow(short),
           " returned).", call. = FALSE)
    }

    # clean_data() renumbers diet_data's stomach_id from the rows it is given,
    # so a table that actually lost rows would come back with the shortened
    # numbering spliced into the full-horizon one. No bundled model stratifies
    # diet by projection year, so this never fires -- but it must not pass
    # silently if one ever does.
    if (identical(nm, "diet_data") && !all(keep)) {
      stop("run_mse(): diet_data is stratified past the shortened operating-model ",
           "horizon; its stomach_id would be renumbered. Widen the horizon or ",
           "drop diet_data from .mse_proj_tables.", call. = FALSE)
    }

    restored <- full
    restored[keep, ] <- short
    fit$data_list[[nm]] <- restored

    if (identical(nm, "catch_data")) {
      for (q in .mse_proj_catch_quantities) {
        x <- fit$quantities[[q]]
        if (is.null(x)) next
        # NA outside the refit's horizon: those years have not been projected,
        # and the loop only ever indexes rows at or before it.
        out <- rep(NA_real_, nrow(full))
        out[keep] <- x
        if (!is.null(names(x))) {
          nms <- rep(NA_character_, nrow(full))
          nms[keep] <- names(x)
          names(out) <- nms
        }
        fit$quantities[[q]] <- out
      }
    }
  }

  fit$data_list$projyr <- data_list_full$projyr
  fit$estimated_params <- .mse_restore_proj_params(fit$estimated_params,
                                                   params_full)
  fit
}

#' Run a management strategy evaluation
#'
#' @description Runs a forward-projecting management strategy evaluation (MSE). Projected selectivity, catchability, foraging days, and weight-at-age are held at the operating model's terminal hindcast year. Survey SD is set to the average over the historical time series, and composition sample size is held at the last year. There is no implementation error and no observation error on catch.
#'
#' @param om CEATTLE model object exported from \code{Rceattle}
#' @param em CEATTLE model object exported from \code{Rceattle}
#' @param nsim Number of simulations to run (default 10)
#' @param start_sim First simulation number to start at. Useful if the code stops at specific seed/sim (default = 1).
#' @param assessment_period Period of years that each assessment is taken
#' @param sampling_period Period of years data sampling is conducted. Single value or vector the same length as the number of fleets.
#' @param simulate_data Include simulated random error proportional to that estimated/provided for the data from the OM.
#' @param regenerate_past Refits the EM to historical/conditioning data prior to the MSE, where the data are generated from the OM with \code{simulate_data = TRUE} or without \code{simulate_data = FALSE} sampling error.
#' @param sample_rec Include resampled recruitment deviations from the hindcast in the OM projection. Resampled deviations are used rather than drawing from N(0, sigmaR) because the initial deviations bias R0 low. If FALSE, uses the mean recruitment deviation.
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier \code{mean rec * (1 + (rec_trend/projection years) * 1:projection years)}. Can be of length 1 or of length nspp. If length 1, all species get the same trend.
#' @param fut_sample future sampling effort relative to last year.  \code{ Log_sd * 1 / fut_sample} for index and \code{ Sample_size * fut_sample} for comps
#' @param cap A cap on the catch in the projection. Can be a single number applied to all species (proportional to recommended catch) or vector of length \code{nspp} applied to each species. Default = NULL
#' @param catch_mult A multiplier for the catch in the projection. Can be a single number or vector of length nspp. Default = NULL
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param file (Optional) Filename where each OM simulation with EMs will be saved. If NULL, no files are saved.
#' @param dir (Optional) Directory where each OM simulation is saved
#' @param seed seed for the simulation
#' @param regenerate_seed seed for regenerating data
#' @param timeout length of time (minutes) estimation will run before stopping a sim (default 999 minutes)
#' @param endyr Terminal year of the MSE projection. Default = NA uses \code{projyr} from the operating model.
#' @param cores Number of cores to use for parallel simulations. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#'
#' @return A list of operating models (differ by simulated recruitment determined by \code{nsim}) and estimation models fit to each operating model (differ by terminal year).
#' @export
run_mse <- function(om, em, nsim = 10, start_sim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, regenerate_past = FALSE, sample_rec = TRUE, rec_trend = 0, fut_sample = 1, cap = NULL, catch_mult = NULL, seed = 666, regenerate_seed = seed, loopnum = 1, file = NULL, dir = NULL, timeout = 999, endyr = NA, cores = NULL){

  # om = om; em = em; nsim = 1; start_sim = 1; assessment_period = 1; sampling_period = 1; simulate_data = TRUE; regenerate_past = FALSE; sample_rec = FALSE; rec_trend = 0; fut_sample = 1; cap = NULL; catch_mult = NULL; seed = 666; regenerate_seed = seed; loopnum = 1; file = NULL; dir = NULL; endyr = NA; timeout = 999

  if (missing(om) || is.null(om) || !inherits(om, "Rceattle")) {
    stop("`om` must be an Rceattle model object (see ?fit_mod).")
  }
  if (missing(em) || is.null(em) || !inherits(em, "Rceattle")) {
    stop("`em` must be an Rceattle model object (see ?fit_mod).")
  }

  # DSEM inside an MSE is a real design question, not plumbing, so refuse rather
  # than produce a quietly wrong projection. Two things are unresolved:
  #
  #  * Under a DSEM the recruitment deviations are DERIVED from the latent
  #    states, so sample_rec()'s draws would have to be written into
  #    x_tj[, rec_dev_col], not into rec_dev. Writing rec_dev (what the code
  #    does today) is simply overwritten on the next evaluation.
  #  * .mse_trim_proj_params() / .mse_restore_proj_params() slice the parameter
  #    blocks that span styr:projyr. Whether x_tj is one of those depends on
  #    build_DSEM(estimate_projection): FALSE (the default) builds it over the
  #    hindcast only, so trimming it would corrupt it; TRUE spans projyr, so NOT
  #    trimming it misaligns the shortened refit. .mse_proj_param_yrdim is a
  #    static table and cannot express that conditional.
  #
  # Neither of the assessments driving this work uses run_mse() with a DSEM, so
  # this is deferred rather than guessed at.
  if (.has_dsem(om) || .has_dsem(em)) {
    stop("run_mse() does not yet support a DSEM. Fit the operating and ",
         "estimation models without one, or use the dev-DSEM branch.",
         call. = FALSE)
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # MSE SETUP ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  '%!in%' <- function(x,y)!('%in%'(x,y))
  set.seed(regenerate_seed)

  Rceattle_OM_list <- list()
  Rceattle_EM_list <- list()

  # * Input checks ----
  # - Set om to project from R0
  om$data_list$proj_mean_rec = TRUE # - Sample rec devs assuming this down the line

  # - Adjust cap
  if(!is.null(cap)){
    if(!length(cap) %in% c(1, om$data_list$nspp)){
      stop("cap is not length 1 or length nspp")
    }
  }

  if(!is.null(catch_mult)){
    if(length(catch_mult) == 1){
      catch_mult = rep(catch_mult, om$data_list$nspp)
    }

    if(length(catch_mult) != om$data_list$nspp){
      stop("catch_mult is not length 1 or length nspp")
    }
  }

  # na.rm: Proj_F_proportion is NA for fleets that never take catch (surveys),
  # which is a legitimate workbook value. Without it the sum is NA and the `if`
  # errors with "missing value where TRUE/FALSE needed" instead of reporting
  # whether any fleet is set to take projected F.
  if(sum(om$data_list$fleet_control$Proj_F_proportion, na.rm = TRUE) == 0){
    stop("F prop per fleet 'Proj_F_proportion' is zero")
  }

  # ** Refit OM if Proj_F_proportion was not activated ----
  if(sum((om$data_list$catch_data$Catch > 0) - (om$quantities$max_catch_hat > 0), na.rm = TRUE) > 0){
    # -- Set estimate mode back to original
    estimate_mode_base <- om$data_list$estimateMode

    # Rerun OM in debug mode to make sure F-prop is set correctly. Reuses the
    # OM's own configuration unchanged (estimateMode = 3 builds without
    # optimizing, so loopnum is inert here).
    om <- .refit_like(
      data_list    = om$data_list,
      inits        = om$estimated_params,
      map          = om$map,
      estimateMode = 3)

    # Adjust back
    om$data_list$estimateMode <- estimate_mode_base
  }

  # - Years for simulations
  hind_yrs <- (em$data_list$styr) : em$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  om_proj_yrs <- (om$data_list$endyr + 1) : om$data_list$projyr
  om_proj_nyrs <- length(om_proj_yrs)

  em_proj_yrs <- (em$data_list$endyr + 1) : em$data_list$projyr
  em_proj_nyrs <- length(em_proj_yrs)
  nflts = nrow(om$data_list$fleet_control)

  # - N sel ages for sel coff dev
  if(all(is.na(om$data_list$fleet_control$N_sel_bins))){
    n_sel_bins_om = dim(om$estimated_params$sel_coff_dev)[3]
  } else {
    n_sel_bins_om <- max(om$data_list$fleet_control$N_sel_bins, na.rm = TRUE)
  }

  if(all(is.na(em$data_list$fleet_control$N_sel_bins))){
    n_sel_bins_em = dim(em$estimated_params$sel_coff_dev)[3]
  } else {
    n_sel_bins_em <- max(em$data_list$fleet_control$N_sel_bins, na.rm = TRUE)
  }

  # - Assessment period
  assess_yrs <- seq(from = om$data_list$endyr + assessment_period, to =  min(c(om$data_list$projyr, em$data_list$projyr, endyr), na.rm = TRUE),  by = assessment_period)

  # - Data sampling period
  if(length(sampling_period)==1){
    sampling_period = rep(sampling_period, nflts)

  }

  if(nflts != nrow(em$data_list$fleet_control)){
    stop("OM and EM fleets do not match or sampling period length is mispecified")
  }
  if(nflts != length(sampling_period)){
    stop("Sampling period length is mispecified, does not match number of fleets")
  }

  # - Set up years of data we are sampling for each fleet
  sample_yrs <- lapply(sampling_period, function(x) seq(from = em$data_list$endyr + x, to = em$data_list$projyr,  by = x))
  fleet_id <- sample_yrs
  # sampling_period is given per fleet in fleet_control row order, but the table
  # built here is matched against the data's Fleet_code column. Carry the fleet's
  # own code rather than its row position: data_check() does require the two to
  # agree, but nothing here should depend on that silently.
  fleet_codes <- em$data_list$fleet_control$Fleet_code
  for(i in 1:length(sample_yrs)){
    fleet_id[[i]] <- replace(fleet_id[[i]], values = fleet_codes[i])
  }
  sample_yrs = data.frame(Fleet_code = unlist(fleet_id), Year = unlist(sample_yrs))


  # * Filter arbitrary "future" data ----
  # -- index_data
  om$data_list$index_data <- om$data_list$index_data |>
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$index_data <- em$data_list$index_data |>
    dplyr::filter(abs(Year) <= em$data_list$endyr)

  # -- comp_data
  om$data_list$comp_data <- om$data_list$comp_data |>
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$comp_data <- em$data_list$comp_data |>
    dplyr::filter(abs(Year) <= em$data_list$endyr)

  # -- caal_data
  om$data_list$caal_data <- om$data_list$caal_data |>
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$caal_data <- em$data_list$caal_data |>
    dplyr::filter(abs(Year) <= em$data_list$endyr)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Regenerate past data from OM and refit EM ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if(regenerate_past){

    # - Simulate index and comp data and updatae EM
    sim_dat <- sim_mod(om, simulate = FALSE)

    em$data_list$index_data <- sim_dat$index_data
    em$data_list$comp_data <- sim_dat$comp_data
    em$data_list$caal_data <- sim_dat$caal_data

    # Re-estimate
    em <- .refit_like(
      data_list    = em$data_list,
      inits        = em$estimated_params,
      estimateMode = ifelse(em$data_list$estimateMode < 3, 0, em$data_list$estimateMode),
      loopnum      = loopnum)

    # Update avg F given model fit to regenerated data
    if(em$data_list$HCR == 2){

      # - Get avg F
      avg_F <- exp(em$estimated_params$log_F) # Average F from last 5 years
      avg_F <- rowMeans(avg_F[,(ncol(avg_F)-4) : ncol(avg_F)])
      avg_F <- data.frame(avg_F = avg_F, spp = em$data_list$fleet_control$Species)
      avg_F <- avg_F |>
        dplyr::group_by(spp) |>
        dplyr::summarise(avg_F = sum(avg_F)) |>
        dplyr::arrange(spp)

      # - Update model: project on the recomputed average F (input-F HCR).
      em <- .refit_like(
        data_list    = em$data_list,
        inits        = em$estimated_params,
        estimateMode = 2,   # Don't estimate
        HCR          = build_hcr(HCR     = 2,   # Input F
                                 Ftarget = avg_F$avg_F,
                                 Ptarget = em$data_list$Ptarget,
                                 Plimit  = em$data_list$Plimit),
        loopnum      = loopnum)
    }
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Expand OM data-dim ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # -- index_data
  proj_srv <- om$data_list$index_data |>
    dplyr::group_by(Fleet_code) |>
    dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
    dplyr::mutate(Year = -om_proj_yrs)
  proj_srv$Log_sd <- proj_srv$Log_sd * 1/fut_sample
  proj_srv$Observation <- NA
  om$data_list$index_data  <- rbind(om$data_list$index_data, proj_srv)
  om$data_list$index_data <- dplyr::arrange(om$data_list$index_data, Fleet_code, abs(Year))

  # -- Nbyage
  if(nrow(om$data_list$NByageFixed) > 0){
    proj_nbyage <- om$data_list$NByageFixed |>
      dplyr::group_by(Species, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = om_proj_yrs)
    proj_nbyage <- proj_nbyage[which(om_proj_yrs %!in% om$data_list$NByageFixed$Year),] # Subset rows already forcasted
    om$data_list$NByageFixed  <- rbind(om$data_list$NByageFixed, proj_nbyage)
    om$data_list$NByageFixed <- dplyr::arrange(om$data_list$NByageFixed, Species, Year)
  }

  # -- comp_data
  if(nrow(om$data_list$comp_data) > 0){
    proj_comp <- om$data_list$comp_data |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = -om_proj_yrs)
    proj_comp$Sample_size <- proj_comp$Sample_size * fut_sample # Adjust future sampling effort
    proj_comp <- proj_comp |>
      dplyr::mutate_at(vars(matches("Comp_")), ~ 1)
    om$data_list$comp_data  <- rbind(om$data_list$comp_data, proj_comp)
    om$data_list$comp_data <- dplyr::arrange(om$data_list$comp_data, Fleet_code, abs(Year))
  }

  # -- caal_data
  if(nrow(om$data_list$caal_data) > 0){
    proj_caal <- om$data_list$caal_data |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = -om_proj_yrs)
    proj_caal$Sample_size <- proj_caal$Sample_size * fut_sample # Adjust future sampling effort
    proj_caal <- proj_caal |>
      dplyr::mutate_at(vars(matches("CAAL_")), ~ 1)
    om$data_list$caal_data  <- rbind(om$data_list$caal_data, proj_caal)
    om$data_list$caal_data <- dplyr::arrange(om$data_list$caal_data, Fleet_code, abs(Year))
  }

  # -- emp_sel - Use terminal year
  if(nrow(om$data_list$emp_sel) > 0){
    proj_emp_sel <- om$data_list$emp_sel |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = om_proj_yrs)
    om$data_list$emp_sel  <- rbind(om$data_list$emp_sel, proj_emp_sel)
    om$data_list$emp_sel <- dplyr::arrange(om$data_list$emp_sel, Fleet_code, Year)
  }

  # -- weight
  #FIXME ignores forecasted growth
  if(nrow(om$data_list$weight) > 0){
    proj_wt <- om$data_list$weight |>
      dplyr::group_by(Wt_index , Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = om_proj_yrs)
    om$data_list$weight  <- rbind(om$data_list$weight, proj_wt)
    om$data_list$weight <- dplyr::arrange(om$data_list$weight, Wt_index, Year)
  }

  # -- ration_data
  if(nrow(om$data_list$ration_data) > 0){
    proj_ration_data <- om$data_list$ration_data |>
      dplyr::group_by(Species, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = om_proj_yrs)
    om$data_list$ration_data  <- rbind(om$data_list$ration_data, proj_ration_data)
    om$data_list$ration_data <- dplyr::arrange(om$data_list$ration_data, Species, Year)
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Expand EM data-dim ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  #FIXME - assuming same as terminal year of hindcast
  # -- EM emp_sel - Use terminal year
  if(nrow(em$data_list$emp_sel) > 0){
    proj_emp_sel <- em$data_list$emp_sel |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  em_proj_nyrs)) |>
      dplyr::mutate(Year = em_proj_yrs)
    em$data_list$emp_sel  <- rbind(em$data_list$emp_sel, proj_emp_sel)
    em$data_list$emp_sel <- dplyr::arrange(em$data_list$emp_sel, Fleet_code, Year)
  }

  # -- EM weight
  if(nrow(em$data_list$weight) > 0){
    proj_wt <- em$data_list$weight |>
      dplyr::group_by(Wt_index , Sex) |>
      dplyr::slice(rep(n(),  em_proj_nyrs)) |>
      dplyr::mutate(Year = em_proj_yrs)
    em$data_list$weight  <- rbind(em$data_list$weight, proj_wt)
    em$data_list$weight <- dplyr::arrange(em$data_list$weight, Wt_index, Year)
  }

  # -- EM ration_data
  if(nrow(em$data_list$ration_data) > 0){
    proj_ration_data <- em$data_list$ration_data |>
      dplyr::group_by(Species, Sex) |>
      dplyr::slice(rep(dplyr::n(),  em_proj_nyrs)) |>
      dplyr::mutate(Year = em_proj_yrs)
    em$data_list$ration_data  <- rbind(em$data_list$ration_data, proj_ration_data)
    em$data_list$ration_data <- dplyr::arrange(em$data_list$ration_data, Species, Year)
  }


  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster.
  # We do *not* use foreach::%dopar% here: under nested test_that
  # backtraces it triggered an 'evaluation nested too deeply: infinite
  # recursion' abort inside rlang's expression deparser (foreach
  # captures call frames that recurse during error formatting). PSOCK
  # clusters work identically on Windows and Unix and avoid that.
  # Respect the CRAN core limit ('_R_CHECK_LIMIT_CORES_' is set during
  # R CMD check; parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- nsim > 1L && cores > 1L

  # TODO: extract run_one_sim as a top-level internal helper with
  # explicit args (om, em, seed, assess_yrs, ...) for testability.
  # Inline closure for now to ship the foreach -> parLapply migration.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-simulation closure (run_one_sim) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_sim <- function(sim) {

    set.seed(seed = seed + sim) # setting unique seed for each simulation
    kill_sim <- list(kill_sim = FALSE, failure = NA)

    # Set models objects
    sim_list <- list(EM = list())# , OM = list())
    sim_list$EM[[1]] <- em
    # sim_list$OM[[1]] <- om

    em_use <- em
    om_use <- om

    # Sample recruitment
    om_use <- Rceattle::sample_rec(om_use, sample_rec = sample_rec, update_model = FALSE, rec_trend = rec_trend)

    # Run through assessment years
    for(k in 1:length(assess_yrs)){

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 1. Get recommended catch from the EM-HCR ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      new_years <- om_proj_yrs[which(om_proj_yrs <= assess_yrs[k] & om_proj_yrs > om_use$data_list$endyr)]

      # - Get projected catch data from EM
      new_catch_data <- em_use$data_list$catch_data
      dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
      new_catch_data$Catch[dat_fill_ind] <- em_use$quantities$catch_hat[dat_fill_ind]

      # * Catch multiplier ----
      if(!is.null(catch_mult)){
        new_catch_data$Catch[dat_fill_ind] <- new_catch_data$Catch[dat_fill_ind] * catch_mult[new_catch_data$Species[dat_fill_ind]]
      }

      # * Apply cap ----
      if(!is.null(cap)){
        # Applied across species
        if(length(cap) == 1){
          new_catch_data$Catch[dat_fill_ind] <- ifelse(sum(new_catch_data$Catch[dat_fill_ind]) > cap,
                                                       cap * new_catch_data$Catch[dat_fill_ind]/sum(new_catch_data$Catch[dat_fill_ind]),
                                                       new_catch_data$Catch[dat_fill_ind]) # FIXME: does not work for assessments that don't occur annually
        } else { # Species-specific
          new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap[new_catch_data$Species[dat_fill_ind]], cap[new_catch_data$Species[dat_fill_ind]], new_catch_data$Catch[dat_fill_ind])
        }
      }

      # * Exploitable biomass limit ----
      # - If projected catch > exploitable biomass in OM, reduce to exploitable biomass
      exploitable_biomass_data <- om_use$data_list$catch_data
      exploitable_biomass_data$Catch[dat_fill_ind] <- om_use$quantities$max_catch_hat[dat_fill_ind]

      new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > exploitable_biomass_data$Catch[dat_fill_ind],
                                                   exploitable_biomass_data$Catch[dat_fill_ind],
                                                   new_catch_data$Catch[dat_fill_ind])

      new_catch_switch <- sum(new_catch_data$Catch[dat_fill_ind]) #Switch to turn off re-running OM if new catch = 0

      # - Update catch data in OM and EM
      om_use$data_list$catch_data <- new_catch_data
      em_use$data_list$catch_data <- new_catch_data


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 2. Update the OM ----
      # - Estimate Fdev and update dynamics
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # - Update endyr of OM
      nyrs_hind <- om_use$data_list$endyr - om_use$data_list$styr + 1
      om_use$data_list$endyr <- assess_yrs[k]

      # * Shorten the projection horizon for this refit ----
      # Reach one assessment step past the new terminal year -- far enough for
      # the next iteration's exploitable-biomass cap -- so the AD tape covers
      # the realized years plus that look-ahead instead of the whole projection.
      # The last assessment keeps the full horizon, so the operating model that
      # is returned (and handed to remove_F) is built over the whole projection,
      # exactly as if the horizon had never been shortened.
      om_dl_full     <- om_use$data_list
      om_params_full <- om_use$estimated_params
      om_projyr_use  <- if (k < length(assess_yrs)) assess_yrs[k + 1] else om_dl_full$projyr
      om_shortened   <- om_projyr_use < om_dl_full$projyr
      if (om_shortened) {
        om_use$data_list$projyr <- om_projyr_use
        om_use$estimated_params <- .mse_trim_proj_params(
          om_use$estimated_params,
          om_projyr_use - om_use$data_list$styr + 1)
      }

      # * Update parameters ----
      # -- log_F
      om_use$estimated_params$log_F <- cbind(om_use$estimated_params$log_F, matrix(0, nrow= nrow(om_use$estimated_params$log_F), ncol = length(new_years)))

      # -- M1_dev
      #FIXME - simulate
      # om_use$estimated_params$log_M1_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- om_use$estimated_params$log_M1_dev[,,,nyrs_hind]

      # -- Time-varying survey catchability - assume last year, filled by columns
      om_use$estimated_params$index_q_dev <- cbind(om_use$estimated_params$index_q_dev, matrix(om_use$estimated_params$index_q_dev[,ncol(om_use$estimated_params$index_q_dev)], nrow= nrow(om_use$estimated_params$index_q_dev), ncol = length(new_years)))

      # -- Time-varying selectivity - assume last year, filled by columns.
      # Use the fitted arrays' own sex extent (max_sex), not a hardcoded 2: a
      # single-sex model has sex-dim 1, and forcing 2 silently recycles the
      # fitted values into a phantom second sex (and, since build_params sizes
      # these by max_sex, mismatches the parameter template on the refit).
      n_sex_om <- dim(om_use$estimated_params$log_sel_slp_dev)[3]
      log_sel_slp_dev = array(0, dim = c(2, nflts, n_sex_om, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, n_sex_om, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_coff_dev = array(0, dim = c(nflts, n_sex_om, n_sel_bins_om, nyrs_hind + length(new_years)))  # selectivity deviation parameters for non-parametric

      log_sel_slp_dev[,,,1:nyrs_hind] <- om_use$estimated_params$log_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_inf_dev
      sel_coff_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_coff_dev

      log_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- log_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

      om_use$estimated_params$log_sel_slp_dev <- log_sel_slp_dev
      om_use$estimated_params$sel_inf_dev <- sel_inf_dev
      om_use$estimated_params$sel_coff_dev <- sel_coff_dev


      # * Update map ----
      # -(Only new parameter we are estimating in OM is the log_F of the new years)
      om_use$map <- build_map(
        data_list = om_use$data_list,
        params = om_use$estimated_params,
        debug = TRUE,
        random_rec = om_use$data_list$random_rec)
      om_use$map$mapFactor$dummy <- as.factor(NA); om_use$map$mapList$dummy <- NA


      # -- Estimate terminal F for catch
      new_f_yrs <- (ncol(om_use$map$mapList$log_F) - length(new_years) + 1) : ncol(om_use$map$mapList$log_F) # - Years of new F
      f_fleets <- om_use$data_list$fleet_control$Fleet_code[which(om_use$data_list$fleet_control$Fleet_type == "Fishery")] # Fleet rows for F
      om_use$map$mapList$log_F[f_fleets,new_f_yrs] <- replace(om_use$map$mapList$log_F[f_fleets,new_f_yrs], values = 1:length(om_use$map$mapList$log_F[f_fleets,new_f_yrs]))

      # -- Map out Fdev for years with 0 catch to very low number
      zero_catch <- om_use$data_list$catch_data |>
        dplyr::filter(Year <= om_use$data_list$endyr &
                        Catch == 0) |>
        dplyr::mutate(Year = Year - om_use$data_list$styr + 1) |>
        dplyr::select(Fleet_code, Year) |>
        as.matrix()
      om_use$estimated_params$log_F[zero_catch] <- -999
      om_use$map$mapList$log_F[zero_catch] <- NA
      om_use$map$mapFactor$log_F <- factor(om_use$map$mapList$log_F)
      rm(zero_catch)

      # -- Set estimate mode
      estimate_mode_base <- om_use$data_list$estimateMode
      estimate_mode_use <- ifelse(
        new_catch_switch == 0, 3, # Run in debug mode if catch is 0 for all species
        ifelse(
          estimate_mode_base < 3, 1, # Estimate hindcast only if estimating
          estimate_mode_base)
      )

      if(new_catch_switch == 0){
        om_use$map = NULL
      }

      # * Fit OM with new catch data ----
      kill_sim <- tryCatch({
        R.utils::withTimeout({
          suppressMessages(
            # Advance the OM with the new catch. The stock-recruit reference
            # period (srr_*) and the suitability window (suit_*) are pinned to
            # the PRISTINE om$ (not the advancing om_use$) so both stay fixed
            # through the projection -- critical for multispecies models, whose
            # predation suitability must not drift over the MSE. The OM projects
            # on mean recruitment across sim iterations (proj_mean_rec = TRUE).
            om_use <- .refit_like(
              data_list        = om_use$data_list,
              inits            = om_use$estimated_params,
              map              = om_use$map,
              estimateMode     = estimate_mode_use,
              loopnum          = loopnum,
              proj_mean_rec    = TRUE,
              srr_mse_switchyr = om$data_list$srr_mse_switchyr,
              srr_hat_styr     = om$data_list$srr_hat_styr,
              srr_hat_endyr    = om$data_list$srr_hat_endyr,
              suit_styr        = om$data_list$suit_styr,
              suit_endyr       = om$data_list$suit_endyr)
          )
          return(list(kill_sim = FALSE, failure = NA))
        },
        timeout = 60*timeout)
      },
      error = function(e){
        return(list(kill_sim = TRUE, failure = "OM"))
      },
      TimeoutException = function(e){
        return(list(kill_sim = TRUE, failure = "OM"))
      })

      if(kill_sim$kill_sim){
        # Nothing to put back: a killed simulation returns only its failure
        # marker, so the partially advanced operating model is discarded.
        break()
      }

      # -- Set estimate mode back to original
      om_use$data_list$estimateMode <- estimate_mode_base

      # sim_mod() reads the operating model's quantities positionally against
      # its own data frames, so it works from the fitted object; the assessment
      # loop below works from the full-horizon one.
      om_fit <- om_use
      if (om_shortened) {
        om_use <- .mse_restore_om_horizon(om_use, om_dl_full, om_params_full)
      }


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 3. Get actual catch from OM ----
      # - Maybe the OM can't support the TAC
      # (should be minimal with the < exploitable biomass check)
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      # - Get realized catch data from OM
      new_catch_data <- om_use$data_list$catch_data
      dat_fill_ind <- which(new_catch_data$Year %in% new_years)
      new_catch_data$Catch[dat_fill_ind] <- om_use$quantities$catch_hat[dat_fill_ind] # Catch from OM

      # - Update catch data in OM and EM
      om_use$data_list$catch_data <- new_catch_data
      em_use$data_list$catch_data <- new_catch_data


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 4. Simulate data from OM ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # - Simulate new survey and comp data
      sim_dat <- Rceattle::sim_mod(om_fit, simulate = simulate_data)

      years_include <- sample_yrs[which(sample_yrs$Year > em_use$data_list$endyr & sample_yrs$Year <= assess_yrs[k]),]

      # sample_yrs pairs each fleet with the years THAT fleet is sampled, so the
      # rows to add are the (fleet, year) pairs it lists -- not every fleet in
      # every year. Matching the year set and the fleet set separately is the
      # same thing only when one year advances per assessment; with a longer
      # assessment period it also admits an every-other-year fleet in its off
      # years, giving the estimation model survey and composition data the
      # sampling design says were never collected.
      sampled_key <- paste(years_include$Fleet_code, years_include$Year)

      # -- Add newly simulated survey data to EM and OM
      # - Get simulated survey data
      new_index_data <- sim_dat$index_data |>
        dplyr::filter(paste(Fleet_code, abs(Year)) %in% sampled_key) |>
        dplyr::mutate(Year = -Year)

      # - Add to EM and OM
      om_use$data_list$index_data <- om_use$data_list$index_data |>
        dplyr::filter(!(paste(Fleet_code, abs(Year)) %in% sampled_key)) |>
        rbind(new_index_data |>
                dplyr::mutate(Year = -abs(Year))) |>
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$index_data <- em_use$data_list$index_data |>
        rbind(new_index_data) |>
        dplyr::arrange(Fleet_code, abs(Year))


      # -- Add newly simulated comp data to EM & OM
      # - Simulated comp data
      new_comp_data <- sim_dat$comp_data |>
        dplyr::filter(paste(Fleet_code, abs(Year)) %in% sampled_key) |>
        dplyr::mutate(Year = -Year)

      new_comp_data$Sample_size <- new_comp_data$Sample_size * as.numeric(rowSums(dplyr::select(new_comp_data, dplyr::contains("Comp_"))) > 0) # Set sample size to 0 if catch is 0
      new_comp_data <- new_comp_data |>
        dplyr::mutate_at(dplyr::vars(dplyr::contains("Comp_")), ~ .x + 1 * (Sample_size == 0)) # Set all values to 1 if catch is 0

      # - Add to EM and OM
      om_use$data_list$comp_data <- om_use$data_list$comp_data |>
        dplyr::filter(!(paste(Fleet_code, abs(Year)) %in% sampled_key)) |>
        rbind(new_comp_data |>
                dplyr::mutate(Year = -abs(Year))) |>
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$comp_data <- em_use$data_list$comp_data |>
        rbind(new_comp_data) |>
        dplyr::arrange(Fleet_code, abs(Year))

      # -- Add newly simulated CAAL to EM & OM
      # - Simulated caal data
      new_caal_data <- sim_dat$caal_data |>
        dplyr::filter(paste(Fleet_code, abs(Year)) %in% sampled_key) |>
        dplyr::mutate(Year = -Year)

      new_caal_data$Sample_size <- new_caal_data$Sample_size * as.numeric(rowSums(dplyr::select(new_caal_data, dplyr::contains("CAAL_"))) > 0) # Set sample size to 0 if catch is 0
      new_caal_data <- new_caal_data |>
        dplyr::mutate_at(dplyr::vars(dplyr::contains("CAAL_")), ~ .x + 1 * (Sample_size == 0)) # Set all values to 1 if catch is 0

      # - Add to EM and OM
      om_use$data_list$caal_data <- om_use$data_list$caal_data |>
        dplyr::filter(!(paste(Fleet_code, abs(Year)) %in% sampled_key)) |>
        rbind(new_caal_data |>
                dplyr::mutate(Year = -abs(Year))) |>
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$caal_data <- em_use$data_list$caal_data |>
        rbind(new_caal_data) |>
        dplyr::arrange(Fleet_code, abs(Year))

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 5. Update EM and HCR ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # Update end year and re-estimate
      em_use$data_list$endyr <- assess_yrs[k]

      # Update parameter size and use previous estimates
      # -- log_F
      em_use$estimated_params$log_F <- cbind(em_use$estimated_params$log_F, matrix(0, nrow= nrow(em_use$estimated_params$log_F), ncol = length(new_years)))

      # # -- log_M1_dev
      # em_use$estimated_params$log_M1_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- em_use$estimated_params$log_M1_dev[,,,nyrs_hind]

      # -- Time-varying survey catchability - assume last year, filled by columns
      em_use$estimated_params$index_q_dev <- cbind(em_use$estimated_params$index_q_dev, matrix(em_use$estimated_params$index_q_dev[,ncol(em_use$estimated_params$index_q_dev)], nrow= nrow(em_use$estimated_params$index_q_dev), ncol = length(new_years)))

      # -- Time-varying selectivity - assume last year, filled by columns.
      # Sex extent from the fitted arrays (max_sex), not a hardcoded 2 -- see
      # the OM extension above.
      n_sex_em <- dim(em_use$estimated_params$log_sel_slp_dev)[3]
      log_sel_slp_dev = array(0, dim = c(2, nflts, n_sex_em, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, n_sex_em, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_coff_dev = array(0, dim = c(nflts, n_sex_em, n_sel_bins_em, nyrs_hind + length(new_years)))  # selectivity deviation parameters for non-parametric

      log_sel_slp_dev[,,,1:nyrs_hind] <- em_use$estimated_params$log_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_inf_dev
      sel_coff_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_coff_dev

      # - Initialize new years with last year
      log_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- log_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

      em_use$estimated_params$log_sel_slp_dev <- log_sel_slp_dev
      em_use$estimated_params$sel_inf_dev <- sel_inf_dev
      em_use$estimated_params$sel_coff_dev <- sel_coff_dev


      # Re-estimate
      kill_sim <- tryCatch({
        R.utils::withTimeout({
          suppressMessages(
            # Re-assess with the newly sampled data. srr_mse_switchyr switches
            # the stock-recruit relationship at the EM's CURRENT assessment
            # endyr (which advances each iteration), not the stored value. The
            # HCR is reconstructed from em_use$, whose HCRorder equals the
            # pristine em$ (copied at `em_use <- em`, never modified in the loop).
            em_use <- .refit_like(
              data_list        = em_use$data_list,
              inits            = em_use$estimated_params,
              estimateMode     = ifelse(em_use$data_list$estimateMode < 3, 0, em_use$data_list$estimateMode),
              loopnum          = loopnum,
              srr_mse_switchyr = em_use$data_list$endyr)
          )
          return(list(kill_sim = FALSE, failure = NA))
        },
        timeout = 60*timeout)
      },
      error = function(ex) {
        return(list(kill_sim = TRUE, failure = "EM"))
      },
      TimeoutException = function(ex) {
        return(list(kill_sim = TRUE, failure = "EM"))
      })

      # A failed re-assessment leaves em_use holding the PREVIOUS year's model,
      # so the simulation cannot carry on: it would store that stale assessment
      # under this year's name and hand its catch advice to the next iteration.
      # Stop the simulation the same way a failed operating-model refit does.
      if(kill_sim$kill_sim){
        break()
      }
      # plot_biomass(list(em_use, om_use), model_names = c("EM", "OM"))
      # End year of assessment

      # - Remove unneeded bits for memory
      em_use$initial_params <- NULL
      em_use$bounds <- NULL
      em_use$map <- NULL
      em_use$phase_params <- NULL
      em_use$obj <- NULL
      em_use$opt <- NULL
      em_use$sdrep <- NULL
      em_use$quantities[names(em_use$quantities) %!in% .mse_keep_quantities] <- NULL

      sim_list$EM[[k+1]] <- em_use
      message(paste0("Sim ",sim, " - EM Year ", assess_yrs[k], " COMPLETE"))
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 6. End year loop ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    }


    # - Rename models
    if (kill_sim$kill_sim) {
      # A simulation that broke off part way through never reached the terminal
      # assessment year, so its operating and estimation models describe a state
      # the MSE did not arrive at. Return only the marker: there is nothing here
      # a summary should average over, and keeping the partial models invites
      # them being read as if the simulation had run to completion.
      sim_list <- list(use_sim = FALSE, failure = kill_sim$failure)
    } else {
      sim_list$use_sim <- TRUE
      sim_list$failure <- NA
      sim_list$OM <- om_use # OM
      # The unfished reference run is a separate numerical problem and can fail
      # on its own (it sdreports a model with F pinned off across the
      # projection). The simulation itself is complete at this point, and its
      # assessments are the catch-advice record, so a failure here must not
      # discard it: dropping these would remove simulations in a
      # stock-state-dependent way and bias every performance metric.
      # Assigned through `[` rather than `$` so a failure leaves OM_no_F present
      # and NULL: `sim_list$OM_no_F <- NULL` would DELETE the element, and the
      # simulation would then be indistinguishable by name from one that never
      # attempted the unfished run.
      sim_list["OM_no_F"] <- list(tryCatch(remove_F(om_use), error = function(e) {
        # Recorded on the object, not just warned about: this runs in a parallel
        # worker, whose warnings are discarded, and the simulation is still
        # usable for everything that does not compare against the unfished run.
        # mse_summary() drops these from the no-F metrics only.
        sim_list$failure <<- paste0("OM_no_F: ", conditionMessage(e))
        NULL
      }))
      names(sim_list$EM) <- c("EM", paste0("OM_Sim_",sim,". EM_yr_", assess_yrs))
    }

    # - Save
    if(!is.null(dir)){
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(sim_list, file = paste0(dir, "/", file, "EMs_from_OM_Sim_",sim, ".rds"))
      # Hand back only the outcome, not the models: run_mse() returns NULL in
      # this mode, and the dispatcher needs to report attrition without reading
      # every saved simulation back off disk.
      list(use_sim = sim_list$use_sim, failure = sim_list$failure)
    } else{
      sim_list # Return simlist
    }
  } # End run_one_sim closure


  # Contain a simulation's failure to that simulation. The refits are already
  # guarded, but the surrounding data reshaping is not, and neither the
  # sequential nor the parallel dispatch has a per-item handler -- so one
  # unguarded error would discard every other simulation of the run. Anything
  # that escapes is reported as a failed simulation, keeping the message so the
  # cause is still visible.
  run_one_sim_guarded <- function(sim) {
    tryCatch(run_one_sim(sim), error = function(e) {
      # Tagged so this stays distinguishable from the enumerated "OM"/"EM"
      # failures the assessment loop reports: those are expected outcomes, this
      # is a bug or a bad input.
      marker <- list(use_sim = FALSE,
                     failure = paste0("unexpected: ", conditionMessage(e)))
      message("Sim ", sim, " FAILED: ", conditionMessage(e))
      if (!is.null(dir)) {
        # Writing the marker must not itself escape, or the containment this
        # handler exists for is lost for every other simulation.
        tryCatch({
          dir.create(dir, showWarnings = FALSE, recursive = TRUE)
          saveRDS(marker,
                  file = paste0(dir, "/", file, "EMs_from_OM_Sim_", sim, ".rds"))
        }, error = function(e2) {
          warning("Sim ", sim, ": could not record the failure (",
                  conditionMessage(e2), ").", call. = FALSE)
        })
      }
      marker
    })
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch sims (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    # FORK where possible (inherits the loaded package + the large OM/EM objects
    # via copy-on-write); PSOCK fallback exports them. See .parallel_lapply().
    sim_list <- .parallel_lapply(start_sim:nsim, run_one_sim_guarded, min(cores, nsim), environment())
  } else {
    sim_list <- lapply(start_sim:nsim, run_one_sim_guarded)
  }

  names(sim_list) <- paste0("Sim_", start_sim:nsim)

  # Report attrition here rather than leaving it to whatever reads the results.
  # With `dir` set the return value is NULL, so a run in which every simulation
  # failed would otherwise be indistinguishable from one in which every
  # simulation succeeded -- the same number of files either way.
  n_failed <- sum(vapply(sim_list, function(x) isFALSE(x$use_sim), logical(1)))
  if (n_failed) {
    warning(n_failed, " of ", length(sim_list),
            " simulations did not complete; they carry use_sim = FALSE and no ",
            "models. Filter on use_sim before summarising.", call. = FALSE)
  }
  # Worker warnings are discarded by the parallel cluster, so a completed
  # simulation whose unfished reference run failed is reported from here.
  n_no_f <- sum(vapply(sim_list, function(x) {
    isTRUE(x$use_sim) && grepl("^OM_no_F: ", paste(x$failure))
  }, logical(1)))
  if (n_no_f) {
    warning(n_no_f, " of ", length(sim_list),
            " simulations completed but their unfished reference run failed; ",
            "these carry use_sim = TRUE and OM_no_F = NULL, and are excluded ",
            "from the no-F metrics only.", call. = FALSE)
  }

  if(is.null(dir)){
    return(sim_list)
  }
}
