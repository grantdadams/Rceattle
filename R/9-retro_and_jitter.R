#' Parallel `lapply` over a cluster, using FORK where the platform allows.
#'
#' @description
#' On non-Windows platforms a FORK cluster inherits the parent process's memory
#' (the already-loaded `Rceattle` namespace and every local, including the large
#' fitted OM/EM objects) via copy-on-write, so it needs neither a per-worker
#' `library(Rceattle)` nor a `clusterExport()` of those objects. The PSOCK path
#' pays both on every worker: a package load and serialization/transfer of the
#' exported bindings. FORK therefore removes the dominant cluster-startup cost
#' for the retrospective / jitter / MSE dispatchers. PSOCK is the cross-platform
#' fallback (Windows), reproducing the previous behavior exactly.
#'
#' A cluster that dies takes the whole call with it, which is worth guarding
#' because the cluster is where the platforms differ: each worker fits complete
#' models, so several large ones at once can exhaust memory and the operating
#' system kills a worker. The parent then reports only what it noticed --
#' `Error in unserialize(node$con) : error reading from connection` -- and every
#' peel or simulation computed so far is lost. Rather than fail, fall back to
#' running the items sequentially, as `osa_residuals()` does when its parallel
#' one-step-ahead loop fails. The run is slower but finishes, and the message
#' names `cores` so the next call can skip the wasted attempt. Note the fallback
#' restarts the work from the beginning: cluster results are not recoverable
#' once a worker has gone.
#'
#' An error raised by `fun` itself propagates out of the sequential retry
#' unchanged, so a real bug still surfaces as itself rather than as a cluster
#' failure.
#'
#' @param items Vector iterated over; each element is passed to `fun`.
#' @param fun Worker closure.
#' @param n_workers Number of cluster workers.
#' @param export_env For the PSOCK fallback only, the environment whose bindings
#'   are exported to the workers (the caller's `environment()`). Ignored for
#'   FORK, where the workers inherit it directly.
#' @return A list of `fun` applied to each element of `items`.
#' @noRd
.parallel_lapply <- function(items, fun, n_workers, export_env) {
  fork <- .Platform$OS.type != "windows"

  run_clustered <- function() {
    cl <- if (fork) {
      parallel::makeCluster(n_workers, type = "FORK")
    } else {
      parallel::makeCluster(n_workers)
    }
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (!fork) {
      parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(Rceattle)))
      parallel::clusterExport(cl, varlist = ls(envir = export_env), envir = export_env)
    }
    parallel::parLapply(cl, items, fun)
  }

  tryCatch(run_clustered(), error = function(e) {
    message("Parallel execution failed (", conditionMessage(e), "); ",
            "running the ", length(items), " tasks sequentially instead. ",
            "A worker that dies this way has usually run out of memory -- ",
            "each one fits a full model. Pass cores = 1 (or a smaller cores) ",
            "to go straight to this path.")
    lapply(items, fun)
  })
}

#' Retrospective peels
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model. The function also evaluates retrospective forecast skill. To evaluate both retrospective bias and forecast skill, the function uses the map functionality of TMB to peel the model:
#' 1. Filters data, filters fixed inputs, and maps out time-varying parameters for the peeled years. All time-varying parameters for the peeled years are set to the terminal year of the model for that peel.
#' 2. Fits the peeled model.
#' 3. Turns off all hindcast parameters, turns on F for the peeled years, and fits to the peeled catch series to update the "forecast" dynamics given projection assumptions and observed catch from the peeled years.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param peels the number of retrospective peels to use in the calculation of rho and for model estimation
#' @param rescale TRUE/FALSE whether to subset and rescale environmental predictors for the range of peel years.
#' @param nyrs_forecast Number of forecast years to calculate Mohn's Rho in addition to terminal year
#' @param cores Number of cores to use for parallel peels. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each peel runs \code{TMB::sdreport} (standard errors).
#'   Mohn's rho uses only point estimates, so \code{FALSE} is faster with no
#'   effect on rho. Default \code{NULL} inherits the input model's setting
#'   (\code{TRUE} if it was fit with \code{getsd = TRUE}, i.e. carries an
#'   \code{sdrep}); the returned peel models then carry standard errors only
#'   when \code{getsd} is \code{TRUE}.
#'
#' @return a list of 1. list of Rceattle models and 2. vector of Mohn's rho for
#'   each species.
#'
#'   A peel that did not converge is dropped, so \code{Rceattle_list} can be
#'   shorter than \code{peels + 1} (a message reports how many). Each entry is
#'   named for its own terminal year (\code{Year_2017}, ...) rather than by
#'   position, so index it by name -- \code{Rceattle_list[[3]]} is not
#'   necessarily the 3-year peel. With no peel left, Mohn's rho is \code{NaN}
#'   and the function warns.
#'
#'   Each peel reports its own terminal year as \code{data_list$endyr}, so plots
#'   draw it only as far as it was fit and the peels fan out.
#'
#'   A peel still estimates the years it dropped: they are its retrospective
#'   forecast, fit to the observed catch with recruitment held at the peel's mean
#'   and the survey and composition data withheld. Three years therefore matter,
#'   and each peel carries all three:
#'   \describe{
#'     \item{\code{endyr}, \code{endyr_peel}}{the peel's terminal year -- what it
#'       was fit through. Equal to each other.}
#'     \item{\code{endyr_full}}{the unpeeled model's terminal year, where the
#'       retrospective forecast ends.}
#'     \item{\code{projyr}}{the end of the harvest-control-rule projection.}
#'   }
#'   So the forecast years are those after \code{endyr_peel} through
#'   \code{endyr_full}, and the projection follows through \code{projyr};
#'   \code{incl_proj = TRUE} plots both. Take the forecast years as
#'   \code{endyr_peel + seq_len(endyr_full - endyr_peel)}, which is empty for the
#'   unpeeled model, rather than \code{(endyr_peel + 1):endyr_full}, which counts
#'   \emph{down} there.
#'
#'   Mohn's rho is computed from \code{endyr_peel} and is unaffected by any of
#'   this.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' retro <- retrospective(ss_run, peels = 10)
#' }
#' @export
retrospective <- function(Rceattle = NULL, peels = 5, rescale = FALSE, nyrs_forecast = 3, cores = NULL, getsd = NULL) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # Peels inherit the input model's sdreport setting unless overridden. Mohn's
  # rho reads only point estimates, so getsd = FALSE is faster and rho-neutral.
  if (is.null(getsd)) getsd <- !is.null(Rceattle$sdrep)

  # Get objects
  Rceattle$data_list$endyr_peel <- Rceattle$data_list$endyr
  # Terminal year of the model being peeled. Each peel reports its OWN terminal
  # year as `endyr` (see run_one_peel), which makes the plots fan out but leaves
  # `endyr` and `endyr_peel` holding the same value -- so without this the
  # unpeeled terminal year is unrecoverable from a peel, and with it the boundary
  # between the retrospective FORECAST years, (endyr_peel + 1):endyr_full, and
  # the true projection, (endyr_full + 1):projyr. Set once here: run_one_peel
  # copies this data_list, and extra fields survive the refits the same way
  # `endyr_peel` already does.
  Rceattle$data_list$endyr_full <- Rceattle$data_list$endyr
  data_list <- Rceattle$data_list # used by Mohn's rho block below
  endyr <- Rceattle$data_list$endyr
  styr <- Rceattle$data_list$styr
  nyrs <- length(styr:endyr)
  projyr <- Rceattle$data_list$projyr
  nyrs_proj <- projyr - styr + 1

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- peels > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-peel closure ----
  # Each peel only reads the original Rceattle, so peels are independent.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_peel <- function(i) {

    # * Get end year of peel ----
    data_list <- Rceattle$data_list
    endyr_peel <- endyr - i
    data_list$endyr_peel <- endyr_peel
    nyrs_peel <- endyr_peel - styr + 1
    peel_prj_yrs <- (endyr_peel+1):endyr
    nyrs_proj_peel <- length(peel_prj_yrs)


    # * Turn off data after endyr_peel ----
    data_list$index_data <- data_list$index_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$comp_data <- data_list$comp_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$caal_data <- data_list$caal_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$diet_data <- data_list$diet_data |>
      dplyr::filter(Year <= endyr_peel)

    peeled_catch_data <- data_list$catch_data |>
      dplyr::filter(Year > endyr_peel)
    data_list$catch_data$Catch[which(data_list$catch_data$Year > endyr_peel)] <- 0 # Set catch data in peeled years to 0 to avoid fitting


    # * Turn off fixed inputs after endyr_peel ----
    #FIXME ignores forecasted growth
    data_list$weight <- data_list$weight |>
      dplyr::filter(Year <= endyr_peel)

    data_list$emp_sel <- data_list$emp_sel |>
      dplyr::filter(Year <= endyr_peel)

    data_list$ration_data <- data_list$ration_data |>
      dplyr::filter(Year <= endyr_peel)

    # * Extend fixed inputs for "projection years"
    # - Assume weight/ration/empirical sel is same as last year of peel

    # - Weight
    #FIXME ignores forecasted growth
    proj_wt <- data_list$weight |>
      dplyr::filter(Year != 0)

    if(nrow(proj_wt) > 0){
      proj_wt <- proj_wt |>
        dplyr::group_by(Wt_index , Sex) |>
        dplyr::slice(rep(dplyr::n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
    }
    data_list$weight  <- rbind(data_list$weight, proj_wt) |>
      dplyr::arrange(Wt_index, Year)

    # - Empirical selectivity
    proj_emp_sel <- data_list$emp_sel |>
      dplyr::filter(Year != 0)

    if(nrow(proj_emp_sel) > 0){
      proj_emp_sel <- proj_emp_sel |>
        dplyr::group_by(Fleet_code, Sex) |>
        dplyr::slice(rep(n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
      data_list$emp_sel  <- rbind(data_list$emp_sel, proj_emp_sel) |>
        dplyr::arrange(Fleet_code, Year)
    }

    # - Ration data
    proj_ration_data <- data_list$ration_data |>
      dplyr::filter(Year != 0)

    if(nrow(proj_ration_data) > 0){
      proj_ration_data <- proj_ration_data |>
        dplyr::group_by(Species, Sex) |>
        dplyr::slice(rep(dplyr::n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
      data_list$ration_data  <- rbind(data_list$ration_data, proj_ration_data) |>
        dplyr::arrange(Species, Year)
    }

    # * Rescale environmental predictors ----
    if(rescale){
      data_list$env_data <- data_list$env_data |>
        dplyr::filter(Year <= endyr_peel)
      data_list$env_data[,2:ncol(data_list$env_data)]<-scale(data_list$env_data[,2:ncol(data_list$env_data)])
    }

    # * Adjust parameters ----
    #FIXME: adjust for forecasting via MVN
    inits <- Rceattle$estimated_params
    inits$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- 0
    inits$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- inits$log_M1_dev[,,,nyrs_peel]
    inits$index_q_dev[,(nyrs_peel+1):nyrs] <- inits$index_q_dev[,nyrs_peel]
    inits$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- inits$log_sel_slp_dev[,,,nyrs_peel]
    inits$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_inf_dev[,,,nyrs_peel]
    inits$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_coff_dev[,,,nyrs_peel]

    # * Adjust map size ----
    # Turn off forecasted parameters
    map <- Rceattle$map
    map$mapList$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- NA
    map$mapFactor$rec_dev <- factor(map$mapList$rec_dev)

    map$mapList$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- NA
    map$mapFactor$log_M1_dev <- factor(map$mapList$log_M1_dev)

    map$mapList$index_q_dev[,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$index_q_dev <- factor(map$mapList$index_q_dev)

    map$mapList$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$log_sel_slp_dev <- factor(map$mapList$log_sel_slp_dev)

    map$mapList$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_inf_dev <- factor(map$mapList$sel_inf_dev)

    map$mapList$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)

    # -- Map out Fdev for years with 0 catch to very low number
    zero_catch <- data_list$catch_data |>
      dplyr::filter(Year <= endyr &
                      Catch == 0) |>
      dplyr::mutate(Year = Year - styr + 1) |>
      dplyr::select(Fleet_code, Year) |>
      as.matrix()
    inits$log_F[zero_catch] <- -999
    map$mapList$log_F[zero_catch] <- NA
    map$mapFactor$log_F <- factor(map$mapList$log_F)

    # * Refit ----
    newmod <- suppressWarnings(
      suppressMessages(
        # Refit the peeled data_list, reusing its HCR / SR / M / growth
        # configuration; clamp the SR-switch, SR-fit-end, and suitability-end
        # years back to this peel's terminal year.
        .refit_like(
          data_list        = data_list,
          inits            = inits,
          map              = map,
          estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
          phase            = TRUE,   # phasing, or the parameters dont wanna move
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
          srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
          suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))
      )
    )

    # Judge the peeled hindcast NOW. `newmod` is overwritten below by the
    # forecast-catch refit, whose map (build_map(debug = TRUE), plus the peeled
    # years' log_F) leaves only a handful of F parameters free -- so its gradient
    # says nothing about whether the hindcast converged. Mohn's rho reads the
    # hindcast quantities, which come from these parameters, so a peel that
    # diverged here has to be dropped even when the F refit lands cleanly.
    hindcast_converged <- .refit_converged(newmod)

    # Forecast ----
    peeled_pars <- newmod$estimated_params

    # - Add in peeled catch to fit to
    data_list$catch_data <- data_list$catch_data |>
      dplyr::filter(Year <= endyr_peel) |>
      rbind(peeled_catch_data) |>
      dplyr::arrange(Fleet_code, Year)

    # - Update map
    # Only new parameter we are estimating in is the log_F of the peeled years
    map <- build_map(
      data_list = data_list,
      params = newmod$estimated_params,
      debug = TRUE,
      random_rec = newmod$data_list$random_rec)
    map$mapFactor$dummy <- as.factor(NA); map$mapList$dummy <- NA

    # - Turn on F for peeled years to fit to catch (matches full model)
    peeled_pars$log_F[,(nyrs_peel+1):nyrs] <- Rceattle$estimated_params$log_F[,(nyrs_peel+1):nyrs]
    map$mapList$log_F[,(nyrs_peel+1):nyrs] <- Rceattle$map$mapList$log_F[,(nyrs_peel+1):nyrs]
    map$mapFactor$log_F <-  factor(map$mapList$log_F)

    # Adjust forecased rec_dev in new mod for bias and refit
    for(sp in 1:newmod$data_list$nspp){

      # -- where SR curve is estimated directly
      if(newmod$data_list$srr_fun == newmod$data_list$srr_pred_fun){
        rec_dev <- log(mean(newmod$quantities$R[sp,1:nyrs_peel]))  - log(newmod$quantities$R0[sp])
      }

      # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
      if(newmod$data_list$srr_fun != newmod$data_list$srr_pred_fun){
        # Already a log-scale deviation, so take the mean directly.
        rec_dev <- mean((log(newmod$quantities$R) - log(newmod$quantities$R_hat))[sp, 1:nyrs_peel])

      }

      # - Update OM with devs
      peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)] <- replace(
        peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)],
        values =  rec_dev)
    }

    newmod <- suppressMessages(
      suppressWarnings(
        # Second refit: same peeled configuration, now started from the
        # bias-adjusted peeled parameters with peeled-year F turned back on.
        .refit_like(
          data_list        = data_list,
          inits            = peeled_pars,
          map              = map,
          estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
          phase            = TRUE,   # phasing, or the parameters dont wanna move
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
          srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
          suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))
      )
    )

    # gc()
    #
    # map$mapFactor <- map$mapFactor[names(newmod$map$mapFactor)]
    # check <- c()
    # check_na <- c()
    # for(j in 1:length(map$mapList)){
    #   check[j] <- sum(map$mapFactor[[j]] != newmod$map$mapFactor[[j]], na.rm = TRUE)
    #   check_na[j] <- sum(is.na(map$mapFactor[[j]]) != is.na(newmod$map$mapFactor[[j]]), na.rm = TRUE)
    # }

    # * Report the peel's own terminal year ----
    # Set here, AFTER both refits, and never before: `endyr` sizes the model,
    # and the forecast refit above turns F back on over `(nyrs_peel+1):nyrs`
    # against the FULL nyrs, so peeling it earlier would index off the end of
    # log_F. At this point it is output metadata only.
    #
    # Every plot builds its year axis per model as `styr:endyr`
    # (`R/7-plot_ceattle.R`), and nothing outside this file reads `endyr_peel`.
    # Without this each peel was drawn to the full model's terminal year, so
    # the peels were indistinguishable -- the opposite of what a retrospective
    # plot is for.
    #
    # Mohn's rho is unaffected: it reads `endyr_peel` off each peel and the
    # full model's `endyr` from this function's enclosing scope, never a peel's
    # `data_list$endyr`.
    #
    # Note this makes the returned peel deliberately inconsistent: its
    # parameters, quantities, and `catch_data` still span the full hindcast,
    # because the peeled years are its retrospective FORECAST. `endyr` marks
    # what was fit, not what was estimated. Plot with `incl_proj = TRUE` to see
    # the forecast years, and read `endyr_full` (carried through from the source
    # model) for where those forecast years end.
    newmod$data_list$endyr <- endyr_peel

    # Return model only if BOTH refits converged, else NULL (dropped
    # post-dispatch)
    if (hindcast_converged && .refit_converged(newmod)) {
      return(newmod)
    }
    return(NULL)
  } # End run_one_peel closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch peels (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    peel_results <- .parallel_lapply(1:peels, run_one_peel, min(cores, peels), environment())
  } else {
    peel_results <- lapply(1:peels, run_one_peel)
  }

  # Drop non-converged peels and prepend the original model
  peel_results <- peel_results[!vapply(peel_results, is.null, logical(1))]
  .report_dropped(peels - length(peel_results), peels, "peel")
  mod_list <- c(list(Rceattle), peel_results)

  # Mohn's rho averages the peels against the full model, so with none left the
  # sums below stay at their initialized zeros and every species column comes
  # back 0/0. Say so rather than return a table of NaN that looks computed.
  if (length(peel_results) == 0L) {
    warning("No peel converged, so Mohn's rho is undefined (every value NaN). ",
            "Inspect a peel with retrospective(..., peels = 1) and read its ",
            "$convergence.", call. = FALSE)
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Calculate Mohs rho ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # * Data frame to save ----
  objects <- c("biomass", "ssb", "R", "F_spp")
  mohns <- data.frame(matrix(0, nrow = length(objects) * (nyrs_forecast+1), ncol = 3 + data_list$nspp))
  colnames(mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  # * Loop through peels ----
  # seq_len, not 1:(n-1): with every peel dropped `mod_list` holds only the
  # original model and `1:0` would count DOWN to mod_list[[2]].
  for (i in seq_len(length(mod_list) - 1)) {
    endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
    ind <- 1

    # * Loop output ----
    for (j in 1:length(objects)) {
      term_quantities <- mod_list[[1]]$quantities[[objects[j]]]
      retro_quantities <- mod_list[[i + 1]]$quantities[[objects[j]]]

      # * Loop forecast years ----
      for(yr in 0:nyrs_forecast){

        # If data exist for forecast (save)
        if(endyr_peel + yr <= endyr){

          # * Get full and peeled models ----
          base <- term_quantities[, nyrs_peel + yr]
          peel <- retro_quantities[, nyrs_peel + yr]
          rel_error <- ((peel - base)/base)

          # * Save and sum relative error ----
          mohns[ind, 1] <- objects[j]         # Object
          mohns[ind, 2] <- yr                 # Year
          mohns[ind, 3] <- mohns[ind, 3] + 1  # N
          mohns[ind, 4:(data_list$nspp + 3) ] <- mohns[ind, 4:(data_list$nspp + 3)] + rel_error # Relative error
        }
        ind = ind+1
      }
    }
  }

  # * Divide N ----
  mohns[, 4:(data_list$nspp + 3) ] <- mohns[, 4:(data_list$nspp + 3)]/mohns[, 3]



  # # * Beta coefficients ----
  # objects <- colnames(Rceattle$estimated_params$beta_rec_pars)
  # beta_mohns <- data.frame(matrix(0, nrow = length(objects), ncol = 3 + data_list$nspp))
  # colnames(beta_mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  #
  #   # * Loop through peels ----
  #   for (i in 1:(length(mod_list) - 1)) {
  #     endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
  #     nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
  #     ind <- 1
  #
  #     # * Loop output ----
  #     for (j in 1:length(objects)) {
  #       base <- mod_list[[1]]$estimated_params$beta_rec_pars[,j]
  #       peel <- mod_list[[i + 1]]$estimated_params$beta_rec_pars[,j]
  #       rel_error <- ((peel - base)/base)
  #
  #       # * Save and sum relative error ----
  #       beta_mohns[j, 1] <- objects[j]        # Object
  #       beta_mohns[j, 2] <- 0                 # Year
  #       beta_mohns[j, 3] <- beta_mohns[j, 3] + 1   # N
  #       beta_mohns[j, 4:(data_list$nspp + 3) ] <- beta_mohns[j, 4:(data_list$nspp + 3)] + rel_error # Relative error
  #     }
  #   }

  # * Divide N ----
  # beta_mohns[, 4:(data_list$nspp + 3) ] <- beta_mohns[, 4:(data_list$nspp + 3)]/beta_mohns[, 3]

  mod_list <- rev(mod_list)
  # Name each entry from its own terminal year rather than assuming all `peels`
  # survived: a dropped peel then leaves a gap in the labels instead of shifting
  # every later one onto the wrong model (and `Year_...` running one short of the
  # vector, which errored). The original model carries endyr_peel = endyr, set
  # at the top of the function.
  names(mod_list) <- paste0(
    "Year_",
    vapply(mod_list, function(x) as.numeric(x$data_list$endyr_peel), numeric(1)))

  return(list(Rceattle_list = mod_list, mohns = rbind(mohns))) #, beta_mohns)))
}




#' Jitter analysis
#'
#' @description Refits the Rceattle model from starting values perturbed by N(0, sd) around the model's initial (pre-fit) parameters, to check convergence robustness.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param njitter the number of jitters to run
#' @param sd standard deviation for jitter (default = 0.2)
#' @param phase as in \code{\link{fit_mod}} default = FALSE. Jitters restart from
#'   perturbed \emph{starting} values, so a model that needed phasing to fit its
#'   real data needs it here too; leave this at \code{FALSE} for such a model and
#'   the jitters end far from any optimum and are dropped as non-converged, which
#'   reads as multimodality rather than as an unphased fit.
#' @param seed random number seed. Each jitter \code{i} uses \code{seed + i}
#'   so results are reproducible under both sequential and parallel execution.
#' @param cores Number of cores to use for parallel jitters. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each jitter runs \code{TMB::sdreport}. Jitter compares
#'   objectives and point estimates across starts, so \code{FALSE} is faster
#'   with no effect on that comparison. Default \code{NULL} inherits the input
#'   model's setting (\code{TRUE} only if it carries an \code{sdrep}).
#' @param timeout elapsed-second limit per jitter, \code{Inf} (default) for none.
#'   A jitter is a deliberately perturbed start and the optimizer runs with no
#'   iteration cap, so this is the diagnostic most likely to send one somewhere
#'   pathological and stall the whole run -- a hang no convergence check can
#'   catch, because the fit never returns. One that exceeds the limit is stopped,
#'   counted as non-converged and reported separately. Approximate: the limit is
#'   checked when control returns to R, so it fires between the optimizer's
#'   function evaluations rather than inside one.
#'
#' @return a list of 1. \code{Rceattle_list}, the converged jitters, and
#'   2. \code{nll}, their objective values. Non-converged (or timed-out) starts
#'   are dropped and reported in a message, so both can be shorter than
#'   \code{njitter} -- and that count is itself the result, since the whole point
#'   is what fraction of random starts reach the same optimum.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' jitters <- jitter(ss_run, njitter = 10)
#' }
#' @export
jitter <- function(Rceattle = NULL, njitter = 50, sd = 0.2, phase = FALSE, seed = 123, cores = NULL, getsd = NULL, timeout = Inf) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # Jitters inherit the input model's sdreport setting unless overridden;
  # multimodality is judged from objectives and point estimates, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(Rceattle$sdrep)

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- njitter > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-jitter closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_jitter <- function(i) {

    set.seed(seed + i) # unique seed per jitter for reproducibility under parallel

    # * Adjust initial values ----
    inits <- Rceattle$initial_params
    mapList <- Rceattle$map$mapList
    data_list <- Rceattle$data_list

    for(j in 1:length(inits)){
      par <- names(inits)[j]
      inits[[j]] <- replace(inits[[j]],
                            values = ifelse(is.na(as.numeric(mapList[[par]])),
                                            as.numeric(inits[[j]]),
                                            as.numeric(inits[[j]]) + stats::rnorm(length(as.numeric(inits[[j]])), 0, sd))
      )
    }


    # * Refit ----
    # Bounded and error-trapped: a jitter is a deliberately perturbed start, so
    # it is the diagnostic most likely to send the optimizer somewhere
    # pathological -- and one such start must not abort the other njitter - 1.
    newmod <- .refit_with_timeout(
      suppressMessages(
        suppressWarnings(
          # Refit from the jittered starting values (map rebuilt from scratch).
          .refit_like(
            data_list        = data_list,
            inits            = inits,
            estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
            phase            = phase,
            getsd            = getsd,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
            suit_endyr       = pmin(data_list$suit_endyr, data_list$endyr))
        )
      ),
      timeout = timeout)

    # Report the verdict beside the model; the dispatcher filters below.
    if (inherits(newmod, "condition")) {
      return(list(model = NULL, converged = FALSE, error = conditionMessage(newmod)))
    }
    list(model = newmod, converged = .refit_converged(newmod))
  } # End run_one_jitter closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch jitters (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    mod_list <- .parallel_lapply(1:njitter, run_one_jitter, min(cores, njitter), environment())
  } else {
    mod_list <- lapply(1:njitter, run_one_jitter)
  }

  # Drop non-converged
  converged <- vapply(mod_list, function(x) isTRUE(x$converged), logical(1))
  errs      <- vapply(mod_list, function(x) x$error %||% NA_character_, character(1))
  mod_list  <- lapply(mod_list, function(x) x$model)[converged]
  .report_dropped(sum(!converged), length(converged), "jitter")
  .report_errors(errs, "jitter")

  # vapply, not sapply: over an empty list sapply returns list(), and the
  # documented `min(jitters$nll)` / `hist(log(nll - min(nll)))` then errors
  # rather than reporting that nothing converged. numeric(0) keeps `nll`
  # numeric at every length.
  jnll <- vapply(mod_list, function(x) as.numeric(x$opt$objective), numeric(1))
  if (length(mod_list) > 0) {
    names(mod_list) <- paste0("Jitter_", seq_along(mod_list))
  }


  # Return ----
  return(list(Rceattle_list = mod_list, nll = unname(jnll)))
}



#' Self-test simulation analysis
#'
#' @description Simulates data from an Rceattle model and refits the model to the simulated data, to check that the fitting procedure recovers the operating-model parameters. TODO add process variation (i.e. random recruitment deviations) to the simulation.
#'
#' @param Rceattle an Rceattle model fit using \code{\link{fit_mod}}
#' @param seed random number seed. Each simulation \code{i} uses \code{seed + i}
#'   so results are reproducible under both sequential and parallel execution.
#' @param nsim number of simulations
#' @param simulate passed to \code{\link{sim_mod}}. If \code{TRUE} (default),
#'   data are simulated with observation error; if \code{FALSE}, expected
#'   values from the model are used.
#' @param cores Number of cores to use for parallel simulations. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each refit runs \code{TMB::sdreport}. Self-test compares
#'   the refit point estimates to the operating model, so \code{FALSE} is faster
#'   with no effect on that comparison. Default \code{NULL} inherits the input
#'   model's setting (\code{TRUE} only if it carries an \code{sdrep}).
#' @param phase as in \code{\link{fit_mod}}. Under the default
#'   \code{start = "initial"} each refit covers the same ground the original fit
#'   did, so a model that needed phasing to fit its real data needs it again for
#'   every simulated one -- without it such a model's refits can end many orders
#'   of magnitude from a zero gradient and be dropped as non-converged. Default
#'   \code{NULL} reads the setting \code{fit_mod()} recorded on the source fit
#'   (\code{fit$run_config$fit_control$phase}), so a model fitted under the
#'   package default of \code{phase = FALSE} is refitted unphased; pass
#'   \code{TRUE} for a model that needs phasing but was not fitted with it.
#' @param start which of the input model's parameter sets each refit starts
#'   from. \code{"initial"} (default) uses \code{initial_params}, the values the
#'   original fit itself started from, so the estimator has to travel the same
#'   distance to the optimum on simulated data that it did on the real data.
#'   \code{"estimated"} starts from \code{estimated_params} instead: much faster
#'   and far more likely to converge, but the fixed effects -- and, with
#'   \code{random_rec = TRUE}, the inner Laplace problem too -- begin at the
#'   generating values, so on a multimodal or weakly identified surface the
#'   optimizer never leaves the basin containing them and recovery is close to
#'   guaranteed by construction. Read it as optimistic about recovery, not
#'   merely less powerful. (Nor is it a complete warm start: \code{fit_mod()}
#'   resets \code{log_Ftarget}, \code{proj_F_prop}, and the stock-recruit
#'   \eqn{\alpha}/\eqn{\beta} from the model's own specification under either
#'   setting.) Non-identifiability shows up in the curvature and so is visible
#'   either way -- via \code{$convergence}'s Hessian conditioning and
#'   estimability checks -- it is \emph{reachability} that a warm start stops
#'   testing.
#' @param debug return every simulation rather than the converged ones. The
#'   dropped runs are the interesting ones when a self-test comes back short, and
#'   each carries its own \code{$convergence} diagnostics. See \strong{Value}.
#' @param timeout elapsed-second limit per simulation, \code{Inf} (default) for
#'   none. The optimizer runs with no iteration cap, so a replicate that wanders
#'   somewhere pathological can stall the whole run -- a hang that no convergence
#'   check can catch, because the fit never returns. One that exceeds the limit
#'   is stopped, counted as non-converged and reported separately. Approximate:
#'   the limit is checked when control returns to R, so it fires between the
#'   optimizer's function evaluations rather than inside one.
#'
#' @return A list of Rceattle models named \code{Sim_1}, \code{Sim_2}, ....
#'   By default only the converged simulations, renumbered contiguously; a
#'   message reports how many were dropped.
#'
#'   With \code{debug = TRUE}, every simulation, with \code{Sim_i} being
#'   simulation \code{i} (so it pairs with the seed \code{seed + i}), and a
#'   logical vector of the convergence verdicts in \code{attr(, "converged")}.
#'   Inspect a failure with \code{sims[[j]]$convergence}. A simulation whose
#'   refit errored outright is returned as the condition object rather than a
#'   model, so it cannot abort the run.
#'
#' @section Interpreting the spread:
#' \code{\link{sim_mod}} redraws the observations only -- indices, catch,
#' compositions and CAAL. It does not redraw recruitment, so with
#' \code{random_rec = TRUE} every replicate shares the operating model's single
#' recruitment realization, and that realization is its shrunk empirical-Bayes
#' modes rather than a draw from N(0, sigmaR). Two consequences: the spread
#' across replicates carries observation error only and is a lower bound on
#' estimation uncertainty in SSB and recruitment (do not read it against the
#' model's own uncertainty bands, which include process error); and sigmaR is
#' re-estimated from deviations that were shrunk toward zero the same way in
#' every replicate, a downward bias that averaging over simulations does not
#' remove.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' sims <- self_test(ss_run, nsim = 10)
#' }
#' @export
self_test <- function(Rceattle = NULL, nsim = 50, simulate = TRUE, seed = 123, cores = NULL, getsd = NULL, phase = NULL, start = c("initial", "estimated"), debug = FALSE, timeout = Inf) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  start <- match.arg(start)
  if (is.null(Rceattle[[paste0(start, "_params")]])) {
    stop("`start = \"", start, "\"` needs the model's ", start,
         "_params, which this fit does not carry.", call. = FALSE)
  }

  # rho/self-test read point estimates, so getsd = FALSE is faster and neutral;
  # default inherits the input model's setting (matches retrospective/jitter).
  if (is.null(getsd)) getsd <- !is.null(Rceattle$sdrep)

  # Phasing likewise inherits the input model's setting. Under the default
  # `start = "initial"` a refit covers the same ground the original fit did, so
  # if that fit needed phasing, so does every refit.
  #
  # Read the value fit_mod() recorded, so a custom phase list carries over as
  # the list rather than collapsing to TRUE and being rebuilt from set_phases()
  # defaults -- that would phase the refit on a different schedule than the fit
  # it is testing. Fall back to whether `phase_params` was attached (fit_mod()
  # does that only when it phased) for models predating `run_config`.
  if (is.null(phase)) {
    phase <- Rceattle$run_config$fit_control$phase
    if (is.null(phase)) phase <- !is.null(Rceattle$phase_params)
  }

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- nsim > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-simulation closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_sim <- function(i) {

    set.seed(seed + i) # unique seed per sim for reproducibility under parallel

    # * Simulate new data
    sim_data <- Rceattle::sim_mod(Rceattle, simulate = simulate)
    data_list <- sim_data

    # * Adjust initial values ----
    inits <- switch(start,
                    initial   = Rceattle$initial_params,
                    estimated = Rceattle$estimated_params)


    # * Refit ----
    # Bounded and error-trapped: a replicate that errors or hangs would otherwise
    # abort self_test() -- and, under a cluster, every other replicate with it --
    # which is the opposite of what `debug` is for.
    newmod <- .refit_with_timeout(
      suppressMessages(
        suppressWarnings(
          # Refit the simulated data set from `start`, under the source model's
          # phasing (see `phase` above).
          .refit_like(
            data_list        = data_list,
            inits            = inits,
            estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
            phase            = phase,
            getsd            = getsd,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
            suit_endyr       = pmin(data_list$suit_endyr, data_list$endyr))
        )
      ),
      timeout = timeout)

    # Report the verdict beside the model rather than dropping it here: the
    # dispatcher filters below, so `debug = TRUE` can hand back the runs that
    # failed with their $convergence diagnostics intact.
    if (inherits(newmod, "condition")) {
      return(list(model = newmod, converged = FALSE, error = conditionMessage(newmod)))
    }
    list(model = newmod, converged = .refit_converged(newmod))
  } # End run_one_sim closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch sims (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    mod_list <- .parallel_lapply(1:nsim, run_one_sim, min(cores, nsim), environment())
  } else {
    mod_list <- lapply(1:nsim, run_one_sim)
  }

  # Split the verdict from the models ----
  # Sim_i is simulation i here, before any filtering, so a debug caller can line
  # a failed run up with the seed (`seed + i`) that produced it.
  converged <- vapply(mod_list, function(x) isTRUE(x$converged), logical(1))
  errs      <- vapply(mod_list, function(x) x$error %||% NA_character_, character(1))
  mod_list  <- lapply(mod_list, function(x) x$model)
  # Unconditional, so `names(which(!attr(sims, "converged")))` works at every
  # length -- including the empty case, where a guarded assignment would have
  # left the attribute unnamed.
  names(mod_list) <- names(converged) <- paste0("Sim_", seq_along(mod_list))
  .report_dropped(sum(!converged), length(converged), "simulation")
  .report_errors(errs, "simulation")


  # Return ----
  # debug: every simulation, `converged` naming which are which, so a run that
  # failed can be read through its own $convergence rather than inferred from a
  # short list. Otherwise the converged runs only, renumbered Sim_1..Sim_k --
  # a self-test is read as a distribution over runs, so the gaps carry no
  # meaning, and plot_biomass(model_names = names(sims)) stays contiguous.
  if (isTRUE(debug)) {
    attr(mod_list, "converged") <- converged
    return(mod_list)
  }

  mod_list <- mod_list[converged]
  if (length(mod_list) > 0) {
    names(mod_list) <- paste0("Sim_", seq_along(mod_list))
  }
  return(mod_list)
}



#' Likelihood profile across one or more parameter cells
#'
#' @description Re-fits an Rceattle model while holding selected cells of a
#'   parameter fixed at user-specified values. Supports profiling a single
#'   cell (e.g. \code{R_log_sd[species = 1]}) and arbitrary N-dimensional
#'   cross-profiles over multiple cells -- e.g. \code{log_M1[1, 1, 1]} and
#'   \code{log_M1[1, 2, 1]} jointly, to profile residual M for males against
#'   females. For each grid point the targeted cells are fixed in the TMB
#'   map and the remaining parameters are re-estimated; the result is a
#'   grid of Rceattle models for downstream NLL surfaces.
#'
#' @param fitted an Rceattle model fit using \code{\link{fit_mod}}
#' @param param Name of the parameter to profile. Two ways to specify it:
#'   \describe{
#'     \item{Raw parameter slot}{any name in
#'       \code{Rceattle$estimated_params}; tested for \code{"R_log_sd"},
#'       \code{"rec_pars"}, and \code{"log_M1"}. \code{slots} must index
#'       into the full array and \code{transform} controls the scale.}
#'     \item{Natural-scale alias}{convenience shortcut for the three
#'       documented parameters. Aliases imply \code{transform = "log"}
#'       (values are taken in natural units and log'd before being
#'       substituted) and, for \code{rec_pars}, fill in the column from
#'       the alias name so \code{slots} only needs the species index:
#'       \itemize{
#'         \item \code{"sigmaR"}, \code{"R_sd"} -> \code{R_log_sd}
#'         \item \code{"M1"} -> \code{log_M1}
#'         \item \code{"R0"} -> \code{rec_pars[, 1]}
#'         \item \code{"alpha"} -> \code{rec_pars[, 2]}
#'         \item \code{"beta"} -> \code{rec_pars[, 3]}
#'       }
#'       If \code{transform} is supplied with an alias it is ignored
#'       (with a warning).}
#'   }
#' @param slots A list whose entries are integer index vectors, one entry
#'   per cell to fix. Each entry's length must equal the number of
#'   dimensions of the resolved parameter -- 1 for vectors
#'   (\code{R_log_sd}), 2 for matrices (\code{rec_pars}), 3 for 3-D arrays
#'   (\code{log_M1}). When using the \code{"R0"}/\code{"alpha"}/\code{"beta"}
#'   aliases, supply only the species index (length 1); the column is
#'   filled in from the alias. E.g. \code{list(c(1, 2, 1))} fixes
#'   \code{log_M1[1, 2, 1]}; \code{list(c(1, 1, 1), c(1, 2, 1))} fixes both
#'   sex cells for a males-vs-females cross-profile of species 1;
#'   \code{list(1, 2)} with \code{param = "sigmaR"} cross-profiles species
#'   1 and 2. If omitted, defaults to a single species-1 slot shaped to
#'   match the resolved parameter (e.g. \code{list(1)} for
#'   \code{R_log_sd}, \code{list(c(1, 1, 1))} for \code{log_M1},
#'   \code{list(1)} for the \code{rec_pars} aliases) and emits a warning;
#'   pass \code{slots} explicitly to silence the warning. Defaulting
#'   requires \code{length(values) == 1L} (otherwise the user must
#'   explicitly say which cell each grid targets).
#' @param values A list of numeric vectors, one per entry of \code{slots}.
#'   The full grid of fits is \code{expand.grid(values)}, so a single slot
#'   gives a 1-D profile and \emph{k} slots give a \emph{k}-D cross-profile.
#' @param transform How to map user values onto the internal parameter scale
#'   before substituting them into \code{inits}. Either \code{"log"}
#'   (default), \code{"identity"}, or a unary function (e.g.
#'   \code{qlogis}). Applied element-wise to every grid value. Aliases
#'   override this with \code{"log"}.
#' @param cores Number of cores to use for parallel fits. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each grid fit runs \code{TMB::sdreport}. The profile
#'   reads only the objective (\code{nll}), so \code{FALSE} is faster with no
#'   effect on the profile. Default \code{NULL} inherits the input model's
#'   setting (\code{TRUE} only if it carries an \code{sdrep}).
#' @param ... Unused; present for consistency with the \code{stats::profile}
#'   generic.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{Rceattle_list}{list of fitted Rceattle models, one per grid
#'       row; entries for non-converged fits are \code{NULL} so positions
#'       stay aligned with \code{grid}.}
#'     \item{grid}{data frame of grid values on the user scale (before
#'       \code{transform}); one column per profiled cell, named
#'       \code{slot_1}, \code{slot_2}, ...}
#'     \item{nll}{numeric vector of joint negative log-likelihoods
#'       (\code{opt$objective}); \code{NA} where the fit did not
#'       converge.}
#'     \item{param}{the profiled parameter name (echoed).}
#'     \item{slots}{the slots list (echoed for downstream plotting).}
#'   }
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#'
#' # 1-D profile of sigmaR for species 1 (alias form -- natural scale)
#' p1 <- profile(ss_run,
#'     param  = "sigmaR",
#'     slots  = list(1),
#'     values = list(seq(0.1, 1.5, by = 0.1)))
#'
#' # Equivalent raw form (log scale -- user does the transform)
#' p1_raw <- profile(ss_run,
#'     param     = "R_log_sd",
#'     slots     = list(1),
#'     values    = list(log(seq(0.1, 1.5, by = 0.1))),
#'     transform = "identity")
#'
#' # 2-D cross-profile of M1 across species 1 and 2 (sex 1, age 1).
#' # BS2017SS is single-sex; with a multi-sex model the same form
#' # (e.g. c(1, 1, 1), c(1, 2, 1)) would cross-profile males vs females.
#' p2 <- profile(ss_run,
#'     param  = "M1",
#'     slots  = list(c(1, 1, 1), c(2, 1, 1)),
#'     values = list(seq(0.1, 0.4, length.out = 3),
#'                   seq(0.1, 0.4, length.out = 3)))
#'
#' # 1-D profile of SRR alpha for species 1 (alias drops the rec_pars column)
#' p3 <- profile(ss_run,
#'     param  = "alpha",
#'     slots  = list(1),
#'     values = list(seq(2, 80, length.out = 20)))
#' }
#' @importFrom stats profile
#' @method profile Rceattle
#' @export
profile.Rceattle <- function(fitted = NULL,
                          param = NULL,
                          slots = NULL,
                          values = NULL,
                          transform = "log",
                          cores = NULL,
                          getsd = NULL,
                          ...) {

  # -- Input validation ----
  if (!inherits(fitted, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }
  # Grid fits inherit the input model's sdreport setting unless overridden;
  # the profile reads only the objective, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(fitted$sdrep)
  if (is.null(param) || !is.character(param) || length(param) != 1L) {
    stop("`param` must be a single character string naming a parameter slot.")
  }
  if (!is.list(values) || length(values) == 0L) {
    stop("`values` must be a non-empty list of numeric grids.")
  }

  # Natural-scale aliases: each maps to a real parameter, implies log()
  # transform, and (for rec_pars aliases) fills in the column index so
  # `slots` only needs the species index.
  alias_table <- list(
    sigmaR = list(param = "R_log_sd", rec_pars_col = NA_integer_),
    R_sd   = list(param = "R_log_sd", rec_pars_col = NA_integer_),
    M1     = list(param = "log_M1",   rec_pars_col = NA_integer_),
    R0     = list(param = "rec_pars", rec_pars_col = 1L),
    alpha  = list(param = "rec_pars", rec_pars_col = 2L),
    beta   = list(param = "rec_pars", rec_pars_col = 3L)
  )

  alias_name   <- NA_character_
  rec_pars_col <- NA_integer_
  if (param %in% names(alias_table)) {
    alias_name <- param
    a <- alias_table[[alias_name]]

    # Aliases force log transform; warn if user passed something else
    if (!identical(transform, "log")) {
      warning(sprintf(
        "`param = \"%s\"` is a natural-scale alias for `%s`; ignoring the supplied `transform` (aliases imply transform = \"log\").",
        alias_name, a$param
      ))
    }
    transform    <- "log"
    rec_pars_col <- a$rec_pars_col
    param        <- a$param   # resolve to real parameter slot
  }

  if (!param %in% names(fitted$estimated_params)) {
    stop("`param` '", param, "' not found in Rceattle$estimated_params.")
  }

  par_array <- fitted$estimated_params[[param]]
  par_ndim  <- if (is.null(dim(par_array))) 1L else length(dim(par_array))

  # Default `slots` to species 1 (a single profile point shaped to match
  # the resolved parameter). For rec_pars aliases the user slot is just
  # the species index; otherwise it's a 1 for every dimension.
  if (is.null(slots)) {
    user_slot_dim <- par_ndim - if (!is.na(rec_pars_col)) 1L else 0L
    default_user_slot <- rep(1L, user_slot_dim)

    if (length(values) != 1L) {
      stop(sprintf(
        "`slots` not supplied but `values` has %d grids -- the species-1 default only covers one cell. Pass `slots` explicitly to profile multiple cells.",
        length(values)
      ))
    }

    pretty_slot <- if (length(default_user_slot) == 1L) {
      as.character(default_user_slot)
    } else {
      paste0("c(", paste(default_user_slot, collapse = ", "), ")")
    }
    warning(sprintf(
      "`slots` not supplied; defaulting to species 1 (slots = list(%s)). Pass `slots` explicitly to silence this warning.",
      pretty_slot
    ))

    slots <- list(default_user_slot)
  }

  if (!is.list(slots) || length(slots) == 0L) {
    stop("`slots` must be a non-empty list of integer index vectors.")
  }
  if (length(values) != length(slots)) {
    stop("`values` must be a list with the same length as `slots`.")
  }

  # Append rec_pars column for rec_pars aliases
  if (!is.na(rec_pars_col)) {
    for (k in seq_along(slots)) {
      if (length(slots[[k]]) != 1L) {
        stop(sprintf(
          "Under alias `\"%s\"`, slots[[%d]] should be a single species index (got length %d). The rec_pars column is filled in from the alias name.",
          alias_name, k, length(slots[[k]])
        ))
      }
      slots[[k]] <- c(as.integer(slots[[k]]), rec_pars_col)
    }
  }

  par_dim <- if (is.null(dim(par_array))) length(par_array) else dim(par_array)

  for (k in seq_along(slots)) {
    if (length(slots[[k]]) != par_ndim) {
      stop(sprintf(
        "slots[[%d]] has length %d but '%s' has %d dimension(s).",
        k, length(slots[[k]]), param, par_ndim
      ))
    }
    if (!all(is.finite(slots[[k]])) || any(slots[[k]] < 1)) {
      stop(sprintf("slots[[%d]] must be a vector of positive integers.", k))
    }
    if (any(slots[[k]] > par_dim)) {
      stop(sprintf(
        "slots[[%d]] = c(%s) is out of bounds for '%s' (dim c(%s)).",
        k,
        paste(slots[[k]], collapse = ", "),
        param,
        paste(par_dim, collapse = ", ")
      ))
    }
  }

  # Build transform fn
  trans_fun <- if (is.function(transform)) {
    transform
  } else if (identical(transform, "log")) {
    log
  } else if (identical(transform, "identity")) {
    function(x) x
  } else {
    stop("`transform` must be \"log\", \"identity\", or a function.")
  }

  # Build grid (user-scale values; transform applied at fit time)
  names(values) <- paste0("slot_", seq_along(values))
  grid <- expand.grid(values, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  ngrid <- nrow(grid)

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (mirrors jitter()/retrospective()). Respect the CRAN core limit.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- ngrid > 1L && cores > 1L

  # Generic [<-]: assign `val` into `arr` at index vector `idx`
  assign_at <- function(arr, idx, val) {
    do.call("[<-", c(list(arr), as.list(idx), list(val)))
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-grid-point closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_point <- function(i) {

    inits     <- fitted$estimated_params
    data_list <- fitted$data_list
    map_obj <- fitted$map

    # Substitute fixed values at each profiled cell
    for (k in seq_along(slots)) {
      inits[[param]] <- assign_at(inits[[param]],
                                  slots[[k]],
                                  trans_fun(grid[i, k]))
    }

    # Force profiled cells to NA
    for (k in seq_along(slots)) {
      map_obj$mapList[[param]] <- assign_at(map_obj$mapList[[param]],
                                            slots[[k]],
                                            NA)
    }
    map_obj$mapFactor <- lapply(map_obj$mapList, factor)

    newmod <-
      suppressMessages(suppressWarnings(
        # Refit with the profiled parameter fixed at its grid value (mapped off
        # in map_obj). estimateMode falls back to 1 -- profile the hindcast fit,
        # not a projection.
        .refit_like(
          data_list        = data_list,
          inits            = inits,
          map              = map_obj,
          estimateMode     = ifelse(data_list$estimateMode < 3, 1, data_list$estimateMode),
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
          suit_endyr       = pmin(data_list$suit_endyr, data_list$endyr))
      ))

    if (.refit_converged(newmod)) {
      return(newmod)
    }
    return(NULL)
  } # End run_one_point closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    mod_list <- .parallel_lapply(seq_len(ngrid), run_one_point, min(cores, ngrid), environment())
  } else {
    mod_list <- lapply(seq_len(ngrid), run_one_point)
  }

  # NLL aligned with grid; NA for non-converged
  nll <- vapply(mod_list,
                function(x) if (is.null(x)) NA_real_ else x$opt$objective,
                numeric(1))

  names(mod_list) <- paste0("Fit_", seq_len(ngrid))

  return(list(
    Rceattle_list = mod_list,
    grid          = grid,
    nll           = nll,
    param         = param,
    slots         = slots
  ))
}

