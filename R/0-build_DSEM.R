
#' Function to fit a dynamic structural equation model related to recruitment
#'
#' @param sem Specification for time-series structural equation model structure including lagged or simultaneous effects. See Details section in \code{dsem::make_dsem_ram} for more description. All variables must be included in and named following variables in \code{env_data}. The default is assumes IID recruitment deviates. NOTE: must include \code{recdevs[spp]} for each species \code{1:nspp} (recdevs1 and recdevs2 for a 2 species model)! If no start value is provided, those model terms are not estimated.
#' @param family Distribution for the measurement error of each \code{env_data}
#'   variable used in the \code{sem}. The latent \code{recdevs<sp>} columns are
#'   never observed and are always treated as \code{fixed} (measured exactly), so
#'   \code{family} describes only the env-data variables retained in the
#'   \code{sem}. Supply a single string (recycled across all env variables) or a
#'   character vector with one element per env variable, each one of: \code{fixed}
#'   (default; no measurement error / identity link), \code{normal} (Gaussian,
#'   identity link), \code{gamma} (gamma, fixed CV and log link), \code{bernoulli}
#'   (logit link), \code{poisson} (log link), \code{lognormal} (log link), or
#'   \code{tweedie}. Translated internally to the named list of \code{family}
#'   objects that \code{dsem::dsem} (>= 3.0.0) expects. You may also pass
#'   \code{dsem} \code{family} objects directly (e.g.
#'   \code{dsem::gaussian_fixed_sd("identity", 0.1)} to fix an observation-error
#'   SD): a single object is recycled across env variables, and a vector or list
#'   mixing strings and objects is matched to env variables by name when named,
#'   otherwise by position.
#' @param sigmaR_prior_sd Optional log-scale SD of a lognormal prior placed on
#'   each species' estimated recruitment SD, centered at the assessment value
#'   \code{data_list$sigma_rec_prior}. Default \code{NA} = no prior (recruitment
#'   SD estimated freely). Supplying a finite value (e.g. 0.5--1) regularizes the
#'   recruitment SD away from the \eqn{1/\sigma_R^2} collapse that occurs when
#'   environmental covariates over-explain the recruitment deviations, while
#'   still estimating it. The prior applies only where the SD is estimated (it is
#'   ignored for species whose SD is fixed in the \code{sem}).
#' @param estimate_projection Whether the SEM extends over the projection period.
#'   Default \code{FALSE}, which builds the latent state space over the hindcast
#'   (\code{styr:endyr}) only, so the projection years are absent from the DSEM
#'   likelihood entirely and projection recruitment comes from the usual
#'   projection method (mean recruitment or the stock-recruit relationship).
#'   \code{TRUE} extends the state space to \code{projyr} and projects
#'   recruitment through the SEM itself; pair it with \code{proj_mean_rec = TRUE}
#'   or the projection states are estimated and then never used.
#'
#'   Note that \code{FALSE} means *absent*, not *fixed*. Retaining the projection
#'   rows and merely pinning them would leave them in the GMRF quadratic form,
#'   where lagged paths couple them to the terminal hindcast years and add a
#'   spurious shrinkage pull on the terminal recruitment deviations -- which
#'   propagate to terminal SSB and the ABC.
#'
#' @description
#' The code links dynamic structural equation models to recruitment within Rceattle. The internals of \code{dsem} were copy and pasted into Rceattle. See \code{??dsem} for more description.
#'
#' @export
#'
build_DSEM <- function(sem = NULL,
                       family = "fixed",
                       sigmaR_prior_sd = NA,
                       estimate_projection = FALSE
){
  return(list(sem = sem, family = family, sigmaR_prior_sd = sigmaR_prior_sd,
              estimate_projection = estimate_projection))
}


#' Parse a single SEM path specification
#'
#' Ported from \code{dsem::parse_path} (originally \code{sem::parse.path}). Used
#' by \code{\link{build_dsem_objects}} to collect the variables referenced in a
#' \code{sem}. dsem is distributed under GPL-3; \code{parse_path} originates from
#' package \code{sem} (GPL >= 2, used with permission from John Fox).
#'
#' @param path text to parse
#' @return tagged list with \code{first}, \code{second}, \code{direction}
#' @keywords internal
#' @noRd
parse_path <- function( path ){
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
}


#' Map an Rceattle DSEM family string to a dsem 3.0.0 family object
#'
#' dsem >= 3.0.0 expects glm-style \code{family} objects rather than the
#' character strings Rceattle's API uses. Returns the input unchanged if it is
#' already a \code{family} object.
#'
#' @param fam a single family string (or a \code{family} object)
#' @return a \code{family} object understood by \code{dsem::dsem}
#' @keywords internal
#' @noRd
dsem_family_object <- function( fam ){
  if(inherits(fam, "family")) return(fam)
  switch(as.character(fam),
         "fixed"     = dsem::fixed(),
         "normal"    = stats::gaussian(),
         "gaussian"  = stats::gaussian(),
         "bernoulli" = stats::binomial(),
         "binomial"  = stats::binomial(),
         "poisson"   = stats::poisson(),
         "gamma"     = stats::Gamma(link = "log"),
         "Gamma"     = stats::Gamma(link = "log"),
         "lognormal" = dsem::lognormal(),
         "tweedie"   = dsem::tweedie(),
         stop("Unsupported DSEM 'family': '", fam, "'. Use one of ",
              "fixed, normal, gamma, bernoulli, poisson, lognormal, tweedie.",
              call. = FALSE))
}


#' Registry of the CEATTLE process deviations DSEM can drive
#'
#' DSEM works by making a process's deviations the latent states of the SEM, so
#' several places need to agree on what those states are *called*: the default
#' sem, the columns injected into the tsdata, the regex that gives them the
#' `fixed` family, the SD self-loop lookup, the column-index map handed to the
#' C++, and the pre-fit spec check. This function is the single place that
#' decides, so adding a process is one edit here rather than a hunt through six
#' `paste0()`s and regexes.
#'
#' Each entry carries the sem variable-name `prefix` users write (e.g.
#' `recdevs` in `recdevs1 <-> recdevs1, 0, sigmaR1, 1`), the number of latent
#' series `n` that process contributes, and the resulting `columns` in model
#' order.
#'
#' Only recruitment is wired today. Growth and M become additional entries when
#' they land -- but note the recruitment-SD and prior machinery further down is
#' still recruitment-specific and needs its own generalization at that point;
#' this registry settles naming, not the per-process statistics.
#'
#' @param data_list a data_list; only `nspp` is used.
#' @noRd
.dsem_latent_registry <- function(data_list) {
  list(
    recruitment = list(
      prefix  = .DSEM_LATENT_PREFIXES[["recruitment"]],
      n       = data_list$nspp,
      columns = paste0(.DSEM_LATENT_PREFIXES[["recruitment"]],
                       seq_len(data_list$nspp))
    )
  )
}

# The sem variable-name stems, one per drivable process. Deliberately separate
# from the registry above: the stems are a property of the grammar, not of any
# particular model, so the pattern below can be built without a data_list.
# Folding them together made check_dsem_spec() -- which needs only the stems --
# start requiring data_list$nspp, and error on inputs it used to accept.
.DSEM_LATENT_PREFIXES <- c(recruitment = "recdevs")

# Latent column names across the whole registry, in model order.
#' @noRd
.dsem_latent_columns <- function(data_list) {
  unlist(lapply(.dsem_latent_registry(data_list), `[[`, "columns"),
         use.names = FALSE)
}

# Regex matching any latent column name. These are states, never observations,
# so they take dsem's `fixed` family regardless of what the user asked for.
# Takes no arguments by design -- see the note on .DSEM_LATENT_PREFIXES.
#' @noRd
.dsem_latent_pattern <- function() {
  paste0("^(", paste(unique(.DSEM_LATENT_PREFIXES), collapse = "|"), ")[0-9]+$")
}


#' Function to build the map and parameter objects for DSEM recruitment linkages
#'
#' @param dsem_settings dsem specifications from \code{\link{build_DSEM}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA
#' @param data_list a data_list read in via \code{\link{read_data}} or built directly in R.
#'
#' @export
build_dsem_objects <- function(dsem_settings = NULL, debug = FALSE, data_list = NULL){

  latent_cols <- .dsem_latent_columns(data_list)
  rec_cols    <- .dsem_latent_registry(data_list)$recruitment$columns

  # Build IID sem if NULL
  if(is.null(dsem_settings$sem)){
    sem = c()
    for(sp in seq_len(data_list$nspp)){
      nm <- rec_cols[sp]
      sem <- c(sem, paste0(nm, " <-> ", nm, ", 0, sigmaR", sp,", 1\n")) # No space after
    }
    sem <- paste0(sem, collapse = " ")
    dsem_settings$sem <- sem
  }

  # DSEM time span. With estimate_projection = FALSE the state space is built
  # over the HINDCAST ONLY, so projection years are absent from the GMRF
  # entirely rather than present-but-pinned.
  #
  # This matters: the density is a quadratic form over the stacked [t, j] state
  # space, and lagged RAM paths couple each row to its predecessors. Merely
  # mapping the projection rows off (the previous behavior) fixes those states
  # at their initial value but leaves them in the quadratic form, so their
  # cross-terms with the terminal hindcast years are still not constant in the
  # estimated states -- the terminal recruitment deviations pick up an extra
  # shrinkage pull toward whatever makes the pinned projection states likely.
  # Terminal recdevs feed terminal SSB and hence the ABC, so this is not a
  # cosmetic difference. Truncating the span removes the coupling at the source.
  dsem_endyr <- if (isTRUE(dsem_settings$estimate_projection)) {
    data_list$projyr
  } else {
    data_list$endyr
  }

  # DSEM data
  dsem_data <- data_list$env_data %>%
    # Adding NA in missing years (match assessment begining)
    dplyr::full_join(data.frame(Year=c(data_list$styr:dsem_endyr)), by = dplyr::join_by(Year)) %>%
    dplyr::filter(Year >= data_list$styr & Year <= dsem_endyr) %>% # FIXME: if including init devs, adjust
    dplyr::arrange(Year)

  # - Prepend one all-NA latent column per registry entry, in model order. All-NA
  #   is what marks them as states rather than observations downstream.
  for(nm in rev(latent_cols)){
    dsem_data <- dsem_data %>%
      dplyr::mutate(!!nm := NA_real_) %>%
      dplyr::relocate(dplyr::all_of(nm))
  }

  # - Drop the Year column once all recdev columns are added
  dsem_data <- dsem_data %>%
    dplyr::select(-Year)

  # - Keep only variables referenced in the sem.
  # Mirror make_dsem_ram(): scan the sem to a model table, then parse_path each
  # path to collect the variables on both heads of every arrow.
  sem_model <- scan(
    text         = dsem_settings$sem,
    what         = list(path = "", lag = 1, par = "", start = 1, dump = ""),
    sep          = ",",
    strip.white  = TRUE,
    comment.char = "#",
    fill         = TRUE,
    quiet        = TRUE
  )
  sem_vars <- unique(unlist(lapply(sem_model$path, function(p){
    pp <- parse_path(p)
    c(pp$first, pp$second)
  })))

  # Preserve dsem_data's column order (recdevs<sp> first, then env vars)
  dsem_data <- dsem_data %>%
    dplyr::select(dplyr::any_of(intersect(colnames(dsem_data), sem_vars)))


  # DSEM family. dsem >= 3.0.0 wants a named list of `family` objects keyed by
  # tsdata column. The latent recdevs<sp> columns are never observed, so they
  # are always `fixed`; the user `family` applies only to the env-data variables
  # retained in the sem. `family` may be:
  #   * a single string or `family` object  -> recycled across all env variables
  #   * a vector/list of strings and/or `family` objects, one per env variable,
  #     matched to env columns by name when named, otherwise by position.
  var_names <- colnames(dsem_data)
  is_recdev <- grepl(.dsem_latent_pattern(), var_names)
  env_names <- var_names[!is_recdev]

  fam_in <- dsem_settings$family
  # A bare `family` object is itself a list; wrap it so it counts as one entry.
  if(inherits(fam_in, "family")) fam_in <- list(fam_in)
  # Normalize so character vectors and lists share one code path.
  fam_in  <- as.list(fam_in)
  fam_nms <- names(fam_in)

  if(length(fam_in) == 1L && is.null(fam_nms)){
    # Single unnamed family -> recycle across every env variable.
    env_family <- stats::setNames(rep(fam_in, length(env_names)), env_names)
  } else if(!is.null(fam_nms) && all(nzchar(fam_nms))){
    # Fully named -> match by env-column name (order-independent).
    unknown <- setdiff(fam_nms, env_names)
    if(length(unknown) > 0){
      stop("'family' names are not env-data variables in the sem (",
           paste(env_names, collapse = ", "), "): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    missing_fam <- setdiff(env_names, fam_nms)
    if(length(missing_fam) > 0){
      stop("'family' is missing an entry for env-data variable(s): ",
           paste(missing_fam, collapse = ", "), ".", call. = FALSE)
    }
    env_family <- fam_in[env_names]
  } else {
    # Unnamed vector/list -> matched to env columns by position.
    if(length(fam_in) != length(env_names)){
      stop("Length of 'family' must be 1 or the number of env-data variables in ",
           "the sem (", length(env_names), "), but got ", length(fam_in), ".",
           call. = FALSE)
    }
    env_family <- stats::setNames(fam_in, env_names)
  }

  dsem_family <- stats::setNames(vector("list", length(var_names)), var_names)
  dsem_family[is_recdev] <- list(dsem::fixed())
  dsem_family[env_names] <- lapply(env_family[env_names], dsem_family_object)

  # The compiled TMB model (src/TMB/dsem.hpp) are
  # matched to dsem 3.0.0. Other dsem versions can emit a different
  # `tmb_inputs$data` layout (e.g. missing `obs_idx`/`unobs_idx` under the
  # gmrf_project parameterization), which surfaces downstream as the cryptic TMB
  # error "Error when reading the variable: 'obs_idx'". Fail fast instead.
  if (!requireNamespace("dsem", quietly = TRUE)) {
    stop("The 'dsem' package (version 3.0.0) is required to build DSEM inputs. ",
         "Install the matching version with:\n",
         "  remotes::install_version(\"dsem\", version = \"3.0.0\")",
         call. = FALSE)
  }
  if (utils::packageVersion("dsem") != "3.0.0") {
    stop("Rceattle's DSEM module requires dsem 3.0.0 (the version src/TMB/dsem.hpp ",
         "is matched to), but dsem ", utils::packageVersion("dsem"),
         " is installed.\n",
         "Other versions can produce an incompatible data layout and a downstream ",
         "TMB error such as \"Error when reading the variable: 'obs_idx'\".\n",
         "Install the matching version with:\n",
         "  remotes::install_version(\"dsem\", version = \"3.0.0\")",
         call. = FALSE)
  }

  fit_dsem <- dsem::dsem(sem = dsem_settings$sem,
                         tsdata = stats::ts(dsem_data),
                         family = dsem_family,
                         control = dsem::dsem_control(use_REML = FALSE,
                                                 quiet = TRUE,
                                                 run_model = FALSE))

  # Extract dsem map and parameter objects
  # - Create mapList object
  mapList <- sapply(fit_dsem$tmb_inputs$parameters, function(x) replace(x, values = c(1:length(x))))

  # - Copy dsem map-factor to map-list
  for(i in 1:length(fit_dsem$tmb_inputs$map)){
    parname <- names(fit_dsem$tmb_inputs$map)[i]
    mapList[[parname]] <- replace(fit_dsem$tmb_inputs$parameters[[parname]], values = as.numeric(fit_dsem$tmb_inputs$map[[i]]))
  }

  # Debug mode
  if(debug){
    mapList <- sapply(mapList, function(x) replace(x, values = rep(NA, length=length(x))))
  }

  # Recruitment variance.
  # Locate the beta_z index of each species' recruitment SD: the two-headed
  # self-loop on recdevs[sp] (recdevs[sp] <-> recdevs[sp]).
  sf       <- fit_dsem$sem_full
  sf_start <- suppressWarnings(as.numeric(sf$start))
  sf_par   <- suppressWarnings(as.numeric(sf$parameter))
  sf_dir   <- abs(suppressWarnings(as.numeric(sf$direction)))
  sigma_rec_prior <- rep_len(data_list$sigma_rec_prior, data_list$nspp)

  rec_sd_idx   <- integer(data_list$nspp)   # 1-based beta_z index; 0 if SD is fixed in the sem
  rec_sd_fixed <- numeric(data_list$nspp)   # fixed SD value used when rec_sd_idx == 0
  for(sp in seq_len(data_list$nspp)){
    nm   <- rec_cols[sp]
    rows <- which(sf$first == nm & sf$second == nm & sf_dir == 2)
    pn   <- unique(sf_par[rows]); pn <- pn[!is.na(pn) & pn > 0]
    if(length(pn) > 0){
      rec_sd_idx[sp] <- pn[1]
      # Initialize the recruitment SD at sigma_rec_prior; fix it if !random_rec
      fit_dsem$tmb_inputs$parameters$beta_z[rec_sd_idx[sp]] <- sigma_rec_prior[sp]
      if(!data_list$random_rec){
        mapList$beta_z[rec_sd_idx[sp]] <- NA
      }
    } else {
      # Recruitment SD fixed in the sem (NA parameter name): use its start value
      val <- sf_start[rows]; val <- val[!is.na(val)]
      rec_sd_fixed[sp] <- if(length(val) > 0) val[1] else sigma_rec_prior[sp]
    }
  }
  fit_dsem$tmb_inputs$data$rec_sd_idx   <- as.integer(rec_sd_idx)
  fit_dsem$tmb_inputs$data$rec_sd_fixed <- as.numeric(rec_sd_fixed)

  # Optional lognormal prior on the estimated recruitment SD, centered at the
  # assessment value (sigma_rec_prior). Regularizes R_sd away from the
  # 1/R_sd^2 collapse when env covariates over-explain recruitment. Enabled per
  # species only when the SD is estimated (rec_sd_idx >= 1) and a finite
  # sigmaR_prior_sd is supplied via build_DSEM(); off (NA) by default.
  prior_sd_in <- dsem_settings$sigmaR_prior_sd
  if(is.null(prior_sd_in)) prior_sd_in <- NA
  prior_sd <- rep_len(prior_sd_in, data_list$nspp)
  rec_sd_use_prior <- as.integer((rec_sd_idx >= 1) & is.finite(prior_sd))
  fit_dsem$tmb_inputs$data$rec_sd_prior     <- as.numeric(sigma_rec_prior)
  fit_dsem$tmb_inputs$data$rec_sd_prior_sd  <- as.numeric(ifelse(is.finite(prior_sd), prior_sd, 0))
  fit_dsem$tmb_inputs$data$rec_sd_use_prior <- as.integer(rec_sd_use_prior)

  # x_tj column of each species' recdevs (0-based for the cpp).
  rec_dev_col <- match(rec_cols, colnames(dsem_data)) - 1L
  if(any(is.na(rec_dev_col))){
    stop("Could not locate a latent column (",
         paste(rec_cols, collapse = ", "),
         ") for every species in the DSEM data")
  }
  fit_dsem$tmb_inputs$data$rec_dev_col <- as.integer(rec_dev_col)

  # estimate_projection (DSEM latent projection states) overlaps with
  # proj_mean_rec (SRR projection method). It is contradictory to estimate
  # projection-period DSEM states while projecting recruitment from mean
  # recruitment (proj_mean_rec == 0 / FALSE), because those states are then
  # never used.
  if( isTRUE(dsem_settings$estimate_projection) &&
      !is.null(data_list$proj_mean_rec) &&
      !isTRUE(as.logical(data_list$proj_mean_rec)) ){
    warning("`estimate_projection = TRUE` with `proj_mean_rec = FALSE` is ",
            "inconsistent: projection recruitment uses mean recruitment, so the ",
            "estimated projection DSEM states are never used. Set ",
            "`proj_mean_rec = TRUE` to project via the SEM, or ",
            "`estimate_projection = FALSE`.")
  }

  # Number of DSEM time steps, and whether they cover the projection. The C++
  # call site needs this: with estimate_projection = FALSE the latent state
  # matrix is SHORTER than nyrs, so the recruitment-deviation copy must be
  # rec_dev.row(sp).head(nyrs_dsem), not a whole-row assignment.
  fit_dsem$tmb_inputs$data$nyrs_dsem <- as.integer(length(data_list$styr:dsem_endyr))
  fit_dsem$covers_projection <- isTRUE(dsem_settings$estimate_projection)

  # Return
  fit_dsem$tmb_inputs$map <- sapply(mapList, function(x) factor(x))
  fit_dsem$mapList <- mapList
  fit_dsem$sem <- dsem_settings$sem
  return(fit_dsem = fit_dsem)
}


# -----------------------------------------------------------------------------
# DSEM <-> CEATTLE integration helpers
# -----------------------------------------------------------------------------
# build_dsem_objects() produces TMB inputs that fit_mod() merges into the
# CEATTLE parameter list, map, random-effects vector, and data list at four
# points in its pipeline. The contract of what DSEM contributes lives here so
# those merges stay consistent:
#   * parameters: beta_z, lnsigma_z, mu_j, delta0_j, x_tj
#   * map:        mapList + mapFactor entries for the above
#   * random:     x_tj (when random_rec = TRUE)
#   * data:       options, RAM, RAMstart, familycode_j, linkcode_j, sigmastart_j,
#                 eps_tj, y_tj, obs_idx, unobs_idx

#' Names of the DSEM parameters contributed to the CEATTLE parameter list
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
dsem_param_names <- function(dsem) {
  names(dsem$tmb_inputs$parameters)
}

#' Merge DSEM parameters into a CEATTLE parameter list
#'
#' Adds any DSEM parameters absent from \code{param_list}, pulling defaults from
#' the built DSEM objects. A no-op for parameters already present, so it is safe
#' for both a fresh \code{build_params()} list (none present -> all added) and a
#' prior fit's \code{estimated_params} (all present -> unchanged). Nothing is
#' overwritten or duplicated.
#'
#' @param param_list CEATTLE parameter list (from build_params or inits)
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
merge_dsem_params <- function(param_list, dsem) {
  dsem_par <- dsem$tmb_inputs$parameters
  missing_par <- setdiff(names(dsem_par), names(param_list))
  if (length(missing_par) > 0) {
    param_list <- c(param_list, dsem_par[missing_par])
  }
  param_list
}

#' Merge DSEM map entries into a CEATTLE map object
#'
#' DSEM contributes both \code{mapList} (named numeric/NA vectors) and
#' \code{mapFactor} (the TMB factor map). Merge them together so the two halves
#' never drift apart.
#'
#' @param map CEATTLE map object (from \code{build_map}) with \code{mapList} and
#'   \code{mapFactor}
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
merge_dsem_map <- function(map, dsem) {
  map$mapList   <- c(map$mapList, dsem$mapList)
  map$mapFactor <- c(map$mapFactor, dsem$tmb_inputs$map)
  map
}


# ---- Pre-fit DSEM spec screen ----
# Moved from R/0-convergence.R during the dev-data-workflow (v5.0) merge so it
# lives with the other DSEM builders. Depends on the .conv_record/.conv_overall
# convergence internals -- reconcile with v5.0 R/0-convergence.R in Tier 3.
#' Pre-fit screen of a DSEM recruitment specification
#'
#' Runs the Tier-1 (\code{"spec"}) convergence checks on a built DSEM and the
#' model's \code{env_data}, returning the same \code{"Rceattle_convergence"}
#' object as [convergence_diagnostics()]. \code{fit_mod()} runs this
#' automatically after building the DSEM (results merged into
#' \code{fit$convergence}); call it directly to screen a spec before fitting.
#'
#' Checks (tier \code{"spec"}):
#' \itemize{
#'   \item \code{rec_predictor_observability} -- recruitment predictors observed
#'     in a low fraction of hindcast years (free latents may absorb recruitment
#'     variance -> sigmaR collapse).
#'   \item \code{rec_design_conditioning} -- condition number / max pairwise
#'     correlation of the lag-aligned recruitment design matrix (collinear
#'     predictors).
#'   \item \code{covariate_scale} -- recruitment covariates spanning orders of
#'     magnitude in SD (ill-conditioning; consider standardizing).
#' }
#'
#' @param data_list A cleaned data list (must carry \code{env_data},
#'   \code{styr}, \code{endyr}).
#' @param dsem A built DSEM object (from [build_dsem_objects()]) carrying
#'   \code{sem_full}.
#'
#' @return An object of class \code{"Rceattle_convergence"}.
#' @export
check_dsem_spec <- function(data_list, dsem) {
  empty <- structure(list(status = "OK", checks = list()),
                     class = "Rceattle_convergence")
  sf <- tryCatch(dsem$sem_full, error = function(e) NULL)
  ed <- data_list$env_data
  if (is.null(sf) || is.null(ed) || is.null(ed$Year)) return(empty)

  # Recruitment paths: one-headed (direction 1) arrows into a latent node.
  is_rec <- grepl(.dsem_latent_pattern(), as.character(sf$second)) &
            as.numeric(sf$direction) == 1
  rec <- sf[is_rec, , drop = FALSE]
  if (nrow(rec) == 0) return(empty)

  preds <- intersect(unique(as.character(rec$first)), colnames(ed))
  if (length(preds) == 0) return(empty)

  styr <- data_list$styr; endyr <- data_list$endyr
  yrs  <- styr:endyr
  nyr  <- length(yrs)
  edh  <- ed[ed$Year %in% yrs, , drop = FALSE]
  out  <- list()

  # T1.1 -- latent observability of recruitment predictors
  cov_frac <- vapply(preds, function(p) sum(!is.na(edh[[p]])) / nyr, numeric(1))
  low <- preds[cov_frac < 0.5]
  if (length(low) > 0) {
    sev <- if (any(cov_frac[low] < 0.25)) "WARN" else "NOTE"
    out$rec_predictor_observability <- .conv_record(
      "rec_predictor_observability", "spec", sev,
      sprintf(paste0(
        "Recruitment predictor(s) observed in <50%% of hindcast years: %s. ",
        "Sparse latents can absorb recruitment variance and drive sigmaR to 0."),
        paste(sprintf("%s (%.0f%%)", low, 100 * cov_frac[low]),
              collapse = ", ")),
      data.frame(predictor = preds, coverage = round(cov_frac, 3),
                 row.names = NULL))
  }

  # T1.2 -- lag-aligned recruitment design conditioning
  cols <- lapply(seq_len(nrow(rec)), function(i) {
    p   <- as.character(rec$first[i]); lag <- as.numeric(rec$lag[i])
    if (!p %in% colnames(ed)) return(NULL)
    ed[[p]][match(yrs - lag, ed$Year)]   # predictor at the lag the SEM uses
  })
  names(cols) <- make.unique(paste0(rec$first, "_lag", rec$lag))
  cols <- cols[!vapply(cols, is.null, logical(1))]
  if (length(cols) >= 2) {
    X  <- do.call(cbind, cols)
    Xc <- X[stats::complete.cases(X), , drop = FALSE]
    if (nrow(Xc) > ncol(Xc)) {
      keep <- apply(Xc, 2, function(z) stats::sd(z) > 0)
      Xc   <- Xc[, keep, drop = FALSE]
      if (ncol(Xc) >= 2) {
        R     <- stats::cor(Xc)
        kap   <- tryCatch(kappa(R, exact = TRUE), error = function(e) NA_real_)
        maxr  <- max(abs(R[upper.tri(R)]))
        sev   <- if ((is.finite(kap) && kap > 100) || maxr > 0.9) "WARN"
                 else if (maxr > 0.8) "NOTE" else "OK"
        if (sev != "OK") {
          out$rec_design_conditioning <- .conv_record(
            "rec_design_conditioning", "spec", sev,
            sprintf(paste0("Collinear lag-aligned recruitment design ",
              "(condition number = %.3g, max |r| = %.2f over %d years)."),
              kap, maxr, nrow(Xc)),
            list(condition_number = kap, max_abs_cor = maxr, R = R))
        }
      }
    }
  }

  # T1.5 -- covariate scale heterogeneity
  sds <- vapply(preds, function(p) stats::sd(ed[[p]], na.rm = TRUE), numeric(1))
  sds <- sds[is.finite(sds) & sds > 0]
  if (length(sds) >= 2) {
    ratio <- max(sds) / min(sds)
    if (ratio > 20) {
      out$covariate_scale <- .conv_record(
        "covariate_scale", "spec", "NOTE",
        sprintf(paste0("Recruitment covariates span %.0fx in SD (%s: %.3g; ",
          "%s: %.3g); consider standardizing env_data."),
          ratio, names(which.max(sds)), max(sds),
          names(which.min(sds)), min(sds)),
        data.frame(predictor = names(sds), sd = signif(sds, 3),
                   row.names = NULL))
    }
  }

  structure(list(status = .conv_overall(out), checks = out),
            class = "Rceattle_convergence")
}
