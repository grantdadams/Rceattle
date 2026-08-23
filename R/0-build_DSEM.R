
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
#'   \code{data_list$sigma_rec}. Default \code{NA} = no prior (recruitment
#'   SD estimated freely). Supplying a finite value (e.g. 0.5--1) regularizes the
#'   recruitment SD away from the \eqn{1/\sigma_R^2} collapse that occurs when
#'   environmental covariates over-explain the recruitment deviations, while
#'   still estimating it. The prior applies only where the SD is estimated (it is
#'   ignored for species whose SD is fixed in the \code{sem}).
#' @param estimate_projection Whether the SEM covers the projection period.
#'   Default \code{FALSE} builds the latent states over the hindcast
#'   (\code{styr:endyr}) only, leaving the projection years out of the DSEM
#'   likelihood; projection recruitment then comes from the usual method (mean
#'   recruitment or the stock-recruit relationship). \code{TRUE} extends them to
#'   \code{projyr} and projects recruitment through the SEM. This overrides
#'   \code{build_srr(proj_mean_rec = )}: under mean-recruitment projection the
#'   projected recruitment is \code{avg_R} and the SEM's projection states
#'   would reach only the dynamic B0/BF reference series, so \code{fit_mod()}
#'   turns mean-recruitment projection off and says so.
#'   \code{FALSE} leaves the projection years out, rather than holding them
#'   fixed: fixed states still sit in the GMRF, where lagged paths tie them to
#'   the last hindcast years and pull on the terminal recruitment deviations.
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

# TMB requires every DATA_ and PARAMETER_ macro in the template to be supplied,
# whether or not the model uses it. These are the stand-ins for a fit with no
# DSEM: zero-length everywhere, dsem_on = 0. That keeps fit_mod() on one code
# path -- the switch is the only branch -- and keeps the non-DSEM objective
# untouched, since nothing zero-length can contribute to it.
#' @noRd
.dsem_null_inputs <- function() {
  list(
    data = list(
      dsem_on           = 0L,
      nyrs_dsem         = 0L,
      dsem_options      = integer(4),
      dsem_RAM          = matrix(0L, 0, 6),
      dsem_RAMstart     = numeric(0),
      dsem_familycode_j = integer(0),
      dsem_linkcode_j   = integer(0),
      dsem_sigmastart_j = integer(0),
      dsem_eps_tj       = array(0, dim = c(0, 0)),
      dsem_y_tj         = array(0, dim = c(0, 0)),
      dsem_obs_idx      = integer(0),
      dsem_unobs_idx    = integer(0),
      rec_dev_col       = integer(0),
      rec_sd_idx        = integer(0),
      rec_sd_fixed      = numeric(0),
      rec_sd_prior      = numeric(0),
      rec_sd_prior_sd   = numeric(0),
      rec_sd_use_prior  = integer(0)
    ),
    parameters = list(
      dsem_beta_z    = numeric(0),
      dsem_lnsigma_z = numeric(0),
      dsem_mu_j      = numeric(0),
      dsem_delta0_j  = numeric(0),
      dsem_x_tj      = array(0, dim = c(0, 0))
    )
  )
}

# The DSEM parameter blocks in a shape merge_dsem_params() can consume: named
# for the template, with the built starting values.
#' @noRd
.dsem_par_shape <- function(dsem) {
  list(tmb_inputs = list(parameters = .dsem_par_template(dsem)))
}

# A warm start must match the DSEM this fit is building. inits from another fit
# -- a different year span, a different sem -- would otherwise be reinterpreted
# silently against the new dimensions.
#' @noRd
.dsem_check_shape <- function(param_list, dsem) {
  if (is.null(dsem)) return(param_list)
  want <- .dsem_par_template(dsem)
  for (nm in names(want)) {
    got <- param_list[[nm]]
    # Absent, or present-but-empty (a warm start from a non-DSEM fit carries the
    # zero-length blocks): nothing to conflict with, take the template's.
    if (is.null(got) || length(got) == 0L) {
      param_list[[nm]] <- want[[nm]]
      next
    }
    if (!identical(dim(got), dim(want[[nm]])) || length(got) != length(want[[nm]])) {
      stop("`inits` carries a `", nm, "` of ",
           paste(dim(got) %||% length(got), collapse = "x"),
           " but this DSEM needs ",
           paste(dim(want[[nm]]) %||% length(want[[nm]]), collapse = "x"),
           ". Rebuild the starting values for this model (inits = NULL).",
           call. = FALSE)
    }
  }
  param_list
}

# DSEM map entries, renamed onto the template's dsem_-prefixed parameter names.
# build_dsem_objects() names them as dsem::dsem() does (beta_z, x_tj, ...); the
# CEATTLE template prefixes them so they cannot collide with anything already in
# the model. The rename has to happen in one place or the map silently fails to
# match the parameters.
#' @noRd
.dsem_map_entries <- function(dsem) {
  if (is.null(dsem)) return(list(mapList = list(), mapFactor = list()))
  ren <- function(x) { names(x) <- paste0("dsem_", names(x)); x }
  list(mapList   = ren(dsem$mapList),
       mapFactor = ren(dsem$tmb_inputs$map))
}

# The DSEM parameter blocks for a fit, in the template's naming. Zero-length
# when there is no DSEM.
#' @noRd
.dsem_par_template <- function(dsem) {
  if (is.null(dsem)) .dsem_null_inputs()$parameters else .dsem_tmb_inputs(dsem)$parameters
}

# The same objects for a fit that HAS a DSEM, mapped from build_dsem_objects()
# onto the template's dsem_-prefixed names.
#' @noRd
.dsem_tmb_inputs <- function(dsem) {
  d <- dsem$tmb_inputs$data
  p <- dsem$tmb_inputs$parameters
  list(
    data = list(
      dsem_on           = 1L,
      nyrs_dsem         = as.integer(d$nyrs_dsem),
      dsem_options      = as.integer(d$options),
      dsem_RAM          = matrix(as.integer(as.matrix(d$RAM)), ncol = 6),
      dsem_RAMstart     = as.numeric(d$RAMstart),
      dsem_familycode_j = as.integer(d$familycode_j),
      dsem_linkcode_j   = as.integer(d$linkcode_j),
      dsem_sigmastart_j = as.integer(d$sigmastart_j),
      dsem_eps_tj       = array(as.numeric(d$eps_tj), dim = dim(d$eps_tj)),
      dsem_y_tj         = array(as.numeric(d$y_tj),   dim = dim(d$y_tj)),
      dsem_obs_idx      = if (is.null(d$obs_idx))   integer(0) else as.integer(d$obs_idx),
      dsem_unobs_idx    = if (is.null(d$unobs_idx)) integer(0) else as.integer(d$unobs_idx),
      rec_dev_col       = as.integer(d$rec_dev_col),
      rec_sd_idx        = as.integer(d$rec_sd_idx),
      rec_sd_fixed      = as.numeric(d$rec_sd_fixed),
      rec_sd_prior      = as.numeric(d$rec_sd_prior),
      rec_sd_prior_sd   = as.numeric(d$rec_sd_prior_sd),
      rec_sd_use_prior  = as.integer(d$rec_sd_use_prior)
    ),
    parameters = list(
      dsem_beta_z    = as.numeric(p$beta_z),
      dsem_lnsigma_z = as.numeric(p$lnsigma_z),
      dsem_mu_j      = as.numeric(p$mu_j),
      dsem_delta0_j  = as.numeric(p$delta0_j),
      dsem_x_tj      = array(as.numeric(p$x_tj), dim = dim(p$x_tj))
    )
  )
}

# The DSEM parameter block names, as dsem::dsem() emits them. Used by fit_mod()
# to strip a DSEM warm start when the current fit has no DSEM: these are absent
# from build_params(), so the inits shape guard passes them through as unflagged
# extras and MakeADFun then rejects parameters the template never declares.
# Kept as a constant rather than derived from a built object, because the scrub
# has to happen precisely when no DSEM object exists.
.DSEM_PARAM_NAMES <- paste0("dsem_", c("beta_z", "lnsigma_z", "mu_j",
                                       "delta0_j", "x_tj"))

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
  # v5.0 renamed this control element sigma_rec_prior -> sigma_rec (the schema
  # keeps the old spelling as an alias and upgrades it on entry to the pipeline),
  # so reading the old name yields NULL and rep_len() then fails. Accept either,
  # preferring the canonical one.
  sigma_rec_in <- data_list$sigma_rec
  if (is.null(sigma_rec_in)) sigma_rec_in <- data_list$sigma_rec_prior
  if (is.null(sigma_rec_in)) {
    stop("`data_list` carries neither `sigma_rec` nor `sigma_rec_prior`; the ",
         "DSEM needs a recruitment-SD starting value.", call. = FALSE)
  }
  sigma_rec_prior <- rep_len(sigma_rec_in, data_list$nspp)

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
  # The recruitment-SD prior is applied per species, which is what a sem with a
  # sigmaR<sp> per species gives. A sem that SHARES one SD parameter across
  # species is different: the prior then lands on that one parameter once per
  # species (tightening it by sqrt(n)), and the per-species prior centres and
  # initial values collide. Warn rather than let either happen silently.
  .shared <- rec_sd_idx[rec_sd_idx >= 1]
  if (anyDuplicated(.shared)) {
    warning("Species share a recruitment-SD parameter in the sem. The ",
            "recruitment-SD prior is applied per species, so it is applied ",
            "once per species to that one parameter; and only one species' ",
            "`sigma_rec` sets its starting value. Give each species its own ",
            "sigmaR term if that is not what you intend.", call. = FALSE)
  }

  # x_tj column of each species' recdevs (0-based for the cpp).
  rec_dev_col <- match(rec_cols, colnames(dsem_data)) - 1L
  if(any(is.na(rec_dev_col))){
    stop("Could not locate a latent column (",
         paste(rec_cols, collapse = ", "),
         ") for every species in the DSEM data")
  }
  fit_dsem$tmb_inputs$data$rec_dev_col <- as.integer(rec_dev_col)

  # estimate_projection (DSEM latent projection states) overlaps with
  # proj_mean_rec (SRR projection method): under mean-recruitment projection
  # the template sets R = avg_R past endyr, so the SEM's projection states
  # reach only the dynamic B0/BF reference series (which does set the ABC under
  # DynamicHCR = 1) and the states are otherwise fitted and then ignored.
  # fit_mod() resolves this by letting estimate_projection win, so this only
  # fires for a direct caller that built the objects itself.
  if( isTRUE(dsem_settings$estimate_projection) &&
      isTRUE(as.logical(data_list$proj_mean_rec)) ){
    warning("`estimate_projection = TRUE` with `proj_mean_rec = TRUE`: ",
            "projected recruitment is set to avg_R, so the SEM's projection ",
            "states do not drive it. Set `proj_mean_rec = FALSE` to project ",
            "recruitment through the SEM, or `estimate_projection = FALSE` to ",
            "leave the projection years out of the DSEM. fit_mod() applies the ",
            "first of those for you.", call. = FALSE)
  }

  # Number of DSEM time steps, and whether they cover the projection. The C++
  # call site needs this: with estimate_projection = FALSE the latent state
  # matrix is SHORTER than nyrs, so the recruitment-deviation copy must be
  # rec_dev.row(sp).head(nyrs_dsem), not a whole-row assignment.
  # Take the span from the BUILT data, not the year range: a duplicated year in
  # env_data makes them differ, and the C++ then copies latent states into the
  # wrong assessment years with no error.
  if (nrow(dsem_data) != length(data_list$styr:dsem_endyr)) {
    stop("env_data resolves to ", nrow(dsem_data), " rows over ",
         data_list$styr, ":", dsem_endyr, " (", length(data_list$styr:dsem_endyr),
         " years) -- check for duplicated or missing years.", call. = FALSE)
  }
  fit_dsem$tmb_inputs$data$nyrs_dsem <- as.integer(nrow(dsem_data))
  fit_dsem$covers_projection <- isTRUE(dsem_settings$estimate_projection)

  # Return
  fit_dsem$tmb_inputs$map <- sapply(mapList, function(x) factor(x))
  fit_dsem$mapList <- mapList
  fit_dsem$sem <- dsem_settings$sem
  # Carry the specification the objects were built from. fit_mod() records it on
  # data_list$dsem_settings, which is what .refit_like() forwards and what the
  # diagnostic guards look for -- so a fit handed pre-built objects would
  # otherwise look DSEM-free to both.
  fit_dsem$dsem_settings <- dsem_settings
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
#' object as [convergence_diagnostics()]. Call it directly to screen a
#' specification before fitting; \code{fit_mod()} does **not** run it for you.
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
#'   \item \code{covariate_variance} -- an observed recruitment predictor with
#'     no exogenous variance of its own (no two-headed path, or one fixed at
#'     zero). Such a variable is deterministic: the model computes it rather
#'     than reading \code{env_data}, so the series supplied has no effect.
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
  if (is.null(sf)) return(empty)
  out <- list()
  finish <- function(out) structure(list(status = .conv_overall(out),
                                         checks = out),
                                    class = "Rceattle_convergence")

  # T1.4 -- a variable with no exogenous variance
  # dsem's projecting parameterizations compute such a variable from its
  # incoming paths instead of estimating it, so the env_data supplied for it is
  # never read (scaling that column leaves the objective bit-identical), it
  # drops out of the reported precision that sample_rec() draws from, and the
  # recruitment bias correction cannot condition on it.
  #
  # The rule is dsem's own: a lag-0 two-headed self-path whose parameter AND
  # start are both zero. Not "has no two-headed path" -- make_dsem_ram() adds an
  # estimated V[x] row for any variable that lacks one, so a missing line is the
  # ordinary way to ask for a freely estimated variance.
  #
  # Over every variable in the sem rather than the recruitment predictors alone,
  # because a variable two steps upstream of recruitment does the same thing.
  # Runs before the env_data checks below, which it does not need.
  sf_dir <- abs(suppressWarnings(as.numeric(sf$direction)))
  sf_par <- suppressWarnings(as.numeric(sf$parameter))
  sf_st  <- suppressWarnings(as.numeric(sf$start))
  sf_lag <- suppressWarnings(as.numeric(sf$lag))
  all_vars <- unique(c(as.character(sf$first), as.character(sf$second)))
  novar <- all_vars[vapply(all_vars, function(p) {
    rows <- which(sf$first == p & sf$second == p & sf_dir == 2 &
                  (is.na(sf_lag) | sf_lag == 0))
    if (length(rows) == 0) return(FALSE)
    all((is.na(sf_par[rows]) | sf_par[rows] == 0) &
        (is.na(sf_st[rows])  | sf_st[rows]  == 0))
  }, logical(1))]
  if (length(novar) > 0) {
    out$covariate_variance <- .conv_record(
      "covariate_variance", "spec", "WARN",
      sprintf(paste0(
        "Variable(s) with no exogenous variance: %s. Each is deterministic, so ",
        "the env_data supplied for it is not read and it drops out of the ",
        "precision sample_rec() draws from. Give it a non-zero variance (e.g. ",
        "`%s <-> %s, 0, sd%s, 1`)."),
        paste(novar, collapse = ", "), novar[1], novar[1], novar[1]),
      data.frame(variable = novar, row.names = NULL))
  }

  # Everything below reads env_data.
  if (is.null(ed) || is.null(ed$Year)) return(finish(out))

  # Recruitment paths: one-headed (direction 1) arrows into a latent node.
  is_rec <- grepl(.dsem_latent_pattern(), as.character(sf$second)) &
            as.numeric(sf$direction) == 1
  rec <- sf[is_rec, , drop = FALSE]
  if (nrow(rec) == 0) return(finish(out))

  preds <- intersect(unique(as.character(rec$first)), colnames(ed))
  if (length(preds) == 0) return(finish(out))

  styr <- data_list$styr; endyr <- data_list$endyr
  yrs  <- styr:endyr
  nyr  <- length(yrs)
  edh  <- ed[ed$Year %in% yrs, , drop = FALSE]

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

  finish(out)
}

# Draw the DSEM's UNKNOWN projection states from the precision that scores them.
#
# This is the conditional draw of the unknown projection nodes given everything
# the fit already knows, taken from the GMRF the model was fitted with -- the
# mechanism dsem's own simulate() uses, rmvnorm_prec(xhat + delta, Q).
# Partitioning the field into drawn (P) and known (H):
#
#   x_P | x_H  ~  N( mu_P - Q_PP^-1 Q_PH (x_H - mu_H),  Q_PP^-1 )
#
# so a projected recruitment deviation inherits whatever the sem says: a
# self-path makes it persist, covariate paths make it respond to the
# environment, and a cross-species path couples the stocks. Drawing iid instead
# would project a different process from the one estimated.
#
# WHAT GOES IN P IS THE WHOLE POINT. Not every projection node is unknown. Under
# family = "fixed" a covariate column of x_tj IS the environmental data, pinned
# by the map; where env_data extends past endyr those projection nodes are KNOWN
# FUTURE VALUES -- a climate scenario the projection is meant to respond to.
# Drawing them instead of conditioning on them replaces the scenario with noise
# and returns the no-covariate answer, which is the opposite of what a
# climate-linked assessment is asking for. It also inflates the draw variance,
# against a bias correction the template computes CONDITIONAL on the
# environment (that is what dsem_cond_k is for in ceattle.cpp section 5.5b, and
# it is built from the same per-cell rule). So P is read off the MAP: draw what
# is estimated, condition on what is fixed, wherever it sits.
#
# Requires the latent field to SPAN the projection -- build_DSEM(
# estimate_projection = TRUE). With estimate_projection = FALSE the GMRF stops
# at endyr, the likelihood says nothing about projection years, and any draw
# there is an extrapolation rather than a draw from the model.
#' @noRd
.dsem_draw_projection <- function(fit, sample = TRUE) {
  q <- fit$quantities
  Q <- q$dsem_Q
  if (is.null(Q) || !length(Q)) {
    stop("This model reports no DSEM precision, so its projection states ",
         "cannot be drawn from the fitted process. Either it carries no DSEM, ",
         "or it predates the reporting of `dsem_Q` and must be refitted.",
         call. = FALSE)
  }
  X   <- as.matrix(fit$estimated_params$dsem_x_tj)
  mu  <- as.matrix(q$dsem_xhat_tj) + as.matrix(q$dsem_delta_tj)
  n_t <- nrow(X); n_j <- ncol(X)
  if (nrow(Q) != n_t * n_j) {
    # Under the projecting parameterizations the template reports the precision
    # of the OBSERVED block only, indexed by position within obs_idx rather than
    # by the stacked cell k, so it is short whenever the sem carries a
    # deterministic variable. Name that as the usual cause without claiming it
    # is the only one -- a mis-wired precision reaches here too.
    stop("The reported DSEM precision is ", nrow(Q), "x", ncol(Q),
         " but the latent field is ", n_t, "x", n_j,
         "; they cannot be aligned. Usually this means a variable in the sem ",
         "has no exogenous variance -- its two-headed path is pinned at zero -- ",
         "so the model computes that variable instead of estimating it. Give it ",
         "a non-zero variance (e.g. `BT <-> BT, 0, sdBT, 1`); ",
         "check_dsem_spec() names the variable.", call. = FALSE)
  }

  n_hind <- fit$data_list$endyr - fit$data_list$styr + 1L
  if (n_t <= n_hind) {
    stop("The DSEM's latent states stop at endyr, so the model says nothing ",
         "about the projection years and they cannot be drawn from it. ",
         "Rebuild with build_DSEM(estimate_projection = TRUE).", call. = FALSE)
  }

  # k = j * n_t + t in the template's 0-based indexing, which is column-major
  # over x_tj -- so as.numeric(X)[k] is the node the GMRF scores at k.
  k_of <- function(t, j) (j - 1L) * n_t + t
  proj <- as.integer(unlist(lapply(seq_len(n_j), function(j)
    k_of((n_hind + 1L):n_t, j))))

  # NA in the map means the node is FIXED -- a known value to condition on, not
  # something to draw. An absent map means nothing is fixed.
  mp <- fit$map$mapList$dsem_x_tj
  estimated <- if (is.null(mp)) rep(TRUE, n_t * n_j) else !is.na(as.numeric(mp))
  if (length(estimated) != n_t * n_j) {
    stop("The DSEM's parameter map covers ", length(estimated), " states but ",
         "the latent field has ", n_t * n_j, "; they cannot be aligned.",
         call. = FALSE)
  }

  P <- proj[estimated[proj]]
  if (!length(P)) {
    stop("Every projection state of this DSEM is fixed, so there is nothing ",
         "to draw. The recruitment deviations must be estimated over the ",
         "projection for the process to supply them.", call. = FALSE)
  }
  H <- setdiff(seq_len(n_t * n_j), P)

  Qpp <- Q[P, P, drop = FALSE]
  Qph <- Q[P, H, drop = FALSE]
  dev_h <- as.numeric(X)[H] - as.numeric(mu)[H]
  cond_mu <- as.numeric(mu)[P] - solve(Qpp, Qph %*% dev_h)

  out <- X
  if (isTRUE(sample)) {
    # Draw N(0, Qpp^-1) from the precision, as dsem::simulate() does.
    L <- chol(Qpp)              # Qpp = L'L, so backsolve(L, z) ~ N(0, Qpp^-1)
    out[P] <- as.numeric(cond_mu) + backsolve(L, stats::rnorm(length(P)))
  } else {
    # The conditional MEAN -- the deterministic counterpart of the draw, and
    # what the process expects given everything already known.
    out[P] <- as.numeric(cond_mu)
  }
  out
}
