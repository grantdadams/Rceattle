# =============================================================================
# Configuration "switch" lifecycle for Rceattle
# =============================================================================
# CEATTLE configuration switches (selectivity type, catchability type,
# composition likelihood, HCR, init mode, ...) can be supplied by the user
# either as strings ("Logistic") or as the integer codes
# TMB template ultimately uses (1L). The maps below are the single source
# of truth for that string <-> integer correspondence, and the four functions
# in this file implement the lifecycle every switch flows through inside
# fit_mod():
#
#   user input (string or integer)
#        |
#        v
#   switch_check()      Fill any missing switches with defaults and normalise
#        |              to canonical *strings* (delegates to revert_switches()).
#        v
#   revert_switches()   Integer code -> canonical string. Provides backwards
#        |              compatibility for older integer-coded data files.
#        v
#   validate_switches() Error early if any switch is not a known code/string.
#        |              Called by data_check().
#        v
#   convert_switches()  Canonical string -> integer code for TMB. Called by
#        |              rearrange_data(), immediately before the fit.
#        v
#   TMB template (integer codes)
#
# revert_switches() and convert_switches() are inverse operations (int->string
# vs string->int), not duplicates; their per-column logic differs (e.g. the
# Time_varying_q / environmental-index handling) and must be kept in sync by
# eye. Keeping the maps and all four functions here means the set of valid
# values and the order of operations are visible in one place.
# =============================================================================


#' Rceattle configuration switch maps
#'
#' Internal lookup tables that translate the human-readable switch strings a
#' user supplies (e.g. \code{"Logistic"}, \code{"NPFMC"}) into the integer
#' codes the TMB template consumes (and back again, via
#' \code{revert_switches()}). They are the single source of truth for the set
#' of valid switch values; \code{validate_switches()} errors on anything not
#' listed here. See the lifecycle note at the top of \code{R/0-switches.R} for
#' how the maps are applied.
#'
#' @details
#' \describe{
#'   \item{\code{sel_map}}{Selectivity form (\code{fleet_control$Selectivity}).}
#'   \item{\code{tv_sel_map}}{Time-varying selectivity structure
#'     (\code{fleet_control$Time_varying_sel}).}
#'   \item{\code{q_map}}{Catchability form (\code{fleet_control$Catchability}).}
#'   \item{\code{tv_q_map}}{Time-varying catchability structure
#'     (\code{fleet_control$Time_varying_q}).}
#'   \item{\code{comp_loglike_map}}{Composition likelihood
#'     (\code{fleet_control$Comp_loglike} and \code{CAAL_loglike}).}
#'   \item{\code{fleet_map}}{Fleet type (\code{fleet_control$Fleet_type}).}
#'   \item{\code{initMode_map}}{Initial age-structure mode
#'     (\code{data_list$initMode}; see \code{\link{fit_mod}}).}
#'   \item{\code{suitMode_map}}{Predator-prey suitability mode, per predator
#'     (\code{data_list$suitMode}).}
#'   \item{\code{hcr_map}}{Harvest control rule (see \code{\link{build_hcr}}).}
#' }
#'
#' @name rceattle_switch_maps
#' @keywords internal
NULL

# --- Forward Mappings (String -> Integer) ---
# - Maintains backwards compatibility
sel_map <- c(
  "Fixed" = 0,
  "Logistic" = 1,
  "NonParametric" = 2,
  "DoubleLogistic" = 3,
  "DescendingLogistic" = 4,
  "Hake" = 5,
  "2DAR1" = 6,
  "3DAR1" = 7,
  "DoubleNormal" = 8,
  "NonParametricPM" = 9,  # Ianelli non-parametric, ADMB AMAK ("pm") selectivity penalty
  "LogisticPM" = 11       # ADMB AMAK ("pm") BTS: logistic (multiplicative inflection/slope devs) + free age-1 log-selectivity
)

tv_sel_map <-c(
  "Off" = 0,
  "IID" = 1,
  "AR1" = 2,
  "Block" = 3,
  "RandomWalk" = 4,
  "RandomWalkAscending" = 5
)

q_map <- c(
  "Fixed" = 0,
  "Estimated" = 1,
  "Estimated-with-prior" = 2,
  "Analytical" = 3,
  "PowerEquation" = 4,
  "Environmental" = 5,
  "AR1" = 6,
  "AnalyticalArith" = 7
)

tv_q_map <- c(
  "Off" = 0,
  "IID" = 1,
  "AR1" = 2,
  "Block" = 3,
  "RandomWalk" = 4
)

comp_loglike_map <- c(
  "MultinomialAFSC" = -1,
  "Multinomial" = 0,
  "DirichletMultinomial" = 1
)

# Survey/index biomass observation likelihood family (per fleet). MVN/MVNORM use
# a user-supplied full variance-covariance matrix (e.g. a VAST-derived Sigma) on
# the natural-scale residual vector (obs - q*pred); see `index_cov`, and pair with
# Catchability = "AnalyticalArith" for the AMAK arithmetic-mean q.
#   Lognormal = independent lognormal on log(obs) (the historical default).
#   MVN       = the AMAK/ebswp `DoCovBTS = 1` *bare* quadratic form 0.5 * r' Sigma^-1 r
#               (drops the normalizing constant, so the reported value matches ADMB).
#   MVNORM    = the *full* multivariate-normal negative log-density
#               0.5 * (r' Sigma^-1 r + logdet(Sigma) + n*log(2*pi)), i.e. TMB's
#               density::MVNORM(Sigma)(r). Identical fit to "MVN" (the extra term is
#               a fixed constant), but a proper normalized likelihood; the reported
#               value is "MVN" + 0.5*(logdet(Sigma) + n*log(2*pi)).
index_loglike_map <- c(
  "Lognormal" = 0,
  "MVN" = 1,
  "MVNORM" = 2
)

fleet_map <- c(
  "Fishery" = 1,
  "Survey" = 2,
  "Off" = 0
)

# Initial age-structure mode
# 0 = Free parameters for initial age-structure
# 1 = Equilibrium, no init devs, Finit = 0 (unfished)
# 2 = Equilibrium + init devs, Finit = 0  [default]
# 3 = Non-equilibrium: Finit estimated, init devs included
# 4 = Non-equilibrium: Finit scales R0
initMode_map <- c(
  "FreeParams"           = 0,
  "Equilibrium"          = 1,
  "NonEquilibrium"       = 2,
  "FishedNonEquilibrium"       = 3,
  "FishedNonEquilibriumScaled" = 4
)

# Predator-prey suitability mode (per predator species)
# 1 and 3 are blocked in data_check() until growth-model validation is added
suitMode_map <- c(
  "Empirical"       = 0,
  "GammaLength"     = 1,   # NOT YET AVAILABLE
  "GammaWeight"     = 2,
  "LognormalLength" = 3,   # NOT YET AVAILABLE
  "LognormalWeight" = 4,
  "NormalLength"    = 5,
  "NormalWeight"    = 6
)

# Harvest control rule
hcr_map <- c(
  "NoFishing"    = 0,
  "CMSY"         = 1,
  "ConstantF"    = 2,
  "ConstantFSSB" = 3,
  "ConstantFSPR" = 4,
  "NPFMC"        = 5,
  "PFMC"         = 6,
  "SESSF"        = 7
)


#' Function to check for missing switches for map and parameter functions
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
switch_check <- function(data_list){

  # Helper to set defaults and notify. Pass msg = NULL to fill the default
  # silently (for low-level switches that should not add console noise).
  set_default <- function(val, default, msg = NULL) {
    if(is.null(val)) {
      if(!is.null(msg)) message(msg)
      return(default)
    }
    return(val)
  }

  if(is.null(data_list$srr_fun)){
    data_list$srr_fun <- 0
    data_list$srr_pred_fun <- 0
    data_list$srr_est_mode <- 0
    message("'srr_fun' are not included in data, assuming 0")
  }

  # Model and multi-species switches
  data_list$estDynamics <- set_default(data_list$estDynamics, rep(0, data_list$nspp), "'estDynamics' are not included in data, assuming 0")
  data_list$Diet_comp_weights <- set_default(data_list$Diet_comp_weights, rep(1, data_list$nspp), "'Diet_comp_weights' are not included in data, assuming 1")
  data_list$Diet_loglike <- set_default(data_list$Diet_loglike, rep(0, data_list$nspp), "'Diet_loglike' are not included in data, assuming 'Multinomial'")
  data_list$alpha_wt_len <- set_default(data_list$alpha_wt_len, 1e-6, "'alpha_wt_len' not specified in data, assuming 1e-6")
  data_list$beta_wt_len <- set_default(data_list$beta_wt_len, 3, "'beta_wt_len' not specified in data, assuming 3")
  data_list$M1_model <- set_default(data_list$M1_model, rep(0, data_list$nspp), "'M1_model' is not included in data, assuming 0")
  data_list$msmMode <- set_default(data_list$msmMode, 0, "'msmMode' is not included in data, assuming single-species (0)")
  data_list$M1_re <- set_default(data_list$M1_re, rep(0, data_list$nspp), "'M1_re' is not in data, assuming 0 for all species")
  data_list$initMode <- set_default(data_list$initMode, 2, "'initMode' is not in the data, setting to 2 (default)")
  data_list$comp_offset <- set_default(data_list$comp_offset, 1e-5, NULL) # Composition proportion offset (added to comp/caal before the multinomial. Filled silently; fit_control(comp_offset=) can override at fit time.

  # Bioenergetics scalars: TMB declares them as DATA_VECTOR length-nspp, so
  # they must exist even when not used. In single-species mode (msmMode == 0)
  # they are never read by the consumption code, so silently fill any missing
  # entries with safe sentinels. When msmMode > 0 we leave them untouched and
  # let data_check() report which ones are missing or wrong-length.
  if(data_list$msmMode == 0){
    bioenergetics_defaults <- list(
      Ceq    = rep(1L,  data_list$nspp),
      Cindex = rep(0L,  data_list$nspp),
      Pvalue = rep(1,   data_list$nspp),
      fday   = rep(365, data_list$nspp),
      CA     = rep(0,   data_list$nspp),
      CB     = rep(0,   data_list$nspp),
      Qc     = rep(0,   data_list$nspp),
      Tco    = rep(0,   data_list$nspp),
      Tcm    = rep(0,   data_list$nspp),
      Tcl    = rep(0,   data_list$nspp),
      CK1    = rep(0,   data_list$nspp),
      CK4    = rep(0,   data_list$nspp)
    )
    for(nm in names(bioenergetics_defaults)){
      if(is.null(data_list[[nm]])){
        data_list[[nm]] <- bioenergetics_defaults[[nm]]
      }
    }
  }

  # 1. Fleet Control defaults ----
  data_list$fleet_control$Sel_norm_bin1 <- set_default(data_list$fleet_control$Sel_norm_bin1, NA, "'Sel_norm_bin1' not specified in 'fleet_control', assuming 'NA'")
  data_list$fleet_control$Sel_norm_bin2 <- set_default(data_list$fleet_control$Sel_norm_bin2, NA, "'Sel_norm_bin2' not specified in 'fleet_control', assuming 'NA'")
  # Sel_curve_pen1/2/3 only matter for non-parametric (Ianelli/PM), Hake, or
  # LogisticPM (random-walk weights) selectivity; only warn about missing columns
  # when such a fleet is present, otherwise default silently (avoids noise for
  # logistic-only models).
  .np_hake <- any(data_list$fleet_control$Selectivity %in%
                    c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM", 5, "Hake", 11, "LogisticPM"))
  data_list$fleet_control$Sel_curve_pen1 <- set_default(data_list$fleet_control$Sel_curve_pen1, 0, if(.np_hake) "'Sel_curve_pen1' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Sel_curve_pen2 <- set_default(data_list$fleet_control$Sel_curve_pen2, 0, if(.np_hake) "'Sel_curve_pen2' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Sel_curve_pen3 <- set_default(data_list$fleet_control$Sel_curve_pen3, 0, if(.np_hake) "'Sel_curve_pen3' not specified in 'fleet_control', assuming '0'")
  data_list$fleet_control$Sel_start_year <- set_default(data_list$fleet_control$Sel_start_year, NA)  # per-fleet selectivity penalty start year (NA -> styr); used by LogisticPM
  data_list$fleet_control$Sel_pen_first_age <- set_default(data_list$fleet_control$Sel_pen_first_age, NA)  # first age for the non-parametric shape penalty (NA -> bin_first_selected)
  data_list$fleet_control$Sel_pen_last_age <- set_default(data_list$fleet_control$Sel_pen_last_age, NA)  # last (left) age of the shape-penalty pairs (NA -> nages-2)
  data_list$fleet_control$Sel_shape_mode <- set_default(data_list$fleet_control$Sel_shape_mode, NA)  # shape-penalty mode: "Directional" (default) or "Smooth" (two-sided d^2, RTMB)
  data_list$fleet_control$Sel_cap_age <- set_default(data_list$fleet_control$Sel_cap_age, NA)  # NonParametricRPM age cap (NA -> no cap)
  data_list$fleet_control$Selectivity_dimension <- set_default(data_list$fleet_control$Selectivity_dimension, "Age", "'Selectivity_dimension' not specified in 'fleet_control', assuming 'Age'")
  data_list$fleet_control$Comp_loglike <- set_default(data_list$fleet_control$Comp_loglike, "MultinomialAFSC", "'Comp_loglike' not specified in 'fleet_control', assuming 'MultinomialAFSC'")
  data_list$fleet_control$CAAL_loglike <- set_default(data_list$fleet_control$CAAL_loglike, "Multinomial", "'CAAL_loglike' not specified in 'fleet_control', assuming 'Multinomial'")
  data_list$fleet_control$Index_loglike <- set_default(data_list$fleet_control$Index_loglike, "Lognormal")  # survey index likelihood family; default preserves the historical lognormal fit
  # Also fill per-element NAs: setting Index_loglike for one fleet (e.g. only the
  # covariance survey) leaves the other rows NA, which should default to Lognormal.
  data_list$fleet_control$Index_loglike[is.na(data_list$fleet_control$Index_loglike)] <- "Lognormal"
  data_list$fleet_control$CAAL_weights <- set_default(data_list$fleet_control$CAAL_weights, 1, "'CAAL_weights' not specified in 'fleet_control', assuming 1")
  data_list$fleet_control$Month <- set_default(data_list$fleet_control$Month, 0, "'Month' not specified in 'fleet_control', assuming 0")

  # Format adjustment for NonParametric
  np_idx <- data_list$fleet_control$Selectivity %in% c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM")
  if(any(np_idx & !is.na(data_list$fleet_control$Time_varying_sel) & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")))){
    data_list$fleet_control <- data_list$fleet_control |>
      dplyr::mutate(
        Sel_curve_pen1 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), Time_varying_sel, Sel_curve_pen1),
        Sel_curve_pen2 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), Time_varying_sel_sd_prior, Sel_curve_pen2),
        Time_varying_sel = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), 0, Time_varying_sel),
        Time_varying_sel_sd_prior = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), 0, Time_varying_sel_sd_prior)
      )
    message("Updating format where 'Selectivity == Non-parametric'. Moving non-parametric penalties to 'Sel_curve_pen1' and 'Sel_curve_pen2'.")
  }

  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen1))) stop("'Sel_curve_pen1' is NA in 'fleet_control' for fleet with non-parametric selectivity")
  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen2))) stop("'Sel_curve_pen2' is NA in 'fleet_control' for fleet with non-parametric selectivity")


  # 2. Sel bins ----
  for(flt in 1:nrow(data_list$fleet_control)){

    sp_idx <- data_list$fleet_control$Species[flt]
    age_selex = data_list$fleet_control$Selectivity_dimension[flt] == "Age"
    selex_text <- ifelse(age_selex, "nages", "nlengths")
    max_bin <- ifelse(age_selex,
                      data_list$nages[sp_idx],
                      data_list$nlengths[sp_idx])

    # - Sel normalization bin
    if(any(data_list$fleet_control$Sel_norm_bin1[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$Sel_norm_bin1[flt] <- max_bin
      message(paste0("'Sel_norm_bin1' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }

    # - Upper sel normalization bin
    if(any(data_list$fleet_control$Sel_norm_bin2[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$Sel_norm_bin2[flt] <- max_bin
      message(paste0("'Sel_norm_bin2' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }

    # - N bins
    if(any(data_list$fleet_control$N_sel_bins[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$N_sel_bins[flt] <- max_bin
      message(paste0("'N_sel_bins' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }
  }

  # Update switches to text-based if necessary for older files
  data_list <- revert_switches(data_list)

  # Auto-Off fleets with no observations contributing to the likelihood.
  # catch_data, index_data, comp_data,or caal_data should have
  # at least one row with Year > 0 referencing its
  # Fleet_code.
  active_codes <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(integer(0))
    if (!all(c("Fleet_code", "Year") %in% colnames(df))) return(integer(0))
    df <- df[!is.na(df$Year) & df$Year > 0 & !is.na(df$Fleet_code), , drop = FALSE]
    if (nrow(df) == 0) return(integer(0))
    unique(as.integer(df$Fleet_code))
  }
  active <- unique(c(
    active_codes(data_list$catch_data),
    active_codes(data_list$index_data),
    active_codes(data_list$comp_data),
    active_codes(data_list$caal_data)
  ))
  fc <- data_list$fleet_control
  inactive_idx <- which(!fc$Fleet_code %in% active & fc$Fleet_type != "Off")
  if (length(inactive_idx) > 0) {
    message(paste0(
      "Auto-Off fleet(s) with no observations in the likelihood (all Year < 0 or missing): ",
      paste(fc$Fleet_name[inactive_idx], collapse = ", "),
      ". Setting Fleet_type = 'Off'."
    ))
    data_list$fleet_control$Fleet_type[inactive_idx]   <- "Off"
  }

  # Default Sel_start_year to each fleet's FIRST YEAR OF DATA (not styr).
  # Selectivity deviations before a fleet's first observation have neither data nor
  # a penalty -- every selectivity penalty in the cpp is anchored at start_yr -- so
  # they are unidentified and leave flat directions the optimiser cannot resolve.
  # Deriving the default from the data drops them automatically instead of relying
  # on the user to know the switch exists. Only fleets with time-varying selectivity
  # are affected (build_map is what consumes this); an explicit per-fleet value wins.
  fleet_obs_yrs <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!all(c("Fleet_code", "Year") %in% colnames(df))) return(NULL)
    df <- df[!is.na(df$Year) & df$Year > 0 & !is.na(df$Fleet_code), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    data.frame(Fleet_code = as.integer(df$Fleet_code), Year = as.integer(df$Year))
  }
  obs_yrs <- do.call(rbind, list(
    fleet_obs_yrs(data_list$catch_data),
    fleet_obs_yrs(data_list$index_data),
    fleet_obs_yrs(data_list$comp_data),
    fleet_obs_yrs(data_list$caal_data)
  ))
  if (!is.null(obs_yrs) && nrow(obs_yrs) > 0) {
    first_obs <- tapply(obs_yrs$Year, obs_yrs$Fleet_code, min)
    fc2 <- data_list$fleet_control
    # Fleets sharing a Selectivity_index share ONE selectivity curve, so the start
    # year must be the earliest first-observation across the whole group -- not the
    # fleet's own. E.g. AVO's data starts 2006 but it mirrors the ATS curve, whose
    # data starts 1994; keying off AVO alone would drop the 1994-2005 deviations
    # that the ATS data actually informs.
    sel_grp <- if (!is.null(fc2$Selectivity_index)) as.character(fc2$Selectivity_index) else
      as.character(fc2$Fleet_code)
    fleet_first <- first_obs[as.character(fc2$Fleet_code)]
    grp_first <- tapply(fleet_first, sel_grp, function(v) if (all(is.na(v))) NA_integer_ else min(v, na.rm = TRUE))
    unset_idx <- which(is.na(fc2$Sel_start_year))
    for (i in unset_idx) {
      fy <- grp_first[sel_grp[i]]
      if (!is.na(fy)) {
        data_list$fleet_control$Sel_start_year[i] <- max(as.integer(data_list$styr), as.integer(fy))
      }
    }
  }

  return(data_list)
}


#' Convert integer switches to intuitive text strings. Maintains backwards compatability.
#'
#' @param data_list Rceattle data list
#' @importFrom rlang .data
#' @keywords internal
#' @noRd
revert_switches <- function(data_list) {

  # Helper: convert integer codes to map names. Numeric-looking strings ("0",
  # " 0", "0.000000000") are canonicalized via as.numeric so they match the
  # integer keys in `map`. Strings that are already a map name (or anything
  # else) pass through unchanged.
  .conv <- function(x, map) {
    x_char <- gsub(" ", "", as.character(x))
    x_num <- suppressWarnings(as.numeric(x_char))
    is_num <- !is.na(x_num)
    x_char[is_num] <- as.character(x_num[is_num])
    idx <- match(x_char, as.character(map))
    matched <- !is.na(idx)
    if (any(matched)) {
      x_char[matched] <- names(map)[idx[matched]]
    }
    return(x_char)
  }

  # - Fleet switches
  # Guard: Index_loglike is a newer column; a hand-built fleet_control (or a
  # data_list passed straight to rearrange_data(), bypassing switch_check) may
  # lack it. Default to Lognormal so the conversion below works without it.
  if(is.null(data_list$fleet_control$Index_loglike)){
    data_list$fleet_control$Index_loglike <- "Lognormal"
  }
  data_list$fleet_control <- data_list$fleet_control |>
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Catchability = .conv(.data$Catchability, q_map),
      Comp_loglike = .conv(.data$Comp_loglike, comp_loglike_map),
      CAAL_loglike = .conv(.data$CAAL_loglike, comp_loglike_map),
      Index_loglike = .conv(.data$Index_loglike, index_loglike_map)
    )

  # Time_varying_q doubles as an environmental-index column when Catchability
  # is "AR1" or "Environmental", so only convert the rows that hold a switch.
  non_env_idx <- !data_list$fleet_control$Catchability %in% c("AR1", "Environmental")
  if (any(non_env_idx)) {
    data_list$fleet_control$Time_varying_q[non_env_idx] <-
      .conv(data_list$fleet_control$Time_varying_q[non_env_idx], tv_q_map)
  }

  # - Population dynamics switches
  data_list$initMode <- .conv(data_list$initMode, initMode_map)
  data_list$HCR <- .conv(data_list$HCR, hcr_map)

  return(data_list)
}


#' Validates switches are correct
#'
#' @param data_list Rceattle data list
#'
#' @keywords internal
#' @noRd
validate_switches <- function(data_list = NULL){
  errors <- character(0)

  # Validate fleet_control inputs ----
  invalid_flt_type <- data_list$fleet_control |>
    dplyr::filter(!.data$Fleet_type %in% c(fleet_map, names(fleet_map)))

  invalid_sel <- data_list$fleet_control |>
    dplyr::filter(.data$Fleet_type != "Off" & !.data$Selectivity %in% c(sel_map, names(sel_map)))

  invalid_tv_sel <- data_list$fleet_control |>
    dplyr::filter(.data$Fleet_type != "Off" & !.data$Time_varying_sel %in% c(tv_sel_map, names(tv_sel_map)))

  invalid_q <- data_list$fleet_control |>
    dplyr::filter(!.data$Catchability %in% c(NA, q_map, names(q_map)))

  invalid_tv_q <- data_list$fleet_control |>
    dplyr::filter(.data$Fleet_type != "Off" & .data$Catchability != "Environmental" & !.data$Time_varying_q %in% c(NA, tv_q_map, names(tv_q_map)))

  invalid_comp_ll <- data_list$fleet_control |>
    dplyr::filter(.data$Fleet_type != "Off" & !.data$Comp_loglike %in% c(comp_loglike_map, names(comp_loglike_map)))

  invalid_caal_ll <- data_list$fleet_control |>
    dplyr::filter(.data$Fleet_type != "Off" & !.data$CAAL_loglike %in% c(comp_loglike_map, names(comp_loglike_map)))

  invalid_index_ll <- data_list$fleet_control |>
    dplyr::filter(.data$Fleet_type != "Off" & !.data$Index_loglike %in% c(index_loglike_map, names(index_loglike_map)))

  # Throw clear errors to guide the user
  if(nrow(invalid_flt_type) > 0) {
    errors <- c(errors, paste("Invalid 'Fleet_type' specified for fleets:",
                              paste(invalid_flt_type$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ",paste(range(fleet_map), collapse = ":")," or one of:",
                              paste(names(fleet_map), collapse = ", ")))
  }

  if(nrow(invalid_sel) > 0) {
    errors <- c(errors, paste("Invalid 'Selectivity' specified for fleets:",
                              paste(invalid_sel$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ",paste(range(sel_map), collapse = ":")," or one of:",
                              paste(names(sel_map), collapse = ", ")))
  }

  if(nrow(invalid_tv_sel) > 0) {
    errors <- c(errors, paste("Invalid 'Time_varying_sel' specified for fleets:",
                              paste(invalid_tv_sel$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ",paste(range(tv_sel_map), collapse = ":")," or one of:",
                              paste(names(tv_sel_map), collapse = ", ")))
  }

  if(nrow(invalid_q) > 0) {
    errors <- c(errors, paste("Invalid 'Catchability' specified for fleets:",
                              paste(invalid_q$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(q_map), collapse = ":")," or one of:",
                              paste(names(q_map), collapse = ", ")))
  }

  if(nrow(invalid_tv_q) > 0) {
    errors <- c(errors, paste("Invalid 'Time_varying_q' specified for fleets:",
                              paste(invalid_tv_q$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(tv_q_map), collapse = ":")," or one of:",
                              paste(names(tv_q_map), collapse = ", ")))
  }

  if(nrow(invalid_comp_ll) > 0) {
    errors <- c(errors, paste("Invalid 'Comp_loglike' specified for fleets:",
                              paste(invalid_comp_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(comp_loglike_map), collapse = ":")," or one of:",
                              paste(names(comp_loglike_map), collapse = ", ")))
  }

  if(nrow(invalid_caal_ll) > 0) {
    errors <- c(errors, paste("Invalid 'CAAL_loglike' specified for fleets:",
                              paste(invalid_caal_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(comp_loglike_map), collapse = ":")," or one of:",
                              paste(names(comp_loglike_map), collapse = ", ")))
  }

  if(nrow(invalid_index_ll) > 0) {
    errors <- c(errors, paste("Invalid 'Index_loglike' specified for fleets:",
                              paste(invalid_index_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(index_loglike_map), collapse = ":")," or one of:",
                              paste(names(index_loglike_map), collapse = ", ")))
  }

  # Validate pop dy controls ----
  # * initMode ----
  invalid_initMode <- (!data_list$initMode %in% c(initMode_map, names(initMode_map)))

  if(sum(invalid_initMode) > 0) {
    errors <- c(errors, paste("Invalid 'initMode' specified:",
                              ".\nPlease use an integer code ", paste(range(initMode_map), collapse = ":")," or one of:",
                              paste(names(initMode_map), collapse = ", ")))
  }

  # * HCR ----
  invalid_hcr <- (!data_list$HCR %in% c(hcr_map, names(hcr_map)))

  if(sum(invalid_hcr) > 0) {
    errors <- c(errors, paste("Invalid 'HCR' specified:",
                              ".\nPlease use an integer code ", paste(range(hcr_map), collapse = ":")," or one of:",
                              paste(names(hcr_map), collapse = ", ")))
  }

  if (data_list$msmMode > 0 & !data_list$HCR %in% c("NoFishing", "CMSY", "ConstantF", "ConstantFSSB", "PFMC")) {
    errors <- c(errors, 'Only HCRs "NoFishing" (0), "CMSY" (1), "ConstantF" (2), "ConstantFSSB" (3), or "PFMC" (6) work in multi-species mode currently')
  }

  return(errors)
}


#' Convert intuitive text strings to integer switches for TMB
#'
#' @param data_list Rceattle data list
#'
#' @importFrom rlang .data
#' @keywords internal
#' @noRd
convert_switches <- function(data_list) {

  # Helper: convert a single string value using a map, pass integers through unchanged
  .conv_single <- function(x, map) {
    if (is.character(x) && x %in% names(map)) unname(map[[x]]) else x
  }
  .conv <- Vectorize(.conv_single, vectorize.args = "x", USE.NAMES = FALSE)

  # Fleet controls ----
  # Guard: default the newer Index_loglike column when a caller (e.g. a direct
  # rearrange_data() unit test) supplies a fleet_control that never went through
  # switch_check.
  if(is.null(data_list$fleet_control$Index_loglike)){
    data_list$fleet_control$Index_loglike <- "Lognormal"
  }
  # If vector is a string that exists in our map, replace it with the integer.
  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Catchability = .conv(.data$Catchability, q_map),
      Time_varying_q = .conv(.data$Time_varying_q, tv_q_map),
      Comp_loglike = .conv(.data$Comp_loglike, comp_loglike_map),
      CAAL_loglike = .conv(.data$CAAL_loglike, comp_loglike_map),
      Index_loglike = .conv(.data$Index_loglike, index_loglike_map)
    ) %>%
    # CRITICAL: Force columns back to integers so TMB doesn't crash expecting ints but getting chars
    dplyr::mutate(
      Fleet_type = as.integer(.data$Fleet_type),
      Selectivity = as.integer(.data$Selectivity),
      Time_varying_sel = as.integer(.data$Time_varying_sel),
      Catchability = as.integer(.data$Catchability),
      Time_varying_q = as.integer(.data$Time_varying_q),
      Comp_loglike = as.integer(.data$Comp_loglike),
      CAAL_loglike = as.integer(.data$CAAL_loglike),
      Index_loglike = as.integer(.data$Index_loglike)
    )

  # Pop dy controls ----
  data_list$initMode <- as.integer(.conv(data_list$initMode, initMode_map))
  data_list$HCR <- as.integer(.conv(data_list$HCR, hcr_map))

  return(data_list)
}
