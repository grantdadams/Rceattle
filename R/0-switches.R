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
#'   \item{\code{sel_norm_scope_map}}{Whether selectivity normalization pools its
#'     reference across sexes (\code{fleet_control$Sel_norm_scope}).}
#'   \item{\code{q_map}}{Catchability form (\code{fleet_control$Catchability}).}
#'   \item{\code{tv_q_map}}{Time-varying catchability structure
#'     (\code{fleet_control$Time_varying_q}).}
#'   \item{\code{comp_loglike_map}}{Composition likelihood
#'     (\code{fleet_control$Comp_distribution} and \code{CAAL_distribution}).}
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

# Whether selectivity normalization pools its reference across sexes. Orthogonal
# to WHERE the reference is taken (Sel_norm_bin: a named bin, or the max), so the
# two combine rather than multiply into more columns.
sel_norm_scope_map <- c(
  "WithinSex"   = 0,  # each sex divided by its own reference; both reach 1
  "AcrossSexes" = 1   # one pooled reference; relative sex selectivity retained
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

# Whether a fleet's selectivity is indexed by age or by length bin. Read in
# rearrange_data(), build_params() and data_check(); a value outside this set
# silently became NA and entered the fit as a missing dimension.
sel_dimension_map <- c(
  "Age"    = 0,
  "Length" = 1
)

# Direction of the (directional-mode) NonParametric shape penalty. The integer
# side is the ADMB sign-of-the-weight convention that .is_incr reads below:
# -1 means increasing, because the penalty weight is negated. Writing 1 here
# would be read as "not increasing" and apply the penalty in the opposite
# direction, so the map must carry the convention the code implements.
sel_shape_dir_map <- c(
  "Decreasing" =  0,
  "Increasing" = -1
)

# Shape-penalty mode: one-sided directional (ADMB/AMAK) or two-sided smooth.
sel_shape_mode_map <- c(
  "Directional" = 0,
  "Smooth"      = 1
)


comp_loglike_map <- c(
  "MultinomialAFSC" = -1,
  "Multinomial" = 0,
  "DirichletMultinomial" = 1
)

# Diet (stomach-content) composition likelihood family, per predator species. The
# diet likelihood supports the multinomial and Dirichlet-multinomial only -- there
# is no AFSC variant -- so this is comp_loglike_map without the -1 case.
diet_loglike_map <- c(
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
#   Normal    = natural-scale normal on the residual (obs - q*pred) with an
#               ABSOLUTE observation sd (the Log_sd column is read as the natural-
#               scale sd, not a log-scale CV), matching the AMAK/ebswp avo_like /
#               cpue_like. No lognormal bias correction; pair with a solved q
#               (Analytical / AnalyticalArith) or an estimated q as needed.
#   TruncatedNormal = the same natural-scale normal, left-truncated at zero:
#               log f(x) = log phi(x; mu, sd) - log Phi(mu/sd). An index cannot
#               be negative and data_check() will not accept one, so this is the
#               only natural-scale family whose simulator and likelihood are the
#               same distribution. Prefer it over "Normal" unless an exact ADMB
#               comparison is needed.
index_distribution_map <- c(
  "Lognormal" = 0,
  "MVN" = 1,
  "MVNORM" = 2,
  "Normal" = 3,
  "TruncatedNormal" = 4
)

#' Read a switch column by canonical name, whichever spelling it arrived in
#'
#' `switch_check()` -> `revert_switches()` upgrades integer codes to names, and
#' every path through `fit_mod()` and `read_data()` runs it. `data_check()` is
#' callable on a hand-built list that has not, so its checks canonicalize first
#' rather than test one spelling and silently pass the other.
#'
#' @param x Column values, names or integer codes.
#' @param map Named integer vector for this switch (e.g. `q_map`).
#' @return Canonical names; `"<blank>"` where the value is missing.
#' @noRd
.canon_switch <- function(x, map) {
  x <- trimws(as.character(x))
  num <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(num), names(map)[match(num, map)], x)
  out[is.na(out)] <- "<blank>"
  out
}


#' Which index rows are fitted on the NATURAL scale?
#'
#' `Index_distribution` splits into two scales, and almost everything downstream
#' of the likelihood assumed the lognormal one. `Lognormal` (0) is a log-scale
#' family whose `Log_sd` is a CV / log-sd; `MVN` (1), `MVNORM` (2), `Normal` (3)
#' and `TruncatedNormal` (4) are natural-scale families whose sd is ABSOLUTE, in
#' the units of the index. Applying a log-scale formula to the second group does
#' not error -- it silently returns nonsense, because `sigma^2 / 2` is then a
#' number the size of the index squared.
#'
#' A new natural-scale family has to be added to the vector below as well as to
#' `index_distribution_map`, or every fleet using it silently reverts to the
#' log-scale treatment this function exists to prevent.
#'
#' @param data_list A `data_list` carrying `fleet_control` and `index_data`.
#' @return Logical, one per `index_data` row; `FALSE` where the fleet is
#'   lognormal or cannot be resolved.
#' @keywords internal
#' @noRd
.index_rows_natural_scale <- function(data_list) {
  idx <- data_list$index_data
  fc  <- data_list$fleet_control
  if (is.null(idx) || !nrow(idx) || is.null(fc)) return(logical(0))
  fam <- fc$Index_distribution
  if (is.null(fam)) return(rep(FALSE, nrow(idx)))
  chr <- trimws(as.character(fam))
  num <- suppressWarnings(as.integer(chr))
  code <- ifelse(!is.na(num), num, as.integer(index_distribution_map[chr]))
  code[is.na(code)] <- 0L
  nat <- code %in% c(1L, 2L, 3L, 4L)   # MVN, MVNORM, Normal, TruncatedNormal
  out <- nat[match(idx$Fleet_code, fc$Fleet_code)]
  out[is.na(out)] <- FALSE
  out
}


#' Fleet codes that carry survey-index observations the model fits
#'
#' An index is a property of the data, not of the fleet type: the template scores
#' an `index_data` row for any fleet that is not `Off`, so a fishery with a CPUE
#' series is fitted like a survey, on that fleet's own selectivity. Selecting on
#' `Fleet_type == "Survey"` instead is what left such a fleet with its
#' catchability frozen and its index absent from `plot_index()`.
#'
#' @param data_list A `data_list` carrying `fleet_control` and `index_data`.
#' @param fitted_only Keep only rows the likelihood uses (positive `Year`, at or
#'   before `endyr`). A prediction-only row is not an observation and must not,
#'   for instance, make catchability estimable.
#' @return Integer vector of `Fleet_code`s, possibly empty.
#' @keywords internal
#' @noRd
.fleets_with_index <- function(data_list, fitted_only = TRUE) {
  idx <- data_list$index_data
  fc  <- data_list$fleet_control
  # nrow() is NULL for anything without dimensions -- a data_list carrying a
  # scalar NA where a table is expected reaches here, and `!NULL` errors.
  if (is.null(idx) || is.null(nrow(idx)) || !nrow(idx) || is.null(fc)) {
    return(integer(0))
  }
  keep <- rep(TRUE, nrow(idx))
  if (fitted_only) {
    keep <- idx$Year > 0
    if (!is.null(data_list$endyr)) keep <- keep & idx$Year <= data_list$endyr
    if (!is.null(idx$Observation)) keep <- keep & !is.na(idx$Observation)
  }
  live <- fc$Fleet_code[!(fc$Fleet_type %in% c("Off", 0))]
  as.integer(sort(intersect(unique(idx$Fleet_code[keep]), live)))
}


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
# 5 = OffsetEquilibrium: F = 0 equilibrium seeded by first-year recruitment
#     (R_init * exp(rec_dev[year 1])), init devs off, no init-dev penalty
#     (Cole Monnahan / AFSC GOA pollock convention). Named for the recruitment
#     offset that seeds it: modes 1 and 5 both start from the initial
#     equilibrium recruitment R_init, but 5 displaces it by the year-1
#     recruitment deviation (the ONLY term separating them in the cpp -- see
#     init_log_scalar in ceattle.cpp). It is an *unfished* (Finit = 0)
#     equilibrium, so "Fished*" would misdescribe it.
initMode_map <- c(
  "FreeParams"                 = 0,
  "Equilibrium"                = 1,
  "NonEquilibrium"             = 2,
  "FishedNonEquilibrium"       = 3,
  "FishedNonEquilibriumScaled" = 4,
  "OffsetEquilibrium"          = 5
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

# Numbers-at-age estimation mode, per species (data_list$estDynamics).
# 0 = estimate the population dynamics; 1 = use fixed input numbers-at-age from
# NByageFixed; 2 = scale the fixed input by one estimated coefficient; 3 = scale
# it by an age-specific estimated coefficient. Named per the Fixed/Estimated
# convention so "0 = not estimated" reads plainly.
estDynamics_map <- c(
  "Estimated"        = 0,
  "Fixed"            = 1,
  "FixedScaled"      = 2,
  "FixedScaledByAge" = 3
)

# Validate a switch value against its map WITHOUT converting it, so a typo is
# reported where it was written rather than several functions later.
#
# model_config() stores what the caller wrote, and fit_mod() rebuilds a
# model_config() from already-resolved values after the fit completes. This must
# therefore accept every legal form: the canonical string, the integer code,
# either as a per-species vector, and NULL/NA for an unset switch. Anything
# stricter would throw on a finished fit and discard it.
.check_switch <- function(x, map, name) {
  if (is.null(x) || !length(x)) return(invisible(x))
  vals <- as.character(x[!is.na(x)])
  if (!length(vals)) return(invisible(x))
  # Treat 2, 2L and "2" alike: normalise anything numeric-looking to its
  # canonical numeric string before matching.
  num <- suppressWarnings(as.numeric(vals))
  vals[!is.na(num)] <- as.character(num[!is.na(num)])
  bad <- unique(vals[!vals %in% c(names(map), as.character(as.numeric(map)))])
  if (length(bad)) {
    stop(sprintf(
      "Invalid `%s` value(s): %s.\nUse one of: %s (or the integer codes %s).",
      name, paste(bad, collapse = ", "),
      paste(names(map), collapse = ", "),
      paste(unname(map), collapse = ", ")), call. = FALSE)
  }
  invisible(x)
}

# Apply a string->integer switch map to a (possibly character) switch value,
# passing integers through unchanged and erroring clearly on an unknown string.
# Used for switches consumed as bare integers before convert_switches() runs
# (e.g. estDynamics is read numerically by build_map()), so the string must be
# resolved early, in switch_check().
.map_switch <- function(x, map, name) {
  if (is.null(x) || !is.character(x)) return(x)
  bad <- setdiff(x[!is.na(x)], names(map))
  if (length(bad) > 0) {
    stop(sprintf("Invalid '%s' value(s): %s. Options: %s (or the integer codes %s).",
                 name, paste(unique(bad), collapse = ", "),
                 paste(names(map), collapse = ", "),
                 paste(unname(map), collapse = ", ")), call. = FALSE)
  }
  unname(map[x])   # map[NA] is NA, so off-fleet NAs pass through
}

# `est_M1` was renamed to `M1_model`. Fold the deprecated name into `M1_model`
# so it is honoured everywhere `M1_model` is (fit_mod's M1Fun reconciliation,
# switch_check's default, combine_data). Only fires when `M1_model` is unset and
# `est_M1` carries a real (non-all-NA) value -- an all-NA `est_M1` (e.g. the
# GOA2018SS placeholder) means "not specified" and must fall through to the
# `M1_model` default, not pin it to NA.
.alias_est_M1 <- function(data_list) {
  if (is.null(data_list$M1_model) && !is.null(data_list$est_M1) &&
      !all(is.na(data_list$est_M1))) {
    data_list$M1_model <- data_list$est_M1
    message("'est_M1' is deprecated; use 'M1_model'.")
  }
  data_list
}

# Observation-SD estimation mode for a survey index or catch series
# (fleet_control$Estimate_index_sd / Estimate_catch_sd). 0 = use the fixed SD
# implied by the data CV (not estimated); 1 = estimate; 2 = analytical
# (Ludwig-Walters 1994). "Fixed" reads plainly for the 0 = not-estimated case.
estimate_sd_map <- c(
  "Fixed"      = 0,
  "Estimated"  = 1,
  "Analytical" = 2
)

# Stock-recruit hyperparameter estimation mode (data_list$srr_est_mode; set via
# build_srr()). 0 = fix alpha/steepness at its prior mean (not estimated);
# 1 = freely estimate; 2 = lognormal prior; 3 = beta prior.
srr_est_mode_map <- c(
  "Fixed"          = 0,
  "Estimated"      = 1,
  "LognormalPrior" = 2,
  "BetaPrior"      = 3
)

# What fit_mod() does with the model (the `estimateMode` argument).
# 0 = fit the hindcast and the HCR projection; 1 = fit the hindcast only;
# 2 = project only, from the supplied inits; 3 = build MakeADFun without
# optimizing (returns the real objective, for pre-fit diagnostics); 4 = optimize
# with every hindcast parameter mapped out (a plumbing smoke test).
estimateMode_map <- c(
  "Estimate"      = 0,
  "Hindcast"      = 1,
  "Projection"    = 2,
  "DebugBuild"    = 3,
  "DebugOptimize" = 4
)

# Predation-mortality mode (the `msmMode` argument). 0 = single-species (no
# predation); 1 = Type II MSVPA (Holsman et al. 2015); 2 = Type III MSVPA. Higher
# integer codes (Kinzey-Punt, Holling forms) are not yet implemented.
msmMode_map <- c(
  "SingleSpecies" = 0,
  "MSVPA"         = 1,
  "TypeIIIMSVPA"  = 2
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

  # Upgrade any deprecated fleet_control column names to their canonical
  # spellings from the single schema-driven migration (`aliases` field), here
  # at the top of switch_check() so the rename lands before build_params()
  # reads the columns and before the non-parametric penalty migration below.
  # Legacy names accepted: Q_prior, Index_sd_prior/Survey_sd_prior,
  # Catch_sd_prior, Time_varying_{q,sel}_sd_prior, Sel_sd_prior, Nselages,
  # Estimate_q, Estimate_survey_sd, Age_first_selected, Age_max_selected(_upper).
  data_list$fleet_control <-
    .rce_upgrade_fleet_control_aliases(data_list$fleet_control)
  # ...and the deprecated control / bioenergetics element names (e.g.
  # `sigma_rec_prior` -> `sigma_rec`, `Diet_loglike` -> `Diet_distribution`).
  data_list <- .rce_upgrade_data_list_aliases(data_list)

  # `Accumulation_age_lower` / `Accumulation_age_upper` are removed. They were
  # only ever range-validated and never applied (no composition-age grouping
  # exists anywhere in the R pipeline or the C++ template), so the columns did
  # nothing. Drop them -- once, with a soft-deprecation message -- if an old
  # workbook / data list still carries either spelling, so a stale column does
  # not silently propagate.
  drop_deprecated_col <- function(fc, old) {
    if (!is.null(fc[[old]])) {
      fc[[old]] <- NULL
      message(sprintf(
        "'%s' is deprecated and ignored; it was never applied to composition data and has been removed.",
        old))
    }
    fc
  }
  data_list$fleet_control <- drop_deprecated_col(
    data_list$fleet_control, "Accumulation_age_lower")
  data_list$fleet_control <- drop_deprecated_col(
    data_list$fleet_control, "Accumulation_age_upper")

  # `growth_re` and `growth_indices` are removed. `growth_re` was documented as
  # the way to put random effects on growth, but nothing consumed it: the
  # deviation array was mapped off in every configuration and the template gave
  # it no density, so setting it changed no fit. Drop them with a message rather
  # than silently, so a data list still carrying `growth_re = 1` says where the
  # feature went instead of appearing to work.
  for (.old in c("growth_re", "growth_indices")) {
    if (!is.null(data_list[[.old]])) {
      data_list[[.old]] <- NULL
      message(sprintf(
        paste0("'%s' is deprecated and ignored; it never had an effect on any ",
               "fit. Specify time-varying growth with build_growth(linkages = ) ",
               "-- see vignette('environmental-linkages-and-priors')."),
        .old))
    }
  }

  if(is.null(data_list$srr_fun)){
    data_list$srr_fun <- 0
    data_list$srr_pred_fun <- 0
    data_list$srr_est_mode <- 0
    message("'srr_fun' are not included in data, assuming 0")
  }

  # Gate optional-input default messages so they fire only when the model actually
  # uses the input (mirrors the PR's conditional data-requirement reporting), rather
  # than nagging about inputs the configuration never consumes. Read here, before the
  # defaults below overwrite the raw fleet_control values:
  #   - growth_estimated: weight-length (alpha_wt_len / beta_wt_len) and the
  #     selectivity dimension are only consumed when growth is estimated (growth_model > 0);
  #   - has_caal: the CAAL distribution / weights defaults only matter with CAAL data;
  #   - sel_norm_upper: the selectivity-normalization upper bin only matters when a fleet
  #     normalizes at a specific bin (Sel_norm_bin >= 0), not max-normalized / off.
  .dflt_when <- list(
    growth_estimated = isTRUE(any(data_list$growth_model > 0)),
    has_caal         = isTRUE(nrow(data_list$caal_data) > 0),
    sel_norm_upper   = isTRUE(any(data_list$fleet_control$Sel_norm_bin >= 0, na.rm = TRUE)),
    #   - sel_norm_scope_flip: the one configuration the "AcrossSexes" default
    #     changes -- a two-sex fleet at a named bin, which used to imply a per-sex
    #     reference. Max-normalized and one-sex fleets are unaffected. Restricted
    #     to fleets the normalization block actually runs on, mirroring the gate
    #     in selectivity.hpp: an "Off" fleet is skipped, Hake normalizes in its own
    #     year/sex block, and LogisticPM reuses Sel_norm_bin1/2 as a penalty
    #     age-range rather than a normalization reference. Without this the
    #     message cries wolf on any AMAK-style model that sets a penalty range.
    sel_norm_scope_flip = isTRUE(any(
      data_list$nsex[data_list$fleet_control$Species] == 2 &
        data_list$fleet_control$Sel_norm_bin >= 0 &
        !data_list$fleet_control$Fleet_type %in% c(0, "Off") &
        !data_list$fleet_control$Selectivity %in%
          c(5, "Hake", 11, "LogisticPM"), na.rm = TRUE)),
    multispecies     = isTRUE(any(data_list$msmMode > 0))
  )

  # Model and multi-species switches
  data_list$estDynamics <- set_default(data_list$estDynamics, rep(0, data_list$nspp), "'estDynamics' are not included in data, assuming 0")
  # Resolve readable strings ("Fixed"/"Estimated"/...) to integer codes now --
  # build_map()/build_params() read estDynamics numerically, before
  # convert_switches() runs.
  data_list$estDynamics <- .map_switch(data_list$estDynamics, estDynamics_map, "estDynamics")
  data_list$suitMode <- .map_switch(data_list$suitMode, suitMode_map, "suitMode")
  if (!is.null(data_list$fleet_control$Estimate_index_sd)) {
    data_list$fleet_control$Estimate_index_sd <-
      .map_switch(data_list$fleet_control$Estimate_index_sd, estimate_sd_map, "Estimate_index_sd")
  }
  if (!is.null(data_list$fleet_control$Estimate_catch_sd)) {
    data_list$fleet_control$Estimate_catch_sd <-
      .map_switch(data_list$fleet_control$Estimate_catch_sd, estimate_sd_map, "Estimate_catch_sd")
  }
  # Diet inputs are consumed only under predation; fill silently otherwise.
  data_list$Diet_comp_weights <- set_default(data_list$Diet_comp_weights, rep(1, data_list$nspp),
    if (.dflt_when$multispecies) "'Diet_comp_weights' are not included in data, assuming 1")
  data_list$Diet_distribution <- set_default(data_list$Diet_distribution, rep(0, data_list$nspp),
    if (.dflt_when$multispecies) "'Diet_distribution' are not included in data, assuming 'Multinomial'")
  # Resolve readable strings ("Multinomial" / "DirichletMultinomial") to the integer
  # codes the C++ diet likelihood reads; integer input passes through unchanged.
  data_list$Diet_distribution <- .map_switch(data_list$Diet_distribution, diet_loglike_map, "Diet_distribution")
  # alpha_wt_len / beta_wt_len (the length -> weight conversion) are only used when
  # growth is estimated, so only announce the default then; otherwise fill silently.
  data_list$alpha_wt_len <- set_default(data_list$alpha_wt_len, 1e-6,
    if (.dflt_when$growth_estimated) "'alpha_wt_len' not specified in data, assuming 1e-6")
  data_list$beta_wt_len <- set_default(data_list$beta_wt_len, 3,
    if (.dflt_when$growth_estimated) "'beta_wt_len' not specified in data, assuming 3")
  data_list <- .alias_est_M1(data_list)
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
  # Per-column defaults + messages come from the canonical schema
  # (R/0-column_schema.R) via .rce_apply_default(), so they live in one place
  # instead of being hand-copied here. The order of the calls (and thus the
  # message order) is unchanged.
  .sch <- .rce_column_schema()
  data_list$fleet_control$Sel_norm_bin <- .rce_apply_default(data_list$fleet_control$Sel_norm_bin, "Sel_norm_bin", .sch)
  data_list$fleet_control$Sel_norm_bin_upper <- .rce_apply_default(data_list$fleet_control$Sel_norm_bin_upper, "Sel_norm_bin_upper", .sch, conditions = .dflt_when)
  # Sel_norm_scope defaults per-cell, not just per-column: .rce_apply_default()
  # returns early once the column exists, but a blank cell is the same silent
  # behaviour flip as a missing column -- it would reach the TMB integer vector
  # as NA, which the C++ reads as WithinSex. Announce whichever case applies, so
  # the "your two-sex named-bin fleet now pools" notice cannot be skipped by
  # supplying the column and leaving one row empty.
  data_list$fleet_control$Sel_norm_scope <- .rce_apply_default(
    data_list$fleet_control$Sel_norm_scope, "Sel_norm_scope", .sch,
    conditions = .dflt_when)
  if (anyNA(data_list$fleet_control$Sel_norm_scope)) {
    data_list$fleet_control$Sel_norm_scope[is.na(data_list$fleet_control$Sel_norm_scope)] <-
      .sch[["Sel_norm_scope"]]$default
    if (isTRUE(.dflt_when$sel_norm_scope_flip))
      message(.sch[["Sel_norm_scope"]]$default_msg)
  }
  # Sel_curve_pen1/2/3 only matter for non-parametric (Ianelli/PM), Hake, or
  # LogisticPM (random-walk weights) selectivity; only warn about missing columns
  # when such a fleet is present, otherwise default silently (avoids noise for
  # logistic-only models).
  .np_hake <- any(data_list$fleet_control$Selectivity %in%
                    c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM", 5, "Hake", 11, "LogisticPM"))
  # Intuitive alternative to the cryptic selectivity penalty WEIGHTS: express each
  # as a standard deviation. Every such penalty is a Gaussian SSQ
  # `weight * x^2 = x^2 / (2*sd^2)`, so `weight = 1/(2*sd^2)`. A fleet may supply
  # the SD column instead of the raw `Sel_curve_pen` weight; convert it here (only
  # where the legacy weight is not already supplied, so legacy models are untouched
  # / bit-identical). Each SD column applies to the selectivity FORMS that use its
  # `Sel_curve_pen` slot as a penalty weight -- on any other form the slot is a
  # logit-scale correlation (2D/3D-AR1) or unused, so setting the SD there is
  # rejected rather than silently converted. A non-positive/non-finite SD (sd = 0
  # would give an Inf penalty; a negative SD squares to a spurious weight) is also
  # rejected.
  .fc  <- data_list$fleet_control
  .col <- function(nm) if (is.null(.fc[[nm]])) rep(NA_real_, nrow(.fc)) else suppressWarnings(as.numeric(.fc[[nm]]))
  .np  <- c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM")
  .lpm <- c(11, "LogisticPM")
  # SD column -> (target Sel_curve_pen slot, forms that use it as a weight). Both
  # NonParametric (2/9) and LogisticPM (11) use pen1 (shape) and pen3 (dev-mag) as
  # Gaussian weights; only NonParametric uses pen2 (curvature) -- LogisticPM leaves
  # it unused.
  .sd_specs <- list(
    Sel_shape_sd     = list(pen = "Sel_curve_pen1", forms = c(.np, .lpm)),
    Sel_curvature_sd = list(pen = "Sel_curve_pen2", forms = .np),
    Sel_devmag_sd    = list(pen = "Sel_curve_pen3", forms = c(.np, .lpm))
  )
  .pen <- list(Sel_curve_pen1 = .col("Sel_curve_pen1"),
               Sel_curve_pen2 = .col("Sel_curve_pen2"),
               Sel_curve_pen3 = .col("Sel_curve_pen3"))
  .dir <- if (is.null(.fc$Sel_shape_dir)) rep("Decreasing", nrow(.fc)) else as.character(.fc$Sel_shape_dir)
  .is_incr <- .dir %in% c("Increasing", "increasing", "-1")
  for (nm in names(.sd_specs)) {
    v <- .col(nm)
    if (all(is.na(v))) next
    spec <- .sd_specs[[nm]]
    bad_form <- which(!is.na(v) & !(.fc$Selectivity %in% spec$forms))
    if (length(bad_form) > 0) {
      stop(sprintf(paste0(
        "'%s' (a penalty-SD column) is set on fleet(s) %s whose Selectivity does ",
        "not use that penalty slot as a weight; it applies to %s only (use ",
        "Sel_curve_pen for other forms)."), nm, paste(bad_form, collapse = ", "),
        if (identical(spec$forms, .np)) "NonParametric/NonParametricPM"
        else "NonParametric/NonParametricPM/LogisticPM"), call. = FALSE)
    }
    bad_val <- which(!is.na(v) & !(is.finite(v) & v > 0))
    if (length(bad_val) > 0) {
      stop(sprintf(paste0(
        "'%s' must be a positive standard deviation; fleet(s) %s are non-positive ",
        "or non-finite."), nm, paste(bad_val, collapse = ", ")), call. = FALSE)
    }
    w <- 1 / (2 * v^2)
    # The shape penalty is DIRECTIONAL for the non-parametric forms only (the sign
    # of pen1 sets penalize-decreasing vs -increasing); LogisticPM's shape term is
    # two-sided (d^2), so its weight must stay positive.
    if (nm == "Sel_shape_sd") {
      np_here <- .fc$Selectivity %in% .np
      bad_dir <- which(!is.na(v) & !np_here & .is_incr)
      if (length(bad_dir) > 0) {
        stop(sprintf(paste0(
          "Sel_shape_dir = 'Increasing' on fleet(s) %s, but their (LogisticPM) ",
          "shape penalty is two-sided -- the weight must be positive."),
          paste(bad_dir, collapse = ", ")), call. = FALSE)
      }
      w <- ifelse(np_here & .is_incr, -w, w)
    }
    u <- !is.na(v) & is.na(.pen[[spec$pen]])
    .pen[[spec$pen]][u] <- w[u]
  }
  data_list$fleet_control$Sel_curve_pen1 <- .pen$Sel_curve_pen1
  data_list$fleet_control$Sel_curve_pen2 <- .pen$Sel_curve_pen2
  data_list$fleet_control$Sel_curve_pen3 <- .pen$Sel_curve_pen3
  data_list$fleet_control$Sel_curve_pen1 <- .rce_apply_default(data_list$fleet_control$Sel_curve_pen1, "Sel_curve_pen1", .sch, .np_hake)
  data_list$fleet_control$Sel_curve_pen2 <- .rce_apply_default(data_list$fleet_control$Sel_curve_pen2, "Sel_curve_pen2", .sch, .np_hake)
  data_list$fleet_control$Sel_curve_pen3 <- .rce_apply_default(data_list$fleet_control$Sel_curve_pen3, "Sel_curve_pen3", .sch, .np_hake)
  data_list$fleet_control$Sel_start_year <- .rce_apply_default(data_list$fleet_control$Sel_start_year, "Sel_start_year", .sch)  # per-fleet selectivity penalty start year (NA -> styr); used by LogisticPM
  # The `Sel_pen_first_age` / `Sel_pen_last_age` / `Sel_cap_age` back-compat
  # renames (these bins were named *_age before they were generalised to work on
  # either dimension) are handled by .rce_upgrade_fleet_control_aliases() at the
  # top of switch_check(), from the schema `aliases` field.
  data_list$fleet_control$Sel_pen_first_bin <- .rce_apply_default(data_list$fleet_control$Sel_pen_first_bin, "Sel_pen_first_bin", .sch)  # first bin (age or length) for the non-parametric shape penalty (NA -> bin_first_selected)
  data_list$fleet_control$Sel_pen_last_bin <- .rce_apply_default(data_list$fleet_control$Sel_pen_last_bin, "Sel_pen_last_bin", .sch)  # last (left) bin of the shape-penalty pairs (NA -> nbins-2)
  data_list$fleet_control$Sel_shape_mode <- .rce_apply_default(data_list$fleet_control$Sel_shape_mode, "Sel_shape_mode", .sch)  # shape-penalty mode: "Directional" (default) or "Smooth" (two-sided d^2, RTMB)
  data_list$fleet_control$Sel_avgsel_pen <- .rce_apply_default(data_list$fleet_control$Sel_avgsel_pen, "Sel_avgsel_pen", .sch)  # weight on the AMAK avgsel base-level penalty (type 9): weight * (log(mean(exp(base coffs))))^2; 0 = off (default), 10 matches AMAK
  data_list$fleet_control$Sel_cap_bin <- .rce_apply_default(data_list$fleet_control$Sel_cap_bin, "Sel_cap_bin", .sch)  # NonParametricRPM bin cap (NA -> no cap)
  data_list$fleet_control$Selectivity_dimension <- .rce_apply_default(data_list$fleet_control$Selectivity_dimension, "Selectivity_dimension", .sch, conditions = .dflt_when)
  data_list$fleet_control$Comp_distribution <- .rce_apply_default(data_list$fleet_control$Comp_distribution, "Comp_distribution", .sch)
  data_list$fleet_control$CAAL_distribution <- .rce_apply_default(data_list$fleet_control$CAAL_distribution, "CAAL_distribution", .sch, conditions = .dflt_when)
  data_list$fleet_control$Index_distribution <- .rce_apply_default(data_list$fleet_control$Index_distribution, "Index_distribution", .sch)  # survey index likelihood family; default preserves the historical lognormal fit
  # Also fill per-element NAs: setting Index_distribution for one fleet (e.g. only
  # the covariance survey) leaves the other rows NA, which should take the same
  # default. Read it from the schema rather than repeating the literal -- the
  # value belongs in one place, and this line sat three lines below the
  # schema-driven call that already knows it.
  data_list$fleet_control$Index_distribution[is.na(data_list$fleet_control$Index_distribution)] <-
    .sch[["Index_distribution"]]$default
  # Same for Selectivity_dimension: .rce_apply_default() fills a MISSING column,
  # never a blank cell, and a partially-assigned column
  # (fleet_control$Selectivity_dimension[i] <- "Length") is a live idiom in the
  # assessment scripts. Without this the blank rows would now be a hard error.
  #
  # Announced under the same gate as the missing-column default (growth
  # estimated, so a length-based selectivity is on the table), and naming the
  # fleets, because only some rows are being filled: a blank left on a
  # growth-estimated model is where the author meant "Length", and an age-based
  # curve on a length-based fleet changes the fit without saying so.
  .sel_dim_blank <- which(is.na(data_list$fleet_control$Selectivity_dimension))
  if (length(.sel_dim_blank) > 0) {
    if (isTRUE(.dflt_when$growth_estimated)) {
      .lbl <- data_list$fleet_control$Fleet_name
      .lbl <- if (is.null(.lbl)) .sel_dim_blank else .lbl[.sel_dim_blank]
      message(sprintf(
        "'Selectivity_dimension' is blank for fleet(s) %s in 'fleet_control', assuming '%s'",
        paste(.lbl, collapse = ", "), .sch[["Selectivity_dimension"]]$default))
    }
    data_list$fleet_control$Selectivity_dimension[.sel_dim_blank] <-
      .sch[["Selectivity_dimension"]]$default
  }
  data_list$fleet_control$CAAL_weights <- .rce_apply_default(data_list$fleet_control$CAAL_weights, "CAAL_weights", .sch, conditions = .dflt_when)
  data_list$fleet_control$Comp_accum_young <- .rce_apply_default(data_list$fleet_control$Comp_accum_young, "Comp_accum_young", .sch)  # young-tail composition accumulation bin (NA -> no accumulation)
  data_list$fleet_control$Comp_accum_old <- .rce_apply_default(data_list$fleet_control$Comp_accum_old, "Comp_accum_old", .sch)  # old-tail composition accumulation bin (NA -> no accumulation)
  data_list$fleet_control$Month <- .rce_apply_default(data_list$fleet_control$Month, "Month", .sch)

  # Format adjustment for NonParametric
  np_idx <- data_list$fleet_control$Selectivity %in% c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM")
  if(any(np_idx & !is.na(data_list$fleet_control$Time_varying_sel) & (!data_list$fleet_control$Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")))){
    data_list$fleet_control <- data_list$fleet_control |>
      dplyr::mutate(
        Sel_curve_pen1 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), Time_varying_sel, Sel_curve_pen1),
        Sel_curve_pen2 = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), Time_varying_sel_sd, Sel_curve_pen2),
        Time_varying_sel = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), 0, Time_varying_sel),
        Time_varying_sel_sd = ifelse(np_idx & (!Time_varying_sel %in% c(NA, 0, 1, "Off", "IID", "RandomWalk")), 0, Time_varying_sel_sd)
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
    if(any(data_list$fleet_control$Sel_norm_bin[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$Sel_norm_bin[flt] <- max_bin
      message(paste0("'Sel_norm_bin' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }

    # - Upper sel normalization bin
    if(any(data_list$fleet_control$Sel_norm_bin_upper[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$Sel_norm_bin_upper[flt] <- max_bin
      message(paste0("'Sel_norm_bin_upper' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
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

  # A fishery's index is predicted from the year-averaged numbers, which has no
  # instant to be taken at, so Month is inert on its index rows (a survey's is
  # read). Say so rather than let a workbook carry a month that does nothing.
  .fsh_idx <- .fleets_with_index(data_list, fitted_only = TRUE)
  if (length(.fsh_idx)) {
    .fc <- data_list$fleet_control
    .fsh_idx <- .fsh_idx[.fsh_idx %in% .fc$Fleet_code[.fc$Fleet_type %in% c(1, "1", "Fishery")]]
    .im <- data_list$index_data
    .dead_month <- .fsh_idx[vapply(.fsh_idx, function(fl) {
      m <- suppressWarnings(as.numeric(.im$Month[.im$Fleet_code == fl & .im$Year > 0]))
      any(!is.na(m) & m != 0)
    }, logical(1))]
    if (length(.dead_month)) {
      message(paste0(
        "'Month' is not read for the index of fishery fleet(s): ",
        paste(.fc$Fleet_name[match(.dead_month, .fc$Fleet_code)], collapse = ", "),
        ". A fishery index is predicted from the year-averaged numbers ",
        "N*(1-exp(-Z))/Z. Put the index on its own Survey fleet, mirroring the ",
        "fishery's Selectivity_index, if it should be timed to a month."))
    }
  }

  # Default Sel_start_year to the fleet's first year of data (not styr). Earlier
  # deviations have neither data nor a penalty (every cpp selectivity penalty is
  # anchored at start_yr), so they are unidentified. Only affects time-varying
  # selectivity; an explicit per-fleet value wins.
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
# Legacy spellings of a canonical switch VALUE (as distinct from the schema's
# column aliases, which rename a whole column). Upgraded on the way in, because
# switch_check() accepted these while validate_switches() rejected them -- so a
# model written with one loaded, ran switch_check() unchanged, and then failed
# its own data check with "Invalid 'Selectivity'".
.rce_switch_value_aliases <- c(
  "Non-parametric" = "NonParametric"
)


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
    # Upgrade a legacy spelling before matching, so it lands on the canonical
    # name rather than passing through to be rejected later.
    alias_hit <- match(x_char, names(.rce_switch_value_aliases))
    if (any(!is.na(alias_hit))) {
      x_char[!is.na(alias_hit)] <-
        unname(.rce_switch_value_aliases[alias_hit[!is.na(alias_hit)]])
    }
    idx <- match(x_char, as.character(map))
    matched <- !is.na(idx)
    if (any(matched)) {
      x_char[matched] <- names(map)[idx[matched]]
    }
    return(x_char)
  }

  # - Fleet switches
  # Guard: Index_distribution and Sel_norm_scope are newer columns; a hand-built
  # fleet_control (or a data_list passed straight to rearrange_data(), bypassing
  # switch_check) may lack them. Fill from the schema so the conversion below
  # works without them, and so the value cannot drift from switch_check()'s.
  .sch_defaults <- .rce_column_schema()
  for (.col in c("Index_distribution", "Sel_norm_scope")) {
    if (is.null(data_list$fleet_control[[.col]]))
      data_list$fleet_control[[.col]] <- .sch_defaults[[.col]]$default
  }
  data_list$fleet_control <- data_list$fleet_control |>
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Sel_norm_scope = .conv(.data$Sel_norm_scope, sel_norm_scope_map),
      Catchability = .conv(.data$Catchability, q_map),
      Comp_distribution = .conv(.data$Comp_distribution, comp_loglike_map),
      CAAL_distribution = .conv(.data$CAAL_distribution, comp_loglike_map),
      Index_distribution = .conv(.data$Index_distribution, index_distribution_map)
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
#' The value map a column is validated against, from the schema
#'
#' @description
#' `.rce_column_schema()` names the governing map on each `type = "switch"` row
#' (`allowed = "sel_map"`). Resolving it here is what makes the schema
#' authoritative for *which* values a column may take: adding a switch column
#' means adding a row, not another hardcoded map reference in this file.
#'
#' The per-column subset predicate and the wording of the error stay at the call
#' site. They are not uniform -- `Time_varying_q` is exempt while `Catchability`
#' is `"Environmental"`, `Catchability` itself allows `NA`, and the newer
#' columns may be absent entirely on a list `switch_check()` has not yet seen --
#' and flattening that into one generic loop would lose real behaviour.
#'
#' @param col Canonical column name.
#' @return The named map, or `NULL` if the schema declares none.
#' @keywords internal
#' @noRd
.rce_allowed_map <- function(col) {
  row <- .rce_column_schema()[[col]]
  if (is.null(row)) {
    stop("Internal: '", col, "' is not in the column schema.", call. = FALSE)
  }
  nm <- row$allowed
  if (is.null(nm) || all(is.na(nm))) {
    # Returning NULL here would make every value invalid and print
    # "Please use one of:" with nothing after it -- loud, but pointing the user
    # at their data instead of at the schema row that is actually missing.
    stop("Internal: the schema declares no `allowed` map for '", col, "'.",
         call. = FALSE)
  }
  get(nm, envir = asNamespace("Rceattle"))
}


validate_switches <- function(data_list = NULL){
  errors <- character(0)

  # Governing maps, resolved from the schema's `allowed` field.
  fleet_map              <- .rce_allowed_map("Fleet_type")
  sel_map                <- .rce_allowed_map("Selectivity")
  tv_sel_map             <- .rce_allowed_map("Time_varying_sel")
  sel_norm_scope_map     <- .rce_allowed_map("Sel_norm_scope")
  q_map                  <- .rce_allowed_map("Catchability")
  tv_q_map               <- .rce_allowed_map("Time_varying_q")
  comp_loglike_map       <- .rce_allowed_map("Comp_distribution")
  index_distribution_map <- .rce_allowed_map("Index_distribution")
  sel_dimension_map      <- .rce_allowed_map("Selectivity_dimension")
  sel_shape_dir_map      <- .rce_allowed_map("Sel_shape_dir")
  sel_shape_mode_map     <- .rce_allowed_map("Sel_shape_mode")

  # Validate fleet_control inputs ----
  # The newer switch columns can be absent here: data_check() is callable on a
  # data list read straight from a workbook, before switch_check() has filled
  # its schema defaults. A column that is not there will be defaulted to a valid
  # value, so there is nothing to validate -- skip the check rather than fail on
  # a missing column.
  .fc_none <- data_list$fleet_control[0, , drop = FALSE]
  .fc_has <- function(col) !is.null(data_list$fleet_control[[col]])

  # Every predicate below exempts an Off fleet. data_check() is callable on a
  # list read straight from a workbook, where Fleet_type may still be the
  # integer 0 -- and `0 != "Off"` coerces to `"0" != "Off"`, which is TRUE, so a
  # bare string comparison would validate the fleets it means to skip.
  # Canonicalize once, and compare against that.
  .fleet_type_canon <- .canon_switch(data_list$fleet_control$Fleet_type, fleet_map)
  data_list$fleet_control$.rce_is_on <- .fleet_type_canon != "Off"

  invalid_flt_type <- data_list$fleet_control |>
    dplyr::filter(!.data$Fleet_type %in% c(fleet_map, names(fleet_map)))

  invalid_sel <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$Selectivity %in% c(sel_map, names(sel_map)))

  invalid_tv_sel <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$Time_varying_sel %in% c(tv_sel_map, names(tv_sel_map)))

  invalid_q <- data_list$fleet_control |>
    dplyr::filter(!.data$Catchability %in% c(NA, q_map, names(q_map)))

  invalid_tv_q <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & .data$Catchability != "Environmental" & !.data$Time_varying_q %in% c(NA, tv_q_map, names(tv_q_map)))

  invalid_sel_norm_scope <- if (.fc_has("Sel_norm_scope")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on &
                      !.data$Sel_norm_scope %in% c(sel_norm_scope_map, names(sel_norm_scope_map)))
  } else .fc_none

  # Three columns that reached the fit unvalidated until 5.12.0. A typo in any
  # of them resolved to NA rather than erroring: `Selectivity_dimension` became
  # a missing selectivity dimension, and the two Sel_shape_* columns a missing
  # penalty mode. Nothing downstream re-checked them.
  # Each of these is validated against exactly what its CONSUMER implements,
  # not against the map alone. rearrange_data() matches Selectivity_dimension on
  # the exact strings and yields NA for anything else -- including an integer --
  # so the integer side of its map must not validate. The two Sel_shape_*
  # columns are read case-insensitively in one spelling each, and Sel_shape_dir
  # additionally accepts "-1" (the ADMB sign convention), so refusing those
  # would reject input the model implements.
  invalid_sel_dim <- if (.fc_has("Selectivity_dimension")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on & !is.na(.data$Selectivity_dimension) &
                      !.data$Selectivity_dimension %in% names(sel_dimension_map))
  } else .fc_none

  .dir_ok <- c(names(sel_shape_dir_map), "increasing", "decreasing", "-1", "0")
  invalid_sel_shape_dir <- if (.fc_has("Sel_shape_dir")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on & !is.na(.data$Sel_shape_dir) &
                      !.data$Sel_shape_dir %in% .dir_ok)
  } else .fc_none

  .mode_ok <- c(names(sel_shape_mode_map), "smooth", "directional",
                unname(sel_shape_mode_map))
  invalid_sel_shape_mode <- if (.fc_has("Sel_shape_mode")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on & !is.na(.data$Sel_shape_mode) &
                      !.data$Sel_shape_mode %in% .mode_ok)
  } else .fc_none

  invalid_comp_ll <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$Comp_distribution %in% c(comp_loglike_map, names(comp_loglike_map)))

  invalid_caal_ll <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$CAAL_distribution %in% c(comp_loglike_map, names(comp_loglike_map)))

  invalid_index_ll <- if (.fc_has("Index_distribution")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on &
                      !.data$Index_distribution %in% c(index_distribution_map, names(index_distribution_map)))
  } else .fc_none

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

  if(nrow(invalid_sel_norm_scope) > 0) {
    errors <- c(errors, paste("Invalid 'Sel_norm_scope' specified for fleets:",
                              paste(invalid_sel_norm_scope$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_norm_scope_map), collapse = ", ")))
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

  if(nrow(invalid_sel_dim) > 0) {
    errors <- c(errors, paste("Invalid 'Selectivity_dimension' specified for fleets:",
                              paste(invalid_sel_dim$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_dimension_map), collapse = ", ")))
  }

  if(nrow(invalid_sel_shape_dir) > 0) {
    errors <- c(errors, paste("Invalid 'Sel_shape_dir' specified for fleets:",
                              paste(invalid_sel_shape_dir$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_shape_dir_map), collapse = ", ")))
  }

  if(nrow(invalid_sel_shape_mode) > 0) {
    errors <- c(errors, paste("Invalid 'Sel_shape_mode' specified for fleets:",
                              paste(invalid_sel_shape_mode$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_shape_mode_map), collapse = ", ")))
  }

  if(nrow(invalid_comp_ll) > 0) {
    errors <- c(errors, paste("Invalid 'Comp_distribution' specified for fleets:",
                              paste(invalid_comp_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(comp_loglike_map), collapse = ":")," or one of:",
                              paste(names(comp_loglike_map), collapse = ", ")))
  }

  if(nrow(invalid_caal_ll) > 0) {
    errors <- c(errors, paste("Invalid 'CAAL_distribution' specified for fleets:",
                              paste(invalid_caal_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(comp_loglike_map), collapse = ":")," or one of:",
                              paste(names(comp_loglike_map), collapse = ", ")))
  }

  if(nrow(invalid_index_ll) > 0) {
    errors <- c(errors, paste("Invalid 'Index_distribution' specified for fleets:",
                              paste(invalid_index_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(index_distribution_map), collapse = ":")," or one of:",
                              paste(names(index_distribution_map), collapse = ", ")))
  }

  # Validate pop dy controls ----
  # * initMode ----
  invalid_initMode <- (!data_list$initMode %in% c(initMode_map, names(initMode_map)))

  if(sum(invalid_initMode) > 0) {
    errors <- c(errors, paste0("Invalid 'initMode' specified: ",
                              paste(unique(data_list$initMode[invalid_initMode]), collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(initMode_map), collapse = ":"), " or one of: ",
                              paste(names(initMode_map), collapse = ", ")))
  }

  # * HCR ----
  # No HCR when the list came from a workbook -- `NULL %in% ...` is logical(0),
  # which propagates through `&` and leaves the `if` below with nothing to test.
  if (length(data_list$HCR)) {
    invalid_hcr <- (!data_list$HCR %in% c(hcr_map, names(hcr_map)))

    if(sum(invalid_hcr) > 0) {
      errors <- c(errors, paste0("Invalid 'HCR' specified: ",
                                paste(unique(data_list$HCR[invalid_hcr]), collapse = ", "),
                                ".\nPlease use an integer code ", paste(range(hcr_map), collapse = ":"), " or one of: ",
                                paste(names(hcr_map), collapse = ", ")))
    }

    if (isTRUE(data_list$msmMode > 0) & !data_list$HCR %in% c("NoFishing", "CMSY", "ConstantF", "ConstantFSSB", "PFMC")) {
      errors <- c(errors, 'Only HCRs "NoFishing" (0), "CMSY" (1), "ConstantF" (2), "ConstantFSSB" (3), or "PFMC" (6) work in multi-species mode currently')
    }
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
  # Guard: default the newer switch columns when a caller supplies a
  # fleet_control that never went through switch_check -- a direct
  # rearrange_data() call (it is exported and documented to work on a data list
  # read straight from a workbook), or a data_check() run before fitting. Every
  # pre-5.8.0 data list, including all the bundled ones, is missing
  # Sel_norm_scope; without this the .data pronoun below fails with a cryptic
  # "column not found". Defaults come from the schema so they cannot drift from
  # the ones switch_check() applies.
  .sch_defaults <- .rce_column_schema()
  for (.col in c("Index_distribution", "Sel_norm_scope")) {
    if (is.null(data_list$fleet_control[[.col]]))
      data_list$fleet_control[[.col]] <- .sch_defaults[[.col]]$default
  }
  # A blank cell in a supplied Sel_norm_scope column would reach TMB as NA, which
  # the C++ reads as WithinSex -- the opposite of the documented default. Fill it
  # here as switch_check() does, so both entry points agree.
  data_list$fleet_control$Sel_norm_scope[is.na(data_list$fleet_control$Sel_norm_scope)] <-
    .sch_defaults[["Sel_norm_scope"]]$default
  # If vector is a string that exists in our map, replace it with the integer.
  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Sel_norm_scope = .conv(.data$Sel_norm_scope, sel_norm_scope_map),
      Catchability = .conv(.data$Catchability, q_map),
      Time_varying_q = .conv(.data$Time_varying_q, tv_q_map),
      Comp_distribution = .conv(.data$Comp_distribution, comp_loglike_map),
      CAAL_distribution = .conv(.data$CAAL_distribution, comp_loglike_map),
      Index_distribution = .conv(.data$Index_distribution, index_distribution_map)
    ) %>%
    # CRITICAL: Force columns back to integers so TMB doesn't crash expecting ints but getting chars
    dplyr::mutate(
      Fleet_type = as.integer(.data$Fleet_type),
      Selectivity = as.integer(.data$Selectivity),
      Time_varying_sel = as.integer(.data$Time_varying_sel),
      Sel_norm_scope = as.integer(.data$Sel_norm_scope),
      Catchability = as.integer(.data$Catchability),
      Comp_distribution = as.integer(.data$Comp_distribution),
      CAAL_distribution = as.integer(.data$CAAL_distribution),
      Index_distribution = as.integer(.data$Index_distribution)
    )

  # `Time_varying_q` is excluded from the blanket coercion above: for "AR1" (6)
  # and "Environmental" (5) catchability it holds 1-based `env_data` column
  # indices rather than a tv_q_map mode, and "Environmental" may carry a
  # comma-separated list ("1,3"). Coerce only the rows holding a mode code.
  .q_env_rows <- data_list$fleet_control$Catchability %in% c(5L, 6L)
  if (!any(.q_env_rows)) {
    data_list$fleet_control$Time_varying_q <-
      as.integer(data_list$fleet_control$Time_varying_q)
  } else if (any(!.q_env_rows)) {
    data_list$fleet_control$Time_varying_q[!.q_env_rows] <-
      as.integer(data_list$fleet_control$Time_varying_q[!.q_env_rows])
  }

  # Pop dy controls ----
  data_list$initMode <- as.integer(.conv(data_list$initMode, initMode_map))
  data_list$HCR <- as.integer(.conv(data_list$HCR, hcr_map))

  return(data_list)
}


#' Report composition weights that are on the log scale
#'
#' @description
#' A composition weight is read on a scale set by its distribution. Under a
#' multinomial it multiplies the input sample size directly; under a
#' Dirichlet-multinomial the model uses `exp(weight)` (`ceattle.cpp`,
#' `DM_pars_comp = exp(comp_weights)`), so the column holds the LOG of the
#' starting weight. A value of 1 is therefore a starting weight of e, and a
#' weight of 1 is written as 0.
#'
#' `write_template()` seeds these columns with 1, so a model built from the
#' template and switched to a Dirichlet-multinomial starts at e without anything
#' saying so. This reports it once per fit, naming the fleets, so the value is a
#' choice rather than an inherited default.
#'
#' Called from `fit_mod()`, not from `switch_check()`: `switch_check()` is a
#' re-entrant normalizer that `build_params()` and `build_map()` also call, so
#' reporting there printed the same message three times per fit and roughly
#' twenty times per `retrospective()`.
#'
#' Only an exact 1 is reported -- the template value. Any other number was typed
#' deliberately and needs no comment. Off fleets and fleets carrying no data for
#' the composition in question are skipped: their weight is never read.
#'
#' @param data_list a data list, after `switch_check()` has resolved the
#'   distribution columns.
#' @return `data_list`, unchanged. This reports; it never edits.
#' @keywords internal
#' @noRd
.rce_flag_dm_weight_scale <- function(data_list) {
  is_dm <- function(x) !is.na(x) & as.character(x) %in% c("1", "DirichletMultinomial")

  # Fleet codes with at least one row actually fitted (Year > 0 marks a fitted
  # row; a negative year is carried but not fitted).
  fitted_codes <- function(df) {
    if (is.null(df) || !nrow(df)) return(integer(0))
    if (!all(c("Fleet_code", "Year") %in% names(df))) return(integer(0))
    unique(as.integer(df$Fleet_code[!is.na(df$Year) & df$Year > 0]))
  }

  fc <- data_list$fleet_control
  if (!is.null(fc)) {
    on <- if (is.null(fc$Fleet_type)) rep(TRUE, nrow(fc)) else
      as.character(fc$Fleet_type) != "Off"

    for (spec in list(list("Comp_distribution", "Comp_weights", data_list$comp_data),
                      list("CAAL_distribution", "CAAL_weights", data_list$caal_data))) {
      dist_col <- spec[[1]]; wt_col <- spec[[2]]
      if (is.null(fc[[dist_col]]) || is.null(fc[[wt_col]])) next
      has_data <- fc$Fleet_code %in% fitted_codes(spec[[3]])
      hit <- which(on & has_data & is_dm(fc[[dist_col]]) &
                     !is.na(fc[[wt_col]]) & fc[[wt_col]] == 1)
      if (length(hit)) {
        who <- if (!is.null(fc$Fleet_name)) fc$Fleet_name[hit] else fc$Fleet_code[hit]
        message("'", wt_col, "' is 1 on Dirichlet-multinomial fleet(s) ",
                paste(who, collapse = ", "),
                ". That likelihood reads the column as a log, so this is a starting ",
                "weight of e (", signif(exp(1), 4), "). Use 0 for a starting weight of 1.")
      }
    }
  }

  dw <- data_list$Diet_comp_weights
  dd <- data_list$Diet_distribution
  # Length-guarded: this runs before data_check(), and recycling a ragged pair
  # here would raise a warning ahead of the check that names the real problem.
  if (!is.null(dw) && !is.null(dd) && length(dw) == length(dd)) {
    hit <- which(is_dm(dd) & !is.na(dw) & dw == 1)
    if (length(hit)) {
      who <- if (!is.null(data_list$spnames)) data_list$spnames[hit] else hit
      message("'Diet_comp_weights' is 1 for Dirichlet-multinomial predator(s) ",
              paste(who, collapse = ", "),
              ". That likelihood reads the value as a log, so this is a starting ",
              "weight of e (", signif(exp(1), 4), "). Use 0 for a starting weight of 1.")
    }
  }

  data_list
}


#' Coerce a string/integer switch argument to its canonical integer code.
#'
#' Shared validation for the integer-returning `build_*()` switch arguments
#' (`srr_fun`, `M1_model`, `sd_plus_group`). Accepts a canonical string key from `map` or a
#' legacy integer code (a value of `map`, or a `deprecated` code); errors on
#' anything else and calls `warn_fn(int)` when a deprecated code is supplied.
#'
#' @param x user input (length-1+ character or numeric).
#' @param map named integer vector mapping canonical string -> integer code.
#' @param what argument name, used in messages.
#' @param deprecated integer codes that are accepted but soft-deprecated.
#' @param warn_fn `function(int)` emitting the deprecation warning, or `NULL`.
#' @param length_exact_one require length 1 (`srr_fun`) vs length >= 1 (`M1_model`).
#' @param legacy_note text appended inside the integer-range error message.
#' @keywords internal
#' @noRd
.coerce_switch_arg <- function(x, map, what,
                               deprecated = integer(0),
                               warn_fn = NULL,
                               length_exact_one = FALSE,
                               legacy_note = "") {
  if (length_exact_one) {
    if (length(x) != 1L) {
      stop(sprintf("`%s` must be length 1", what), call. = FALSE)
    }
  } else if (length(x) == 0L) {
    stop(sprintf("`%s` must have length >= 1", what), call. = FALSE)
  }
  if (is.numeric(x)) {
    int <- as.integer(x)
    allowed <- c(unname(map), deprecated)
    if (anyNA(int) || any(!int %in% allowed)) {
      stop(sprintf(
        "integer `%s` must be one of: %s (= %s%s)",
        what,
        paste(map, collapse = ", "),
        paste(names(map), collapse = "/"),
        legacy_note), call. = FALSE)
    }
    if (length(deprecated) > 0L && any(int %in% deprecated) && !is.null(warn_fn)) {
      warn_fn(int)
    }
    return(int)
  }
  if (is.character(x)) {
    bad <- setdiff(x, names(map))
    if (length(bad) > 0L) {
      stop(sprintf(
        "unknown `%s` value(s): %s; allowed: %s",
        what,
        paste(unique(bad), collapse = ", "),
        paste(names(map), collapse = ", ")), call. = FALSE)
    }
    return(unname(map[x]))
  }
  stop(sprintf("`%s` must be a string or integer; got %s",
               what, class(x)[1]), call. = FALSE)
}
