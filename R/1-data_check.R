#' Check an Rceattle data list for errors. Does not modify the data.
#'
#' @param data_list Rceattle data list
#'
#' @keywords internal
data_check <- function(data_list) {
  errors <- character(0)

  # Upgrade any deprecated fleet_control column names to canonical first, so the
  # checks below (and a standalone data_check() on a legacy data list) see the
  # canonical names. In the fit_mod() pipeline switch_check() has already run, so
  # this is a silent no-op there.
  data_list$fleet_control <-
    .rce_upgrade_fleet_control_aliases(data_list$fleet_control)
  data_list <- .rce_upgrade_data_list_aliases(data_list)

  # ---- Helpers ----
  has_data <- function(df) !is.null(df) && nrow(df) > 0
  fc_num <- function(fc, col, flt){
    if(!col %in% colnames(fc)) return(NA_real_)
    val <- suppressWarnings(as.numeric(fc[[col]][flt]))
    if(length(val) == 0) NA_real_ else val
  }

  # =======================================================================
  # 1. Top-level scalars and run-time switches ----
  # =======================================================================

  # msmMode: Kinzey & Punt (2009) functional responses are deprecated.
  if(!is.null(data_list$msmMode) && data_list$msmMode %in% 3:9){
    errors <- c(errors, paste0(
      "msmMode = ", data_list$msmMode, " selects a Kinzey & Punt (2009) ",
      "predation formulation (Holling I/II/III, predator interference, ",
      "preemption, Hassell-Varley, or Ecosim). These are deprecated and ",
      "have not been validated in the current code base. Use msmMode = 0 ",
      "(single-species), 1 (Holsman et al. 2015 MSVPA), or 2 (Holling ",
      "Type III MSVPA)."
    ))
  }

  # suitMode: length-based modes are not yet implemented.
  # 1 = GammaLength, 3 = LognormalLength, 5 = NormalLength.
  if(any(data_list$suitMode %in% c(1, 3, 5))){
    errors <- c(errors, "Length-based suitability (suitMode 1, 3, or 5) is not yet implemented; use a weight-based mode (2, 4, or 6) or empirical suitability (0)")
  }

  # Catchability = "PowerEquation" is not yet implemented: the power coefficient
  # (index_q_pow) is not built as a parameter and the template does not apply it,
  # so the fleet would silently get a plain estimated q instead.
  if(!is.null(data_list$fleet_control$Catchability) &&
     any(data_list$fleet_control$Catchability %in% c("PowerEquation", 4), na.rm = TRUE)){
    errors <- c(errors, "'PowerEquation' catchability not yet implemented")
  }

  # Catchability = "AR1" is the QAR1 form (Rogers et al. 2024):
  # q = exp(log_q + beta * dev_y), with `index_q_dev` a latent AR1 process and
  # the environmental index an OBSERVATION of it.
  #
  # It does not work. build_map() gates the deviates on
  # `Time_varying_q %in% c("IID", "AR1", "RandomWalk")`, but under this form
  # `Time_varying_q` holds an `env_data` COLUMN INDEX rather than a mode -- so a
  # QAR1 fleet never matches, `index_q_dev` stays mapped out, and q comes back
  # constant. Nothing errors.
  #
  # The damage is not confined to q, so the warning does not stop at "q is
  # constant". Measured on BS2017SS fleet 7, estimateMode = 3: the
  # "Catchability deviates" likelihood row accumulates 54.8 from deviates that
  # are identically zero (the AR1 normalizing constant, plus the environmental
  # index fitted as noise about zero), so the reported objective is not
  # comparable with any other model's; and `index_q_dev_log_sd` is left free
  # with a gradient of nyrs_hind and nothing opposing it, so its sigma is
  # driven to zero. Divergent, not merely flat.
  #
  # Stop rather than warn. A warned fit still returns a summary() that looks
  # ordinary, and nothing downstream can tell that its q is constant or its
  # objective inflated -- which is the failure this package cannot ship. It is
  # also the severity 'PowerEquation' already carries a few lines above, and
  # that switch is merely unimplemented rather than actively divergent.
  #
  # The cost of stopping is now small: GOA pollock 2025 runs the linkage form
  # (../Rceattle-models: GOA pollock/2025/04-fit-and-diagnostics.R), and the
  # remaining Catchability = 6 call sites are Rceattle 3.3.1-era scripts.
  #
  # Note this is a DIFFERENT switch from `Time_varying_q = "AR1"`, which is an
  # AR1 structure on an ordinary "Estimated" q and works correctly.
  if(!is.null(data_list$fleet_control$Catchability) &&
     any(data_list$fleet_control$Catchability %in% c("AR1", 6), na.rm = TRUE)){
    qar1 <- which(data_list$fleet_control$Catchability %in% c("AR1", 6))

    # Name each fleet with the environmental series it points at, so the
    # replacement below can be filled in without decoding the column index.
    # `Time_varying_q` indexes the covariates, i.e. env_data without its Year
    # column. Anything outside 1..ncovs is named generically and left to the
    # range check further down; force it to NA first, because a 0 or a negative
    # index would drop or invert elements rather than return one per fleet.
    covs <- setdiff(colnames(data_list$env_data), "Year")
    # `[[` not `$`: on a data.frame `$` partial-matches, so a fleet_control
    # without `Time_varying_q` would silently hand back `Time_varying_q_sd` --
    # a starting SD read as a covariate index. The NULL branch is defensive
    # only: validate_switches() requires the column, so a data_list without it
    # dies there before this message is ever printed.
    tvq <- data_list$fleet_control[["Time_varying_q"]]
    idx <- if (is.null(tvq)) rep(NA_integer_, length(qar1))
           else suppressWarnings(as.integer(tvq[qar1]))
    idx[is.na(idx) | idx < 1L | idx > length(covs)] <- NA_integer_
    cov_of <- covs[idx]
    cov_of[is.na(cov_of)] <- "<env_data column>"

    codes <- data_list$fleet_control$Fleet_code
    codes <- if (is.null(codes)) qar1 else codes[qar1]
    qar1_flts <- paste0(data_list$fleet_control$Fleet_name[qar1],
                        " (Fleet_code ", codes,
                        ", env_data column '", cov_of, "')")

    errors <- c(errors, paste0(
      "Catchability = 'AR1' (QAR1) is removed: it never worked, and is now an ",
      "error rather than a silently constant q. Fleet(s) ",
      paste(qar1_flts, collapse = ", "), ". The AR1 deviates on log-q were ",
      "never estimated, so q came back constant, the objective carried the ",
      "AR1 normalizing constant and fitted the environmental index as noise ",
      "about zero, and the deviation sd was left free and divergent. Express ",
      "it as a linkage, which implements the Rogers et al. (2024) form ",
      "correctly:\n",
      "  build_catchability(linkages = list(q = linkage_spec(\n",
      "    ~ ar1(1 | Year), by = ~ fleet, fleet = <Fleet_code>,\n",
      "    observe = \"<that fleet's env_data column>\", obs_sd = <its ",
      "measurement SD>)))\n",
      "and pass it to fit_mod(qFun = ). `observe` is what makes it the QAR1 ",
      "form rather than a free AR1 on q: it names the series the deviates are ",
      "observed against. `obs_sd` is that series' measurement SD, which the ",
      "old switch never asked for and you must supply. GOA pollock 2025 is a ",
      "worked example.\n",
      "This is not the same switch as Time_varying_q = 'AR1', which works."))
  }

  # Catchability = "Environmental" (Estimate_q = 5) is superseded by a q
  # linkage. The old form still fits -- it keeps its own C++ path and exact
  # numerics -- but names covariates by position in Time_varying_q rather than
  # by a formula, and cannot carry priors, bounds or a phase. Redirect the
  # user without changing their result.
  if(!is.null(data_list$fleet_control$Catchability) &&
     any(data_list$fleet_control$Catchability %in% c("Environmental", 5),
         na.rm = TRUE)){
    env_flts <- data_list$fleet_control$Fleet_name[
      data_list$fleet_control$Catchability %in% c("Environmental", 5)]
    warning(paste0(
      "Catchability = 'Environmental' is soft-deprecated. Fleet(s) ",
      paste(env_flts, collapse = ", "), " still fit with their current ",
      "numerics, but the environmental effect on q is better expressed as a ",
      "linkage:\n  build_catchability(linkages = list(q = ",
      "linkage_spec(~ your_covariate, by = ~ fleet)))\n",
      "which names the covariate by a formula and can carry priors, bounds ",
      "and an estimation phase."), call. = FALSE)
  }

  # Time_varying_q / Time_varying_sel / M1_re are superseded by random-effect
  # linkages. The legacy switches still fit with their exact numerics (they keep
  # their own C++ path); the formula grammar expresses the same IID / random-walk
  # / AR1 deviations through build_*() with (1 | Year) / rw(1 | Year) /
  # ar1(1 | Year), and additionally allows a prior on -- or free estimation of --
  # the deviation SD. Warn only where a grammar equivalent exists: the
  # environmental / Rogers-AR1 catchability modes overload Time_varying_q to name
  # env columns (not a time-varying mode), the non-parametric selectivity forms
  # have no additive slot, and the separable M1_re (age x year, code 6) has no
  # 1-D grammar structure -- those stay on the legacy path without a nudge.
  fc <- data_list$fleet_control
  if (!is.null(fc$Time_varying_q)) {
    tvq_dep <- fc$Time_varying_q %in% c(1, 2, 4, "IID", "AR1", "RandomWalk")
    q_est   <- fc$Catchability %in% c(1, 2, "Estimated", "Estimated-with-prior")
    tvq_flts <- fc$Fleet_name[tvq_dep & q_est]
    if (length(tvq_flts) > 0) {
      warning(paste0(
        "Time_varying_q is soft-deprecated for fleet(s) ",
        paste(tvq_flts, collapse = ", "), ". They still fit with their current ",
        "numerics, but a time-varying catchability is better expressed as a ",
        "random-effect linkage:\n  build_catchability(linkages = list(q = ",
        "linkage_spec(~ (1 | Year), by = ~ fleet)))\n",
        "(use rw(1 | Year) for a random walk, ar1(1 | Year) for AR1); see ",
        "vignette('environmental-linkages-and-priors')."), call. = FALSE)
    }
  }
  if (!is.null(fc$Time_varying_sel)) {
    tvs_dep  <- fc$Time_varying_sel %in%
      c(1, 2, 4, 5, "IID", "AR1", "RandomWalk", "RandomWalkAscending")
    sel_para <- fc$Selectivity %in%
      c(1, 3, 4, 8, 11, "Logistic", "DoubleLogistic", "DescendingLogistic",
        "DoubleNormal", "LogisticPM")
    tvs_flts <- fc$Fleet_name[tvs_dep & sel_para]
    if (length(tvs_flts) > 0) {
      warning(paste0(
        "Time_varying_sel is soft-deprecated for fleet(s) ",
        paste(tvs_flts, collapse = ", "), ". They still fit with their current ",
        "numerics, but time-varying parametric selectivity is better expressed ",
        "as a random-effect linkage on the relevant parameter, e.g.\n",
        "  build_selectivity(linkages = list(inf_asc = ",
        "linkage_spec(~ (1 | Year), by = ~ fleet)))\n",
        "(use rw(1 | Year) for a random walk); see ",
        "vignette('environmental-linkages-and-priors')."), call. = FALSE)
    }
  }
  if (!is.null(data_list$M1_re)) {
    # modes 1-5 (iid_age / iid_year / iid_age_year / ar1_age / ar1_year) have a
    # grammar equivalent; mode 6 (separable ar1_age_year) does not -- exclude it.
    m1_dep <- data_list$M1_re %in%
      c(1, 2, 3, 4, 5, "iid_age", "iid_year", "iid_age_year", "ar1_age", "ar1_year")
    if (any(m1_dep, na.rm = TRUE)) {
      m1_sp <- if (!is.null(data_list$spnames)) data_list$spnames[which(m1_dep)] else which(m1_dep)
      warning(paste0(
        "M1_re is soft-deprecated for species ", paste(m1_sp, collapse = ", "),
        ". It still fits with its current numerics, but time-varying natural ",
        "mortality is better expressed as a random-effect linkage:\n",
        "  build_M1(linkages = list(M1 = linkage_spec(~ (1 | Year))))\n",
        "(use ar1(Year) / ar1(age_bin) for correlated deviations; note the ",
        "grammar AR1 uses the marginal SD); see ",
        "vignette('environmental-linkages-and-priors')."), call. = FALSE)
    }

    # The age-by-year modes carry one deviation per age-year cell -- nages *
    # nyrs_hind per species, which is 882 on GOA2018SS where a mapping defect
    # gave 42 before 5.9.0. That field is flexible enough to absorb a trend that
    # belongs to selectivity or to recruitment, and nothing else in a
    # single-species model pins the level of M. Say so where the user can act
    # on it; this is not an error, and a well-informed model is free to run it.
    m1_2d <- data_list$M1_re %in% c(3, 6, "iid_age_year", "ar1_age_year")
    m1_pr <- data_list$M1_use_prior
    if (is.null(m1_pr)) m1_pr <- rep(0, length(m1_2d))
    m1_2d <- m1_2d & !(as.numeric(m1_pr) > 0)
    if (any(m1_2d, na.rm = TRUE)) {
      m1_sp <- if (!is.null(data_list$spnames)) data_list$spnames[which(m1_2d)] else which(m1_2d)
      warning(paste0(
        "M1_re estimates an age-by-year deviation field for species ",
        paste(m1_sp, collapse = ", "), " with no prior on M ",
        "(M1_use_prior = FALSE). That is one deviation per age and year, which ",
        "is confounded with selectivity and recruitment in a single-species ",
        "model. Set build_M1(M1_use_prior = TRUE), or check that M is ",
        "identified before reading the estimates."), call. = FALSE)
    }
  }

  # `Q_block` is vestigial: it is never read (q time-blocking reuses the
  # `Selectivity_block` column, R/3-build_map.R), so a supplied `Q_block` is
  # silently ignored. Warn once so users don't rely on it.
  if (any(vapply(c("index_data", "catch_data"),
                 function(b) "Q_block" %in% names(data_list[[b]]), logical(1)))) {
    warning("'Q_block' in index_data/catch_data is deprecated and ignored: q ",
            "time-blocking (Time_varying_q = 'Block') uses the 'Selectivity_block' ",
            "column, not 'Q_block'.", call. = FALSE)
  }

  # minage: < 0 error
  if(any(data_list$minage < 0)){
    errors <- c(errors, "Minimum age is < 0. Check 'minage'.")
  }

  # Year ordering
  if(!is.null(data_list$styr) && !is.null(data_list$endyr) && data_list$styr > data_list$endyr){
    errors <- c(errors, paste0("styr (", data_list$styr, ") must be <= endyr (", data_list$endyr, ")"))
  }
  if(!is.null(data_list$projyr) && !is.null(data_list$endyr) && data_list$projyr < data_list$endyr){
    errors <- c(errors, paste0("projyr (", data_list$projyr, ") must be >= endyr (", data_list$endyr, ")"))
  }
  # suit_styr/suit_endyr are per-predator (length nspp) or a recycled scalar.
  # The suitability window must sit inside the hindcast [styr, endyr]: the C++
  # loops index the year arrays by the offset window, so a start below styr
  # reads out of bounds (segfault in MakeADFun) and an end past endyr silently
  # pads the nyrs_suit divisor with empty years, under-normalising predation
  # mortality on modelled prey.
  if(!is.null(data_list$suit_styr) && !is.null(data_list$suit_endyr) && any(data_list$suit_styr > data_list$suit_endyr)){
    bad <- which(data_list$suit_styr > data_list$suit_endyr)
    errors <- c(errors, paste0("suit_styr must be <= suit_endyr for every predator; violated for species ",
                               paste(bad, collapse = ", "),
                               " (suit_styr = ", paste(data_list$suit_styr[bad], collapse = ", "),
                               " > suit_endyr = ", paste(data_list$suit_endyr[bad], collapse = ", "), ")"))
  }
  if(!is.null(data_list$suit_styr) && !is.null(data_list$styr) && any(data_list$suit_styr < data_list$styr)){
    bad <- which(data_list$suit_styr < data_list$styr)
    errors <- c(errors, paste0("suit_styr must be >= styr (", data_list$styr,
                               ") for every predator; violated for species ", paste(bad, collapse = ", "),
                               " (suit_styr = ", paste(data_list$suit_styr[bad], collapse = ", "), ")"))
  }
  if(!is.null(data_list$suit_endyr) && !is.null(data_list$endyr) && any(data_list$suit_endyr > data_list$endyr)){
    bad <- which(data_list$suit_endyr > data_list$endyr)
    errors <- c(errors, paste0("suit_endyr must be <= endyr (", data_list$endyr,
                               ") for every predator; violated for species ", paste(bad, collapse = ", "),
                               " (suit_endyr = ", paste(data_list$suit_endyr[bad], collapse = ", "), ")"))
  }

  # =======================================================================
  # 2. Per-species dimensions ----
  # =======================================================================

  if(data_list$nspp != max(data_list$weight$Species)){
    errors <- c(errors, "`nspp` does not match the number of species in the weight data. Check `nspp` or `weight`")
  }
  if(length(data_list$spnames)     != data_list$nspp) errors <- c(errors, "'spnames' (species names) not included for all species (must have length nspp)")
  if(length(data_list$spawn_month) != data_list$nspp) errors <- c(errors, "'spawn_month' not included for all species")
  if(length(data_list$nages)       != data_list$nspp) errors <- c(errors, "'nages' not included for all species")
  if(length(data_list$nlengths)    != data_list$nspp) errors <- c(errors, "'nlengths' not included for all species")
  if(length(data_list$other_food)  != data_list$nspp) errors <- c(errors, "'other_food' not included for all species")

  # spawn_month is used as a fraction of year (spawn_month/12 in TMB), so 0 is
  # a valid sentinel meaning "spawn at start of year".
  if(!is.null(data_list$spawn_month) &&
     any(data_list$spawn_month < 0 | data_list$spawn_month > 12, na.rm = TRUE)){
    errors <- c(errors, "spawn_month values must be in 0:12")
  }

  # =======================================================================
  # 3. Biology: weight, maturity, sex_ratio, M1, ration, ALK, age error, ----
  #             NByageFixed, bioenergetics
  # =======================================================================

  if(length(data_list$M1_base) == 1)    errors <- c(errors, "M1 is a single value, please make it age/species specific")
  if(sum(data_list$other_food < 0) > 0) errors <- c(errors, paste0("'other_food' is negative for species ", paste(which(data_list$other_food < 0), collapse = ", "), "; it must be >= 0"))

  # Weight: year coverage (only when time-varying)
  wt_yr <- data_list$weight |>
    dplyr::group_by(Wt_index, Sex) |>
    dplyr::distinct(Year) |>
    dplyr::mutate(Tmp_ind = paste0("index = ", Wt_index, " & sex = ", Sex))
  for(ind in unique(wt_yr$Tmp_ind)){
    tmp_wt <- wt_yr |> dplyr::filter(Tmp_ind == ind) |> dplyr::distinct(Year) |> dplyr::pull(Year)
    if(length(tmp_wt) > 1 && any(!(data_list$styr:data_list$endyr) %in% tmp_wt)){
      errors <- c(errors, paste0("Weight data for ", ind, " does not span all hindcast years"))
    }
  }

  # Weight: pop_wt_index / ssb_wt_index reference valid Wt_index values
  wt_index <- data_list$weight |> dplyr::distinct(Wt_index, Species, Sex)
  if(any(!data_list$pop_wt_index %in% wt_index$Wt_index)){
    errors <- c(errors, "'pop_wt_index' references a Wt_index that is not present in the 'weight' data")
  }
  if(any(!data_list$ssb_wt_index %in% wt_index$Wt_index)){
    errors <- c(errors, "'ssb_wt_index' references a Wt_index that is not present in the 'weight' data")
  }

  # Weight: Wt_index must be unique per species
  for(sp in 1:data_list$nspp){
    wt_no  <- data_list$weight |> dplyr::filter(Species != sp) |> dplyr::distinct(Wt_index) |> dplyr::pull(Wt_index)
    wt_yes <- data_list$weight |> dplyr::filter(Species == sp) |> dplyr::distinct(Wt_index) |> dplyr::pull(Wt_index)
    if(any(wt_no %in% wt_yes)){
      errors <- c(errors, "Check weight indices (Wt_index), the same weight index was used for multiple species")
    }
  }

  # Weight / maturity / sex_ratio age coverage
  if(any(data_list$weight |>
         dplyr::select(-c(Wt_name, Wt_index, Species, Sex, Year)) |>
         ncol() < data_list$nages)){
    errors <- c(errors, "Weight data does not span range of ages")
  }
  if(ncol(data_list$maturity)  <= max(data_list$nages)) errors <- c(errors, "Maturity-at-age (maturity) does not span all ages")
  if(ncol(data_list$sex_ratio) <= max(data_list$nages)) errors <- c(errors, "Sex ratio does not span all ages")

  # Maturity / sex_ratio value ranges
  mat_vals <- data_list$maturity[, grep("^Age", colnames(data_list$maturity)), drop = FALSE]
  if(any(mat_vals < 0, na.rm = TRUE)){ # | mat_vals > 1
    errors <- c(errors, "maturity values must be > 0")
  }
  sr_vals <- data_list$sex_ratio[, grep("^Age", colnames(data_list$sex_ratio)), drop = FALSE]
  if(any(sr_vals < 0 | sr_vals > 1, na.rm = TRUE)){
    errors <- c(errors, "sex_ratio values must be in [0, 1]")
  }

  # Sex consistency: max Sex in M1_base / weight / ration_data must be <= nsex
  m1_sex <- data_list$M1_base |> dplyr::group_by(Species) |>
    dplyr::summarise(max_sex = max(Sex)) |> dplyr::arrange(Species)
  if(any(m1_sex$max_sex > data_list$nsex)){
    errors <- c(errors, "'M1_base' has more sexes than specified in 'nsex'")
  }
  wt_sex <- data_list$weight |> dplyr::group_by(Species) |>
    dplyr::summarise(max_sex = max(Sex)) |> dplyr::arrange(Species)
  if(any(wt_sex$max_sex > data_list$nsex)){
    errors <- c(errors, "'weight' has more sexes than specified in 'nsex'")
  }

  # ration_data
  if(has_data(data_list$ration_data)){
    if(any(data_list$ration_data |>
           dplyr::select(-c(Species, Sex, Year)) |>
           ncol() < data_list$nages)){
      errors <- c(errors, "'ration_data' data does not span range of ages")
    }
    ration_sex <- data_list$ration_data |> dplyr::group_by(Species) |>
      dplyr::summarise(max_sex = max(Sex)) |> dplyr::arrange(Species)
    if(any(ration_sex$max_sex > data_list$nsex)){
      errors <- c(errors, "'ration_data' has more sexes than specified in 'nsex'")
    }
  }
  # ration_data presence (msmMode > 0): message-only requirement (see the
  # declarative requirement table). Not an error -- missing ration data is a
  # notice, not a blocker.
  .rce_notify_presence(data_list, "ration_data")

  # Age transition matrix: length coverage when growth is fixed and length data are used
  if(any(data_list$growth_model == 0) &
     has_data(data_list$comp_data |> dplyr::filter(Age0_Length1 == 1)) &
     any(data_list$age_trans_matrix |>
         dplyr::select(-c(Age_transition_name, Age_transition_index, Species, Sex, Age)) |>
         ncol() < data_list$nlengths)){
    errors <- c(errors, "`age_trans_matrix` data does not span range of lengths")
  }

  # Age error matrix: observed-age column count
  if(any(data_list$age_error |>
         as.data.frame() |>
         dplyr::select(-c(Species, True_age)) |>
         ncol() < data_list$nages)){
    errors <- c(errors, "`age_error` observed ages do not span range of ages")
  }
  # ALK & age_error: per-species age coverage (fillable with 0s downstream -- message-level)
  for(sp in 1:data_list$nspp){
    expected_ages <- data_list$minage[sp]:(data_list$minage[sp] + data_list$nages[sp] - 1)

    atm_ages <- data_list$age_trans_matrix |> as.data.frame() |>
      dplyr::filter(Species == sp) |> dplyr::pull(Age)
    if(!all(expected_ages %in% atm_ages)){
      message(paste("`age_trans_matrix` data does not span range of age for species", sp, "will fill with 0s"))
    }

    ae_ages <- data_list$age_error |> as.data.frame() |>
      dplyr::filter(Species == sp) |> dplyr::pull(True_age)
    if(!all(expected_ages %in% ae_ages)){
      message(paste("`age_error` data does not span range of true ages for species", sp, "will fill with 0s"))
    }
  }
  # ALK / age_error: row sums (warning -- rearrange may renormalize)
  if(has_data(data_list$age_error)){
    ae_cols <- setdiff(colnames(data_list$age_error), c("Species", "True_age"))
    if(length(ae_cols) > 0){
      ae_sums <- rowSums(data_list$age_error[, ae_cols, drop = FALSE], na.rm = TRUE)
      if(any(ae_sums > 0 & abs(ae_sums - 1) > 1e-3)){
        warning("Some `age_error` rows do not sum to 1")
      }
    }
  }
  if(has_data(data_list$age_trans_matrix) &
     has_data(data_list$comp_data |> dplyr::filter(Age0_Length1 == 1))
  ){
    atm_cols <- setdiff(colnames(data_list$age_trans_matrix),
                        c("Age_transition_name", "Age_transition_index", "Species", "Sex", "Age"))
    if(length(atm_cols) > 0){
      atm_sums <- rowSums(data_list$age_trans_matrix[, atm_cols, drop = FALSE], na.rm = TRUE)
      if(any(atm_sums > 0 & abs(atm_sums - 1) > 1e-3)){
        warning("Some `age_trans_matrix` rows do not sum to 1")
      }
    }
  }

  # NByageFixed: presence required when estDynamics > 0 (declarative requirement
  # table); the column-count adequacy check stays imperative below.
  errors <- c(errors, .rce_check_presence(data_list, "NByageFixed"))
  if(any(data_list$estDynamics > 0) && has_data(data_list$NByageFixed)){
    expected_cols <- 4 + max(data_list$nages)
    if(ncol(data_list$NByageFixed) != expected_cols){
      errors <- c(errors, paste0("NByageFixed should have ", expected_cols,
                                 " columns (Species_name, Species, Sex, Year, Age1...Age",
                                 max(data_list$nages), "), but has ", ncol(data_list$NByageFixed)))
    }
  }

  # Bioenergetics: temperature-dependent consumption requires environmental driver
  if(!is.null(data_list$Ceq)){
    # NA reaches here from a workbook whose bioenergetics columns were left
    # blank. Report it rather than letting `if (NA > 1)` throw R's bare
    # "missing value where TRUE/FALSE needed", which names nothing.
    .ceq_na <- which(is.na(data_list$Ceq[1:data_list$nspp]))
    if(length(.ceq_na)){
      errors <- c(errors, paste0("`Ceq` is missing for species ",
                                 paste(.ceq_na, collapse = ", "),
                                 ". Set the consumption equation (1 = temperature-independent) for every species."))
    }
    for(sp in setdiff(1:data_list$nspp, .ceq_na)){
      if(data_list$Ceq[sp] > 1){
        if(is.null(data_list$env_data) || ncol(data_list$env_data) < (data_list$Cindex[sp] + 1)){
          errors <- c(errors, paste0("Species ", sp, " uses temperature-dependent consumption (Ceq = ",
                                     data_list$Ceq[sp], ") but Cindex (", data_list$Cindex[sp],
                                     ") exceeds available environmental indices"))
        }
      }
    }
  }

  # =======================================================================
  # 4. Fleet control (structure, references, per-fleet bin/q checks) ----
  # =======================================================================
  if(!has_data(data_list$fleet_control)) {
    stop("Missing 'fleet_control' from data.")
  } else{
    fc <- data_list$fleet_control

    # Fleet_code uniqueness; Species in 1:nspp
    fcodes <- if("Fleet_code" %in% colnames(fc)) suppressWarnings(as.numeric(fc$Fleet_code)) else integer(0)
    fsp    <- if("Species"    %in% colnames(fc)) suppressWarnings(as.numeric(fc$Species))    else integer(0)
    if(any(duplicated(fcodes[!is.na(fcodes)]))){
      errors <- c(errors, "fleet_control$Fleet_code values must be unique")
    }
    # Fleet_code is used directly as the fleet slot of the per-fleet parameter
    # and map arrays, which are built in fleet_control row order. The two must
    # therefore agree, or parameters are silently attached to the wrong fleet.
    if(length(fcodes) == nrow(fc) && !identical(as.integer(fcodes), seq_len(nrow(fc)))){
      errors <- c(errors, paste0(
        "fleet_control$Fleet_code must equal the row number (1, 2, ... ", nrow(fc),
        "); got ", paste(fcodes, collapse = ", "),
        ". Fleet_code indexes the per-fleet parameter/map arrays, which are built in row order."))
    }
    if(any(!is.na(fsp) & (fsp < 1 | fsp > data_list$nspp))){
      errors <- c(errors, paste0("fleet_control$Species values must be in 1:", data_list$nspp))
    }

    # Month range (0 = unspecified sentinel)
    if(!is.null(fc$Month)){
      fm <- suppressWarnings(as.numeric(fc$Month))
      if(any(!is.na(fm) & (fm < 0 | fm > 12))){
        errors <- c(errors, "fleet_control$Month values must be in 0:12")
      }
    }

    # Proj_F_proportion: must be 0 (NoFishing) or sum to 1 across fleets
    flt_spp <- unique(fc$Species)
    for(sp in flt_spp){
      fc_spp <- fc |> dplyr::filter(Species == sp &
                                      Fleet_type %in% c(1, "Fishery"))
      total_proj_f <- sum(fc_spp$Proj_F_proportion, na.rm = TRUE)

      # HCR is a fit_mod() argument, not a data field, so a list read from a
      # workbook carries none and `NULL %in% ...` is logical(0). Nothing to
      # contradict, so nothing to report.
      hcr <- data_list$HCR
      fishing_planned <- length(hcr) == 1L && !is.na(hcr) &&
        !(hcr %in% c(0, "NoFishing"))
      if(total_proj_f == 0 && fishing_planned){
        errors <- c(errors, "HCR is > 0 and 'Proj_F_proportion' is 0")
      }

      if(total_proj_f > 0 && abs(total_proj_f - 1) > 1e-6){
        errors <- c(errors, paste0("Proj_F_proportion values should sum to 1 for species ", sp, ", but sum to ", total_proj_f))
      }
    }

    # Per-fleet checks
    for(flt in 1:nrow(fc)){
      sp_idx <- fc_num(fc, "Species", flt)
      if(is.na(sp_idx) || sp_idx < 1 || sp_idx > data_list$nspp) next

      flt_name   <- fc$Fleet_name[flt]
      dim_is_age <- !isTRUE(fc$Selectivity_dimension[flt] == "Length")
      max_bin    <- if(dim_is_age) data_list$nages[sp_idx] else data_list$nlengths[sp_idx]

      # Bin_first_selected
      bfs <- fc_num(fc, "Bin_first_selected", flt)
      if(!is.na(bfs) && (bfs < 1 || bfs > max_bin)){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Bin_first_selected (", bfs, ") must be in 1:", max_bin))
      }

      # N_sel_bins
      nsb <- fc_num(fc, "N_sel_bins", flt)
      if(!is.na(nsb) && (nsb < 1 || nsb > max_bin)){
        errors <- c(errors, paste0("Fleet '", flt_name, "': N_sel_bins (", nsb, ") must be in 1:", max_bin))
      }

      # Non-parametric shape-penalty range and cap, given on the fleet's own
      # selectivity dimension: an age (from minage) for age-based fleets, a
      # 1-based length-bin ordinal for length-based. Out-of-range values would
      # index past the selectivity array in the template.
      bin_lo <- if(dim_is_age) data_list$minage[sp_idx] else 1L
      bin_hi <- bin_lo + max_bin - 1L
      for(col in c("Sel_pen_first_bin", "Sel_pen_last_bin", "Sel_cap_bin")){
        val <- fc_num(fc, col, flt)
        if(!is.na(val) && (val < bin_lo || val > bin_hi)){
          errors <- c(errors, paste0("Fleet '", flt_name, "': ", col, " (", val, ") must be in ",
                                     bin_lo, ":", bin_hi,
                                     if(dim_is_age) " (age-based selectivity)" else " (length-based selectivity)"))
        }
      }
      pf <- fc_num(fc, "Sel_pen_first_bin", flt); pl <- fc_num(fc, "Sel_pen_last_bin", flt)
      if(!is.na(pf) && !is.na(pl) && pf > pl){
        errors <- c(errors, paste0("Fleet '", flt_name, "': Sel_pen_first_bin (", pf,
                                   ") must be <= Sel_pen_last_bin (", pl, ")"))
      }

      # Composition young/old tail-accumulation bins (Comp_accum_young/old) are
      # 1-based ordinals on the fleet's COMPOSITION dimension -- age or length,
      # from comp_data$Age0_Length1 -- and are PER SEX BLOCK for joint-sex (Sex 3)
      # comps (so the bound is nages/nlengths, not the doubled joint row). An
      # out-of-range value or young >= old would fold into a nonexistent bin,
      # build a negative-length vector in the template, or collapse the whole
      # composition into a single (zero-information) bin, so reject them here.
      # A single per-fleet column drives every composition row on the fleet, so
      # the bound is the MOST restrictive dimension present (min): a value that
      # is out of range for one dimension would otherwise silently no-op on it.
      if(!is.null(data_list$comp_data) && nrow(data_list$comp_data) > 0){
        a0l1 <- unique(data_list$comp_data$Age0_Length1[
          data_list$comp_data$Fleet_code == fc$Fleet_code[flt]])
        a0l1 <- a0l1[!is.na(a0l1)]
        if(length(a0l1) > 0){
          comp_max_bin <- min(ifelse(a0l1 == 0, data_list$nages[sp_idx],
                                                 data_list$nlengths[sp_idx]))
          ay <- fc_num(fc, "Comp_accum_young", flt)
          ao <- fc_num(fc, "Comp_accum_old",   flt)
          # Effective old bin: the NA / 0 sentinels both mean "no old accumulation",
          # which the template reads as the last bin.
          ao_eff <- if(is.na(ao) || ao == 0) comp_max_bin else ao
          if(!is.na(ay) && (ay < 1 || ay > comp_max_bin)){
            errors <- c(errors, paste0("Fleet '", flt_name, "': Comp_accum_young (", ay,
                                       ") must be in 1:", comp_max_bin,
                                       " (a per-sex-block bin on the fleet's composition dimension)"))
          } else if(!is.na(ay) && ay >= ao_eff){
            errors <- c(errors, paste0("Fleet '", flt_name, "': Comp_accum_young (", ay,
                                       ") must be < Comp_accum_old (", ao_eff,
                                       "); folding into a single bin discards the composition"))
          }
          if(!is.na(ao) && ao != 0 && (ao < 1 || ao > comp_max_bin)){
            errors <- c(errors, paste0("Fleet '", flt_name, "': Comp_accum_old (", ao,
                                       ") must be in 1:", comp_max_bin, " (or 0/NA for no old accumulation)"))
          }
        }
      }

      # AR1 catchability: Time_varying_q must be a valid env_data column index (1..ncol-1)
      if(!is.na(fc$Catchability[flt]) && fc$Catchability[flt] == "AR1"){
        if((fc$Time_varying_q[flt] > (ncol(data_list$env_data) - 1) ||
            fc$Time_varying_q[flt] < 1)){
          errors <- c(errors, "For 'AR1' catchability, environmental index specified in 'Time_varying_q' is greater than number of indices in 'env_data'")
        }
      }

      # Time-varying form is selectivity-type specific:
      #  - NonParametric (Ianelli, type 2) / NonParametricPM (type 9): the cpp
      #    penalizes year-to-year log selectivity-at-age, i.e. a RANDOM WALK ->
      #    allow only "Off"/"RandomWalk".
      #  - Hake (Taylor, type 5): IID coefficient deviates -> allow only "Off"/"IID".
      if(fc$Selectivity[flt] %in% c("NonParametric", "NonParametricPM") &&
         !fc$Time_varying_sel[flt] %in% c("Off", "RandomWalk")){
        errors <- c(errors, "For 'NonParametric'/'NonParametricPM' selectivity, 'Time_varying_sel' must be 'Off' or 'RandomWalk'")
      }
      if(fc$Selectivity[flt] == "Hake" &&
         !fc$Time_varying_sel[flt] %in% c("Off", "IID")){
        errors <- c(errors, "For 'Hake' selectivity, 'Time_varying_sel' must be 'Off' or 'IID'")
      }
      #  - LogisticPM (ADMB AMAK "pm" BTS, type 11): random-walk deviates on
      #    slope/inflection/age-1 -> allow only "Off"/"RandomWalk".
      if(fc$Selectivity[flt] == "LogisticPM" &&
         !fc$Time_varying_sel[flt] %in% c("Off", "RandomWalk")){
        errors <- c(errors, "For 'LogisticPM' selectivity, 'Time_varying_sel' must be 'Off' or 'RandomWalk'")
      }

      # Non-parametric (Ianelli) selectivity penalties (Sel_curve_pen1 =
      # decreasing penalty, Sel_curve_pen2 = curvature) must be present and
      # numeric to identify the free coefficients. Catch the case where they are
      # missing / non-numeric (e.g. a Time_varying_sel mode string accidentally
      # written into Sel_curve_pen) before it surfaces as a cryptic
      # "inits not within bounds" error in build_bounds.
      if(fc$Selectivity[flt] %in% c("NonParametric", "NonParametricPM")){
        cp1 <- suppressWarnings(as.numeric(fc$Sel_curve_pen1[flt]))
        cp2 <- suppressWarnings(as.numeric(fc$Sel_curve_pen2[flt]))
        if(is.na(cp1) || is.na(cp2)){
          errors <- c(errors, paste0(
            "Fleet '", fc$Fleet_name[flt], "' has Selectivity = '", fc$Selectivity[flt], "' but ",
            "'Sel_curve_pen1'/'Sel_curve_pen2' are missing or non-numeric (got '",
            fc$Sel_curve_pen1[flt], "' / '", fc$Sel_curve_pen2[flt],
            "'). Non-parametric selectivity requires numeric curvature/smoothness penalties."))
        }
      }

      # LogisticPM: Sel_curve_pen1/2/3 are the random-walk weights on the
      # slope / inflection / age-1 deviates (ADMB 50 / 50 / 8). Require numeric
      # when time-varying so a stray mode string is caught early.
      if(fc$Selectivity[flt] == "LogisticPM" && fc$Time_varying_sel[flt] == "RandomWalk"){
        cps <- suppressWarnings(as.numeric(c(fc$Sel_curve_pen1[flt], fc$Sel_curve_pen2[flt], fc$Sel_curve_pen3[flt])))
        if(any(is.na(cps))){
          errors <- c(errors, paste0(
            "Fleet '", fc$Fleet_name[flt], "' has Selectivity = 'LogisticPM' with time-varying ",
            "selectivity but 'Sel_curve_pen1'/'Sel_curve_pen2'/'Sel_curve_pen3' (slope/inflection/age-1 ",
            "random-walk weights) are missing or non-numeric."))
        }
      }

      # 2DAR1/3DAR1: 'Sel_curve_pen1'/'Sel_curve_pen2' are reused as logit-scale AR1
      # correlation parameters via rho_trans(x) = 2/(1+exp(-2x)) - 1. |x| > ~10 saturates
      # rho at +-1 and the AR1 log-density evaluates to NaN.
      if(fc$Selectivity[flt] %in% c("2DAR1", "3DAR1")){
        for(col in c("Sel_curve_pen1", "Sel_curve_pen2")){
          val <- suppressWarnings(as.numeric(fc[[col]][flt]))
          if(!is.na(val) && abs(val) > 10){
            errors <- c(errors, sprintf(
              "Fleet '%s' has Selectivity = '%s'. For 2DAR1/3DAR1, '%s' sets the AR1 correlation between bins on a transformed scale that maps to (-1, 1); a magnitude above 10 pushes the correlation to +-1 and the fit fails to evaluate. Got %s = %s. Suggested replacement: %s = 0 (no correlation). Use a value in roughly [-5, 5] to set a non-trivial correlation.",
              flt_name, fc$Selectivity[flt], col, col, val, col
            ))
          }
        }
      }
    }

    # emp_sel presence required when any fleet has Selectivity = "Fixed"
    # (declarative requirement table).
    errors <- c(errors, .rce_check_presence(data_list, "emp_sel"))

    # Estimated selectivity (Selectivity != "Fixed" and Fleet_type != "Off")
    # requires comp or CAAL data with Year > 0 to be identifiable. Otherwise
    # the selectivity parameters are unconstrained and the optimizer wanders.
    # EXCEPTION: a fleet whose Selectivity_index is shared (mirrored) with
    # another fleet that DOES have active comp/CAAL data is identifiable through
    # the master fleet's data, so it is not flagged.
    has_active_age_data <- function(flt_code, df) {
      if (!has_data(df) || !all(c("Fleet_code", "Year") %in% colnames(df))) return(FALSE)
      any(df$Fleet_code == flt_code & !is.na(df$Year) & df$Year > 0 & df$Sample_size > 0)
    }
    # Fleets sharing a Selectivity_index share one deviation block, so a
    # differing Sel_start_year within the group is resolved to the group minimum
    # (the earliest year any sharing fleet has data). Surface it so the
    # resolution is deliberate rather than silent.
    if (all(c("Selectivity_index", "Sel_start_year") %in% colnames(fc))) {
      for (si in unique(fc$Selectivity_index[!is.na(fc$Selectivity_index)])) {
        rows <- which(!is.na(fc$Selectivity_index) & fc$Selectivity_index == si)
        if (length(rows) > 1) {
          ys <- suppressWarnings(as.integer(fc$Sel_start_year[rows]))
          ys[is.na(ys)] <- as.integer(data_list$styr)
          if (length(unique(ys)) > 1) {
            warning(paste0("Fleets sharing Selectivity_index ", si, " (",
                           paste(fc$Fleet_name[rows], collapse = ", "),
                           ") have different Sel_start_year (", paste(ys, collapse = ", "),
                           "); the shared selectivity deviations use the earliest (", min(ys), ")."))
          }
        }
      }
    }

    # Per-fleet settings a shared Selectivity_index does not reconcile. Checked
    # here rather than in build_map(), whose warnings fit_mod() suppresses.
    #
    # SHAPING columns are read per fleet by the cpp when it builds the curve, so
    # a difference means the group does not share one selectivity. NA counts as
    # its own value among them -- a blank Sel_norm_bin means "do not normalize",
    # so blank against 2 is two different curves, not "inherit the lead's".
    #
    # Time_varying_sel is resolved by build_map(), which copies the lead fleet's
    # deviation map over the group: the curves match and the other fleets'
    # settings are discarded. NA there really is "unset", so it is skipped.
    # Sel_norm_scope and Sel_cap_bin belong here too: both are per-fleet
    # DATA_IVECTORs read inside the curve builder (selectivity.hpp: the
    # across-sex normalization reference, and the NonParametricRPM bin cap),
    # not behind a flt_sel_lead gate.
    .sel_shaping_cols <- c("Selectivity", "Selectivity_dimension",
                           "Bin_first_selected", "N_sel_bins",
                           "Sel_norm_bin", "Sel_norm_bin_upper",
                           "Sel_norm_scope", "Sel_cap_bin")
    if ("Selectivity_index" %in% colnames(fc)) {
      .live <- if ("Fleet_type" %in% colnames(fc)) {
        !(fc$Fleet_type %in% c("Off", 0, "0"))
      } else rep(TRUE, nrow(fc))
      for (si in unique(fc$Selectivity_index[!is.na(fc$Selectivity_index) & .live])) {
        rows <- which(!is.na(fc$Selectivity_index) & fc$Selectivity_index == si & .live)
        if (length(rows) < 2) next
        .differs <- function(col, na_is_a_value) {
          if (!col %in% colnames(fc)) return(FALSE)
          v <- as.character(fc[[col]][rows])
          if (na_is_a_value) v[is.na(v)] <- "<blank>" else v <- v[!is.na(v)]
          length(unique(v)) > 1
        }
        shaping <- Filter(function(cl) .differs(cl, na_is_a_value = TRUE),
                          .sel_shaping_cols)
        if (length(shaping)) {
          warning(paste0(
            "Fleets sharing Selectivity_index ", si, " (",
            paste(fc$Fleet_name[rows], collapse = ", "), ") differ in ",
            paste(shaping, collapse = ", "),
            ". These are read per fleet when the curve is built, so the fleets ",
            "will not share one selectivity. To mirror a fleet, copy its ",
            "fleet_control row and change only the identity and catchability ",
            "columns."))
        }
        if (.differs("Time_varying_sel", na_is_a_value = FALSE)) {
          warning(paste0(
            "Fleets sharing Selectivity_index ", si, " (",
            paste(fc$Fleet_name[rows], collapse = ", "),
            ") have different Time_varying_sel; the shared deviation block uses ",
            "the first estimated fleet's setting and the others are ignored."))
        }
      }
    }

    # The same for a shared Catchability_index, and for the same reason.
    #
    # Fixed / Estimated / Estimated-with-prior share the one index_log_q the map
    # wires up, so a difference resolves to the lead fleet's answer.
    #
    # Analytical and AnalyticalArith do not: they solve q from the fleet's own
    # OBSERVATIONS (ceattle.cpp 8.2, 8.2b), bypassing the shared parameter, so a
    # group containing one shares no catchability -- two Analytical fleets still
    # solve separately. Reported on the form, not only on a disagreement.
    #
    # Environmental and AR1 also overwrite index_q per fleet (ceattle.cpp 6.4),
    # but from index_log_q, index_q_beta and index_q_dev, all of which build_map()
    # maps to the lead fleet's, against an env_index row that is not fleet
    # specific. Those groups DO share, including when the fleets name different
    # env series -- the lead's is used. Verified by fitting, not by inspection.
    .q_solved <- c("Analytical", "AnalyticalArith")
    # Accept either spelling: a data list reaching data_check() straight from a
    # workbook may still carry the integer codes.
    .q_canon <- function(x) .canon_switch(x, q_map)
    if (all(c("Catchability_index", "Catchability") %in% colnames(fc))) {
      .live <- if ("Fleet_type" %in% colnames(fc)) {
        !(fc$Fleet_type %in% c("Off", 0, "0"))
      } else rep(TRUE, nrow(fc))
      for (qi in unique(fc$Catchability_index[!is.na(fc$Catchability_index) & .live])) {
        rows <- which(!is.na(fc$Catchability_index) & fc$Catchability_index == qi & .live)
        if (length(rows) < 2) next
        qv <- .q_canon(fc$Catchability[rows])

        # The lead settles the form for the group: build_map_catchability()
        # reads it, and adjust_map_shared_params() copies that slice over the
        # rest. .group_lead() picks the same row.
        lead_form <- qv[1]
        init <- if ("Catchability_init" %in% colnames(fc)) {
          suppressWarnings(as.numeric(fc$Catchability_init[rows]))
        } else rep(NA_real_, length(rows))

        solved <- intersect(unique(qv), .q_solved)
        if (identical(lead_form, "Fixed")) {
          # Nothing is estimated, so there is no shared parameter to hold: the
          # group's index_log_q is mapped out and each fleet sits at its own
          # Catchability_init. Reported whenever that is not what the columns
          # ask for -- a member wanting an estimated q, or inits that differ.
          if (length(unique(qv)) > 1 ||
              length(unique(init[!is.na(init)])) > 1) {
            warning(paste0(
              "Fleets sharing Catchability_index ", qi, " (",
              paste(fc$Fleet_name[rows], collapse = ", "),
              ") have a Fixed lead fleet, so no catchability is estimated for ",
              "the group and none is shared: each fleet uses its own ",
              "Catchability_init (",
              paste(ifelse(is.na(init), "<blank>", format(init)), collapse = ", "),
              ")."))
          }
        } else if (length(solved)) {
          warning(paste0(
            "Fleets sharing Catchability_index ", qi, " (",
            paste(fc$Fleet_name[rows], collapse = ", "), ") include ",
            paste(solved, collapse = ", "),
            ", which computes catchability per fleet. The group does not share ",
            "one catchability despite sharing the index; give each fleet its own ",
            "Catchability_index, or use an estimated q for the whole group."))
        } else if (length(unique(qv)) > 1) {
          warning(paste0(
            "Fleets sharing Catchability_index ", qi, " (",
            paste(fc$Fleet_name[rows], collapse = ", "),
            ") have different Catchability (", paste(unique(qv), collapse = ", "),
            "); the shared catchability uses the first estimated fleet's setting ",
            "and the others are ignored."))
        }

        # Time_varying_q is overloaded -- under Environmental and AR1 it names
        # env_data columns rather than a mode -- but either way the group takes
        # the lead fleet's, so a difference is worth reporting in both readings.
        if ("Time_varying_q" %in% colnames(fc)) {
          tv <- as.character(fc$Time_varying_q[rows])
          tv <- tv[!is.na(tv)]
          if (length(unique(tv)) > 1) {
            warning(paste0(
              "Fleets sharing Catchability_index ", qi, " (",
              paste(fc$Fleet_name[rows], collapse = ", "),
              ") have different Time_varying_q; the group uses the first ",
              "estimated fleet's setting and the others are ignored."))
          }
        }
      }
    }

    est_sel_flts <- fc[!is.na(fc$Selectivity) &
                         fc$Selectivity != "Fixed" &
                         (!"Fleet_type" %in% colnames(fc) | fc$Fleet_type != "Off"),
                       , drop = FALSE]
    # Selectivity_index values that have active age data in ANY sharing fleet
    sel_idx_has_data <- if ("Selectivity_index" %in% colnames(fc)) {
      vapply(unique(fc$Selectivity_index), function(si) {
        codes <- fc$Fleet_code[!is.na(fc$Selectivity_index) & fc$Selectivity_index == si]
        any(vapply(codes, function(cc)
          has_active_age_data(cc, data_list$comp_data) ||
            has_active_age_data(cc, data_list$caal_data), logical(1)))
      }, logical(1))
    } else NULL
    if (!is.null(sel_idx_has_data)) names(sel_idx_has_data) <- as.character(unique(fc$Selectivity_index))
    if (nrow(est_sel_flts) > 0) {
      missing_age_data <- vapply(seq_len(nrow(est_sel_flts)), function(i) {
        fc_code <- est_sel_flts$Fleet_code[i]
        own_data <- has_active_age_data(fc_code, data_list$comp_data) ||
          has_active_age_data(fc_code, data_list$caal_data)
        # mirrored identifiability: another fleet sharing this Selectivity_index
        # supplies the comp/CAAL data
        si <- est_sel_flts$Selectivity_index[i]
        mirror_data <- !is.null(sel_idx_has_data) && !is.na(si) &&
          isTRUE(sel_idx_has_data[[as.character(si)]])
        !own_data && !mirror_data
      }, logical(1))
      if (any(missing_age_data)) {
        errors <- c(errors, paste0(
          "Fleet(s) with estimated Selectivity but no comp_data or caal_data ",
          "rows in the likelihood (all Year/Sample_size == 0 or missing): ",
          paste(est_sel_flts$Fleet_name[missing_age_data], collapse = ", "),
          ". Either provide composition / CAAL data, mark Selectivity = 'Fixed' ",
          "with emp_sel, or set Fleet_type = 'Off'."
        ))
      }
    }


    # Mirroring (informational) NOW in configuration
    # mirror_sel <- fc |> dplyr::group_by(Selectivity_index) |>
    #   dplyr::filter(dplyr::n() > 1) |> dplyr::ungroup()
    # if(nrow(mirror_sel) > 0){
    #   message(paste0("Selectivity for ", paste(mirror_sel$Fleet_name, collapse = ", "),
    #                  " is mirrored with another fleet"))
    # }
    # mirror_q <- fc |> dplyr::filter(!is.na(Catchability)) |>
    #   dplyr::group_by(Catchability_index) |> dplyr::filter(dplyr::n() > 1) |> dplyr::ungroup()
    # if(nrow(mirror_q) > 0){
    #   message(paste0("Catchability for ", paste(mirror_q$Fleet_name, collapse = ", "),
    #                  " is mirrored with another fleet"))
    # }
  }

  # =======================================================================
  # 5. Observation tables: catch_data, index_data, comp_data, caal_data ----
  # =======================================================================

  # Required columns
  required_cols_by_table <- list(
    index_data = c("Fleet_code", "Year", "Observation"),
    catch_data = c("Fleet_code", "Year", "Catch"),
    comp_data  = c("Fleet_code", "Species", "Sex", "Year")
  )
  for(df_name in names(required_cols_by_table)){
    df <- data_list[[df_name]]
    if(has_data(df)){
      missing_cols <- setdiff(required_cols_by_table[[df_name]], colnames(df))
      if(length(missing_cols) > 0){
        errors <- c(errors, paste(df_name, "is missing required columns:",
                                  paste(missing_cols, collapse = ", ")))
      }
    }
  }

  # Fleet_code referential integrity
  if(has_data(data_list$fleet_control)){
    fcodes <- suppressWarnings(as.numeric(data_list$fleet_control$Fleet_code))
    valid_codes <- fcodes[!is.na(fcodes)]
    for(df_name in c("catch_data", "index_data", "comp_data", "caal_data")){
      df <- data_list[[df_name]]
      if(has_data(df) && "Fleet_code" %in% colnames(df)){
        bad <- setdiff(unique(df$Fleet_code), valid_codes)
        if(length(bad) > 0){
          errors <- c(errors, paste0(df_name, " references Fleet_code(s) not in fleet_control: ",
                                     paste(bad, collapse = ", ")))
        }
      }
    }
  }

  # Lognormal SDs must be > 0
  for(df_name in c("index_data", "catch_data")){
    df <- data_list[[df_name]]
    if(has_data(df) && "Log_sd" %in% colnames(df) && any(!is.na(df$Log_sd) & df$Log_sd <= 0 & df$Year > 0)){
      errors <- c(errors, paste0(df_name, " has 'Log_sd' <= 0; lognormal likelihood requires positive SD"))
    }
  }

  # Sample_size >= 0
  for(df_name in c("comp_data", "caal_data", "diet_data")){
    df <- data_list[[df_name]]
    if(has_data(df) && "Sample_size" %in% colnames(df) && any(df$Sample_size < 0 & df$Year > 0, na.rm = TRUE)){
      errors <- c(errors, paste0(df_name, " contains 'Sample_size' values < 0"))
    }
  }

  # Observation values
  if(has_data(data_list$index_data) && "Observation" %in% colnames(data_list$index_data) &&
     any(!is.na(data_list$index_data$Observation) & data_list$index_data$Observation <= 0 & data_list$index_data$Year > 0)){
    errors <- c(errors, "index_data$Observation must be > 0 (lognormal likelihood breaks at 0)")
  }
  # NA Catch is allowed (clean_data sets NA for projection rows)
  catch_df <- data_list$catch_data |>
    dplyr::filter(Year <= data_list$endyr)
  if(has_data(catch_df) && "Catch" %in% colnames(catch_df) &&
     any(is.na(catch_df$Catch) | catch_df$Catch < 0)){
    errors <- c(errors, "catch_data$Catch must be >= 0")
  }

  # MVN survey covariance requirement ----
  # A fleet using Index_distribution == "MVN" or "MVNORM" must supply a square,
  # symmetric variance-covariance matrix in data_list$index_cov (keyed by
  # Fleet_name or Fleet_code) whose dimension equals the number of fitted survey
  # observations for that fleet. Validated here so the requirement surfaces with a
  # clear message at the flag rather than as a cryptic error in rearrange_data().
  fc <- data_list$fleet_control
  # Initialised here, not inside the branch: with no Index_distribution column at
  # all there are no covariance fleets, which is precisely the case the stray
  # index_cov warning below is meant to catch.
  mvn_flts <- integer(0)
  # The analytical sd (Ludwig and Walters 1994) is accumulated from squared LOG
  # residuals, so it is a log-scale sd. What that costs depends on whether the
  # family actually reads it, and the two groups differ:
  #
  #   Normal / TruncatedNormal read the sd as an ABSOLUTE value in index units,
  #     so the likelihood itself is evaluated on the wrong scale. Refuse it.
  #   MVN / MVNORM score through index_cov_mat and never read the scalar sd, so
  #     the FIT is unaffected. But index_sd is still reported from it, and that
  #     is what residuals(type = "pearson") and plot_index()'s interval divide
  #     by, so the diagnostics are on the wrong scale. Warn rather than refuse a
  #     model that fits correctly.
  if(has_data(fc) && all(c("Index_distribution", "Estimate_index_sd") %in% colnames(fc))){
    is_analytical <- fc$Estimate_index_sd %in% c("Analytical", 2, "2")
    is_on <- !(fc$Fleet_type %in% c("Off", 0, "0"))
    fitted_scale <- fc$Index_distribution %in% c("Normal", "TruncatedNormal", 3, 4, "3", "4")
    cov_scale    <- fc$Index_distribution %in% c("MVN", "MVNORM", 1, 2, "1", "2")

    bad <- which(is_analytical & is_on & fitted_scale)
    if(length(bad)){
      errors <- c(errors, paste0(
        "Fleet(s) ", paste(fc$Fleet_name[bad], collapse = ", "),
        " combine Estimate_index_sd = 'Analytical' with Index_distribution = '",
        paste(unique(as.character(fc$Index_distribution[bad])), collapse = "', '"),
        "'. The analytical sd is computed from log residuals, so it is a ",
        "log-scale sd, while these families read the sd as an absolute value in ",
        "the units of the index -- the likelihood would be evaluated on the ",
        "wrong scale. Use Estimate_index_sd = 'Fixed' with an absolute Log_sd, ",
        "or 'Estimated', or switch the fleet to Lognormal."))
    }

    noisy <- which(is_analytical & is_on & cov_scale)
    if(length(noisy)){
      warning(paste0(
        "Fleet(s) ", paste(fc$Fleet_name[noisy], collapse = ", "),
        " combine Estimate_index_sd = 'Analytical' with a covariance index ",
        "family. The fit is unaffected -- MVN/MVNORM score through index_cov ",
        "and never read the scalar sd -- but the reported index_sd is then a ",
        "log-scale number, and residuals(type = 'pearson') and plot_index()'s ",
        "observation interval divide by it. Read those two on this fleet with ",
        "care, or set Estimate_index_sd = 'Fixed'."), call. = FALSE)
    }
  }

  if(has_data(fc) && "Index_distribution" %in% colnames(fc)){
    mvn_flts <- which(fc$Index_distribution %in% c("MVN", "MVNORM", 1, 2, "1", "2"))
    for(flt in mvn_flts){
      flt_name <- fc$Fleet_name[flt]
      flt_code <- fc$Fleet_code[flt]
      Sigma <- NULL
      if(!is.null(data_list$index_cov)){
        Sigma <- data_list$index_cov[[as.character(flt_name)]]
        if(is.null(Sigma)) Sigma <- data_list$index_cov[[as.character(flt_code)]]
      }
      if(is.null(Sigma)){
        errors <- c(errors, paste0("Fleet '", flt_name, "' has Index_distribution == 'MVN' but no covariance matrix ",
                                   "was found in data_list$index_cov (expected an element named '", flt_name,
                                   "' or '", flt_code, "')."))
        next
      }
      Sigma <- as.matrix(Sigma)
      n_fit <- 0
      if(has_data(data_list$index_data)){
        n_fit <- sum(data_list$index_data$Fleet_code == flt_code &
                       data_list$index_data$Year > 0 &
                       data_list$index_data$Year <= data_list$endyr &
                       data_list$index_data$Observation > 0, na.rm = TRUE)
      }
      if(nrow(Sigma) != ncol(Sigma)){
        errors <- c(errors, paste0("index_cov matrix for fleet '", flt_name, "' must be square (got ",
                                   nrow(Sigma), " x ", ncol(Sigma), ")."))
      } else if(nrow(Sigma) != n_fit){
        errors <- c(errors, paste0("index_cov matrix for fleet '", flt_name, "' is ", nrow(Sigma), " x ",
                                   ncol(Sigma), " but the fleet has ", n_fit, " fitted survey observations ",
                                   "(Year in [styr, endyr], Observation > 0). Sigma must match those rows in index_data order."))
      } else if(!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-4, check.attributes = FALSE))){
        errors <- c(errors, paste0("index_cov matrix for fleet '", flt_name, "' is not symmetric."))
      } else if(!.is_positive_definite(Sigma)){
        # The covariance likelihood factorizes Sigma (MVNORM / Eigen::LLT), so a
        # symmetric but indefinite or singular matrix clears every check above
        # and then fails inside the TMB objective, where the message names
        # neither the fleet nor Sigma.
        ev <- suppressWarnings(tryCatch(
          min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values),
          error = function(e) NA_real_))
        errors <- c(errors, paste0(
          "index_cov matrix for fleet '", flt_name, "' is not positive definite",
          if(is.finite(ev)) paste0(" (smallest eigenvalue ", signif(ev, 3), ")") else "",
          ". Sigma is factorized by the covariance index likelihood, so it must ",
          "be positive definite, not only symmetric."))
      }
    }
  }

  # An index_cov entry for a fleet that is not using the covariance likelihood is
  # ignored both here and in .align_index_cov(), so the fleet quietly fits its own
  # index likelihood (lognormal by default) and the covariance has no effect.
  # Almost always this means Index_distribution was never set to 'MVN'/'MVNORM'.
  if(has_data(fc) && !is.null(data_list$index_cov) && length(data_list$index_cov)){
    mvn_keys <- character(0)
    if(length(mvn_flts)){
      mvn_keys <- unique(c(as.character(fc$Fleet_name[mvn_flts]),
                           as.character(fc$Fleet_code[mvn_flts])))
    }
    stray <- setdiff(names(data_list$index_cov), mvn_keys)
    if(length(stray)){
      warning(paste0(
        "data_list$index_cov has ", if(length(stray) > 1) "entries" else "an entry",
        " for ", paste0("'", stray, "'", collapse = ", "),
        " but ", if(length(stray) > 1) "those fleets are" else "that fleet is",
        " not using Index_distribution 'MVN' or 'MVNORM'. The covariance is ",
        "ignored and the fleet fits its own index likelihood (lognormal by ",
        "default). Set Index_distribution to use it."))
    }
  }

  # Duplicates per (Fleet_code, Year, Month)
  for(df_name in c("catch_data", "index_data")){
    df <- data_list[[df_name]]
    if(has_data(df) && all(c("Fleet_code", "Year", "Month") %in% colnames(df))){
      n_dup <- sum(duplicated(df[, c("Fleet_code", "Year", "Month")]))
      if(n_dup > 0){
        errors <- c(errors, paste0(df_name, " has ", n_dup,
                                   " duplicated row(s) for (Fleet_code, Year, Month)"))
      }
    }
  }

  # catch_data must span hindcast years (use 0 where no catch occurred);
  # A fleet carrying fitted index observations needs its catchability columns,
  # whatever its Fleet_type. The template scores an index row for any non-Off
  # fleet, so a fishery with a CPUE series is fitted like a survey -- but these
  # columns have no schema default, so on a fishery they arrive NA and the index
  # would be fitted at an undefined q with an undefined sd. Required rather than
  # defaulted: guessing a catchability form for someone's CPUE series is the
  # kind of silent default this check exists to prevent.
  #
  # These two are read whatever the catchability form. The rest are conditional
  # and are handled below -- Catchability_init in particular is unread under
  # Analytical, which several working GOA hake configurations rely on.
  .idx_fleets <- .fleets_with_index(data_list)
  if (length(.idx_fleets)) {
    .fc  <- data_list$fleet_control
    .need <- c("Catchability", "Estimate_index_sd")
    .have <- intersect(.need, names(.fc))
    for (.cl in .have) {
      .bad <- .fc$Fleet_code %in% .idx_fleets & is.na(.fc[[.cl]])
      if (any(.bad)) {
        stop("Fleet(s) ", paste(unique(.fc$Fleet_name[.bad]), collapse = ", "),
             " carry index_data but have no `", .cl, "`. An index is fitted for ",
             "any fleet that is not Off -- a fishery with a CPUE series included ",
             "-- so it needs the same catchability settings a survey does. Set `",
             .cl, "` for those fleets. Setting their index rows to Year < 0 ",
             "stops them being fitted but does not remove the fleet's index: ",
             "its catchability stays unmapped, so index_hat comes back NA and ",
             "reaches sdreport. Switch the fleet Off instead if it should carry ",
             "no index at all.", call. = FALSE)
      }
    }
  }

  # The remaining q columns are each required by the one switch that reads them.
  # All are logged when the parameter list is built, so a blank or non-positive
  # entry becomes a NaN or -Inf starting value and the objective is not finite
  # at the first evaluation -- loudly, but from inside MakeADFun, where the
  # message names neither the fleet nor the column.
  #
  # Every condition mirrors build_map()'s own gate, so a setting that reads no
  # starting value is not asked for one. `Block` is deliberately absent from the
  # time-varying set (a time block carries no penalty), and the analytical
  # catchability forms are exempted below.
  if (length(.idx_fleets)) {
    .fc   <- data_list$fleet_control
    .isq  <- .fc$Fleet_code %in% .idx_fleets
    .col  <- function(nm) if (nm %in% names(.fc)) .fc[[nm]] else rep(NA, nrow(.fc))
    # Read the switches by canonical name whichever spelling they arrived in.
    # fit_mod() and read_data() both run switch_check() -> revert_switches()
    # first, but data_check() is callable on a hand-built list that has not, and
    # the shared-block checks above already canonicalize for the same reason.
    .qform <- .canon_switch(.col("Catchability"), q_map)
    .qtv   <- .canon_switch(.col("Time_varying_q"), tv_q_map)
    .qest  <- .qform %in% c("Estimated", "Estimated-with-prior")

    .req <- list(
      list(col = "Index_sd",
           when = .isq & .col("Estimate_index_sd") %in% c(1, "1"),
           why  = "`Estimate_index_sd` estimates the observation sd, and its starting value is log(Index_sd)"),
      list(col = "Catchability_prior_sd",
           when = .isq & .qform %in% c("Estimated-with-prior", "AR1"),
           why  = "the catchability prior is scored at it, and AR1 starts its estimated sd from it"),
      list(col = "Time_varying_q_sd",
           when = .isq & (.qest | .qform == "AR1") &
                  .qtv %in% c("IID", "AR1", "RandomWalk"),
           why  = "the time-varying catchability deviations are penalized at this standard deviation"),
      # Analytical / AnalyticalArith solve q in closed form and overwrite
      # index_q (ceattle.cpp section 8.2), so they never read this column --
      # several GOA hake configurations leave it at 0 and fit correctly.
      list(col = "Catchability_init",
           when = .isq & !(.qform %in% c("Analytical", "AnalyticalArith")),
           why  = "catchability is held at, or starts from, log(Catchability_init)")
    )

    for (.r in .req) {
      .w <- .r$when
      .w[is.na(.w)] <- FALSE
      if (!any(.w)) next
      .v   <- suppressWarnings(as.numeric(.col(.r$col)))
      .bad <- .w & (is.na(.v) | !is.finite(.v) | .v <= 0)
      if (any(.bad)) {
        stop("Fleet(s) ", paste(unique(.fc$Fleet_name[.bad]), collapse = ", "),
             " need a positive `", .r$col, "`: ", .r$why, ". Without it the ",
             "objective is not finite at the first evaluation, which surfaces ",
             "as a TMB error naming none of this.", call. = FALSE)
      }
    }
  }

  # index_data gaps are normal (biennial / triennial surveys, missed years).
  if(!is.null(data_list$styr) && !is.null(data_list$endyr) && has_data(data_list$catch_data)){
    missing_years <- setdiff(data_list$styr:data_list$endyr, unique(data_list$catch_data$Year))
    if(length(missing_years) > 0){
      errors <- c(errors, paste0("catch_data is missing data for years: ",
                                 paste(missing_years, collapse = ", ")))
    }
  }

  # Empirical growth (growth_model == 0) does NOT populate growth_matrix from
  # age_trans_matrix; the C++ pred_CAAL = N * sel * growth_matrix then sees
  # growth_matrix = 0 and pred_CAAL collapses to ~0. The multinomial NLL
  # evaluates to a constant offset with no useful gradient, and the optimizer
  # finds spurious minima. Check per-species since growth_model is per-species.
  if(has_data(data_list$caal_data) && "Species" %in% colnames(data_list$caal_data)){
    bad_sp <- integer(0)
    for(sp in seq_along(data_list$growth_model)){
      if(data_list$growth_model[sp] == 0 &&
         any(data_list$caal_data$Species == sp & data_list$caal_data$Year > 0  & data_list$caal_data$Sample_size > 0)){
        bad_sp <- c(bad_sp, sp)
      }
    }
    if(length(bad_sp) > 0){
      sp_names <- if(!is.null(data_list$spnames)) data_list$spnames[bad_sp] else as.character(bad_sp)
      errors <- c(errors, paste0(
        "Empirical growth (growth_model == 0) is incompatible with CAAL data ",
        "for species: ", paste(sp_names, collapse = ", "), ". The C++ ",
        "growth_matrix is not populated from age_trans_matrix in the empirical ",
        "branch, so pred_CAAL = 0 and the CAAL likelihood gradient is ",
        "uninformative. Either (a) switch to parametric growth for these ",
        "species via build_growth(fun = 'vonBertalanffy'), or (b) drop ",
        "caal_data rows for these species (set to empty, Sample_size = 0, or Year < 0)."
      ))
    }
  }

  # The sibling of the rule above, for the other reason pred_CAAL comes back
  # zero. CAAL is a composition of ages WITHIN a length bin, so the prediction is
  # selectivity-at-length convolved with the growth matrix (ceattle.cpp, section
  # 10.2). Selectivity at length only exists for a length-dimensioned fleet:
  # selectivity.hpp writes sel_at_length only under `is_length_based`, leaving it
  # zero otherwise, so an age-dimensioned fleet predicts nothing.
  #
  # It does not fail quietly in the harmless sense -- the CAAL likelihood is
  # still evaluated, against a prediction that is uniform once comp_offset is
  # added, so the observations are scored against a flat composition and the
  # objective carries a term that no parameter can move. Selectivity_dimension
  # defaults to "Age", so this is the default outcome rather than an unusual one.
  if(has_data(data_list$caal_data) &&
     all(c("Fleet_code", "Year", "Sample_size") %in% colnames(data_list$caal_data)) &&
     !is.null(data_list$fleet_control$Selectivity_dimension)){
    fc <- data_list$fleet_control
    active <- data_list$caal_data$Year > 0 & data_list$caal_data$Sample_size > 0
    bad_flt <- character(0)
    for(flt in unique(data_list$caal_data$Fleet_code[active])){
      i <- which(fc$Fleet_code == flt)
      if(!length(i)) next
      if(!identical(as.character(fc$Selectivity_dimension[i[1]]), "Length")){
        bad_flt <- c(bad_flt, as.character(fc$Fleet_name[i[1]]))
      }
    }
    # A warning rather than an error, unlike the empirical-growth rule above,
    # only because this is a released package and such a model currently fits.
    # It should become an error once existing configurations have been checked:
    # the CAAL data in one of these fits inform nothing at all.
    if(length(bad_flt) > 0){
      warning(paste0(
        "CAAL data on fleet(s) ", paste(bad_flt, collapse = ", "),
        " whose Selectivity_dimension is not 'Length'. CAAL is the age ",
        "composition within a length bin, so it is predicted from ",
        "selectivity-at-length; an age-dimensioned fleet has none, so pred_CAAL ",
        "is zero and the CAAL likelihood scores these observations against a ",
        "flat composition it cannot move -- they add a term to the objective ",
        "but inform nothing. Either set Selectivity_dimension = 'Length' for ",
        "these fleets, or drop their caal_data rows (set to empty, ",
        "Sample_size = 0, or Year < 0)."), call. = FALSE)
    }
  }

  # CAAL: presence required when growth is being estimated (declarative
  # requirement table); the column / length adequacy checks stay imperative.
  errors <- c(errors, .rce_check_presence(data_list, "caal_data"))
  if(any(data_list$growth_model > 0) && has_data(data_list$caal_data)){
    missing_cols <- setdiff(c("Fleet_code", "Species", "Year", "Length", "Sample_size"),
                            colnames(data_list$caal_data))
    if(length(missing_cols) > 0){
      errors <- c(errors, paste("caal_data is missing required columns:",
                                paste(missing_cols, collapse = ", ")))
    }
    for(sp in 1:data_list$nspp){
      sp_lengths <- unique(data_list$caal_data$Length[data_list$caal_data$Species == sp])
      if(length(sp_lengths) != data_list$nlengths[sp]){
        errors <- c(errors, paste0("Species ", sp, " has ", length(sp_lengths),
                                   " unique lengths in caal_data, but nlengths[", sp, "] = ",
                                   data_list$nlengths[sp]))
      }
    }
    caal_cols <- grep("^CAAL_", colnames(data_list$caal_data), value = TRUE)
    if(length(caal_cols) == 0){
      errors <- c(errors, "caal_data is missing CAAL_ columns (CAAL_1, CAAL_2, etc.)")
    } else {
      missing_caal <- setdiff(paste0("CAAL_", 1:max(data_list$nages)), caal_cols)
      if(length(missing_caal) > 0){
        errors <- c(errors, paste("caal_data is missing CAAL columns:",
                                  paste(missing_caal, collapse = ", ")))
      }
    }
  }

  # =======================================================================
  # 6. Diet & predation ----
  # =======================================================================

  if(has_data(data_list$diet_data)){
    dd <- data_list$diet_data

    # Pred / prey age bounds, against each species' own oldest age. `nages` is a
    # COUNT of age bins, so that is minage + nages - 1. Index by the species the
    # group is for: group_by() drops species absent from the column, so lining
    # the maxima up by position pairs species with the wrong limit. match()
    # sends an out-of-range species code to NA rather than to a short vector or
    # an error.
    maxage_for <- function(sp){
      i <- match(sp, seq_along(data_list$nages))
      data_list$minage[i] + data_list$nages[i] - 1L
    }
    pred_max <- dd |> dplyr::group_by(Pred) |>
      dplyr::summarise(Max_age = max(Pred_age))
    if(any(pred_max$Max_age > maxage_for(pred_max$Pred), na.rm = TRUE)){
      errors <- c(errors, "Pred ages in 'diet_data' > oldest age ('minage' + 'nages' - 1)")
    }
    prey_max <- dd |> dplyr::group_by(Prey) |>
      dplyr::summarise(Max_age = max(Prey_age))
    if(any(prey_max$Max_age > maxage_for(prey_max$Prey), na.rm = TRUE)){
      errors <- c(errors, "Prey ages in 'diet_data' > oldest age ('minage' + 'nages' - 1)")
    }
    # Duplicates
    if(sum(duplicated(dd)) > 0){
      errors <- c(errors, "Diet data includes duplicated rows")
    }
    # Stomach proportion ranges
    if("Stomach_proportion_by_weight" %in% colnames(dd)){
      if(any(dd$Stomach_proportion_by_weight < 0, na.rm = TRUE)){
        errors <- c(errors, "diet_data contains negative Stomach_proportion_by_weight values")
      }
      diet_sum <- dd |> dplyr::group_by(Pred, Pred_age, Pred_sex, Year) |>
        dplyr::summarise(diet_sum = sum(Stomach_proportion_by_weight))
      # A stomach whose prey account for the whole diet sums to exactly 1, and a
      # simulated one reaches that through a division that can land a bit above.
      # Reject a real excess, not floating-point noise: the worst case there is
      # a few ulps per prey bin, so 1e-12 is ample and stays strict enough to
      # catch a genuinely malformed stomach.
      if(any(diet_sum$diet_sum > 1 + 1e-12)){
        errors <- c(errors, "Stomach proportion in `diet_data` for some predators-at-age/sex/year is > 1")
      }
    }
    # Stomach grouping. The TMB diet likelihood walks diet_ctl with one forward
    # cursor, taking stomach i's prey as the run of rows where stomach_id == i
    # (ceattle.cpp, section 13.2). That needs the ids sorted: 0, 1, 2, ... with
    # no gaps. Sorted order is what makes each stomach's rows consecutive AND
    # puts them where the cursor looks for them, so testing only that the rows
    # are grouped is not enough -- a table whose blocks are each intact but out
    # of order (say re-sorted by predator age) passes that test while the cursor
    # runs past nearly all of them. Every stomach it misses drops out of the
    # likelihood silently, with a lower jnll. clean_data() sorts by stomach_id,
    # so anything that came through it is fine; this catches a hand-built or
    # re-sorted diet table.
    if("stomach_id" %in% colnames(dd)){
      sid <- as.integer(dd$stomach_id)
      if(any(is.na(sid))){
        errors <- c(errors, "'stomach_id' in 'diet_data' contains NA")
      } else if(length(sid) > 0){
        if(!identical(sort(unique(sid)), 0:max(sid))){
          errors <- c(errors, paste0(
            "'stomach_id' in 'diet_data' must be numbered 0, 1, ... with no gaps ",
            "(the TMB diet likelihood indexes stomachs by that number). ",
            "Run clean_data() to renumber."))
        } else if(is.unsorted(sid)){
          errors <- c(errors, paste0(
            "'diet_data' must be sorted by 'stomach_id' (the TMB diet likelihood ",
            "scans forward for each stomach in turn, so rows out of that order ",
            "silently drop from the likelihood). Run clean_data() to reorder."))
        }
      }
    }

    # Age coverage under empirical (MSVPA) suitability.
    .check_diet_age_coverage(data_list, dd)
  }
  # diet_data presence required when msmMode > 0 (declarative requirement table).
  errors <- c(errors, .rce_check_presence(data_list, "diet_data"))

  # Diet composition likelihood weights & types (only enforced when suitMode > 0)
  if(any(is.na(data_list$Diet_comp_weights) & data_list$suitMode > 0)){
    errors <- c(errors, "Diet composition likelihood weight for a species with estimated suitability is NA")
  }
  if(any(!(data_list$Diet_distribution %in% c(0, 1)) & data_list$suitMode > 0)){
    errors <- c(errors, "Diet composition likelihood for a species with estimated suitability is not 0 or 1")
  }

  # other_food = 0 causes divide-by-zero in suitability
  if(!is.null(data_list$msmMode) && data_list$msmMode > 0 &&
     !is.null(data_list$other_food) && any(data_list$other_food == 0, na.rm = TRUE)){
    errors <- c(errors, "msmMode > 0 requires other_food > 0 for all species; zero values cause divide-by-zero in suitability")
  }

  # Bioenergetics scalars: required (length nspp) when msmMode > 0. switch_check
  # fills these with safe sentinels in single-species mode, so the only way
  # they're missing or wrong-length here is in multispecies mode. The 12 scalars
  # form one grouped requirement (declarative requirement table).
  errors <- c(errors, .rce_check_presence(data_list, "bioenergetics"))

  # =======================================================================
  # 7. Environmental data ----
  # =======================================================================

  # env_data may have just a Year column (no indices) -- clean_data fills that
  # default. Downstream checks (Cindex when Ceq > 1, env-q catchability,
  # srr_indices, M1_indices) report when an index is actually needed but
  # missing.
  if(has_data(data_list$env_data)){
    if(!"Year" %in% colnames(data_list$env_data)){
      errors <- c(errors, "env_data is missing required 'Year' column")
    }
  }
  if(any(data_list$M1_indices > ncol(data_list$env_index))){
    errors <- c(errors, "'M1_indices' greater than the number of indices included")
  }

  # =======================================================================
  # 8. Switches (delegated to validate_switches) ----
  # =======================================================================

  errors <- c(errors, validate_switches(data_list))

  if(length(errors) > 0){
    stop(paste(errors, collapse = "\n"))
  }
}


#' Warn about ages that empirical suitability cannot see
#'
#' Under `suitMode = 0` suitability comes straight from `diet_data`, so an age
#' with no diet row is switched off rather than estimated: a predator age with
#' no rows gets `suit_other = 1` from `calculate_msvpa_suitability()`
#' (`predation.hpp`) and exerts no predation, and a prey age with no rows is
#' never eaten. Neither raises an error or moves the likelihood.
#'
#' Only prey-at-age-in-predator-at-age rows count. `organize_diet_obs()`
#' (`diet_data.hpp`) forms `Pred_age - minage(rsp)` / `Prey_age - minage(ksp)`
#' and skips the row when either is negative, so the aggregated diet formats,
#' which sit below `minage` (see the `diet_data` entry in
#' `.rce_column_schema()`), never reach the suitability array and cannot close
#' a gap.
#'
#' Coverage ignores `Year`: a diet sampled in only some years is a design, not
#' a gap.
#'
#' @param data_list An Rceattle data list, post-`switch_check()`.
#' @param dd Its `diet_data` table.
#' @return `NULL`, invisibly. Called for the warnings.
#' @noRd
.check_diet_age_coverage <- function(data_list, dd){

  if(!isTRUE(any(data_list$msmMode > 0))) return(invisible(NULL))
  if(is.null(data_list$suitMode) || is.null(data_list$minage) ||
     is.null(data_list$nages)   || is.null(data_list$nsex)) return(invisible(NULL))
  if(!all(c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age") %in%
          colnames(dd))) return(invisible(NULL))

  nspp    <- data_list$nspp
  suit    <- rep(data_list$suitMode, length.out = nspp)
  minage  <- data_list$minage
  nages   <- data_list$nages
  nsex    <- data_list$nsex
  spnames <- data_list$spnames
  if(is.null(spnames) || length(spnames) != nspp) spnames <- paste("Species", seq_len(nspp))

  empirical <- which(suit == 0)
  if(length(empirical) == 0) return(invisible(NULL))

  # Rows the suitability array reads: both species in range, both ages at or
  # above minage.
  in_range <- dd$Pred %in% seq_len(nspp) & dd$Prey %in% seq_len(nspp)
  dd <- dd[in_range, , drop = FALSE]
  if(nrow(dd) == 0) return(invisible(NULL))
  typed <- dd[dd$Pred_age >= minage[dd$Pred] & dd$Prey_age >= minage[dd$Prey], ,
              drop = FALSE]
  # Only an empirical-suitability predator's rows build suitability, so only
  # they can cover a prey age.
  typed <- typed[typed$Pred %in% empirical, , drop = FALSE]

  # "11, 12, 13, 20" -> "11-13, 20"
  runs <- function(x){
    if(length(x) == 0) return("")
    brk <- c(0, which(diff(x) != 1), length(x))
    paste(vapply(seq_len(length(brk) - 1), function(i){
      seg <- x[(brk[i] + 1):brk[i + 1]]
      if(length(seg) == 1) as.character(seg) else paste0(seg[1], "-", seg[length(seg)])
    }, character(1)), collapse = ", ")
  }

  # A sex-combined row (sex 0) covers both sexes; a sexed row covers only its
  # own. `organize_diet_obs()` reads the sex column only for a two-sex species,
  # so for a one-sex species every row counts whatever it says.
  covered <- function(tab, sp, sex, nsex_sp, age_col, sex_col, sp_col){
    hit <- tab[[sp_col]] == sp
    if(nsex_sp == 2){
      hit <- hit & (tab[[sex_col]] == 0 | tab[[sex_col]] == sex)
    }
    unique(tab[[age_col]][hit])
  }

  # Predator coverage is owed only by a species that derives its own
  # suitability from diet data; prey coverage is owed by every species that is
  # eaten at all, since any of them can be eaten by an empirical predator.
  # `suitMode` describes how a species feeds, not how it is fed on.
  gaps <- character(0)
  for(sp in seq_len(nspp)){
    # `nages` counts age bins; the ages run minage .. minage+nages-1, as
    # everywhere else in the package.
    ages <- seq_len(nages[sp]) - 1L + minage[sp]
    for(sex in seq_len(nsex[sp])){
      sex_lab <- if(nsex[sp] == 1) "" else paste0(" (", c("female", "male")[sex], ")")

      # A species with no rows in a role at all is reported as such rather
      # than as "ages 1-N".
      miss_pred <- if(sp %in% empirical){
        setdiff(ages, covered(typed, sp, sex, nsex[sp],
                              "Pred_age", "Pred_sex", "Pred"))
      } else integer(0)
      if(length(miss_pred)){
        gaps <- c(gaps, paste0(
          "  ", spnames[sp], sex_lab, " as PREDATOR: ",
          if(length(miss_pred) == length(ages)){
            "no diet data at any age -- it exerts no predation."
          } else {
            paste0("no diet data at age ", runs(sort(miss_pred)),
                   " -- those ages exert no predation.")
          }))
      }

      # A species with NO prey rows at any age is not a truncated diet table --
      # it is a species nothing in the model eats, which is a modelling choice
      # and a common one (an apex predator in a two-species run). Warning about
      # it fires on every fit of a correctly specified model and says nothing
      # the author did not intend. Only a PARTIAL gap is evidence of truncation,
      # so that is what is reported.
      #
      # The predator role keeps its all-ages case: a species that asked for
      # empirical suitability and supplied no diet data at all did not choose
      # to exert no predation, it just gets that.
      prey_seen <- covered(typed, sp, sex, nsex[sp],
                           "Prey_age", "Prey_sex", "Prey")
      miss_prey <- if(length(prey_seen)) setdiff(ages, prey_seen) else integer(0)
      if(length(miss_prey)){
        gaps <- c(gaps, paste0(
          "  ", spnames[sp], sex_lab, " as PREY: no diet data at age ",
          runs(sort(miss_prey)), " -- those ages are never eaten."))
      }
    }
  }

  if(length(gaps)){
    warning(paste0(
      "'diet_data' does not cover every age of a species using empirical ",
      "suitability (suitMode = 0), and empirical suitability is read straight ",
      "from the diet data, so an uncovered age is switched out of the ",
      "predation calculation rather than estimated:\n",
      paste(gaps, collapse = "\n"),
      "\nSupply diet rows for those ages, or pool them into the plus group. ",
      "Rows with an age below the species' 'minage' are the aggregated diet ",
      "formats and do not count here -- only prey-at-age-in-predator-at-age ",
      "rows build MSVPA suitability."), call. = FALSE)
  }

  invisible(NULL)
}


#' Is a symmetric matrix positive definite?
#'
#' Tested by attempting a Cholesky factorization rather than by inspecting
#' eigenvalues, because that is the operation the covariance index likelihood
#' actually performs (`MVNORM()` / `Eigen::LLT` in `ceattle.cpp`): a matrix that
#' factorizes here is one TMB can use. Assumes symmetry has already been
#' checked -- `chol()` reads only the upper triangle, so it would accept an
#' asymmetric matrix whose upper triangle happens to be positive definite.
#'
#' @param x A numeric matrix, assumed square and symmetric.
#' @return `TRUE` if `x` admits a Cholesky factorization, otherwise `FALSE`.
#' @noRd
.is_positive_definite <- function(x) {
  if(!is.matrix(x) || nrow(x) != ncol(x) || nrow(x) == 0) return(FALSE)
  if(anyNA(x) || any(!is.finite(x))) return(FALSE)
  isTRUE(tryCatch({
    chol(x)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE))
}


#' Warn when a q linkage breaks a shared `Catchability_index` group
#'
#' `adjust_map_shared_params()` ties the group's `index_log_q` to the lead
#' fleet's, but `ceattle.cpp` adds `q_linkage_offset(flt, yr)` per fleet
#' afterwards, and nothing reconciles that. A linkage naming only some of a
#' group's fleets therefore gives them different catchabilities while the
#' `fleet_control` still says they share one. Measured on `BS2017SS` with fleets
#' 4 and 7 in one group and the linkage on fleet 7: fleet 4 flat at 0.035, fleet
#' 7 running 0.087-0.537. With every fleet in the group named, and equal
#' coefficients, they stay together -- so this fires on a strict subset only.
#'
#' Separate from `data_check()` because the linkage table does not exist yet
#' when that runs: `fit_mod()` pools it after the check.
#'
#' @param data_list A `data_list` carrying `linkage_table` and `fleet_control`.
#' @keywords internal
#' @noRd
.warn_q_linkage_shared_group <- function(data_list) {
  lt <- data_list$linkage_table
  fc <- data_list$fleet_control
  if (is.null(lt) || is.null(nrow(lt)) || !nrow(lt)) return(invisible())
  if (is.null(fc) || is.null(fc$Catchability_index)) return(invisible())

  qrows <- lt[as.character(lt$process) %in% c("q", "4"), , drop = FALSE]
  if (!nrow(qrows)) return(invisible())
  # `fleet` on the linkage table is a Fleet_code (linkage_spec() documents it as
  # a 1-based Fleet_code), so compare against Fleet_code, not the row number.
  linked <- unique(stats::na.omit(as.integer(qrows$fleet)))
  if (!length(linked)) return(invisible())

  live <- !(fc$Fleet_type %in% c("Off", 0, "0"))
  for (qi in unique(fc$Catchability_index[!is.na(fc$Catchability_index) & live])) {
    rows <- which(!is.na(fc$Catchability_index) & fc$Catchability_index == qi & live)
    if (length(rows) < 2) next
    has <- fc$Fleet_code[rows] %in% linked
    if (any(has) && !all(has)) {
      warning(paste0(
        "A catchability linkage names ", paste(fc$Fleet_name[rows][has], collapse = ", "),
        " but not ", paste(fc$Fleet_name[rows][!has], collapse = ", "),
        ", which share Catchability_index ", qi,
        ". The linkage offset is added per fleet and is not shared, so the group ",
        "will not have one catchability. Name every fleet in the group, or give ",
        "the linked fleet its own Catchability_index."))
    }
  }
  invisible()
}
