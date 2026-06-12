#' Build the flat observation vector and metadata for OSA residuals
#'
#' @description
#' Internal helper used by [rearrange_data()] to assemble the inputs the TMB
#' template needs to compute one-step-ahead (OSA) residuals via
#' [TMB::oneStepPredict()].
#'
#' OSA residuals require every observation that enters the likelihood to live in
#' a single flat vector (`obsvec`) with a companion `keep` indicator. This
#' function walks the already-built control/observation matrices and produces:
#'
#' \itemize{
#'   \item `obsvec` -- the flat vector of observations. Aggregate (catch and
#'     index) observations are stored as log(observation) because their
#'     likelihood is lognormal; composition (comp) and conditional-age-at-length
#'     (caal) observations are stored as bin counts, `(proportion + 1e-5) * N`,
#'     matching exactly what the TMB likelihood forms during fitting.
#'   \item `obs_ctl` -- a data frame with one row per element of `obsvec`,
#'     mapping each position back to its data type, fleet, species, year,
#'     age/length bin, etc., so residuals stay interpretable. R-side metadata
#'     only; removed before the data list is passed to TMB.
#'   \item `*_obsvec_idx` -- per observation-row index vectors. For aggregate
#'     series this is the 0-based `obsvec` position of each row; for composition
#'     and caal it is the 0-based `obsvec` position of the row's FIRST bin (the
#'     template reads the row's bins as `obsvec.segment(start, n_bins)`). `-1`
#'     marks rows excluded from the likelihood.
#' }
#'
#' The inclusion rules and per-row bin counts below must match the guards in the
#' TMB template exactly, so that every observation the template evaluates has a
#' valid `obsvec` position. All fitted observation types are handled: aggregate
#' catch and index, comp and caal composition, and predator diet composition.
#'
#' @param data_list A `data_list` whose `*_ctl`/`*_obs`/`*_n` matrices have
#'   already been built (and coerced to matrices) by [rearrange_data()].
#' @param build_osa Logical. When `TRUE`, build the full OSA observation data
#'   for every type (aggregate, composition, caal, and diet) so
#'   [osa_residuals()] can be computed. When `FALSE` (the default), only the
#'   aggregate index/catch entries the TMB template always reads are built and
#'   the (much larger) composition/caal/diet metadata is skipped -- this is the
#'   fast path for simulation testing (e.g. [run_mse()]), where the fitted
#'   objective is identical but OSA composition residuals are not produced.
#'
#' @return The input `data_list` with `obsvec`, `obs_ctl`, `osa_mode`, and the
#'   per-type `*_obsvec_idx` vectors added.
#'
#' @keywords internal
build_osa_data <- function(data_list, build_osa = FALSE) {

  endyr    <- data_list$endyr
  flt_type <- data_list$flt_type   # indexed by Fleet_code, matching the TMB cpp
  nages    <- data_list$nages
  nlengths <- data_list$nlengths

  obsvec_parts <- list()           # numeric blocks; concatenated once at the end
  ctl_parts    <- list()           # obs_ctl row blocks; rbind'd once at the end
  next_pos     <- 0L               # 0-based position counter (TMB convention)

  # Append a block of observations to `obsvec`/`obs_ctl` and return their
  # 0-based positions. `value` is what the likelihood reads. When
  # `one_group = TRUE` all elements share a single group id (one decomposition
  # unit, e.g. the bins of one composition); otherwise each element is its own
  # group (e.g. one aggregate observation). Blocks are accumulated via `<<-` and
  # joined once after the loops (see the assembly step at the end).
  append_obs <- function(value, source, data_row, fleet_code, species, year,
                         sex = NA_integer_, age_length_bin = NA_integer_,
                         length = NA_integer_, bin_index = NA_integer_,
                         comp_type = NA_integer_, is_last_bin = FALSE,
                         stomach_id = NA_integer_, one_group = FALSE) {
    n <- length(value)
    if (n == 0) return(integer(0))
    pos      <- next_pos + seq_len(n) - 1L
    group_id <- if (one_group) rep(pos[1], n) else pos
    k <- length(obsvec_parts) + 1L
    obsvec_parts[[k]] <<- as.numeric(value)
    ctl_parts[[k]]    <<- data.frame(
      obs_pos        = as.integer(pos),
      source         = source,
      data_row       = as.integer(data_row),
      fleet_code     = as.integer(fleet_code),
      species        = as.integer(species),
      sex            = as.integer(sex),
      year           = as.integer(year),
      age_length_bin = as.integer(age_length_bin),
      length         = as.integer(length),
      bin_index     = as.integer(bin_index),
      comp_type     = as.integer(comp_type),
      is_last_bin   = as.logical(is_last_bin),
      stomach_id    = as.integer(stomach_id),
      group_id      = as.integer(group_id),
      stringsAsFactors = FALSE)
    next_pos <<- next_pos + n
    pos
  }

  # Append one composition row (comp or caal): its bin counts =
  # (proportion + 1e-5) * Neff, one decomposition group, final bin flagged.
  # Returns the obsvec start position of the row's first bin.
  append_composition <- function(source, obs_row, n_bins, Neff, fleet, sp, sex,
                                 yr, data_row, comp_type = NA_integer_,
                                 length = NA_integer_) {
    counts <- (as.numeric(obs_row[seq_len(n_bins)]) + 0.00001) * Neff
    append_obs(value = counts, source = source, data_row = data_row,
               fleet_code = fleet, species = sp, sex = sex, year = yr,
               age_length_bin = seq_len(n_bins), length = length,
               bin_index = seq_len(n_bins) - 1L, comp_type = comp_type,
               is_last_bin = seq_len(n_bins) == n_bins, one_group = TRUE)[1]
  }

  # ---- Index (survey) observations: lognormal, stored as log(obs) ----
  # TMB guard: Year in (0, endyr], fleet on (flt_type > 0), observation > 0.
  index_ctl <- data_list$index_ctl
  index_obs <- data_list$index_obs
  index_obsvec_idx <- rep(-1L, nrow(index_obs))
  if (nrow(index_obs) > 0) {
    inc <- which(index_ctl[, 3] > 0 & index_ctl[, 3] <= endyr &
                   flt_type[index_ctl[, 1]] > 0 & index_obs[, 1] > 0)
    if (length(inc) > 0) {
      index_obsvec_idx[inc] <- append_obs(
        value = log(index_obs[inc, 1]), source = "index", data_row = inc,
        fleet_code = index_ctl[inc, 1], species = index_ctl[inc, 2],
        year = index_ctl[inc, 3])
    }
  }

  # ---- Catch (fishery) observations: lognormal, stored as log(obs) ----
  # TMB guard: Year in (0, endyr], fishery fleet (flt_type == 1), catch > 0.
  catch_ctl <- data_list$catch_ctl
  catch_obs <- data_list$catch_obs
  catch_obsvec_idx <- rep(-1L, nrow(catch_obs))
  if (nrow(catch_obs) > 0) {
    inc <- which(catch_ctl[, 3] > 0 & catch_ctl[, 3] <= endyr &
                   flt_type[catch_ctl[, 1]] == 1 & catch_obs[, 1] > 0)
    if (length(inc) > 0) {
      catch_obsvec_idx[inc] <- append_obs(
        value = log(catch_obs[inc, 1]), source = "catch", data_row = inc,
        fleet_code = catch_ctl[inc, 1], species = catch_ctl[inc, 2],
        year = catch_ctl[inc, 3])
    }
  }

  # ---- Age/length composition: bin counts (proportion + 1e-5) * N ----
  # TMB guard: Year in (0, endyr], fleet on, Sample_size > 0. Bins per row:
  # nages (age comp) or nlengths (length comp), doubled for joint-sex (Sex==3).
  # One decomposition unit per row; the last bin is fixed by the sum-to-N
  # constraint (is_last_bin -> residual NA).
  comp_ctl <- data_list$comp_ctl
  comp_obs <- data_list$comp_obs
  comp_n   <- data_list$comp_n
  comp_obsvec_idx <- rep(-1L, nrow(comp_obs))
  if (build_osa && nrow(comp_obs) > 0) {
    for (r in seq_len(nrow(comp_obs))) {
      fleet     <- comp_ctl[r, 1]; sp <- comp_ctl[r, 2]; sex <- comp_ctl[r, 3]
      comp_type <- comp_ctl[r, 4]; yr <- comp_ctl[r, 5]; Neff <- comp_n[r, 2]
      if (!(yr > 0 && yr <= endyr && flt_type[fleet] > 0 && Neff > 0)) next
      joint_adjust <- if (sex == 3) 2L else 1L          # joint-sex doubles the bins
      n_comp <- (if (comp_type == 0) nages[sp] else nlengths[sp]) * joint_adjust
      comp_obsvec_idx[r] <- append_composition(
        "comp", comp_obs[r, ], n_comp, Neff, fleet, sp, sex, yr, r,
        comp_type = comp_type)
    }
  }

  # ---- Conditional age-at-length (caal): age-bin counts per length bin ----
  # TMB guard: Year in (0, endyr], fleet on, Sample_size > 0. n_caal = nages.
  caal_ctl <- data_list$caal_ctl
  caal_obs <- data_list$caal_obs
  caal_n   <- data_list$caal_n
  caal_obsvec_idx <- rep(-1L, nrow(caal_obs))
  if (build_osa && nrow(caal_obs) > 0) {
    for (r in seq_len(nrow(caal_obs))) {
      fleet <- caal_ctl[r, 1]; sp <- caal_ctl[r, 2]; sex <- caal_ctl[r, 3]
      yr    <- caal_ctl[r, 4]; len_bin <- caal_ctl[r, 5]; Neff <- caal_n[r, 1]
      if (!(yr > 0 && yr <= endyr && flt_type[fleet] > 0 && Neff > 0)) next
      caal_obsvec_idx[r] <- append_composition(
        "caal", caal_obs[r, ], nages[sp], Neff, fleet, sp, sex, yr, r,
        comp_type = 0L, length = len_bin)            # age bins conditioned on length
    }
  }

  # ---- Diet composition: per-stomach prey counts + "other prey" ----
  # TMB guard: the predator's suitability is estimated (suitMode[pred] > 0) and
  # the stomach has >= 1 prey item. The composition for a stomach is its prey
  # items plus an "other prey" category (the last bin, fixed by the sum-to-1
  # constraint and dropped). One decomposition unit per stomach; diet_obsvec_idx
  # is indexed by 0-based stomach id, so it has length n_stomach_obs.
  n_stomach <- data_list$n_stomach_obs
  diet_obsvec_idx <- if (is.null(n_stomach)) integer(0) else rep(-1L, n_stomach)
  if (build_osa && !is.null(n_stomach) && n_stomach > 0 &&
      !is.null(data_list$diet_ctl) && nrow(data_list$diet_ctl) > 0) {
    diet_ctl   <- data_list$diet_ctl
    diet_obs   <- data_list$diet_obs
    stomach_id <- data_list$stomach_id
    suitMode   <- data_list$suitMode
    for (i in seq_len(n_stomach) - 1L) {        # stomach ids are 0-based
      rows   <- which(stomach_id == i)
      n_prey <- length(rows)
      if (n_prey == 0) next
      rsp <- diet_ctl[rows[1], 1]               # predator species (1-based)
      if (is.null(suitMode) || suitMode[rsp] <= 0) next
      N_s  <- diet_obs[rows[1], 1]              # sample size (number of stomachs)
      # observed prey proportions plus the residual "other prey" category, then
      # the same offset/normalize the TMB likelihood applies, scaled to counts.
      obs_p  <- as.numeric(diet_obs[rows, 2])
      vec    <- c(obs_p, 1 - min(sum(obs_p), 1)) + 0.00001
      counts <- vec / sum(vec) * N_s
      diet_obsvec_idx[i + 1L] <- append_obs(
        value          = counts, source = "diet",
        data_row       = c(rows, NA_integer_),
        fleet_code     = rsp, species = rsp,
        year           = diet_ctl[rows[1], 7],
        age_length_bin = c(diet_ctl[rows, 2], NA_integer_),   # prey id; NA = other
        bin_index      = seq_len(n_prey + 1L) - 1L,
        is_last_bin   = seq_len(n_prey + 1L) == (n_prey + 1L),
        stomach_id    = i, one_group = TRUE)[1]
    }
  }

  # Assemble the flat vector and metadata in a single pass from the blocks
  # accumulated above, instead of growing them with c() / rbind() inside the
  # loops (which is quadratic and runs on every fit_mod() call).
  obsvec  <- if (length(obsvec_parts)) unlist(obsvec_parts, use.names = FALSE) else numeric(0)
  obs_ctl <- if (length(ctl_parts)) do.call(rbind, ctl_parts) else .new_obs_ctl()
  rownames(obs_ctl) <- NULL

  # TMB needs a non-empty DATA_VECTOR; use a length-1 sentinel when no
  # observations are included (the *_obsvec_idx vectors are then all -1 and the
  # sentinel value is never read by the template).
  if (length(obsvec) == 0) obsvec <- 0

  data_list$obsvec           <- as.numeric(obsvec)
  data_list$obs_ctl          <- obs_ctl
  data_list$index_obsvec_idx <- as.integer(index_obsvec_idx)
  data_list$catch_obsvec_idx <- as.integer(catch_obsvec_idx)
  data_list$comp_obsvec_idx  <- as.integer(comp_obsvec_idx)
  data_list$caal_obsvec_idx  <- as.integer(caal_obsvec_idx)
  data_list$diet_obsvec_idx  <- as.integer(diet_obsvec_idx)
  data_list$osa_mode         <- 0L

  data_list
}


#' Empty `obs_ctl` metadata frame with the correct column types
#'
#' @return A 0-row data frame with the columns used to map each `obsvec` element
#'   back to its observation. See [build_osa_data()].
#' @keywords internal
.new_obs_ctl <- function() {
  data.frame(
    obs_pos        = integer(0),
    source         = character(0),
    data_row       = integer(0),
    fleet_code     = integer(0),
    species        = integer(0),
    sex            = integer(0),
    year           = integer(0),
    age_length_bin = integer(0),
    length         = integer(0),
    bin_index     = integer(0),
    comp_type     = integer(0),
    is_last_bin   = logical(0),
    stomach_id    = integer(0),
    group_id      = integer(0),
    stringsAsFactors = FALSE)
}
