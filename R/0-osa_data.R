#' Build the flat observation vector and metadata for OSA residuals
#'
#' @description
#' Internal helper used by [rearrange_data()] to assemble the inputs the TMB
#' template needs to compute one-step-ahead (OSA) residuals via
#' [TMB::oneStepPredict()].
#'
#' OSA residuals require every observation that enters the likelihood to live in
#' a single flat vector (`obsvec`) with a companion `keep` indicator. This
#' function walks the already-built control/observation matrices (`index_ctl` /
#' `index_obs`, `catch_ctl` / `catch_obs`, ...) and produces:
#'
#' \itemize{
#'   \item `obsvec` -- the flat vector of observations. Aggregate (catch and
#'     index) observations are stored on the log scale because their likelihood
#'     is lognormal, matching what the TMB template reads.
#'   \item `obs_ctl` -- a data frame with one row per element of `obsvec`,
#'     mapping each position back to its data type, fleet, species, year, etc.,
#'     so residuals stay interpretable. This is R-side metadata only and is
#'     removed before the data list is passed to TMB.
#'   \item `*_obsvec_idx` -- for each observation type, an integer vector giving
#'     the 0-based position of each `*_obs` row within `obsvec`, or `-1` when the
#'     row is excluded from the likelihood (projection years, non-positive
#'     observations). The TMB template reads these to map each observation to its
#'     `obsvec` element.
#' }
#'
#' The inclusion rules below must match the guards in the TMB template exactly,
#' so that every observation the template evaluates has a valid `obsvec`
#' position. Phase 1 covers the aggregate catch and index series; composition,
#' CAAL, and diet observations are added in later development phases.
#'
#' @param data_list A `data_list` whose `*_ctl`/`*_obs` matrices have already
#'   been built by [rearrange_data()].
#'
#' @return The input `data_list` with `obsvec`, `obs_ctl`, `osa_mode`, and the
#'   per-type `*_obsvec_idx` vectors added.
#'
#' @keywords internal
build_osa_data <- function(data_list) {

  endyr    <- data_list$endyr
  flt_type <- data_list$flt_type   # indexed by Fleet_code, matching the TMB cpp

  obsvec   <- numeric(0)
  obs_ctl  <- .new_obs_ctl()
  next_pos <- 0L                   # 0-based position counter (TMB convention)

  # Append a block of observations to `obsvec`/`obs_ctl` and return their
  # 0-based positions. `value` is what the likelihood reads (e.g. log(obs) for
  # the lognormal aggregate series). Uses `<<-` to update the running state in
  # the enclosing function environment.
  append_obs <- function(value, type, data_row, fleet_code, species, year,
                         sex = NA_integer_, age_or_length = NA_integer_,
                         bin_index = NA_integer_, comp_type = NA_integer_,
                         is_last_bin = FALSE, stomach_id = NA_integer_,
                         group_id = NA_integer_) {
    n <- length(value)
    if (n == 0) return(integer(0))
    pos <- next_pos + seq_len(n) - 1L
    if (all(is.na(group_id))) group_id <- pos    # default: one group per obs
    obsvec  <<- c(obsvec, as.numeric(value))
    obs_ctl <<- rbind(obs_ctl, data.frame(
      obs_pos       = as.integer(pos),
      type          = type,
      data_row      = as.integer(data_row),
      fleet_code    = as.integer(fleet_code),
      species       = as.integer(species),
      sex           = as.integer(sex),
      year          = as.integer(year),
      age_or_length = as.integer(age_or_length),
      bin_index     = as.integer(bin_index),
      comp_type     = as.integer(comp_type),
      is_last_bin   = as.logical(is_last_bin),
      stomach_id    = as.integer(stomach_id),
      group_id      = as.integer(group_id),
      stringsAsFactors = FALSE))
    next_pos <<- next_pos + n
    pos
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
        value      = log(index_obs[inc, 1]),
        type       = "index",
        data_row   = inc,
        fleet_code = index_ctl[inc, 1],
        species    = index_ctl[inc, 2],
        year       = index_ctl[inc, 3])
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
        value      = log(catch_obs[inc, 1]),
        type       = "catch",
        data_row   = inc,
        fleet_code = catch_ctl[inc, 1],
        species    = catch_ctl[inc, 2],
        year       = catch_ctl[inc, 3])
    }
  }

  # TMB needs a non-empty DATA_VECTOR; use a length-1 sentinel when no
  # observations are included (the *_obsvec_idx vectors are then all -1 and the
  # sentinel value is never read by the template).
  if (length(obsvec) == 0) obsvec <- 0

  data_list$obsvec           <- as.numeric(obsvec)
  data_list$obs_ctl          <- obs_ctl
  data_list$index_obsvec_idx <- as.integer(index_obsvec_idx)
  data_list$catch_obsvec_idx <- as.integer(catch_obsvec_idx)
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
    obs_pos       = integer(0),
    type          = character(0),
    data_row      = integer(0),
    fleet_code    = integer(0),
    species       = integer(0),
    sex           = integer(0),
    year          = integer(0),
    age_or_length = integer(0),
    bin_index     = integer(0),
    comp_type     = integer(0),
    is_last_bin   = logical(0),
    stomach_id    = integer(0),
    group_id      = integer(0),
    stringsAsFactors = FALSE)
}
