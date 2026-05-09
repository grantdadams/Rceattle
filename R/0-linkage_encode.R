#' Encode an Rceattle linkage table into TMB-friendly inputs
#'
#' The TMB template consumes the pooled linkage table as a set of
#' parallel `IVECTOR` / `VECTOR` inputs plus a dense design matrix.
#' R encodes string-valued columns (`process`, `param`, `link`,
#' `prior_family`) as 0-based integer codes (TMB-friendly); converts
#' `NA` stratum ids (`species`, `sex`, `age_bin`) to a sentinel `0`
#' meaning "applies to all levels" (TMB-side dispatch expands the
#' shared rows over the relevant 1-based levels); converts the 1-based
#' `X_col` to 0-based.
#'
#' This is mechanism only; no model behavior depends on the encoding
#' until the TMB template wires up the corresponding inputs.
#'
#' @keywords internal
#' @name linkage_encode
NULL


#' Integer codes for the `process` column
#'
#' Stable across versions. New processes append at the end so existing
#' codes never shift.
#'
#' @keywords internal
LINKAGE_PROCESS_CODES <- c(
  recruitment = 0L,
  M           = 1L,
  growth      = 2L,
  q           = 3L,
  sel         = 4L
)


#' Integer codes for the `link` column
#' @keywords internal
LINKAGE_LINK_CODES <- c(
  identity = 0L,
  log      = 1L,
  logit    = 2L
)


#' Integer codes for the `prior_family` column
#' @keywords internal
LINKAGE_PRIOR_CODES <- c(
  none      = 0L,
  normal    = 1L,
  lognormal = 2L,
  gamma     = 3L,
  beta      = 4L
)


#' Integer codes for the `param` column, per process
#'
#' Each process has its own independent parameter namespace (`log_K`
#' means something different to growth than `log_alpha` does to
#' recruitment), so the param column is encoded against a per-process
#' table. Unknown parameter names error.
#'
#' @keywords internal
LINKAGE_PARAM_CODES <- list(
  growth      = c(log_K = 0L, log_L1 = 1L, log_Linf = 2L, log_m = 3L),
  M           = c(log_M1 = 0L),
  recruitment = c(log_R0 = 0L, log_alpha = 1L, log_beta = 2L),
  q           = c(log_q = 0L),
  sel         = character(0)   # not yet wired
)


#' Sentinel meaning "applies to every level of this stratum"
#'
#' Stored where the user passed `NA` for `species`, `sex`, or
#' `age_bin`. The TMB template expands these by iterating over
#' `1:nspp` / `1:nsex` / `1:nages` rather than indexing a single cell.
#'
#' @keywords internal
LINKAGE_STRATUM_ALL <- 0L


#' Encode a linkage table into TMB-friendly parallel vectors
#'
#' @param table an `Rceattle_linkage_table` (typically the output of
#'   [pool_linkages()]).
#' @param X the global design matrix from [pool_linkages()] (passed
#'   through unchanged so callers can stash it alongside).
#'
#' @return A named list with components:
#' \describe{
#'   \item{n_linkage}{`integer(1)`. Number of coefficients
#'     (`nrow(table)`).}
#'   \item{linkage_process, linkage_param, linkage_species,
#'     linkage_sex, linkage_age_bin, linkage_X_col, linkage_link,
#'     linkage_prior_family}{Parallel `integer(n_linkage)` vectors.}
#'   \item{linkage_prior_p1, linkage_prior_p2}{Parallel
#'     `numeric(n_linkage)` vectors.}
#'   \item{linkage_X}{The design matrix as passed in.}
#'   \item{linkage_design_names}{`character` colnames of `linkage_X`.}
#' }
#' @keywords internal
encode_linkage_for_tmb <- function(table, X) {
  if (!is_linkage_table(table)) {
    stop("`table` must be an Rceattle_linkage_table", call. = FALSE)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix", call. = FALSE)
  }
  n <- nrow(table)
  if (n == 0L) {
    return(.empty_tmb_encoding(X))
  }

  # process: known-good (validated upstream); fail loudly on surprises
  proc_int <- unname(LINKAGE_PROCESS_CODES[table$process])
  if (anyNA(proc_int)) {
    stop("encode_linkage_for_tmb: unknown process(es) in table",
         call. = FALSE)
  }

  # param: encoded against per-process namespace
  param_int <- vapply(seq_len(n), function(i) {
    codes <- LINKAGE_PARAM_CODES[[table$process[i]]]
    code  <- codes[table$param[i]]
    if (is.na(code)) {
      stop(sprintf(
        "encode_linkage_for_tmb: unknown param '%s' for process '%s'",
        table$param[i], table$process[i]), call. = FALSE)
    }
    unname(code)
  }, integer(1))

  link_int  <- unname(LINKAGE_LINK_CODES[table$link])
  prior_int <- unname(LINKAGE_PRIOR_CODES[table$prior_family])

  # NA stratum ids => sentinel 0 ("applies to all"); else 1-based
  to_stratum <- function(v) {
    out <- as.integer(v)
    out[is.na(out)] <- LINKAGE_STRATUM_ALL
    out
  }

  list(
    n_linkage             = n,
    linkage_process       = proc_int,
    linkage_param         = param_int,
    linkage_species       = to_stratum(table$species),
    linkage_sex           = to_stratum(table$sex),
    linkage_age_bin       = to_stratum(table$age_bin),
    linkage_X_col         = as.integer(table$X_col) - 1L,   # 0-based for TMB
    linkage_link          = link_int,
    linkage_is_intercept  = as.integer(table$design_col == "(Intercept)"),
    linkage_prior_family  = prior_int,
    linkage_prior_p1      = as.numeric(table$prior_p1),
    linkage_prior_p2      = as.numeric(table$prior_p2),
    linkage_X             = X,
    linkage_design_names  = colnames(X) %||% character(ncol(X))
  )
}


#' @keywords internal
#' @noRd
.empty_tmb_encoding <- function(X) {
  list(
    n_linkage             = 0L,
    linkage_process       = integer(0),
    linkage_param         = integer(0),
    linkage_species       = integer(0),
    linkage_sex           = integer(0),
    linkage_age_bin       = integer(0),
    linkage_X_col         = integer(0),
    linkage_link          = integer(0),
    linkage_is_intercept  = integer(0),
    linkage_prior_family  = integer(0),
    linkage_prior_p1      = numeric(0),
    linkage_prior_p2      = numeric(0),
    linkage_X             = X,
    linkage_design_names  = colnames(X) %||% character(ncol(X))
  )
}
