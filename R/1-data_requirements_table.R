# =============================================================================
# Declarative data-requirement table
# =============================================================================
#
# `data_check()` historically hard-coded, at scattered sites, the *conditional
# presence requirements* of the ~40-element Rceattle data list: e.g. "diet_data
# is required when msmMode > 0", "NByageFixed is required when estDynamics > 0".
# Those conditions are exactly what a user needs to know up front, and what
# `data_requirements()` reports and `build_data()` pre-checks. Rather than let
# the same condition live in three places (drifting apart over time), it lives
# once here, and `data_check()` *consumes* this table.
#
# Scope, deliberately narrow (see dev/PLAN-data-workflow-and-linkage-grammar.md,
# PR 4 step 1): this table drives only the **pure presence-requirement** gates.
# Dimension / value / referential / structural checks and the two
# mirroring-dependent *adequacy* gates (index_cov-MVN, comp/caal-vs-estimated-
# selectivity) stay imperative in `data_check()` -- they depend on per-fleet row
# counts and cross-fleet lookups that do not reduce to a declarative row. Rows
# marked `driven = FALSE` are introspection-only (read by `data_requirements()`)
# and are NOT consulted by `data_check()`; the authoritative check for those
# stays in `data_check()`.
#
# Each row is a list with fields:
#   element        chr  data_list element name (the row's identity / lookup key).
#   category       chr  grouping for display ("dimensions", "biology", "fishery",
#                       "composition", "predation", "environment").
#   required_when  fn   function(dl) -> logical; TRUE when the element is
#                       required for this configuration. Evaluated with isTRUE()
#                       by the callers, so a length-0 / NA result reads as "not
#                       required" (robust to a not-yet-defaulted switch).
#   condition_label chr human-readable form of `required_when`, for reports.
#   driven         lgl  TRUE if `data_check()` emits this row's requirement via
#                       the evaluators below; FALSE = introspection-only.
#   severity       chr  ("error" | "message") how data_check surfaces an unmet
#                       requirement. Only meaningful when driven = TRUE.
#   adequate       fn   function(dl) -> logical; TRUE when the element is present
#                       and adequate enough to satisfy the *presence* gate (the
#                       separate dimension/value checks stay in data_check).
#                       Only required when driven = TRUE.
#   message        fn   function(dl) -> chr; the exact text data_check emits when
#                       required_when && !adequate. Only required when driven.
#   optional_status chr when NOT required: "defaulted" (clean_data fills a safe
#                       empty default) or "none". Used by data_requirements().
#   default_label  chr  human description of the clean_data default, for reports.
#
# NOTE: `required_when` reproduces each gate's original guard. Where the original
# guard would error on a NULL switch (e.g. the diet gate's bare `msmMode > 0`),
# the isTRUE() wrapping here reads NULL as "not required" instead -- which is
# unreachable in the real pipeline (switch_check defaults msmMode to 0 before
# data_check runs) and strictly more robust. Flagged for review.

# Canonical presence predicate: a data.frame element counts as "present" only
# when it is non-NULL and has at least one row. Mirrors the `has_data()` helper
# local to data_check(); defined here so the table and data_check share one
# definition.
.rce_has_data <- function(df) !is.null(df) && nrow(df) > 0

# The 12 bioenergetics scalars are required together (one gate, one message
# listing whichever are missing), so they share a single group row keyed on
# "bioenergetics".
.RCE_BIOENERGETICS_SCALARS <- c(
  "Ceq", "Cindex", "Pvalue", "fday", "CA", "CB",
  "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4"
)

#' Conditional data-requirement rows consumed by data_check()/data_requirements()
#'
#' @return A named list of requirement rows (see file header for the schema),
#'   keyed by the row `element`.
#' @keywords internal
#' @noRd
.rce_requirement_table <- function() {
  rows <- list(

    # ---- Predation / diet (msmMode > 0) -----------------------------------
    list(
      element        = "diet_data",
      category       = "predation",
      required_when  = function(dl) dl$msmMode > 0,
      condition_label = "msmMode > 0 (multispecies)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$diet_data),
      message        = function(dl) "No diet data included",
      optional_status = "defaulted",
      default_label  = "empty diet_data (clean_data)"
    ),

    list(
      element        = "ration_data",
      category       = "predation",
      required_when  = function(dl) dl$msmMode > 0,
      condition_label = "msmMode > 0 (multispecies)",
      driven         = TRUE,
      severity       = "message",
      adequate       = function(dl) .rce_has_data(dl$ration_data),
      message        = function(dl) "No ration data",
      optional_status = "defaulted",
      default_label  = "empty ration_data (clean_data)"
    ),

    list(
      element        = "bioenergetics",
      category       = "predation",
      required_when  = function(dl) !is.null(dl$msmMode) && dl$msmMode > 0,
      condition_label = "msmMode > 0 (multispecies)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) {
        all(vapply(.RCE_BIOENERGETICS_SCALARS, function(nm) {
          !is.null(dl[[nm]]) && length(dl[[nm]]) == dl$nspp
        }, logical(1)))
      },
      message        = function(dl) {
        bad <- .RCE_BIOENERGETICS_SCALARS[vapply(.RCE_BIOENERGETICS_SCALARS,
          function(nm) is.null(dl[[nm]]) || length(dl[[nm]]) != dl$nspp,
          logical(1))]
        paste0(
          "msmMode > 0 requires bioenergetics scalars of length nspp = ",
          dl$nspp, "; missing or wrong length: ", paste(bad, collapse = ", ")
        )
      },
      optional_status = "defaulted",
      default_label  = "safe sentinels in single-species mode (switch_check)"
    ),

    # ---- Fixed selectivity requires empirical selectivity ------------------
    list(
      element        = "emp_sel",
      category       = "fishery",
      required_when  = function(dl) {
        fc <- dl$fleet_control
        !is.null(fc) && any(!is.na(fc$Selectivity) & fc$Selectivity == "Fixed")
      },
      condition_label = "any fleet Selectivity == 'Fixed'",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$emp_sel),
      message        = function(dl) {
        fc <- dl$fleet_control
        fixed_flts <- fc[!is.na(fc$Selectivity) & fc$Selectivity == "Fixed",
                         , drop = FALSE]
        paste0(
          "Fleet(s) with Selectivity = 'Fixed' (",
          paste(fixed_flts$Fleet_name, collapse = ", "),
          ") require emp_sel data"
        )
      },
      optional_status = "defaulted",
      default_label  = "empty emp_sel (clean_data)"
    ),

    # ---- Fixed numbers-at-age requires NByageFixed -------------------------
    list(
      element        = "NByageFixed",
      category       = "dimensions",
      required_when  = function(dl) any(dl$estDynamics > 0),
      condition_label = "any estDynamics > 0 (fixed numbers-at-age)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$NByageFixed),
      message        = function(dl)
        "estDynamics > 0 requires NByageFixed data to be provided",
      optional_status = "defaulted",
      default_label  = "empty NByageFixed (clean_data)"
    ),

    # ---- Growth estimation requires conditional age-at-length --------------
    # caal_data is ALSO required (as an *adequacy* gate, kept imperative in
    # data_check) for fleets with estimated selectivity and no comp data; that
    # cross-fleet, mirroring-aware check does not reduce to a declarative row.
    list(
      element        = "caal_data",
      category       = "composition",
      required_when  = function(dl) any(dl$growth_model > 0),
      condition_label = "any growth_model > 0 (estimated growth)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$caal_data),
      message        = function(dl)
        "Growth estimation (growth_model > 0) requires caal_data to be provided",
      optional_status = "defaulted",
      default_label  = "empty caal_data (clean_data)"
    ),

    # ---- Environmental index required by temperature-dependent consumption -
    # Introspection-only: the real check (per-species Cindex column count) is a
    # dimension gate and stays imperative in data_check().
    list(
      element        = "env_data",
      category       = "environment",
      required_when  = function(dl) !is.null(dl$Ceq) && any(dl$Ceq > 1),
      condition_label = "any Ceq > 1 (temperature-dependent consumption)",
      driven         = FALSE,
      optional_status = "defaulted",
      default_label  = "Year-only env_data (clean_data)"
    ),

    # ---- MVN survey covariance --------------------------------------------
    # Introspection-only: presence + square + dimension + symmetry all live in
    # data_check() (they need the per-fleet fitted-observation count).
    list(
      element        = "index_cov",
      category       = "fishery",
      required_when  = function(dl) {
        fc <- dl$fleet_control
        !is.null(fc) && "Index_loglike" %in% colnames(fc) &&
          any(fc$Index_loglike %in% c("MVN", "MVNORM", 1, 2, "1", "2"))
      },
      condition_label = "any fleet Index_loglike == 'MVN'",
      driven         = FALSE,
      optional_status = "defaulted",
      default_label  = "empty list (clean_data)"
    )
  )

  names(rows) <- vapply(rows, function(r) r$element, character(1))
  rows
}

# -----------------------------------------------------------------------------
# Evaluators used by data_check() at each driven gate site. Placing the call at
# the original site preserves the accumulation order of `errors`; the table
# supplies the condition, adequacy test, and message text.
# -----------------------------------------------------------------------------

#' Presence-requirement check for a `severity = "error"` row
#'
#' @param dl Rceattle data list.
#' @param element Row key (data_list element name).
#' @return `character(0)` if the requirement is not active or is satisfied;
#'   otherwise the single error string for the caller to accumulate.
#' @keywords internal
#' @noRd
.rce_check_presence <- function(dl, element) {
  row <- .rce_requirement_table()[[element]]
  if (is.null(row) || !isTRUE(row$driven) || row$severity != "error") {
    stop("Internal: no driven error requirement '", element, "'")
  }
  if (isTRUE(row$required_when(dl)) && !isTRUE(row$adequate(dl))) {
    row$message(dl)
  } else {
    character(0)
  }
}

#' Presence-requirement notice for a `severity = "message"` row
#'
#' Emits the row message via [message()] when the requirement is active and
#' unmet (matching the historical `message("No ration data")` behavior).
#'
#' @param dl Rceattle data list.
#' @param element Row key.
#' @return Invisibly `NULL`.
#' @keywords internal
#' @noRd
.rce_notify_presence <- function(dl, element) {
  row <- .rce_requirement_table()[[element]]
  if (is.null(row) || !isTRUE(row$driven) || row$severity != "message") {
    stop("Internal: no driven message requirement '", element, "'")
  }
  if (isTRUE(row$required_when(dl)) && !isTRUE(row$adequate(dl))) {
    message(row$message(dl))
  }
  invisible(NULL)
}
