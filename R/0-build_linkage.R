#' Formula-driven linkage specifications for Rceattle processes
#'
#' These helpers let users describe how a process parameter (e.g.
#' `log_alpha` for the Beverton-Holt SRR, `log_M` for natural mortality,
#' `log_K` for von Bertalanffy growth) depends on environmental
#' covariates and on stratifying factors (species, sex). They produce:
#'
#'   1. an `Rceattle_linkage_spec` object that captures the user's
#'      intent (formula + grouping) without committing to a global
#'      column index, and
#'   2. a [materialize_linkage()] step that, given the env data and
#'      stratum levels, expands the spec into the canonical long-format
#'      linkage-table rows consumed by TMB.
#'
#' Splitting capture from materialization lets `fit_mod()` pool specs
#' from every process, build a single shared design matrix, and assign
#' globally consistent `X_col` indices in one place.
#'
#' @keywords internal
#' @name build_linkage
NULL


#' Capture a linkage specification
#'
#' @param formula one-sided R formula whose RHS describes the linear
#'   predictor for `param` (e.g. `~ 1`, `~ temp`, `~ temp + PDO`).
#' @param param target parameter name on the linear predictor scale
#'   (e.g. `"log_alpha"`, `"log_M"`, `"log_K"`).
#' @param by one-sided formula naming stratifying factors that should
#'   each get their own coefficients (e.g. `~species`,
#'   `~species + sex`). `NULL` means the same coefficients apply to
#'   every species/sex.
#' @param link link function applied to the predictor when assembling
#'   process values; one of [LINKAGE_LINKS].
#' @param init optional named numeric vector of initial values keyed by
#'   the design-matrix column name (e.g.
#'   `c(`(Intercept)` = -1, temp = 0)`). Missing entries default to `0`.
#' @param bounds optional named list of `c(lower, upper)` keyed the same
#'   way as `init`.
#' @param priors optional named list whose entries are
#'   [Rceattle_priors] objects, keyed by design-matrix column name.
#'   Inside this argument the unprefixed shorthand `normal()`,
#'   `lognormal()`, `gamma()`, and `beta()` resolves to the
#'   corresponding `prior_*` constructors via a private data mask, so
#'   `priors = list(temp = normal(0, 1))` works without masking
#'   [base::gamma()] or [base::beta()] at the package level. Equivalent
#'   to `priors = list(temp = prior_normal(0, 1))`.
#' @param re_group optional character: name of a random-effect grouping
#'   for these coefficients. `NA` (default) means fixed.
#' @param est_phase optional integer estimation phase. Default `1L`.
#'
#' @return An `Rceattle_linkage_spec` object.
#' @keywords internal
#' @importFrom rlang enquo eval_tidy
linkage_spec <- function(formula,
                         param,
                         by        = NULL,
                         link      = "identity",
                         init      = NULL,
                         bounds    = NULL,
                         priors    = NULL,
                         re_group  = NA_character_,
                         est_phase = 1L) {
  priors_quo <- rlang::enquo(priors)
  priors_obj <- rlang::eval_tidy(priors_quo, data = .prior_dispatch_mask())
  priors_obj <- .validate_priors_arg(priors_obj)

  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula (e.g. ~ temp + PDO)")
  }
  if (length(formula) != 2L) {
    stop("`formula` must be one-sided (RHS only); got: ",
         deparse(formula))
  }
  if (!is.null(by) && !inherits(by, "formula")) {
    stop("`by` must be a one-sided formula (e.g. ~species + sex) or NULL")
  }
  link <- match.arg(link, LINKAGE_LINKS)

  structure(
    list(
      formula   = formula,
      param     = as.character(param),
      by        = by,
      link      = link,
      init      = init   %||% list(),
      bounds    = bounds %||% list(),
      priors    = priors_obj,
      re_group  = as.character(re_group),
      est_phase = as.integer(est_phase)
    ),
    class = "Rceattle_linkage_spec"
  )
}


#' Validate the `priors` argument once it has been NSE-evaluated.
#'
#' Accepts `NULL` or a named list of `Rceattle_prior` objects. Returns
#' the canonicalized list (always a named list, possibly empty).
#'
#' @keywords internal
#' @noRd
.validate_priors_arg <- function(priors) {
  if (is.null(priors)) return(list())
  if (!is.list(priors) || (length(priors) > 0 && is.null(names(priors)))) {
    stop("`priors` must be a named list keyed by design-matrix ",
         "column name (e.g. list(temp = normal(0, 1)))",
         call. = FALSE)
  }
  bad <- !vapply(priors, is_rceattle_prior, logical(1))
  if (any(bad)) {
    nm <- names(priors)[bad][1]
    stop(sprintf(
      paste0("priors$%s is not an Rceattle_prior. Use ",
             "normal()/lognormal()/gamma()/beta() inside ",
             "`priors = list(...)`, or prior_normal() etc. ",
             "when constructing programmatically."),
      nm), call. = FALSE)
  }
  priors
}


#' Null-coalesce
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' @export
print.Rceattle_linkage_spec <- function(x, ...) {
  cat("<Rceattle linkage spec>\n")
  cat("  param:   ", x$param, "\n", sep = "")
  cat("  formula: ", deparse(x$formula), "\n", sep = "")
  cat("  by:      ",
      if (is.null(x$by)) "(shared)" else deparse(x$by),
      "\n", sep = "")
  cat("  link:    ", x$link, "\n", sep = "")
  cat("  phase:   ", x$est_phase, "\n", sep = "")
  if (!is.na(x$re_group)) {
    cat("  re_group:", x$re_group, "\n", sep = "")
  }
  invisible(x)
}


#' Materialize a linkage spec into linkage-table rows
#'
#' Expands an [linkage_spec()] into rows of the canonical linkage table.
#' One row is produced for every combination of:
#'
#'   * design-matrix column from `model.matrix(spec$formula, env_data)`,
#'   * stratum implied by `spec$by` (e.g. one per species when
#'     `by = ~species`).
#'
#' The `X_col` column is initially set to a *local* column index into
#' the per-spec design matrix; it is the caller's responsibility (in
#' `fit_mod()`'s pooling step) to remap these to global column indices
#' once all specs have been combined into a single shared `X` matrix.
#'
#' @param spec an `Rceattle_linkage_spec`.
#' @param process one of [LINKAGE_PROCESSES].
#' @param env_data data frame of environmental covariates (one row per
#'   model year). Must contain every variable named on the RHS of
#'   `spec$formula`.
#' @param strata named list giving the discrete levels for each
#'   stratifying factor named in `spec$by`. For example, for
#'   `by = ~species` the user must supply `strata = list(species = 1:3)`.
#'   Each element should be a 1-based integer vector of stratum ids.
#'   Allowed names are `"species"`, `"sex"`, and `"age_bin"`.
#'
#' @return An `Rceattle_linkage_table` with one row per coefficient.
#' @keywords internal
materialize_linkage <- function(spec, process, env_data, strata = list()) {
  if (!inherits(spec, "Rceattle_linkage_spec")) {
    stop("`spec` must be an Rceattle_linkage_spec")
  }
  process <- match.arg(process, LINKAGE_PROCESSES)
  if (!is.data.frame(env_data)) {
    stop("`env_data` must be a data.frame")
  }
  X <- stats::model.matrix(spec$formula, data = env_data)
  X_names <- colnames(X)
  n_cols  <- ncol(X)

  by_vars <- if (is.null(spec$by)) character(0) else all.vars(spec$by)
  unknown_by <- setdiff(by_vars, c("species", "sex", "age_bin"))
  if (length(unknown_by) > 0) {
    stop("unknown grouping variable(s) in `by`: ",
         paste(unknown_by, collapse = ", "),
         "; allowed: species, sex, age_bin")
  }
  for (v in by_vars) {
    if (is.null(strata[[v]])) {
      stop(sprintf("`strata$%s` must be supplied because `by` references it",
                   v))
    }
  }
  level_grid <- if (length(by_vars) == 0L) {
    data.frame(.dummy = 1L)
  } else {
    do.call(expand.grid,
            c(strata[by_vars], list(KEEP.OUT.ATTRS = FALSE,
                                    stringsAsFactors = FALSE)))
  }

  rows <- vector("list", n_cols * nrow(level_grid))
  k <- 0L
  for (col_idx in seq_len(n_cols)) {
    cn <- X_names[col_idx]
    init_val   <- spec$init[[cn]]   %||% 0
    bounds_val <- spec$bounds[[cn]] %||% c(-Inf, Inf)
    prior <- spec$priors[[cn]]
    if (is.null(prior)) {
      pf <- "none"; pp1 <- NA_real_; pp2 <- NA_real_
    } else {
      pf <- prior$family; pp1 <- prior$p1; pp2 <- prior$p2
    }
    for (g in seq_len(nrow(level_grid))) {
      sp_id <- if ("species" %in% by_vars) level_grid$species[g] else NA_integer_
      sx_id <- if ("sex"     %in% by_vars) level_grid$sex[g]     else NA_integer_
      ab_id <- if ("age_bin" %in% by_vars) level_grid$age_bin[g] else NA_integer_
      k <- k + 1L
      rows[[k]] <- linkage_row(
        process      = process,
        param        = spec$param,
        X_col        = col_idx,
        species      = sp_id,
        sex          = sx_id,
        age_bin      = ab_id,
        link         = spec$link,
        init         = init_val,
        lower        = bounds_val[1],
        upper        = bounds_val[2],
        prior_family = pf,
        prior_p1     = pp1,
        prior_p2     = pp2,
        re_group     = spec$re_group,
        est_phase    = spec$est_phase
      )
    }
  }
  out <- bind_linkage(rows)
  attr(out, "design_colnames") <- X_names
  out
}
