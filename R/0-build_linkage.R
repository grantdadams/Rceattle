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
#'   (e.g. `"log_alpha"`, `"log_M1"`, `"log_K"`). May be `NULL` when the
#'   spec is built inside a `build_*()` call that infers the parameter
#'   name from the enclosing list key (see [build_growth()]).
#' @param by one-sided formula naming stratifying factors that should
#'   each get their own coefficients. Allowed names are `species`,
#'   `sex`, and `age_bin`. The default `~species` produces one
#'   coefficient set per species (the typical multispecies
#'   assessment use case); pass `~species + sex` for per-(species,
#'   sex) coefficients, or `NULL` to share a single coefficient
#'   set across every species/sex.
#' @param species optional integer vector of 1-based species ids that
#'   this spec applies to. `NULL` (default) means every species in
#'   `strata$species` at materialization time. Use this to give
#'   different species different formulas, e.g. by registering
#'   multiple specs against the same parameter -- see
#'   [build_growth()] for the multi-spec syntax.
#' @param link link function applied to the predictor when assembling
#'   process values; one of [LINKAGE_LINKS]. TODO
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
#' @export
#' @importFrom rlang enquo eval_tidy
linkage_spec <- function(formula,
                         param     = NULL,
                         by        = ~ species,
                         species   = NULL,
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
  if (!is.null(species)) {
    species <- as.integer(species)
    if (anyNA(species) || any(species < 1L)) {
      stop("`species` must be a vector of positive 1-based species ids",
           call. = FALSE)
    }
  }
  link <- match.arg(link, LINKAGE_LINKS)
  param_str <- if (is.null(param)) NA_character_ else as.character(param)

  structure(
    list(
      formula   = formula,
      param     = param_str,
      by        = by,
      species   = species,
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


#' Generic linkage-list validator used by `build_*()` helpers
#'
#' Each `build_*()` function has near-identical validation:
#'
#'   1. The argument must be either `NULL` or a non-empty named list.
#'   2. Each list key must be one of the per-process allowed
#'      parameter names (e.g. `c("log_K", "log_L1", ...)` for growth).
#'   3. Each list value must be either an `Rceattle_linkage_spec` or
#'      a list of them (the per-species-formula form).
#'
#' This helper centralizes those three checks; per-process wrappers
#' (`.validate_growth_linkages`, `.validate_M_linkages`,
#' `.validate_recruitment_linkages`) call this and add their own
#' domain-specific warnings on top, then run [`.stamp_param`] to
#' fill in `param` from the list keys.
#'
#' @param linkages the user-supplied `linkages` argument.
#' @param allowed_params character vector of valid parameter names
#'   for this process (e.g. [GROWTH_LINKAGE_PARAMS]).
#' @param process_label one of `"growth"`, `"M"`, `"recruitment"`,
#'   used only to make error messages actionable.
#'
#' @return `NULL` (when `linkages = NULL`) or the canonicalized list
#'   with `param` filled in on every spec. Errors loudly on shape
#'   problems.
#' @keywords internal
#' @noRd
.validate_process_linkages <- function(linkages, allowed_params,
                                       process_label) {
  if (is.null(linkages)) return(NULL)
  if (!is.list(linkages) || length(linkages) == 0L ||
      is.null(names(linkages)) || any(!nzchar(names(linkages)))) {
    stop(sprintf(
      "`linkages` must be a named list keyed by %s parameter (one of: %s)",
      process_label,
      paste(allowed_params, collapse = ", ")), call. = FALSE)
  }
  bad_names <- setdiff(names(linkages), allowed_params)
  if (length(bad_names) > 0) {
    stop(sprintf(
      "unknown %s linkage parameter(s): %s; allowed: %s",
      process_label,
      paste(bad_names, collapse = ", "),
      paste(allowed_params, collapse = ", ")), call. = FALSE)
  }
  for (nm in names(linkages)) {
    val <- linkages[[nm]]
    if (inherits(val, "Rceattle_linkage_spec")) next
    if (is.list(val) &&
        all(vapply(val, inherits, logical(1),
                   what = "Rceattle_linkage_spec"))) next
    stop(sprintf(
      "linkages$%s must be a linkage_spec() or a list of linkage_spec() objects.",
      nm), call. = FALSE)
  }
  Map(.stamp_param, linkages, names(linkages))
}


#' Set or override the target parameter name on a linkage spec
#'
#' Used by `build_*()` helpers that infer the parameter name from the
#' list key under which a spec is registered (e.g.
#' `linkages = list(log_K = linkage_spec(~temp))` -> set `param =
#' "log_K"`). If the spec already names a different parameter the
#' function errors to surface user mistakes.
#'
#' @param spec an `Rceattle_linkage_spec`.
#' @param param target parameter name (character scalar).
#' @return The spec, with `param` set.
#' @keywords internal
.set_linkage_param <- function(spec, param) {
  param <- as.character(param)
  cur <- spec$param
  if (!is.na(cur) && nzchar(cur) && cur != param) {
    stop(sprintf(
      "linkage spec param '%s' conflicts with the list key '%s'",
      cur, param), call. = FALSE)
  }
  spec$param <- param
  spec
}


#' Validate the `priors` argument once it has been NSE-evaluated.
#'
#' Accepts `NULL` or a named list keyed by design-matrix column name.
#' Each entry may be:
#'
#'   * an `Rceattle_prior` object -- the prior applies to every
#'     species / sex / age row that uses this column, or
#'   * a named list keyed by species id (character or integer); each
#'     value is itself either
#'       - an `Rceattle_prior` (applies to every sex of that species),
#'         or
#'       - a named list keyed by sex id, mapping to `Rceattle_prior`
#'         objects (per-(species, sex) priors).
#'
#' Rows for species/sex combinations not in the keyset get no prior.
#' Per-species and per-sex keying only meaningfully apply when the
#' spec's `by` includes the corresponding stratum.
#'
#' Returns the canonicalized list (always a named list, possibly
#' empty).
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
  for (nm in names(priors)) {
    val <- priors[[nm]]
    if (is_rceattle_prior(val)) next
    if (is.list(val)) {
      if (is.null(names(val)) || any(!nzchar(names(val)))) {
        stop(sprintf(
          paste0("priors$%s must be a named list keyed by species ",
                 "id when supplying per-species priors (e.g. ",
                 "list(`1` = normal(0, 1), `2` = normal(0, 0.5)))."),
          nm), call. = FALSE)
      }
      for (sp_key in names(val)) {
        sp_branch <- val[[sp_key]]
        if (is_rceattle_prior(sp_branch)) next
        if (is.list(sp_branch)) {
          if (is.null(names(sp_branch)) ||
              any(!nzchar(names(sp_branch)))) {
            stop(sprintf(
              paste0("priors$%s$`%s` must be a named list keyed by ",
                     "sex id (e.g. list(`1` = normal(0, 0.1), ",
                     "`2` = normal(0, 0.2)))."),
              nm, sp_key), call. = FALSE)
          }
          bad_sex <- !vapply(sp_branch, is_rceattle_prior, logical(1))
          if (any(bad_sex)) {
            stop(sprintf(
              "priors$%s$`%s`$`%s` is not an Rceattle_prior",
              nm, sp_key, names(sp_branch)[bad_sex][1]
            ), call. = FALSE)
          }
          next
        }
        stop(sprintf(
          paste0("priors$%s$`%s` must be an Rceattle_prior or a ",
                 "named list of priors keyed by sex id."),
          nm, sp_key), call. = FALSE)
      }
      next
    }
    stop(sprintf(
      paste0("priors$%s must be an Rceattle_prior or a named list ",
             "of priors keyed by species id. Use ",
             "normal()/lognormal()/gamma()/beta() inside ",
             "`priors = list(...)`, or prior_normal() etc. when ",
             "constructing programmatically."),
      nm), call. = FALSE)
  }
  priors
}


#' Resolve the prior to apply on a specific (column, species, sex) row.
#'
#' Walks the nested `priors` list: scalar -> applies to all,
#' species-keyed -> applies to that species' rows, species-then-sex
#' -> applies to that exact (species, sex) cell. A missing key at
#' either level returns `NULL` (no prior).
#'
#' @param prior_spec value at `spec$priors[[colname]]`.
#' @param sp_id 1-based species id, or `NA` for shared rows.
#' @param sx_id 1-based sex id, or `NA` for sex-shared rows.
#' @return an `Rceattle_prior` or `NULL`.
#' @keywords internal
#' @noRd
.resolve_prior <- function(prior_spec, sp_id, sx_id = NA_integer_) {
  if (is.null(prior_spec)) return(NULL)
  if (is_rceattle_prior(prior_spec)) return(prior_spec)
  if (is.list(prior_spec)) {
    if (is.na(sp_id)) return(NULL)
    sp_branch <- prior_spec[[as.character(sp_id)]]
    if (is.null(sp_branch)) return(NULL)
    if (is_rceattle_prior(sp_branch)) return(sp_branch)
    if (is.list(sp_branch)) {
      if (is.na(sx_id)) return(NULL)
      return(sp_branch[[as.character(sx_id)]])
    }
  }
  NULL
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
  if (!is.null(x$species)) {
    cat("  species: ", paste(x$species, collapse = ", "), "\n", sep = "")
  }
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
  if (is.na(spec$param) || !nzchar(spec$param)) {
    stop("linkage spec is missing a `param`; supply it via ",
         "`linkage_spec(param = ...)` or via the list key in a ",
         "build_*() linkages argument.", call. = FALSE)
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
  if ("sex" %in% by_vars) {
    sex_levels_one <- FALSE
    if (is.list(strata$sex)) {
      if (!"species" %in% by_vars) {
        stop("species-specific `strata$sex` values are only supported when `by` includes species",
             call. = FALSE)
      }
      sex_levels_one <- all(vapply(strata$sex, length, integer(1)) == 1L)
    } else if (length(strata$sex) == 1L) {
      sex_levels_one <- TRUE
    }
    if (sex_levels_one) {
      warning("`by = ~ ... + sex` was requested but only one sex level is available in `strata$sex`; \
              sex-specific coefficients will collapse to a single shared sex level.")
    }
  }

  expand_linkage_strata <- function(strata, by_vars) {
    if (length(by_vars) == 0L) {
      return(data.frame(.dummy = 1L))
    }

    if (!"species" %in% by_vars) {
      if (any(vapply(strata[by_vars], is.list, logical(1)))) {
        stop("species-specific strata values are only supported when `by` includes species",
             call. = FALSE)
      }
      return(do.call(expand.grid,
                     c(strata[by_vars], list(KEEP.OUT.ATTRS = FALSE,
                                             stringsAsFactors = FALSE))))
    }

    species_vals <- as.integer(strata$species)
    rows <- vector("list", length(species_vals))
    for (i in seq_along(species_vals)) {
      sp <- species_vals[i]
      row_vars <- list(species = sp)
      for (v in setdiff(by_vars, "species")) {
        vals <- strata[[v]]
        if (is.list(vals)) {
          if (is.null(names(vals))) {
            if (length(vals) != length(species_vals)) {
              stop(sprintf(
                "species-specific strata$%s must be a named list keyed by species id or a list of length %d",
                v, length(species_vals)), call. = FALSE)
            }
            vals <- vals[[i]]
          } else {
            sp_key <- as.character(sp)
            if (!sp_key %in% names(vals)) {
              stop(sprintf(
                "species-specific strata$%s must include an entry for species %s",
                v, sp_key), call. = FALSE)
            }
            vals <- vals[[sp_key]]
          }
        }
        row_vars[[v]] <- vals
      }
      rows[[i]] <- do.call(expand.grid,
                            c(row_vars, list(KEEP.OUT.ATTRS = FALSE,
                                            stringsAsFactors = FALSE)))
    }
    out <- do.call(rbind, rows)
    row.names(out) <- NULL
    out
  }

  level_grid <- expand_linkage_strata(strata, by_vars)

  # Honor an optional species filter set on the spec. When `species`
  # is supplied, only rows for those species ids are emitted; species
  # not represented in `by` are unaffected (the filter is a no-op
  # against a spec that doesn't stratify by species).
  if (!is.null(spec$species) && "species" %in% names(level_grid)) {
    level_grid <- level_grid[level_grid$species %in% spec$species, ,
                             drop = FALSE]
  }
  if (nrow(level_grid) == 0L) {
    return(.empty_materialized(X, X_names))
  }

  rows <- vector("list", n_cols * nrow(level_grid))
  k <- 0L
  for (col_idx in seq_len(n_cols)) {
    cn <- X_names[col_idx]
    init_val   <- spec$init[[cn]]   %||% 0
    bounds_val <- spec$bounds[[cn]] %||% c(-Inf, Inf)
    prior_spec <- spec$priors[[cn]]
    for (g in seq_len(nrow(level_grid))) {
      sp_id <- if ("species" %in% by_vars) level_grid$species[g] else NA_integer_
      sx_id <- if ("sex"     %in% by_vars) level_grid$sex[g]     else NA_integer_
      ab_id <- if ("age_bin" %in% by_vars) level_grid$age_bin[g] else NA_integer_

      prior <- .resolve_prior(prior_spec, sp_id, sx_id)
      if (is.null(prior)) {
        pf <- "none"; pp1 <- NA_real_; pp2 <- NA_real_
      } else {
        pf <- prior$family; pp1 <- prior$p1; pp2 <- prior$p2
      }

      k <- k + 1L
      rows[[k]] <- linkage_row(
        process      = process,
        param        = spec$param,
        X_col        = col_idx,
        species      = sp_id,
        sex          = sx_id,
        age_bin      = ab_id,
        design_col   = cn,
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
  attr(out, "design_matrix")   <- X
  out
}


#' Empty materialized table that still carries design metadata.
#'
#' Used by `materialize_linkage()` when a `species` filter on the spec
#' eliminates every row of the level grid. We still need to expose the
#' design matrix in the attributes so the pooler can union columns
#' across specs that did emit rows.
#'
#' @keywords internal
#' @noRd
.empty_materialized <- function(X, X_names) {
  out <- new_linkage_table()
  attr(out, "design_colnames") <- X_names
  attr(out, "design_matrix")   <- X
  out
}


#' Pool a collection of linkage specs into a global table + design matrix
#'
#' Collects every spec that's been registered against a process,
#' materializes each against the supplied `env_data` and `strata`,
#' unions their design-matrix columns by name into a single shared
#' design matrix `X`, and remaps each row's local `X_col` index to the
#' global column. Used inside `fit_mod()` once all `build_*()` outputs
#' have populated `data_list`.
#'
#' Two design columns from different specs that share a name are
#' assumed to carry identical numeric values (which holds for
#' deterministic terms generated by [stats::model.matrix()] over the
#' same `env_data`). The first occurrence wins on order; subsequent
#' specs reuse the column.
#'
#' @param spec_groups named list whose names are process labels
#'   (one of [LINKAGE_PROCESSES]) and whose elements are themselves
#'   named lists of [linkage_spec()] objects keyed by parameter name.
#'   Empty entries are skipped.
#' @param env_data data.frame of environmental covariates (one row per
#'   model year, including any factor columns referenced by spec
#'   formulas).
#' @param strata named list of stratum-level integer vectors, e.g.
#'   `list(species = 1:nspp, sex = 1:nsex)`.
#'
#' @return A list with components:
#' \describe{
#'   \item{table}{An `Rceattle_linkage_table` with `X_col` indices
#'     remapped to the global design matrix. May have zero rows.}
#'   \item{X}{A numeric matrix `[nrow(env_data) x n_global_cols]` with
#'     `colnames(X)` = the union of design-column names across specs.}
#'   \item{design_names}{`colnames(X)` (separate copy for caller
#'     convenience).}
#' }
#' @keywords internal
pool_linkages <- function(spec_groups, env_data, strata = list()) {
  has_specs <- !is.null(spec_groups) && length(spec_groups) > 0L &&
    any(vapply(spec_groups, length, integer(1)) > 0L)
  if (!has_specs) return(.empty_pool(env_data))
  if (!is.data.frame(env_data)) {
    stop("`env_data` must be a data.frame when linkages are ",
         "supplied; set `data_list$env_data` before calling ",
         "`fit_mod()`.", call. = FALSE)
  }
  per_spec <- list()
  for (proc in names(spec_groups)) {
    specs <- spec_groups[[proc]]
    if (is.null(specs) || length(specs) == 0L) next
    for (param in names(specs)) {
      # Each value may be a single spec, or a list of specs that all
      # target the same parameter (e.g. species-specific formulas
      # registered together).
      this <- specs[[param]]
      spec_list <- if (inherits(this, "Rceattle_linkage_spec")) {
        list(this)
      } else if (is.list(this)) {
        this
      } else {
        stop(sprintf(
          "spec_groups$%s$%s must be a linkage_spec() or a list of them",
          proc, param), call. = FALSE)
      }
      for (one in spec_list) {
        if (!inherits(one, "Rceattle_linkage_spec")) {
          stop(sprintf(
            "an entry under spec_groups$%s$%s is not a linkage_spec()",
            proc, param), call. = FALSE)
        }
        tbl <- materialize_linkage(.set_linkage_param(one, param),
                                   process = proc,
                                   env_data = env_data,
                                   strata   = strata)
        per_spec[[length(per_spec) + 1L]] <- tbl
      }
    }
  }
  if (length(per_spec) == 0L) return(.empty_pool(env_data))

  # Union of design-column names, in first-seen order.
  global_names <- unique(unlist(
    lapply(per_spec, function(t) attr(t, "design_colnames")),
    use.names = FALSE))

  # Build the global design matrix by overlaying per-spec columns by
  # name. Same-named columns must carry identical numeric values; we
  # take the first occurrence and verify on subsequent ones.
  nyr <- nrow(env_data)
  global_X <- matrix(0, nrow = nyr, ncol = length(global_names))
  colnames(global_X) <- global_names
  filled <- rep(FALSE, length(global_names))
  for (i in seq_along(per_spec)) {
    X_i <- attr(per_spec[[i]], "design_matrix")
    for (cn in colnames(X_i)) {
      gi <- match(cn, global_names)
      col_vals <- as.numeric(X_i[, cn])
      if (!filled[gi]) {
        global_X[, gi] <- col_vals
        filled[gi] <- TRUE
      } else if (!isTRUE(all.equal(as.numeric(global_X[, gi]), col_vals))) {
        stop(sprintf(
          "design column '%s' has inconsistent values across linkage specs",
          cn), call. = FALSE)
      }
    }
  }

  # Remap each per-spec table's local X_col -> global X_col.
  remapped <- vector("list", length(per_spec))
  for (i in seq_along(per_spec)) {
    tbl <- per_spec[[i]]
    local_names <- attr(tbl, "design_colnames")
    tbl$X_col <- match(local_names[tbl$X_col], global_names)
    attr(tbl, "design_colnames") <- NULL
    attr(tbl, "design_matrix")   <- NULL
    remapped[[i]] <- tbl
  }

  list(
    table        = bind_linkage(remapped),
    X            = global_X,
    design_names = global_names
  )
}


#' @keywords internal
#' @noRd
.empty_pool <- function(env_data) {
  nyr <- if (is.data.frame(env_data)) nrow(env_data) else 0L
  list(
    table        = new_linkage_table(),
    X            = matrix(numeric(0), nrow = nyr, ncol = 0),
    design_names = character(0)
  )
}
