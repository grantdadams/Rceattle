#' Formula-driven linkage specifications for Rceattle processes
#'
#' These helpers let users describe how a process parameter (e.g.
#' `alpha` for the Beverton-Holt SRR, `M1` for natural mortality,
#' `K` for von Bertalanffy growth) depends on environmental
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
#' @param param target parameter name on the natural scale
#'   (e.g. `"alpha"`, `"M1"`, `"K"`). May be `NULL` when the
#'   spec is built inside a `build_*()` call that infers the parameter
#'   name from the enclosing list key (see [build_growth()]).
#' @param data (Optional) data frame for formula validation. Currently
#'   validation happens at materialization time inside [fit_mod()].
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
#' @param sex optional vector of sex ids that this spec applies to.
#'   May be supplied as integers (`1L` = female, `2L` = male) or as
#'   character strings (`"Females"`/`"Males"`, case-insensitive;
#'   `"female"`, `"male"`, `"f"`, `"m"` are also accepted). `NULL`
#'   (default) means every sex in `strata$sex` at materialization
#'   time. Only meaningful when `by` includes `sex`; otherwise the
#'   filter is a no-op. Use this to register separate specs per sex
#'   (e.g. one prior on females, another on males) against the same
#'   parameter.
#' @param fleet optional integer vector of 1-based `Fleet_code`s this spec
#'   applies to. `NULL` (default) means every fleet in `strata$fleet` at
#'   materialization time. Only meaningful when `by` includes `fleet`;
#'   otherwise the filter is a no-op. Used by catchability and selectivity
#'   linkages to give different fleets different formulas.
#' @param link link function relating the linear predictor to the
#'   natural-scale target parameter. One of `"log"` (default) or
#'   `"identity"`. With `link = "log"`, `log(param) = X * beta` -- slope
#'   contributions are multiplicative on the natural-scale parameter. With `link = "identity"`,
#'   `param = X * beta` -- slope contributions are additive on the
#'   natural scale. All linkage targets currently expose log-scale TMB
#'   parameters, so `"log"` is the natural default; `"identity"` is
#'   reserved for future processes (e.g. logit for steepness).
#' @param init optional named numeric vector of initial values keyed by
#'   the design-matrix column name (e.g.
#'   `c(`(Intercept)` = 4, temp = 0)`). Missing entries default to `0`.
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
                         data      = NULL,
                         by        = ~ species,
                         species   = NULL,
                         sex       = NULL,
                         fleet     = NULL,
                         link      = "log",
                         init      = NULL,
                         bounds    = NULL,
                         priors    = NULL,
                         re_group  = NA_character_,
                         est_phase = 1L,
                         observe   = NULL,
                         obs_sd    = NULL) {
  priors_quo <- rlang::enquo(priors)
  priors_obj <- rlang::eval_tidy(priors_quo, data = .prior_dispatch_mask())
  priors_obj <- .validate_priors_arg(priors_obj)

  # State-space (Rogers et al. 2024 QAR1) option: an ar1 latent that is also
  # OBSERVED as an env_data column with a fixed measurement SD, and enters the
  # linked parameter through an estimated effect size (beta). Only meaningful
  # with an `ar1(1 | group)` term.
  if (!is.null(observe)) {
    if (!is.character(observe) || length(observe) != 1L || !nzchar(observe)) {
      stop("`observe` must be a single env_data column name.", call. = FALSE)
    }
    if (is.null(obs_sd) || !is.numeric(obs_sd) || length(obs_sd) != 1L ||
        !is.finite(obs_sd) || obs_sd <= 0) {
      stop("`observe` requires a positive fixed `obs_sd` (the measurement SD).",
           call. = FALSE)
    }
  } else if (!is.null(obs_sd)) {
    stop("`obs_sd` is only used with `observe`.", call. = FALSE)
  }

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
  if (!is.null(sex)) {
    sex <- .coerce_sex_arg(sex)
  }
  if (!is.null(fleet)) {
    fleet <- as.integer(fleet)
    if (anyNA(fleet) || any(fleet < 1L)) {
      stop("`fleet` must be a vector of positive 1-based Fleet_code values",
           call. = FALSE)
    }
  }
  link <- match.arg(link, LINKAGE_LINKS)
  .check_link_implemented(link)
  param_str <- if (is.null(param)) NA_character_ else as.character(param)

  init   <- .validate_linkage_init_arg(init)
  bounds <- .validate_linkage_bounds_arg(bounds)

  structure(
    list(
      formula   = formula,
      param     = param_str,
      by        = by,
      species   = species,
      sex       = sex,
      fleet     = fleet,
      link      = link,
      init      = init,
      bounds    = bounds,
      priors    = priors_obj,
      re_group  = as.character(re_group),
      est_phase = as.integer(est_phase),
      observe   = observe,
      obs_sd    = obs_sd
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
#'      parameter names (e.g. `c("K", "L1", ...)` for growth).
#'   3. Each list value must be either an `Rceattle_linkage_spec` or
#'      a list of them (the per-species-formula form).
#'
#' This helper centralizes those three checks; per-process wrappers
#' (\code{.validate_growth_linkages}, \code{.validate_M_linkages},
#' \code{.validate_recruitment_linkages}) call this and add their own
#' domain-specific warnings on top, then run \code{.stamp_param} to
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
#' `linkages = list(K = linkage_spec(~temp))` -> set `param =
#' "K"`). If the spec already names a different parameter the
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

#' Coerce a `sex` argument to 1-based integer ids.
#'
#' Accepts numeric (1, 2) or character ("female"/"females",
#' "male"/"males"; case-insensitive) input and returns a positive
#' integer vector. Errors with a pointed message on unrecognized
#' values.
#'
#' @keywords internal
#' @noRd
.coerce_sex_arg <- function(sex) {
  if (is.character(sex)) {
    key <- tolower(trimws(sex))
    out <- integer(length(key))
    out[] <- NA_integer_
    out[key %in% c("female", "females", "f")] <- 1L
    out[key %in% c("male", "males", "m")]     <- 2L
    if (anyNA(out)) {
      bad <- unique(sex[is.na(out)])
      stop(sprintf(
        "`sex` must be 1/2 or \"Females\"/\"Males\"; got: %s",
        paste(shQuote(bad), collapse = ", ")), call. = FALSE)
    }
    return(out)
  }
  sex <- as.integer(sex)
  if (anyNA(sex) || any(sex < 1L)) {
    stop("`sex` must be a vector of positive 1-based sex ids ",
         "(1 = female, 2 = male) or the strings \"Females\"/\"Males\"",
         call. = FALSE)
  }
  sex
}


.validate_linkage_init_arg <- function(init) {
  if (is.null(init)) return(list())
  if (!is.list(init) || (length(init) > 0 && is.null(names(init)))) {
    stop("`init` must be a named list keyed by design-matrix column name ",
         "(e.g. list(temp = 0)).", call. = FALSE)
  }
  if (length(init) > 0 && any(!nzchar(names(init)))) {
    stop("`init` must be a named list keyed by design-matrix column name.",
         call. = FALSE)
  }
  for (nm in names(init)) {
    val <- init[[nm]]
    if (!is.numeric(val) || length(val) != 1L) {
      stop(sprintf("init$%s must be a numeric scalar", nm),
           call. = FALSE)
    }
  }
  init
}

.validate_linkage_bounds_arg <- function(bounds) {
  if (is.null(bounds)) return(list())
  if (!is.list(bounds) || (length(bounds) > 0 && is.null(names(bounds)))) {
    stop("`bounds` must be a named list keyed by design-matrix column name ",
         "(e.g. list(temp = c(-Inf, Inf))).", call. = FALSE)
  }
  if (length(bounds) > 0 && any(!nzchar(names(bounds)))) {
    stop("`bounds` must be a named list keyed by design-matrix column name.",
         call. = FALSE)
  }
  for (nm in names(bounds)) {
    val <- bounds[[nm]]
    if (!is.numeric(val) || length(val) != 2L) {
      stop(sprintf("bounds$%s must be a numeric vector of length 2", nm),
           call. = FALSE)
    }
  }
  bounds
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
  if (!is.null(x$sex)) {
    cat("  sex:     ", paste(x$sex, collapse = ", "), "\n", sep = "")
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
  .check_process_implemented(process)
  if (!is.data.frame(env_data)) {
    stop("`env_data` must be a data.frame")
  }
  # Split fixed from random before building the design matrix. Passing a
  # formula containing a bar straight to model.matrix() evaluates `1 | Year`
  # as a logical OR and yields a nonsense column ("1 | YearTRUE"), so the
  # split has to happen first.
  parsed <- .parse_linkage_formula(spec$formula)

  # Fixed part: standard design matrix. Random parts (bar terms) are expanded
  # into per-level indicator columns appended after the fixed columns, so the
  # C++ accumulator adds beta[i] * indicator(level) exactly as for a fixed
  # column -- the only differences are which beta vector supplies the
  # coefficient (beta_linkage_re) and that a density is placed on it.
  X_fixed <- stats::model.matrix(parsed$fixed, data = env_data)
  re <- .materialize_re_design(parsed$re_terms, parsed$re_structures, env_data)

  X <- cbind(X_fixed, re$X)
  X_names <- colnames(X)
  n_cols  <- ncol(X)
  # Per-column RE metadata, aligned with the columns of X (NA for fixed).
  col_re_group  <- c(rep(NA_character_, ncol(X_fixed)), re$group)
  col_re_struct <- c(rep(NA_character_, ncol(X_fixed)), re$struct)
  col_re_time   <- c(rep(NA_real_,      ncol(X_fixed)), re$time)

  # State-space (Rogers QAR1) `observe` needs an ar1 term and a real column.
  if (!is.null(spec$observe)) {
    if (!"ar1" %in% re$struct) {
      stop("`observe` (a state-space covariate observation) requires an ",
           "`ar1(1 | group)` term in the formula.", call. = FALSE)
    }
    if (!spec$observe %in% names(env_data)) {
      stop(sprintf("`observe` column `%s` is not in env_data.", spec$observe),
           call. = FALSE)
    }
  }

  by_vars <- if (is.null(spec$by)) character(0) else all.vars(spec$by)
  unknown_by <- setdiff(by_vars, c("species", "sex", "age_bin", "fleet"))
  if (length(unknown_by) > 0) {
    stop("unknown grouping variable(s) in `by`: ",
         paste(unknown_by, collapse = ", "),
         "; allowed: species, sex, age_bin, fleet")
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
      warning(
        "`by = ~ ... + sex` was requested but the model is single-sex ",
        "(`nsex = 1` for every relevant species); sex-specific ",
        "coefficients will collapse to a single shared sex level.",
        call. = FALSE
      )
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
  if (!is.null(spec$sex) && "sex" %in% names(level_grid)) {
    level_grid <- level_grid[level_grid$sex %in% spec$sex, ,
                             drop = FALSE]
  }
  if (!is.null(spec$fleet) && "fleet" %in% names(level_grid)) {
    level_grid <- level_grid[level_grid$fleet %in% spec$fleet, ,
                             drop = FALSE]
  }
  if (nrow(level_grid) == 0L) {
    return(.empty_materialized(X, X_names))
  }

  rows <- vector("list", n_cols * nrow(level_grid))
  k <- 0L
  for (col_idx in seq_len(n_cols)) {
    cn <- X_names[col_idx]
    init_user_supplied <- cn %in% names(spec$init)
    init_val   <- spec$init[[cn]]   %||% 0
    bounds_val <- spec$bounds[[cn]] %||% c(-Inf, Inf)
    prior_spec <- spec$priors[[cn]]
    for (g in seq_len(nrow(level_grid))) {
      sp_id <- if ("species" %in% by_vars) level_grid$species[g] else NA_integer_
      sx_id <- if ("sex"     %in% by_vars) level_grid$sex[g]     else NA_integer_
      ab_id <- if ("age_bin" %in% by_vars) level_grid$age_bin[g] else NA_integer_
      fl_id <- if ("fleet"   %in% by_vars) level_grid$fleet[g]   else NA_integer_

      prior <- .resolve_prior(prior_spec, sp_id, sx_id)
      if (is.null(prior)) {
        pf <- "none"; pp1 <- NA_real_; pp2 <- NA_real_
      } else {
        pf <- prior$family; pp1 <- prior$p1; pp2 <- prior$p2
      }

      # Per-group RE-SD (sigma) routing: only on random-effect columns. The
      # reserved `sigma` key in the spec's init/priors sets the deviation SD's
      # start value (or, with est_phase == 0, its fixed value) and an optional
      # prior. Resolved per (species, sex) so a `by = ~ species` RE can carry a
      # species-specific sigma prior; the fixed-effect columns keep NA.
      if (!is.na(col_re_struct[col_idx])) {
        s_init  <- spec$init[["sigma"]] %||% NA_real_
        s_prior <- .resolve_prior(spec$priors[["sigma"]], sp_id, sx_id)
        if (is.null(s_prior)) {
          spf <- NA_character_; sp1 <- NA_real_; sp2 <- NA_real_
        } else {
          spf <- s_prior$family; sp1 <- s_prior$p1; sp2 <- s_prior$p2
        }
        # rho routing + state-space observation apply only to ar1 columns.
        if (col_re_struct[col_idx] == "ar1") {
          r_init  <- spec$init[["rho"]] %||% NA_real_
          r_prior <- .resolve_prior(spec$priors[["rho"]], sp_id, sx_id)
          if (is.null(r_prior)) {
            rpf <- NA_character_; rp1 <- NA_real_; rp2 <- NA_real_
          } else {
            rpf <- r_prior$family; rp1 <- r_prior$p1; rp2 <- r_prior$p2
          }
          # Rogers QAR1: the latent ar1 deviate is observed as env_data[[observe]]
          # at this level's time, with fixed SD obs_sd.
          if (!is.null(spec$observe)) {
            grp_col <- col_re_group[col_idx]
            o_val <- env_data[[spec$observe]][
              match(col_re_time[col_idx], env_data[[grp_col]])]
            o_val <- if (length(o_val) == 1L) as.numeric(o_val) else NA_real_
            r_obs_v <- o_val; r_obs_sd <- as.numeric(spec$obs_sd)
          } else {
            r_obs_v <- NA_real_; r_obs_sd <- NA_real_
          }
        } else {
          r_init <- NA_real_; rpf <- NA_character_; rp1 <- NA_real_; rp2 <- NA_real_
          r_obs_v <- NA_real_; r_obs_sd <- NA_real_
        }
      } else {
        s_init <- NA_real_; spf <- NA_character_; sp1 <- NA_real_; sp2 <- NA_real_
        r_init <- NA_real_; rpf <- NA_character_; rp1 <- NA_real_; rp2 <- NA_real_
        r_obs_v <- NA_real_; r_obs_sd <- NA_real_
      }

      k <- k + 1L
      rows[[k]] <- linkage_row(
        process       = process,
        param         = spec$param,
        X_col         = col_idx,
        species       = sp_id,
        sex           = sx_id,
        age_bin       = ab_id,
        fleet         = fl_id,
        design_col    = cn,
        link          = spec$link,
        init          = init_val,
        init_supplied = init_user_supplied,
        lower         = bounds_val[1],
        upper         = bounds_val[2],
        prior_family  = pf,
        prior_p1      = pp1,
        prior_p2      = pp2,
        # A random-effect column carries its own group id (made unique per
        # process/param/stratum below in pool_linkages) and structure; a fixed
        # column inherits the spec's re_group (usually NA).
        re_group      = col_re_group[col_idx] %||% spec$re_group,
        re_struct     = col_re_struct[col_idx],
        re_time       = col_re_time[col_idx],
        re_sigma_init = s_init,
        re_sigma_prior_family = spf,
        re_sigma_prior_p1     = sp1,
        re_sigma_prior_p2     = sp2,
        re_rho_init   = r_init,
        re_rho_prior_family = rpf,
        re_rho_prior_p1     = rp1,
        re_rho_prior_p2     = rp2,
        re_obs_value  = r_obs_v,
        re_obs_sd     = r_obs_sd,
        est_phase     = spec$est_phase
      )
    }
  }
  out <- bind_linkage(rows)
  attr(out, "design_colnames") <- X_names
  attr(out, "design_matrix")   <- X
  out
}


#' Build the indicator design for random-effect (bar) terms
#'
#' @description
#' Each bar term `struct(<lhs> | <group>)` becomes one indicator column per
#' level of its grouping variable, so the deviation for that level enters the
#' offset the same way a fixed coefficient does. This increment supports the
#' unstructured / IID case (a plain `(1 | group)` bar, structure `"us"`);
#' `rw()` / `ar1()` and other structures are recognised by the parser but
#' rejected here until their densities are wired.
#'
#' @param re_terms list of bar expressions from [.parse_linkage_formula()].
#' @param re_structures character vector of structure names, one per term.
#' @param env_data the covariate/time table.
#'
#' @return A list with `X` (indicator design, 0 columns if no RE terms),
#'   `group` (character, one per column: the grouping variable name), and
#'   `struct` (character, one per column: the covariance structure).
#' @keywords internal
#' @noRd
.materialize_re_design <- function(re_terms, re_structures, env_data) {
  if (length(re_terms) == 0L) {
    return(list(X = matrix(numeric(0), nrow = nrow(env_data), ncol = 0L),
                group = character(0), struct = character(0),
                time = numeric(0)))
  }

  not_wired <- setdiff(unique(re_structures), c("us", "rw", "ar1"))
  if (length(not_wired) > 0) {
    stop(sprintf(
      paste0("random-effect structure(s) not yet wired: %s.\n",
             "  Supported: `(1 | group)` (IID), `rw(1 | group)` (random walk), ",
             "`ar1(1 | group)` (first-order autoregressive)."),
      paste0(not_wired, "()", collapse = ", ")), call. = FALSE)
  }

  cols <- list(); groups <- character(0); structs <- character(0)
  times <- numeric(0)
  for (i in seq_along(re_terms)) {
    # A bar expression `lhs | group`: the grouping variable is the RHS.
    bar <- re_terms[[i]]
    grp_var <- all.vars(bar[[3L]])
    if (length(grp_var) != 1L) {
      stop(sprintf("random-effect grouping must be a single variable; got `%s`",
                   deparse1(bar)), call. = FALSE)
    }
    if (!grp_var %in% names(env_data)) {
      stop(sprintf("random-effect grouping variable `%s` is not a column of env_data",
                   grp_var), call. = FALSE)
    }
    # Only random intercepts `struct(1 | group)` are supported: the grouping
    # variable indexes the deviations (one per level) and the structure runs
    # over them. A random slope (a non-`1` bar LHS, e.g. `ar1(Year + 0 | fleet)`)
    # is a different, unsupported model.
    if (length(all.vars(bar[[2L]])) > 0L) {
      stop(sprintf(
        paste0("random-effect term `%s` is a random slope; only random ",
               "intercepts `%s(1 | group)` are supported."),
        deparse1(bar), re_structures[i]), call. = FALSE)
    }
    # rw()/ar1() couple deviations by elapsed time, so the grouping variable
    # must be numeric (a real lag). IID is order-invariant, so it accepts any
    # grouping (e.g. a regime label).
    if (re_structures[i] != "us" &&
        !is.numeric(env_data[[grp_var]]) &&
        anyNA(suppressWarnings(as.numeric(as.character(env_data[[grp_var]]))))) {
      stop(sprintf(
        paste0("random-effect structure `%s(1 | %s)` needs a numeric grouping ",
               "variable so deviations can be ordered in real elapsed time; ",
               "`%s` is not numeric."),
        re_structures[i], grp_var, grp_var), call. = FALSE)
    }
    # One indicator column per level (no reference dropped: the density pins
    # them, so all levels are identifiable). model.matrix orders columns by
    # sorted factor level, which is numeric order for a numeric grouping.
    ind <- stats::model.matrix(
      stats::as.formula(sprintf("~ 0 + factor(%s)", grp_var)), data = env_data)
    lev <- sub(sprintf("^factor\\(%s\\)", grp_var), "", colnames(ind))
    colnames(ind) <- sprintf("%s_re::%s", grp_var, lev)
    # Numeric level value drives the real-elapsed-time ordering that rw()/ar1()
    # require (the numeric-Year rule). A non-numeric grouping (e.g. a regime
    # label) yields NA here; that is fine for IID (order-invariant) and is
    # rejected for rw()/ar1() where a true lag is needed.
    lev_num <- suppressWarnings(as.numeric(lev))
    cols[[i]] <- ind
    groups  <- c(groups,  rep(grp_var, ncol(ind)))
    structs <- c(structs, rep(re_structures[i], ncol(ind)))
    times   <- c(times,   lev_num)
  }
  list(X = do.call(cbind, cols), group = groups, struct = structs,
       time = times)
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
    for (i in seq_along(specs)) {
      param <- names(specs)[i]
      # Each value may be a single spec, or a list of specs that all
      # target the same parameter (e.g. species-specific formulas
      # registered together).
      this <- specs[[i]]
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
    table        = .assign_re_registry(bind_linkage(remapped)),
    X            = global_X,
    design_names = global_names
  )
}


#' Assign the random-effect registry on a pooled linkage table
#'
#' @description
#' Fills the `re_index` and `sigma_index` columns for the materialized
#' random-effect rows (those with a non-`NA` `re_struct`); fixed rows keep
#' `NA`. A distinct **sigma group** is one unique
#' `process|param|species|sex|age_bin|fleet|re_group|re_struct` combination --
#' the same key `map_linkage_adjuster()` uses, extended by the RE group and
#' structure so different fleets/params/groups each estimate their own variance.
#' `re_index` (the 0-based slot in `beta_linkage_re`) is assigned in
#' `(sigma group, elapsed time)` order, so each group's deviations are
#' contiguous and time-ordered -- `rw()`/`ar1()` rely on the ordering; IID is
#' order-invariant. The assignment is a bijection: `re_index` takes each value
#' in `0:(n_re - 1)` exactly once.
#'
#' @param tbl a pooled `Rceattle_linkage_table`.
#' @return `tbl` with `re_index`/`sigma_index` filled on RE rows.
#' @keywords internal
#' @noRd
.assign_re_registry <- function(tbl) {
  if (nrow(tbl) == 0L) return(tbl)
  is_re <- !is.na(tbl$re_struct)
  if (!any(is_re)) return(tbl)          # no RE rows: leave the registry NA

  re_rows <- which(is_re)
  key <- paste(tbl$process, tbl$param,
               ifelse(is.na(tbl$species), "*", tbl$species),
               ifelse(is.na(tbl$sex),     "*", tbl$sex),
               ifelse(is.na(tbl$age_bin), "*", tbl$age_bin),
               ifelse(is.na(tbl$fleet),   "*", tbl$fleet),
               tbl$re_group, tbl$re_struct, sep = "|")[re_rows]
  uniq_keys <- unique(key)
  tbl$sigma_index[re_rows] <- match(key, uniq_keys) - 1L

  ord <- order(tbl$sigma_index[re_rows], tbl$re_time[re_rows],
               method = "radix", na.last = TRUE)
  tbl$re_index[re_rows[ord]] <- seq_along(re_rows) - 1L
  tbl
}


#' Per-group random-effect SD (sigma) table
#'
#' @description
#' Collapses the per-row RE registry to one row per sigma group (ordered by
#' `sigma_index`, so row `i` is the 0-based group `i - 1`). The `sigma`-key
#' routing on `linkage_spec()` is identical across a group's rows, so this is a
#' plain de-duplication. Encodes the fix-vs-estimate contract:
#' * `init = list(sigma = v)` **and no** sigma prior -> `sigma_fixed = TRUE`,
#'   the SD is held at `v` (reproduces the reference `Time_varying_*_sd_prior`
#'   fixed-input SD);
#' * a sigma prior -> estimated with that prior (started from `init` if given);
#' * neither -> estimated from a default start.
#'
#' @param tbl a pooled `Rceattle_linkage_table` (post-registry).
#' @return a data.frame with one row per group, or `NULL` if there are none.
#' @keywords internal
#' @noRd
.re_group_table <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0L) return(NULL)
  re <- tbl[!is.na(tbl$sigma_index), , drop = FALSE]
  if (nrow(re) == 0L) return(NULL)
  g <- re[!duplicated(re$sigma_index), , drop = FALSE]
  g <- g[order(g$sigma_index), , drop = FALSE]
  has_prior <- !is.na(g$re_sigma_prior_family) & g$re_sigma_prior_family != "none"
  has_rprior <- !is.na(g$re_rho_prior_family) & g$re_rho_prior_family != "none"
  data.frame(
    sigma_index  = as.integer(g$sigma_index),
    re_struct    = as.character(g$re_struct),
    sigma_start  = ifelse(is.na(g$re_sigma_init), 0.3, as.numeric(g$re_sigma_init)),
    sigma_fixed  = !is.na(g$re_sigma_init) & !has_prior,
    prior_family = ifelse(has_prior, g$re_sigma_prior_family, "none"),
    prior_p1     = as.numeric(g$re_sigma_prior_p1),
    prior_p2     = as.numeric(g$re_sigma_prior_p2),
    # ar1 correlation, same fix-vs-estimate contract on the (-1,1) scale.
    rho_start    = ifelse(is.na(g$re_rho_init), 0, as.numeric(g$re_rho_init)),
    rho_fixed    = !is.na(g$re_rho_init) & !has_rprior,
    rho_prior_family = ifelse(has_rprior, g$re_rho_prior_family, "none"),
    rho_prior_p1 = as.numeric(g$re_rho_prior_p1),
    rho_prior_p2 = as.numeric(g$re_rho_prior_p2),
    # State-space (Rogers QAR1) observation: observed groups carry a fixed obs SD
    # and an estimated effect size (beta) applied to the deviate.
    observed     = !is.na(g$re_obs_sd),
    obs_sd       = as.numeric(g$re_obs_sd),
    stringsAsFactors = FALSE
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
