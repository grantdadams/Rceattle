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
#' @param data Optional data frame for formula validation; validation is
#'   performed at fit time inside [fit_mod()].
#' @param by one-sided formula naming stratifying factors that should
#'   each get their own coefficients. Allowed names are `species`,
#'   `sex`, `age_bin`, and `fleet` (`fleet` for catchability and
#'   selectivity linkages). **When omitted, `by` defaults to the base
#'   stratum of whichever process the spec is attached to** -- `~fleet`
#'   for catchability, selectivity, and the fleet composition weights
#'   (`theta_comp` / `theta_caal`), and `~species` for recruitment, M,
#'   growth, and the diet weight (`theta_diet`) -- so you rarely need to
#'   spell it out for the base case. Pass it explicitly to override:
#'   e.g. `~species + sex` for per-(species, sex) coefficients, or
#'   `NULL` to share a single coefficient across every stratum. An
#'   explicit `by` (including `NULL`) is always kept as given.
#' @param species optional vector of species that this spec applies to,
#'   given either as 1-based species ids (`c(1L, 2L)`) or as species
#'   **names** matching `data_list$spnames` (`c("Pollock", "Cod")`).
#'   Names are matched exactly, after trimming whitespace, when the model
#'   is assembled in [fit_mod()]; an unrecognized name is an error that
#'   lists the model's species. Give ids or names, not a mix -- R coerces
#'   `c(1, "Cod")` to `c("1", "Cod")`. `NULL` (default) means every species in
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
#' @param fleet optional vector of fleets this spec applies to, given
#'   either as 1-based `Fleet_code`s (`c(1L, 3L)`) or as fleet **names**
#'   matching `fleet_control$Fleet_name` (`c("Shelikof", "Summer BT")`).
#'   Names are matched exactly, after trimming whitespace, when the model
#'   is assembled in [fit_mod()]; an unrecognized name -- or one that is
#'   not unique in `fleet_control` -- is an error that lists the model's
#'   fleets. Prefer names: a `Fleet_code` that is wrong but in range
#'   attaches the linkage to a different fleet and still fits, whereas a
#'   misspelled name cannot. Give ids or names, not a mix -- R coerces
#'   `c(7, "Pollock")` to `c("7", "Pollock")`. `NULL` (default) means every fleet in
#'   `strata$fleet` at materialization time. Only meaningful when `by`
#'   includes `fleet`; otherwise the filter is a no-op. Used by
#'   catchability and selectivity linkages to give different fleets
#'   different formulas.
#' @param link link function relating the linear predictor to the
#'   natural-scale target parameter. One of `"log"` (default) or
#'   `"identity"`. With `link = "log"`, `log(param) = X * beta` -- slope
#'   contributions are multiplicative on the natural-scale parameter. With `link = "identity"`,
#'   `param = X * beta` -- slope contributions are additive on the
#'   natural scale. The linkage targets are estimated on the log scale,
#'   so `"log"` is the default.
#' @param init optional named numeric vector of initial values keyed by
#'   the design-matrix column name (e.g.
#'   `c(`(Intercept)` = 4, temp = 0)`). Missing entries default to `0`.
#' @param bounds optional named list of `c(lower, upper)` keyed the same
#'   way as `init`.
#' @param priors optional named list of [Rceattle_priors] objects, keyed by
#'   design-matrix column name. Inside this argument you may write `normal()`,
#'   `lognormal()`, `gamma()`, or `beta()` directly, e.g.
#'   `priors = list(temp = normal(0, 1))` -- equivalent to
#'   `priors = list(temp = prior_normal(0, 1))`.
#' @param re_group optional character: name of a random-effect grouping
#'   for these coefficients. `NA` (default) means fixed.
#' @param est_phase optional integer estimation phase. Default `1L`.
#' @param observe optional character: for an `ar1(1 | group)` term, the name of
#'   an `env_data` column that measures the AR1 latent (a state-space covariate,
#'   sensu Rogers et al. 2024). The latent enters the linked parameter through
#'   an estimated effect size and is observed against this column. `NULL`
#'   (default) leaves the AR1 as a plain random effect. The latent is zero-mean
#'   (no estimated level), so standardize the observed covariate to mean 0; a
#'   non-zero-mean column confounds its level with the intercept and warns.
#' @param obs_sd optional positive numeric: the measurement SD for the `observe`
#'   covariate (one per observed group). Required with `observe`, unused
#'   otherwise. Held **fixed** at this value by default (`obs_sd_est = FALSE`); it
#'   is the *starting* value when `obs_sd_est = TRUE`.
#' @param obs_sd_est optional single `TRUE`/`FALSE` (default `FALSE`): estimate the
#'   `observe` measurement SD instead of holding it fixed, as the state-space
#'   survey-catchability (GOA pollock) model does. **Caveat:** the effect size and
#'   `obs_sd` are only jointly identified when the observed covariate is
#'   informative; on a smooth series the AR1 latent can track it exactly and the
#'   freely-estimated `obs_sd` collapses toward 0. Keep it fixed unless the
#'   covariate is informative. Only used with `observe`.
#'
#' @details The reserved keys `sigma` and `rho` in `init` / `priors` route the
#'   random-effect deviation SD and (for `ar1`) the correlation: e.g.
#'   `init = list(sigma = 0.1)` fixes the SD, `priors = list(rho = normal(0,
#'   0.3))` places a prior on the correlation. `sigma` means different things by
#'   structure: for `rw()` it is the innovation (per-step) SD; for `ar1()` it is
#'   the marginal (stationary) SD. The two are not directly comparable across
#'   structures -- see `vignette("environmental-linkages-and-priors")`.
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
                         obs_sd    = NULL,
                         obs_sd_est = FALSE) {
  # Whether `by` was left to its default. If so, the process that the spec is
  # later attached to (via build_*()) fills in its base stratum -- `~ fleet` for
  # catchability / selectivity / fleet composition weights, `~ species` otherwise
  # -- so users need not spell out `by = ~ fleet` / `by = ~ species` for the base
  # case. An explicitly-passed `by` (including `NULL` for a shared coefficient) is
  # always kept as given. Resolved in `.resolve_auto_by()`.
  by_auto <- missing(by)
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
  if (!is.logical(obs_sd_est) || length(obs_sd_est) != 1L || is.na(obs_sd_est)) {
    stop("`obs_sd_est` must be a single TRUE/FALSE.", call. = FALSE)
  }
  if (obs_sd_est && is.null(observe)) {
    stop("`obs_sd_est` is only used with `observe`.", call. = FALSE)
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
  # `species` / `fleet` take either 1-based ids or the model's own names. A name
  # cannot be resolved here -- linkage_spec() never sees a data_list, by design
  # (capture is split from materialization) -- so it is carried through as a
  # string and resolved in materialize_linkage() against the strata labels
  # fit_mod() attaches. Ids are coerced here.
  if (!is.null(species)) {
    species <- .coerce_stratum_arg(species, "species", "species ids")
  }
  if (!is.null(sex)) {
    sex <- .coerce_sex_arg(sex)
  }
  if (!is.null(fleet)) {
    fleet <- .coerce_stratum_arg(fleet, "fleet", "Fleet_code values")
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
      by_auto   = by_auto,
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
      obs_sd    = obs_sd,
      obs_sd_est = obs_sd_est
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
  Map(function(val, key) .stamp_param(val, key, process_label), linkages, names(linkages))
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


# The design-matrix column name model.matrix() gives the intercept. A constant
# because it is compared against in several places and is easy to mistype.
.INTERCEPT_COL <- "(Intercept)"


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


#' Accept `intercept` as a spelling of the design column `(Intercept)`
#'
#' `init` / `bounds` / `priors` are keyed by design-matrix column name. Every
#' covariate keys by a name the user wrote (`temp`), but the intercept keys by
#' the name `model.matrix()` gives it, which needs backticks. Both spellings are
#' accepted and stored canonically so lookups stay a single exact match.
#'
#' @keywords internal
#' @noRd
.canonical_intercept_names <- function(x, X_names) {
  if (is.null(x) || !length(x) || is.null(names(x))) return(x)
  if (!(.INTERCEPT_COL %in% X_names)) return(x)
  # Never translate a key that already names a design column: a covariate can be
  # called `intercept` (or `Intercept`), and there the key means that covariate.
  # Translating it would move the prior onto a different coefficient -- silently,
  # since both names are then valid. Keying by an existing column always wins.
  alias <- names(x) %in% c("intercept", "Intercept") & !(names(x) %in% X_names)
  if (any(alias)) {
    if (any(names(x) == .INTERCEPT_COL)) {
      stop("give the intercept prior/init/bounds once: both `",
           .INTERCEPT_COL, "` and `", names(x)[alias][1], "` were supplied.",
           call. = FALSE)
    }
    names(x)[alias] <- .INTERCEPT_COL
  }
  x
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

#' Coerce a `species` / `fleet` selector to ids-or-names.
#'
#' Unlike `sex` -- whose "Females"/"Males" mapping is universal and so can be
#' resolved on the spot -- a species or fleet name only means something relative
#' to a particular `data_list`. Character input is therefore validated for shape
#' and returned as-is, to be resolved later by
#' \code{.resolve_spec_strata_names()}. Numeric input is coerced to positive
#' 1-based ids.
#'
#' A factor is converted to its labels first: `as.integer()` on a factor
#' silently returns level codes, which would attach the linkage to whichever
#' fleet happens to sit at that position.
#'
#' @param x the user-supplied `species` / `fleet` argument (non-NULL).
#' @param arg argument name, for error messages.
#' @param id_label how to describe the integer form ("species ids").
#' @return Either a positive integer vector or a character vector of names.
#' @keywords internal
#' @noRd
.coerce_stratum_arg <- function(x, arg, id_label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    if (anyNA(x) || !all(nzchar(x))) {
      stop(sprintf("`%s` names must be non-empty, non-NA strings.", arg),
           call. = FALSE)
    }
    return(x)
  }
  x <- as.integer(x)
  if (anyNA(x) || any(x < 1L)) {
    stop(sprintf(paste0("`%s` must be a vector of positive 1-based %s, or the ",
                        "corresponding name(s)."), arg, id_label),
         call. = FALSE)
  }
  x
}

#' Resolve one character stratum selector against the model's own labels.
#'
#' Matching is exact after whitespace trimming (deliberately case-sensitive:
#' `Fleet_name` is user data, and case-folding could make two genuinely distinct
#' fleets collide). Both failure modes name the offender and print the valid set,
#' because the alternative -- an id that is wrong but in range -- attaches the
#' linkage to a different fleet and still fits.
#'
#' @param x the selector; returned untouched unless it is character.
#' @param labels the stratum labels, or `NULL` when the caller supplied none.
#' @param arg argument name, for error messages.
#' @param source_label where the labels come from, for error messages.
#' @return A positive integer vector of 1-based ids.
#' @keywords internal
#' @noRd
.resolve_stratum_names <- function(x, labels, arg, source_label) {
  if (is.null(x) || !is.character(x)) return(x)
  # A blank or NA label is unmatchable (`NA == k` is NA, which `which()`
  # drops), so only a wholly unusable label set is fatal -- one missing name
  # must not stop the other strata being selected by theirs.
  usable <- !is.null(labels) && any(!is.na(labels) & nzchar(labels))
  if (!usable) {
    stop(sprintf(paste0(
      "`%s` was given as name(s) (%s), but this model carries no %s to match ",
      "against. Supply 1-based ids instead."),
      arg, paste(shQuote(x), collapse = ", "), source_label), call. = FALSE)
  }
  hits <- lapply(x, function(k) which(trimws(labels) == k))

  bad <- x[lengths(hits) == 0L]
  if (length(bad) > 0L) {
    # `c(7L, "Pollock")` is coerced by R to `c("7", "Pollock")` before it ever
    # reaches us, so a whole-number "name" almost always means the caller mixed
    # ids and names in one vector. Say so -- otherwise the error reads as though
    # a fleet really is called "7".
    hint <- if (any(grepl("^[0-9]+$", bad))) sprintf(paste0(
      "\n  %s looks like an id, not a name: `%s` takes ids or names but not ",
      "both, and R coerces a mixed c(7, \"%s\") to c(\"7\", \"%s\"). Use all ",
      "ids or all names."),
      paste(shQuote(unique(bad[grepl("^[0-9]+$", bad)])), collapse = ", "),
      arg, labels[1], labels[1]) else ""
    stop(sprintf("unknown `%s` name(s): %s.\n  This model's %s: %s.%s",
                 arg, paste(shQuote(unique(bad)), collapse = ", "),
                 source_label, paste(shQuote(labels), collapse = ", "), hint),
         call. = FALSE)
  }
  dup <- x[lengths(hits) > 1L]
  if (length(dup) > 0L) {
    stop(sprintf(paste0(
      "`%s` name(s) %s match more than one entry in the model's %s, so they ",
      "cannot identify a single stratum. Make the names unique, or select by ",
      "1-based id."),
      arg, paste(shQuote(unique(dup)), collapse = ", "), source_label),
      call. = FALSE)
  }
  as.integer(unlist(hits))
}

#' Resolve a spec's `species` / `fleet` names against the strata labels.
#'
#' @param spec an `Rceattle_linkage_spec`.
#' @param strata the `strata` list passed to [materialize_linkage()]; its
#'   `species` / `fleet` elements may carry names (see [pool_linkages()]).
#' @return `spec`, with `species` / `fleet` as 1-based integer ids.
#' @keywords internal
#' @noRd
.resolve_spec_strata_names <- function(spec, strata) {
  spec$species <- .resolve_stratum_names(
    spec$species, names(strata$species), "species",
    "species names (data_list$spnames)")
  spec$fleet <- .resolve_stratum_names(
    spec$fleet, names(strata$fleet), "fleet",
    "fleet names (fleet_control$Fleet_name)")
  spec
}

#' Drop the label names from atomic strata vectors.
#'
#' Names on `strata$species` / `strata$fleet` are labels for name resolution
#' only; left in place they ride through `expand.grid()` into the level grid and
#' out onto the materialized table's id columns. Species-keyed *lists*
#' (`strata$sex`, `strata$age_bin`) keep their names, which are load-bearing.
#'
#' @param strata the `strata` list.
#' @return `strata`, with names stripped from its atomic elements.
#' @keywords internal
#' @noRd
.strata_drop_labels <- function(strata) {
  lapply(strata, function(v) if (is.list(v)) v else unname(v))
}

#' Attach the model's labels to a stratum-level vector.
#'
#' Labels are what let a spec select a species / fleet by name. Missing or
#' wrong-length labels are dropped rather than recycled: an unlabelled stratum
#' gives a clear "this model carries no names" error at resolution time, whereas
#' a recycled label would resolve to the wrong id.
#'
#' @param ids 1-based integer ids.
#' @param labels the corresponding names, or `NULL`.
#' @return `ids`, named by `labels` when usable.
#' @keywords internal
#' @noRd
.label_strata <- function(ids, labels) {
  if (is.null(labels) || length(labels) != length(ids)) return(ids)
  stats::setNames(ids, as.character(labels))
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
  by_txt <- if (is.null(x$by)) {
    "(shared)"
  } else if (isTRUE(x$by_auto) && is.na(x$param)) {
    "(process default)"            # unattached: resolves to ~fleet / ~species per process
  } else if (isTRUE(x$by_auto)) {
    paste0(deparse(x$by), " (default)")
  } else {
    deparse(x$by)
  }
  cat("  by:      ", by_txt, "\n", sep = "")
  if (!is.null(x$species)) {
    cat("  species: ", paste(x$species, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$sex)) {
    cat("  sex:     ", paste(x$sex, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$fleet)) {
    cat("  fleet:   ", paste(x$fleet, collapse = ", "), "\n", sep = "")
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
#' The `X_col` column starts as a *local* column index into the
#' per-spec design matrix; the caller (`fit_mod()`'s pooling step)
#' remaps these to global column indices once all specs have been
#' combined into a single shared `X` matrix.
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
#'   Allowed names are `"species"`, `"sex"`, `"age_bin"`, and `"fleet"`.
#'   The `species` and `fleet` vectors may additionally be *named* with the
#'   model's own labels (`data_list$spnames` /
#'   `fleet_control$Fleet_name`); those labels are what a spec built with
#'   `linkage_spec(fleet = "Shelikof")` is resolved against. Without them,
#'   such a spec errors -- ids always work.
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
  # Fixed-effect covariates must be finite over the model range: model.matrix()
  # silently DROPS rows whose covariate is NA (e.g. a fixed-effect covariate with
  # missing years after env_data extension), which would misalign X_fixed against
  # the RE design and error cryptically in the cbind below. Reject a missing year
  # up front with a clear message. A state-space `observe` covariate is NOT a
  # fixed term, so it is not checked here (it is masked in the observation).
  fixed_vars <- intersect(all.vars(parsed$fixed), names(env_data))
  na_fixed <- fixed_vars[vapply(fixed_vars,
                                function(v) anyNA(env_data[[v]]), logical(1))]
  if (length(na_fixed) > 0) {
    stop(sprintf(paste0(
      "fixed-effect linkage covariate(s) have missing (NA) values over the ",
      "model years: %s. Provide the covariate for every year styr:endyr (a ",
      "state-space `observe` covariate may be partial, but a fixed-effect ",
      "covariate may not)."), paste(na_fixed, collapse = ", ")), call. = FALSE)
  }
  X_fixed <- stats::model.matrix(parsed$fixed, data = env_data)
  re <- .materialize_re_design(parsed$re_terms, parsed$re_structures, env_data)

  X <- cbind(X_fixed, re$X)
  X_names <- colnames(X)
  n_cols  <- ncol(X)

  # Translate the backtick-free intercept spelling now that the design columns
  # are known (see .canonical_intercept_names: a covariate named `intercept`
  # keeps its own meaning).
  spec$init   <- .canonical_intercept_names(spec$init,   X_names)
  spec$bounds <- .canonical_intercept_names(spec$bounds, X_names)
  spec$priors <- .canonical_intercept_names(spec$priors, X_names)

  # `init` / `bounds` / `priors` are looked up by exact design-column name, so a
  # key that matches no column is simply never read -- the prior or starting
  # value is dropped with no sign that it was. Reject those here instead: an
  # unapplied prior is indistinguishable from an unpenalised model, and on the
  # models where a prior is load-bearing that is the difference between a fit
  # that converges and one that runs to a boundary.
  # `sigma` is read only for random-effect columns and `rho` only for ar1 ones,
  # so accepting them on a spec without those terms would reinstate exactly the
  # silent drop this check exists to close.
  .reserved <- c(if (length(re$struct)) "sigma",
                 if (any(re$struct == "ar1")) "rho")
  for (.arg in c("init", "bounds", "priors")) {
    .given <- names(spec[[.arg]])
    .bad   <- setdiff(.given, c(X_names, .reserved))
    if (length(.bad)) {
      stop(sprintf(
        paste0("unknown `%s` name(s) for this linkage: %s.\n",
               "  This formula's design columns: %s.%s"),
        .arg, paste(sQuote(.bad), collapse = ", "),
        paste(sQuote(X_names), collapse = ", "),
        if (any(.bad %in% c("intercept", "Intercept")))
          paste0("\n  This formula has no intercept to key to -- `~ 0 + ...` and ",
                 "`~ ... - 1` suppress it. Drop the `0 +`/`- 1` to add one.")
        else ""), call. = FALSE)
    }
  }
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
    # Latent ar1 is zero-mean (no estimated level), so a non-zero-mean covariate
    # confounds its mean with the linked parameter's intercept. Warn; standardize.
    obs_present <- env_data[[spec$observe]][is.finite(env_data[[spec$observe]])]
    if (length(obs_present) > 1L &&
        abs(mean(obs_present)) > 0.1 * stats::sd(obs_present)) {
      warning(sprintf(
        paste0("`observe` covariate `%s` is not ~zero-mean (mean %.3g); the ",
               "latent ar1 state carries no estimated level, so its mean is ",
               "confounded with the linked parameter's intercept. Standardize ",
               "the covariate (mean 0) unless this is intended."),
        spec$observe, mean(obs_present)), call. = FALSE)
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

  # `species` / `fleet` may have been given as names rather than ids. Resolve
  # them here -- this is the first point that knows the model -- so everything
  # below, and every id column of the materialized table, stays integer.
  # Resolution runs even when the corresponding filter is a no-op (`by` does not
  # stratify by that factor), so a misspelled name is always an error rather
  # than a silently ignored one.
  spec   <- .resolve_spec_strata_names(spec, strata)
  strata <- .strata_drop_labels(strata)

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
            r_obs_v <- o_val; r_obs_sd <- as.numeric(spec$obs_sd); r_obs_est <- isTRUE(spec$obs_sd_est)
          } else {
            r_obs_v <- NA_real_; r_obs_sd <- NA_real_; r_obs_est <- NA
          }
        } else {
          r_init <- NA_real_; rpf <- NA_character_; rp1 <- NA_real_; rp2 <- NA_real_
          r_obs_v <- NA_real_; r_obs_sd <- NA_real_; r_obs_est <- NA
        }
      } else {
        s_init <- NA_real_; spf <- NA_character_; sp1 <- NA_real_; sp2 <- NA_real_
        r_init <- NA_real_; rpf <- NA_character_; rp1 <- NA_real_; rp2 <- NA_real_
        r_obs_v <- NA_real_; r_obs_sd <- NA_real_; r_obs_est <- NA
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
        re_obs_est    = r_obs_est,
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
    # Numeric level value drives the elapsed-time ordering that rw()/ar1()
    # require (the numeric-Year rule). A non-numeric grouping (e.g. a regime
    # label) yields NA here; that is fine for IID (order-invariant) and is
    # rejected for rw()/ar1() where a true lag is needed.
    lev_num <- suppressWarnings(as.numeric(lev))
    # rw()/ar1() couple CONSECUTIVE deviations with a unit lag, so the time
    # points must be EQUALLY SPACED. A gap (e.g. a missing year) would be
    # silently treated as a single step -- wrong variance (rw) / wrong
    # correlation (ar1). Require equal spacing rather than mis-specify; a
    # gap-aware (continuous-time) density is future work.
    if (re_structures[i] != "us" && !anyNA(lev_num)) {
      ts <- sort(unique(lev_num))
      if (length(ts) >= 3L) {
        dts <- diff(ts)
        if (any(abs(dts - dts[1]) > 1e-8 * max(1, abs(dts[1])))) {
          stop(sprintf(
            paste0("random-effect structure `%s(1 | %s)` requires equally-spaced ",
                   "time points (it couples consecutive deviations with a unit ",
                   "lag); `%s` has gaps/irregular spacing. Fill the missing ",
                   "levels or use `(1 | %s)` (IID)."),
            re_structures[i], grp_var, grp_var, grp_var), call. = FALSE)
        }
      }
    }
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
#' eliminates every row of the level grid. The design matrix stays in
#' the attributes so the pooler can union columns across specs that did
#' emit rows.
#'
#' @keywords internal
#' @noRd
.empty_materialized <- function(X, X_names) {
  out <- new_linkage_table()
  attr(out, "design_colnames") <- X_names
  attr(out, "design_matrix")   <- X
  out
}


#' Extend env_data to start at the model start year (prepend + gap-fill with NA)
#'
#' Linkages align env_data POSITIONALLY (row r -> model year styr + r - 1). A
#' state-space `observe` covariate (Rogers QAR1) commonly starts later than the
#' model (e.g. a survey that began mid-series), which would misalign every row.
#' Rather than force the user to hand-pad, prepend the missing leading years
#' (styr .. first year - 1) and fill any interior gaps, with `NA` in the
#' covariate columns. The QAR1 observation is applied only where the covariate is
#' present (NA years are masked out), so no observation is fabricated. Years
#' BEYOND the last supplied row are left alone (they keep the existing
#' zero-offset-beyond-env_data behaviour). Malformed year columns (unsorted, NA,
#' or starting before styr) are left untouched for `.check_env_data_years()` to
#' reject.
#'
#' @param env_data the covariate/time table (or `NULL`).
#' @param styr the model start year.
#' @return the (possibly extended) env_data.
#' @keywords internal
#' @noRd
.extend_env_data <- function(env_data, styr) {
  if (is.null(env_data) || !is.data.frame(env_data) ||
      !"Year" %in% names(env_data) || is.null(styr)) {
    return(env_data)
  }
  yrs <- env_data$Year
  if (anyNA(yrs) || is.unsorted(yrs, strictly = TRUE) || yrs[1] < styr) {
    return(env_data)
  }
  full <- styr:max(yrs)
  if (identical(as.integer(yrs), as.integer(full))) return(env_data)  # already contiguous from styr
  out <- merge(data.frame(Year = full), env_data, by = "Year",
               all.x = TRUE, sort = TRUE)
  message(sprintf(paste0(
    "env_data extended to start at styr (%s) with NA for missing years; a ",
    "state-space `observe` covariate is used only where present."), styr))
  out
}


#' Validate that env_data aligns to the model years for linkages
#'
#' @description
#' Linkages consume `env_data` positionally: row `r` supplies the offset for
#' model year `styr + r - 1`, and years beyond the last row get a zero offset
#' (env_data need not span the projection horizon). A misaligned `env_data`
#' therefore applies a covariate or random-effect deviate to the wrong year. If
#' `env_data` carries a `Year` column, require it to be sorted, start at `styr`,
#' and be contiguous (no gaps), erroring loudly otherwise. Without a `Year`
#' column the positional contract cannot be checked and is assumed.
#'
#' @param env_data the covariate/time table (or `NULL`).
#' @param styr the model start year.
#' @return invisibly `NULL`; errors on misalignment.
#' @keywords internal
#' @noRd
.check_env_data_years <- function(env_data, styr) {
  if (is.null(env_data) || !is.data.frame(env_data) ||
      !"Year" %in% names(env_data) || is.null(styr)) {
    return(invisible())
  }
  yrs <- env_data$Year
  if (anyNA(yrs) || is.unsorted(yrs, strictly = TRUE)) {
    stop("env_data$Year must be sorted ascending with no duplicates or NA; ",
         "linkages align each row to model year styr + (row - 1).",
         call. = FALSE)
  }
  if (yrs[1] != styr) {
    stop(sprintf(paste0(
      "env_data$Year must start at the model start year styr (%s); it starts ",
      "at %s. Linkages align by row position, so env_data must begin at styr."),
      styr, yrs[1]), call. = FALSE)
  }
  if (length(yrs) > 1L && any(diff(yrs) != 1L)) {
    bad <- yrs[which(diff(yrs) != 1L) + 1L]
    stop(sprintf(paste0(
      "env_data$Year has gaps (missing year(s) before %s); linkages align by ",
      "row = styr + offset, so the years must be contiguous. Fill the missing ",
      "rows (env_data may still stop short of the projection horizon)."),
      paste(bad, collapse = ", ")), call. = FALSE)
  }
  invisible()
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
#'   `list(species = 1:nspp, sex = 1:nsex)`. `fit_mod()` names the `species`
#'   and `fleet` vectors with `data_list$spnames` /
#'   `fleet_control$Fleet_name`, which is how a spec that selects a species or
#'   fleet by name is resolved to an id.
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
    obs_sd_est   = !is.na(g$re_obs_est) & g$re_obs_est,
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
