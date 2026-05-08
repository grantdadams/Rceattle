# Capture a linkage specification

Capture a linkage specification

## Usage

``` r
linkage_spec(
  formula,
  param = NULL,
  data = NULL,
  by = ~species,
  species = NULL,
  sex = NULL,
  link = "identity",
  init = NULL,
  bounds = NULL,
  priors = NULL,
  re_group = NA_character_,
  est_phase = 1L
)
```

## Arguments

- formula:

  one-sided R formula whose RHS describes the linear predictor for
  `param` (e.g. `~ 1`, `~ temp`, `~ temp + PDO`).

- param:

  target parameter name on the linear predictor scale (e.g.
  `"log_alpha"`, `"log_M1"`, `"log_K"`). May be `NULL` when the spec is
  built inside a `build_*()` call that infers the parameter name from
  the enclosing list key (see
  [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)).

- data:

  (Optional) data frame for formula validation. Currently validation
  happens at materialization time inside
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).

- by:

  one-sided formula naming stratifying factors that should each get
  their own coefficients. Allowed names are `species`, `sex`, and
  `age_bin`. The default `~species` produces one coefficient set per
  species (the typical multispecies assessment use case); pass
  `~species + sex` for per-(species, sex) coefficients, or `NULL` to
  share a single coefficient set across every species/sex.

- species:

  optional integer vector of 1-based species ids that this spec applies
  to. `NULL` (default) means every species in `strata$species` at
  materialization time. Use this to give different species different
  formulas, e.g. by registering multiple specs against the same
  parameter – see
  [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  for the multi-spec syntax.

- sex:

  optional vector of sex ids that this spec applies to. May be supplied
  as integers (`1L` = female, `2L` = male) or as character strings
  (`"Females"`/`"Males"`, case-insensitive; `"female"`, `"male"`, `"f"`,
  `"m"` are also accepted). `NULL` (default) means every sex in
  `strata$sex` at materialization time. Only meaningful when `by`
  includes `sex`; otherwise the filter is a no-op. Use this to register
  separate specs per sex (e.g. one prior on females, another on males)
  against the same parameter.

- link:

  link function applied to the predictor when assembling process values;
  one of
  [LINKAGE_LINKS](https://grantdadams.github.io/Rceattle/reference/LINKAGE_PROCESSES.md).
  TODO

- init:

  optional named numeric vector of initial values keyed by the
  design-matrix column name (e.g. `c(`(Intercept)` = -1, temp = 0)`).
  Missing entries default to `0`.

- bounds:

  optional named list of `c(lower, upper)` keyed the same way as `init`.

- priors:

  optional named list whose entries are
  [Rceattle_priors](https://grantdadams.github.io/Rceattle/reference/Rceattle_priors.md)
  objects, keyed by design-matrix column name. Inside this argument the
  unprefixed shorthand `normal()`, `lognormal()`,
  [`gamma()`](https://rdrr.io/r/base/Special.html), and
  [`beta()`](https://rdrr.io/r/base/Special.html) resolves to the
  corresponding `prior_*` constructors via a private data mask, so
  `priors = list(temp = normal(0, 1))` works without masking
  [`base::gamma()`](https://rdrr.io/r/base/Special.html) or
  [`base::beta()`](https://rdrr.io/r/base/Special.html) at the package
  level. Equivalent to `priors = list(temp = prior_normal(0, 1))`.

- re_group:

  optional character: name of a random-effect grouping for these
  coefficients. `NA` (default) means fixed.

- est_phase:

  optional integer estimation phase. Default `1L`.

## Value

An `Rceattle_linkage_spec` object.
