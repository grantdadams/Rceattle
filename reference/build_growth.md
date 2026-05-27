# Specify the growth model for Rceattle

Specify the growth model for Rceattle

## Usage

``` r
build_growth(fun = "empirical", growth_age_L1 = NA, linkages = NULL)
```

## Arguments

- fun:

  Growth function. Either a string
  ([GROWTH_FUNS](https://grantdadams.github.io/Rceattle/reference/GROWTH_FUNS.md):
  `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
  equivalent integer code (`0`, `1`, `2`). The canonical string form is
  stored on the returned object.

- growth_age_L1:

  Von Bertalanffy / Richards anchor age (the age at which mean length
  equals `L1`). Matches SS3's `Growth_Age_for_L1` control input. Scalar
  (recycled across species) or a length-`nspp` vector for per-species
  values. Default `NA` inherits `data_list$growth_age_L1` if supplied
  (e.g. from the SS3 converter), otherwise falls back to
  `max(0.5, minage[sp])` so `minage >= 1` models stay
  backwards-compatible and `minage = 0` models pick up an SS3-consistent
  half-year anchor.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by parameter name (must be one of
  [GROWTH_LINKAGE_PARAMS](https://grantdadams.github.io/Rceattle/reference/GROWTH_LINKAGE_PARAMS.md)).
  The mean-growth keys (`K`, `L1`, `Linf`, `m`) accept arbitrary
  one-sided formulas and produce year-varying offsets applied inside
  `growth.hpp`. The SD-endpoint keys (`sd_L1`, `sd_Linf`) only honor
  intercept-bearing formulas (typically `~ 1`) – they thread `init`,
  `bounds`, and `priors` onto the underlying `growth_log_sd` parameter,
  giving the SDs the same prior/fix/initial-value contract as the mean
  parameters. Slope rows on SD specs raise a warning and have no effect;
  slope-only formulas (`~ 0 + temp`) error. Materialization into the
  global linkage table happens inside
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  once `data_list$env_data` is in scope.

## Value

A list of switches defining the growth model.

## Examples

``` r
if (FALSE) { # \dontrun{
# Sex-specific von Bertalanffy with temperature on K, by species + sex
build_growth(
  fun = "vonBertalanffy",   # or fun = 1
  linkages = list(
    K = linkage_spec(
      formula = ~ temp,
      by      = ~ species + sex,
      priors  = list(temp = normal(0, 1))
    )
  )
)
} # }
```
