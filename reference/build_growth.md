# Specify the growth model for Rceattle

Specify the growth model for Rceattle

## Usage

``` r
build_growth(fun = "empirical", linkages = NULL)
```

## Arguments

- fun:

  Growth function. Either a string
  ([GROWTH_FUNS](https://grantdadams.github.io/Rceattle/reference/GROWTH_FUNS.md):
  `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
  equivalent integer code (`0`, `1`, `2`). The canonical string form is
  stored on the returned object.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by parameter name (must be one of
  [GROWTH_LINKAGE_PARAMS](https://grantdadams.github.io/Rceattle/reference/GROWTH_LINKAGE_PARAMS.md)).
  Each spec describes how that growth parameter depends on environmental
  covariates and on stratifying factors (species, sex). The parameter
  name on each spec is filled in from the list key. Materialization into
  the global linkage table happens inside
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
    log_K = linkage_spec(
      formula = ~ temp,
      by      = ~ species + sex,
      priors  = list(temp = normal(0, 1))
    )
  )
)
} # }
```
