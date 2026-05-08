# Observed-vs-fitted residuals from an Rceattle fit

Returns a long-format data frame of residuals across one or more of the
four fitted data sources: `"index"` (survey indices), `"catch"` (fishery
catches), `"comp"` (age- or length-composition proportions), and
`"caal"` (conditional age-at-length proportions).

## Usage

``` r
# S3 method for class 'Rceattle'
residuals(object, type = "index", scale = "log", ...)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).

- type:

  One or more of `"index"`, `"catch"`, `"comp"`, `"caal"`, or `"all"`
  (default `"index"`).

- scale:

  `"log"` (default) or `"natural"`. Only affects `"index"` and `"catch"`
  residuals.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `Source`, `Fleet_code`, `Fleet_name`,
`Species`, `Sex`, `Year`, `Length`, `Bin`, `Age0_Length1`,
`Sample_size`, `Observed`, `Fitted`, `Residual`. `Sex`, `Length`, `Bin`,
`Age0_Length1`, and `Sample_size` are `NA` where they do not apply (e.g.
for index/catch rows).

## Details

For `"index"` and `"catch"`, the `Residual` column is on the log scale
by default (matching the lognormal observation likelihood) and can be
switched to the natural scale via `scale = "natural"`. For `"comp"` and
`"caal"`, residuals are Pearson residuals on the fitted proportions:
\$\$r = (p - \hat p)/\sqrt{\hat p (1 - \hat p)/N}\$\$ where N is the
input sample size. Composition rows are returned in long form: one row
per (observation, age/length bin).

Composition rows carry the `Age0_Length1` flag from `comp_data` (`0` for
age comps, `1` for length comps) so age and length comps can be filtered
apart. CAAL rows carry both the conditioning `Length` and the age `Bin`.
