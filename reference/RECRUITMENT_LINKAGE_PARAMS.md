# Allowed recruitment-parameter names for `linkages` in [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)

Natural-scale names of the underlying recruitment parameters that the
linkage system can address. Linkages on `R0` are meaningful for any
`srr_fun` (the offset is added to the log of equilibrium / mean
recruitment when the default log link is used); linkages on `alpha` and
`beta` only do work when the chosen `srr_fun` actually uses alpha / beta
(Beverton-Holt, Ricker).

## Usage

``` r
RECRUITMENT_LINKAGE_PARAMS
```

## Format

An object of class `character` of length 3.
