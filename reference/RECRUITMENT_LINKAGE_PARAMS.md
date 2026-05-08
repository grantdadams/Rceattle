# Allowed recruitment-parameter names for `linkages` in [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)

Linear-predictor names of the underlying recruitment parameters that the
linkage system can address. Linkages on `log_R0` are meaningful for any
`srr_fun` (the offset is added to the log of equilibrium / mean
recruitment); linkages on `log_alpha` and `log_beta` only do work when
the chosen `srr_fun` actually uses alpha / beta (Beverton-Holt, Ricker),
where they enter on the log scale before exponentiation.

## Usage

``` r
RECRUITMENT_LINKAGE_PARAMS
```

## Format

An object of class `character` of length 3.
