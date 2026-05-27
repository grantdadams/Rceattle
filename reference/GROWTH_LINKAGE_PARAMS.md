# Allowed growth-parameter names for `linkages` in [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)

Natural-scale names of the underlying growth-function parameters. Von
Bertalanffy uses `K`, `L1`, `Linf`; Richards adds `m`. `sd_L1` /
`sd_Linf` are the standard deviations of length-at-age anchored at `L1`
and `Linf` (the SD-at-age interpolation endpoints from `growth.hpp`).
Only intercept-only specs (`~ 1`) are honored on the SD endpoints – they
thread through `init` / `bounds` / `priors` onto `growth_log_sd` but do
not vary by year. The empirical weight-at-age model admits no linkages.

## Usage

``` r
GROWTH_LINKAGE_PARAMS
```

## Format

An object of class `character` of length 6.
