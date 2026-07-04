# SDNR and tail statistics with standard-normal null intervals

SDNR and tail statistics with standard-normal null intervals

## Usage

``` r
.osa_sdnr_tails(resid, nsim = 10000, probs = c(0.025, 0.975), seed = 123)
```

## Arguments

- resid:

  Numeric vector of residuals (assumed standard normal under H0).

- nsim, probs, seed:

  See
  [`osa_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/osa_diagnostics.md).

## Value

A one-row data frame of statistics and their null intervals.
