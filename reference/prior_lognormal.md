# Lognormal prior on a linkage coefficient

A normal density on the log of the coefficient. Under
`fit_control(bias_adjust_proc = TRUE)`, the default, it is centred at
`meanlog - sdlog^2/2`, so `exp(meanlog)` is the prior mean; with `FALSE`
it matches [`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html)
and `exp(meanlog)` is the prior median.

## Usage

``` r
prior_lognormal(meanlog, sdlog)
```

## Arguments

- meanlog:

  log of the prior mean (of the median when `bias_adjust_proc = FALSE`).

- sdlog:

  prior standard deviation of the log (must be positive).

## Value

An `Rceattle_prior` of family `"lognormal"`.
