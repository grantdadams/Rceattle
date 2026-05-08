# Lognormal prior on a linkage coefficient

Parameterized on the log scale (mean and sd of the log of the
coefficient), matching
[`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html).

## Usage

``` r
prior_lognormal(meanlog, sdlog)
```

## Arguments

- meanlog:

  prior mean of the log of the coefficient.

- sdlog:

  prior standard deviation of the log (must be positive).

## Value

An `Rceattle_prior` of family `"lognormal"`.
