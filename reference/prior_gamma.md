# Gamma prior on a linkage coefficient

Shape-rate parameterization, matching
[`stats::dgamma()`](https://rdrr.io/r/stats/GammaDist.html).

## Usage

``` r
prior_gamma(shape, rate)
```

## Arguments

- shape:

  positive shape parameter.

- rate:

  positive rate parameter.

## Value

An `Rceattle_prior` of family `"gamma"`.
