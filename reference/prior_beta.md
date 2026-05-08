# Beta prior on a linkage coefficient

Standard Beta(shape1, shape2) on (0, 1), matching
[`stats::dbeta()`](https://rdrr.io/r/stats/Beta.html). For priors on
stock-recruit steepness on the standard (0.2, 1) interval, transform
inside the model and use `prior_beta()` on the rescaled quantity.

## Usage

``` r
prior_beta(shape1, shape2)
```

## Arguments

- shape1, shape2:

  positive shape parameters.

## Value

An `Rceattle_prior` of family `"beta"`.
