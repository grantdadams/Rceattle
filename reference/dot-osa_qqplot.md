# Q-Q plot of OSA residuals with standard-normal null envelope

Q-Q plot of OSA residuals with standard-normal null envelope

## Usage

``` r
.osa_qqplot(osa, nsim = 10000, seed = 123)
```

## Arguments

- osa:

  An `rceattle_osa` data frame with a `source` column.

- nsim, seed:

  Passed to the SDNR / tail-statistic annotation.

## Value

A `ggplot` object.
