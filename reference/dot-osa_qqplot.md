# Q-Q plot of OSA residuals with standard-normal null envelope

SDNR upper left, tail order statistics against their exact nulls lower
right, following `afscOSA`. The annotation names the nominal probability
each order statistic sits at, since it is not exactly 2.5%.

## Usage

``` r
.osa_qqplot(osa, add_sdnr_ci = TRUE, add_qq_quantiles = TRUE)
```

## Arguments

- osa:

  An `rceattle_osa` data frame with a `source` column.

- add_sdnr_ci:

  Show the chi-square null interval beside SDNR.

- add_qq_quantiles:

  Annotate the tail order statistics and their nulls.

## Value

A `ggplot` object.
