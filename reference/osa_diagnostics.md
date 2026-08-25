# Statistical diagnostics for OSA residuals

Computes the Stewart and Monnahan (2025) statistical diagnostics for a
set of OSA residuals: the standard deviation of the normalized residuals
(SDNR) and the lower/upper tail statistics, each with the 95% interval
expected under the standard normal null hypothesis (so departures can be
judged objectively rather than by eye). Computed per data source (type x
fleet) and overall.

Under a correctly specified model OSA residuals are already iid standard
normal, so the SDNR is simply their sample standard deviation. Its null
interval follows the chi-square result for the sample standard deviation
of `n` standard normals (Francis 2014); the tail-statistic null
intervals are obtained by simulation.

## Usage

``` r
osa_diagnostics(osa, nsim = 10000, probs = c(0.025, 0.975), seed = 123)
```

## Arguments

- osa:

  An `rceattle_osa` object from
  [`osa_residuals()`](https://grantdadams.github.io/Rceattle/reference/osa_residuals.md),
  or a data frame with `residual` and (optionally) `type`/`fleet`
  columns.

- nsim:

  Number of simulations for the tail-statistic null intervals. Default
  10000.

- probs:

  Lower/upper tail probabilities. Default `c(0.025, 0.975)`.

- seed:

  Seed for the tail-interval simulation (reproducibility). Default 123.

## Value

A data frame (class `"rceattle_osa_diagnostics"`, so it prints as a
compact severity-tagged summary; every column is still there and `$`
works as before) with one row per data source plus an `"all"` row, with
columns: `group` (the `"<source> fleet <n>"` label), `source`, `fleet`,
`n`, `sdnr`, `sdnr_lo`, `sdnr_hi`, `lower`, `lower_lo`, `lower_hi`,
`upper`, `upper_lo`, `upper_hi`, and the logical flags `sdnr_ok`,
`lower_ok`, `upper_ok` (TRUE when the statistic is inside its null
interval). On the `"all"` row `source` and `fleet` are `NA`.

## References

Francis, R.I.C.C. 2014. Replacing the multinomial in stock assessment
models: a first step. Fish. Res. 151:70-84.

Stewart, I.J., and Monnahan, C.C. 2025. Can. J. Fish. Aquat. Sci.
82:1-13.

## See also

[`osa_residuals()`](https://grantdadams.github.io/Rceattle/reference/osa_residuals.md)
