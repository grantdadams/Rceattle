# Self test simulation analysis analysis

Simulates data from an Rceattle model and refits the model to the
simulated data. TODO add process variation (i.e. random devs) to
simulation.

## Usage

``` r
self_test(
  Rceattle = NULL,
  nsim = 50,
  simulate = TRUE,
  seed = 123,
  cores = NULL
)
```

## Arguments

- Rceattle:

  an Rceattle model fit using
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)

- nsim:

  number of simulations

- simulate:

  passed to
  [`sim_mod`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md).
  If `TRUE` (default), data are simulated with observation error; if
  `FALSE`, expected values from the model are used.

- seed:

  random number seed. Each simulation `i` uses `seed + i` so results are
  reproducible under both sequential and parallel execution.

- cores:

  Number of cores to use for parallel simulations. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

## Value

a list of Rceattle models

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS,
    inits = NULL, file = NULL,
    estimateMode = 0, random_rec = FALSE,
    msmMode = 0, avgnMode = 0,
    phase = FALSE, verbose = 0)
#> Warning: Passing ‘phase’, ‘verbose’ directly to fit_mod() is deprecated and will be removed in a future release. Bundle these into fit_control() instead, e.g. fit_control(phase = ..., verbose = ...). Forwarding for now.
#> 'Diet_loglike' are not included in data, assuming 'Multinomial'
#> 'Sel_curve_pen3' not specified in 'fleet_control', assuming '0'
#> 'Selectivity_dimension' not specified in 'fleet_control', assuming 'Age'
#> 'CAAL_weights' not specified in 'fleet_control', assuming 1
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
sims <- self_test(ss_run, nsim = 10)
# }
```
