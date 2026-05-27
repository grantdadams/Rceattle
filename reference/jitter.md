# Jitter analysis

Run's the Rceattle model at initial values that are +- N(0, sd) from the
initial parameters.

## Usage

``` r
jitter(
  Rceattle = NULL,
  njitter = 50,
  sd = 0.2,
  phase = FALSE,
  seed = 123,
  cores = NULL
)
```

## Arguments

- Rceattle:

  an Rceattle model fit using
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)

- njitter:

  the number of jitters to run

- sd:

  standard deviation for jitter (default = 0.2)

- phase:

  as in
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  default = FALSE

- seed:

  random number seed. Each jitter `i` uses `seed + i` so results are
  reproducible under both sequential and parallel execution.

- cores:

  Number of cores to use for parallel jitters. Default `NULL` picks
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
#> 'Selectivity_dimension' not specified in 'fleet_control', assuming 'Age'
#> 'CAAL_weights' not specified in 'fleet_control', assuming 1
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
jitters <- jitter(ss_run, njitter = 10)
# }
```
