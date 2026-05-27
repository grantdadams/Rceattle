# Retrospective peels

Calculate Mohn's rho and run retrospective peels for an Rceattle model.
The function also evaluates retrospective forecast skill. To evaluate
both retrospective bias and forecast skill, the function uses the map
functionality of TMB to peel the model:

1.  Filters data, filters fixed inputs, and maps out time-varying
    parameters for the peeled years. All time-varying parameters for the
    peeled years are set to the terminal year of the model for that
    peel.

2.  Fits the peeled model.

3.  Turns off all hindcast parameters, turns on F for the peeled years,
    and fits to the peeled catch series to update the "forecast"
    dynamics given projection assumptions and observed catch from the
    peeled years.

## Usage

``` r
retrospective(
  Rceattle = NULL,
  peels = 5,
  rescale = FALSE,
  nyrs_forecast = 3,
  cores = NULL
)
```

## Arguments

- Rceattle:

  an Rceattle model fit using
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)

- peels:

  the number of retrospective peels to use in the calculation of rho and
  for model estimation

- rescale:

  TRUE/FALSE whether to subset and rescale environmental predictors for
  the range of peel years.

- nyrs_forecast:

  Number of forecast years to calculate Mohn's Rho in addition to
  terminal year

- cores:

  Number of cores to use for parallel peels. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

## Value

a list of 1. list of Rceattle models and 2. vector of Mohn's rho for
each species

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
retro <- retrospective(ss_run, peels = 10)
# }
```
