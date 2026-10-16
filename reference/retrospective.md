# Retrospective peels

Calculate Mohn's rho and run retrospective peels for an Rceattle model,
and evaluate retrospective forecast skill. Each peel uses TMB's map to:

1.  filter the data and fixed inputs and map out time-varying parameters
    for the peeled years, holding them at that peel's terminal year;

2.  fit the peeled model;

3.  turn the hindcast parameters off and F for the peeled years on, and
    fit to the peeled catch series, updating the "forecast" dynamics
    given the projection assumptions and the observed catch.

## Usage

``` r
retrospective(
  object = NULL,
  peels = 5,
  rescale = FALSE,
  nyrs_forecast = 3,
  cores = NULL,
  getsd = NULL,
  phase = TRUE,
  fit_control = NULL,
  Rceattle = NULL
)
```

## Arguments

- object:

  an Rceattle model fit using
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- peels:

  number of retrospective peels used for rho and estimation

- rescale:

  TRUE/FALSE whether to subset and rescale environmental predictors over
  the peel years.

- nyrs_forecast:

  forecast years over which to calculate Mohn's rho, in addition to the
  terminal year

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- getsd:

  whether each peel runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html), which
  costs an extra model build per peel (see Details). Mohn's rho uses
  point estimates only, so `FALSE` is faster with no effect on it.
  `NULL` (default) inherits the input model's setting (`TRUE` only if it
  holds an `sdrep`), and the returned peels hold standard errors only
  when it is `TRUE`.

- phase:

  whether each peel is refitted in phases, default `TRUE`. A peel
  restarts from the unpeeled fit's starting values with a year removed,
  so without phasing the parameters barely move, the peels sit on top of
  the full model, and Mohn's rho is biased toward zero. Change it only
  deliberately.

- fit_control:

  optional
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  bundle for the refits. Only `phase` and `getsd` are read; see **What
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  reaches**.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

a list of 1. a list of Rceattle models and 2. Mohn's rho per species.

A peel that did not converge is dropped, so `Rceattle_list` can be
shorter than `peels + 1`, with a message. Entries are named for their
own terminal year (`Year_2017`, ...), so index by name:
`Rceattle_list[[3]]` is not necessarily the 3-year peel. With no peel
left, rho is `NaN` and the function warns. Each peel reports its own
terminal year as `data_list$endyr`, so plots draw it only as far as it
was fit.

A peel still estimates the years it dropped: that is its retrospective
forecast, fit to observed catch with the survey and composition data
withheld. Its recruitment deviation is the one
[`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md)
sets with `sample_rec = FALSE`, computed from the peel's own fit. Three
years matter, and every peel carries all three:

- `endyr`, `endyr_peel`:

  the peel's terminal year, what it was fit through. Equal to each
  other.

- `endyr_full`:

  the unpeeled model's terminal year, where the retrospective forecast
  ends.

- `projyr`:

  the end of the harvest-control-rule projection.

A peel whose standard errors could not be rebuilt carries
`attr(peel, "sd_report_error")`, a one-line reason naming the peel and
its terminal year; it is absent on every other peel.

Take the forecast years as
`endyr_peel + seq_len(endyr_full - endyr_peel)`, which is empty for the
unpeeled model, rather than `(endyr_peel + 1):endyr_full`, which counts
*down* there. Mohn's rho is computed from `endyr_peel` and is unaffected
by any of this.

Catchability is estimated only for a fleet holding fitted index rows
(see
[`build_map`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)),
and a peel moves `endyr`, so a survey whose index observations all fall
in the peeled-off years has no q estimated in that peel: the parameter
count is not constant across peels. That is deliberate, since a q with
no index to inform it is a flat direction in the likelihood. Rho is
computed from SSB and is unaffected, but `npar` and the reported
catchability do differ between a shallow and a deep peel for such a
fleet.

## Details

Each peel is fitted twice: a peeled hindcast, then a forecast refit
estimating only the peeled years' F. The second holds every hindcast
parameter fixed, so on its own it reports a standard error of zero
across the hindcast. Under `getsd = TRUE` the peel is therefore rebuilt
at those same parameters with the hindcast free in the map and reported
from there; nothing is re-estimated, so no point estimate moves. Where
that rebuild fails – it throws, or comes back with a Hessian it could
not factor, or with no finite standard error – the peel keeps its point
estimates and carries no `sdrep` at all, rather than the forecast
refit's zero standard errors. The reason is named in one warning per
call and carried on the peel;
[`print()`](https://rdrr.io/r/base/print.html) reports how many peels
could not report standard errors.

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
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
retro <- retrospective(ss_run, peels = 10)
# }
```
