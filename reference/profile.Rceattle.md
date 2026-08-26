# Likelihood profile across one or more parameter cells

Re-fits an Rceattle model while holding selected cells of a parameter
fixed at user-specified values. Supports profiling a single cell (e.g.
`R_log_sd[species = 1]`) and arbitrary N-dimensional cross-profiles over
multiple cells – e.g. `log_M1[1, 1, 1]` and `log_M1[1, 2, 1]` jointly,
to profile residual M for males against females. For each grid point the
targeted cells are fixed in the TMB map and the remaining parameters are
re-estimated; the result is a grid of Rceattle models for downstream NLL
surfaces.

## Usage

``` r
# S3 method for class 'Rceattle'
profile(
  fitted = NULL,
  param = NULL,
  slots = NULL,
  values = NULL,
  transform = "log",
  cores = NULL,
  getsd = NULL,
  ...
)
```

## Arguments

- fitted:

  an Rceattle model fit using
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)

- param:

  Name of the parameter to profile. Two ways to specify it:

  Raw parameter slot

  :   any name in `Rceattle$estimated_params`; tested for `"R_log_sd"`,
      `"rec_pars"`, and `"log_M1"`. `slots` must index into the full
      array and `transform` controls the scale.

  Natural-scale alias

  :   convenience shortcut for the three documented parameters. Aliases
      imply `transform = "log"` (values are taken in natural units and
      log'd before being substituted) and, for `rec_pars`, fill in the
      column from the alias name so `slots` only needs the species
      index:

      - `"sigmaR"`, `"R_sd"` -\> `R_log_sd`

      - `"M1"` -\> `log_M1`

      - `"R0"` -\> `rec_pars[, 1]`

      - `"alpha"` -\> `rec_pars[, 2]`

      - `"beta"` -\> `rec_pars[, 3]`

      If `transform` is supplied with an alias it is ignored (with a
      warning).

- slots:

  A list whose entries are integer index vectors, one entry per cell to
  fix. Each entry's length must equal the number of dimensions of the
  resolved parameter – 1 for vectors (`R_log_sd`), 2 for matrices
  (`rec_pars`), 3 for 3-D arrays (`log_M1`). When using the
  `"R0"`/`"alpha"`/`"beta"` aliases, supply only the species index
  (length 1); the column is filled in from the alias. E.g.
  `list(c(1, 2, 1))` fixes `log_M1[1, 2, 1]`;
  `list(c(1, 1, 1), c(1, 2, 1))` fixes both sex cells for a
  males-vs-females cross-profile of species 1; `list(1, 2)` with
  `param = "sigmaR"` cross-profiles species 1 and 2. If omitted,
  defaults to a single species-1 slot shaped to match the resolved
  parameter (e.g. `list(1)` for `R_log_sd`, `list(c(1, 1, 1))` for
  `log_M1`, `list(1)` for the `rec_pars` aliases) and emits a warning;
  pass `slots` explicitly to silence the warning. Defaulting requires
  `length(values) == 1L` (otherwise the user must explicitly say which
  cell each grid targets).

- values:

  A list of numeric vectors, one per entry of `slots`. The full grid of
  fits is `expand.grid(values)`, so a single slot gives a 1-D profile
  and *k* slots give a *k*-D cross-profile.

- transform:

  How to map user values onto the internal parameter scale before
  substituting them into `inits`. Either `"log"` (default),
  `"identity"`, or a unary function (e.g. `qlogis`). Applied
  element-wise to every grid value. Aliases override this with `"log"`.

- cores:

  Number of cores to use for parallel fits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

- getsd:

  whether each grid fit runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). The
  profile reads only the objective (`nll`), so `FALSE` is faster with no
  effect on the profile. Default `NULL` inherits the input model's
  setting (`TRUE` only if it carries an `sdrep`).

- ...:

  Unused; present for consistency with the
  [`stats::profile`](https://rdrr.io/r/stats/profile.html) generic.

## Value

A list with elements:

- Rceattle_list:

  list of fitted Rceattle models, one per grid row; entries for
  non-converged fits are `NULL` so positions stay aligned with `grid`.

- grid:

  data frame of grid values on the user scale (before `transform`); one
  column per profiled cell, named `slot_1`, `slot_2`, ...

- nll:

  numeric vector of joint negative log-likelihoods (`opt$objective`);
  `NA` where the fit did not converge.

- param:

  the profiled parameter name (echoed).

- slots:

  the slots list (echoed for downstream plotting).

- alias:

  the name you asked for, when it was one of the natural-scale aliases.
  Profiling `"M1"` returns `param = "log_M1"`, because the model
  estimates M on the log scale, while `grid` holds the M values you
  supplied. `alias` keeps `"M1"`, so a figure can label the axis in the
  units profiled rather than in log units. `NA` if you named the
  parameter slot directly.

Carries class `"Rceattle_profile"`, so printing it reports whether the
grid brackets the minimum; see
[`print.Rceattle_profile`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_profile.md).
Every element indexes exactly as before.

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

# 1-D profile of sigmaR for species 1 (alias form -- natural scale)
p1 <- profile(ss_run,
    param  = "sigmaR",
    slots  = list(1),
    values = list(seq(0.1, 1.5, by = 0.1)))

# Equivalent raw form (log scale -- user does the transform)
p1_raw <- profile(ss_run,
    param     = "R_log_sd",
    slots     = list(1),
    values    = list(log(seq(0.1, 1.5, by = 0.1))),
    transform = "identity")

# 2-D cross-profile of M1 across species 1 and 2 (sex 1, age 1).
# BS2017SS is single-sex; with a multi-sex model the same form
# (e.g. c(1, 1, 1), c(1, 2, 1)) would cross-profile males vs females.
p2 <- profile(ss_run,
    param  = "M1",
    slots  = list(c(1, 1, 1), c(2, 1, 1)),
    values = list(seq(0.1, 0.4, length.out = 3),
                  seq(0.1, 0.4, length.out = 3)))

# 1-D profile of SRR alpha for species 1 (alias drops the rec_pars column)
p3 <- profile(ss_run,
    param  = "alpha",
    slots  = list(1),
    values = list(seq(2, 80, length.out = 20)))
# }
```
