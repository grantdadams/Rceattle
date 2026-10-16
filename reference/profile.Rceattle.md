# Likelihood profile across one or more parameter cells

Re-fits an Rceattle model while holding selected cells of a parameter
fixed. Profiles a single cell (`R_log_sd[species = 1]`) or crosses
several, such as `log_M1[1, 1, 1]` against `log_M1[1, 2, 1]` for
residual M in males versus females. For each grid point the targeted
cells are fixed in the TMB map and everything else is re-estimated.

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
  joint = c("none", "multiply", "add", "value"),
  ...
)
```

## Arguments

- fitted:

  an Rceattle model fit using
  [`fit_mod`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)

- param:

  Name of the parameter to profile: a raw slot in
  `Rceattle$estimated_params`, or one of the natural-scale aliases
  `"sigmaR"`, `"R_sd"`, `"M1"`, `"R0"`, `"alpha"`, `"beta"`, `"q"`. See
  Details.

- slots:

  A list of integer index vectors, one per cell to fix, each as long as
  the resolved parameter has dimensions: 1 for `R_log_sd`, 2 for
  `rec_pars`, 3 for `log_M1`. The `"R0"`/`"alpha"`/ `"beta"` aliases
  take the species index alone. So `list(c(1, 1, 1), c(1, 2, 1))`
  crosses the two sex cells of species 1. If omitted it defaults to a
  single species-1 cell shaped to the parameter and warns, which
  requires `length(values) == 1L`.

- values:

  A list of numeric vectors, one per entry of `slots`. The grid is
  `expand.grid(values)`, so *k* slots give a *k*-D cross-profile.

- transform:

  How user values map onto the internal scale before substitution:
  `"log"` (default), `"identity"`, or a unary function such as `qlogis`,
  applied element-wise.

- cores:

  Cores for the parallel fits. `NULL` (default) picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check`; 1
  forces sequential.

- getsd:

  Whether each grid fit runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). The
  profile reads only `nll`, so `FALSE` is faster and changes nothing.
  `NULL` (default) inherits the input model's setting (`TRUE` only if it
  holds an `sdrep`).

- joint:

  How the grid moves the cells in `slots`: `"none"` (default),
  `"multiply"`, `"add"` or `"value"`. See Details.

- ...:

  Unused; present for consistency with the
  [`stats::profile`](https://rdrr.io/r/stats/profile.html) generic.

## Value

A list of class `"Rceattle_profile"` with elements:

- Rceattle_list:

  the fitted models, one per grid row, `NULL` where a fit did not
  converge so positions stay aligned with `grid`.

- grid:

  grid values on the user scale, one column per profiled cell, named
  `slot_1`, `slot_2`, ...

- nll:

  joint negative log-likelihoods, `NA` where the fit did not converge.

- param:

  the profiled parameter name.

- slots:

  the slots list, echoed for plotting.

- alias:

  the natural-scale name asked for, or `NA`. Profiling `"M1"` returns
  `param = "log_M1"` because the model estimates M on the log scale,
  while `grid` holds the M values supplied, so a figure can label the
  axis in the units profiled.

- joint:

  the mode used, so a figure can say whether `grid` holds a multiplier,
  an offset or a value.

Printing it reports whether the grid brackets the minimum; see
[`print.Rceattle_profile`](https://afsc-assessments.github.io/Rceattle/reference/print.Rceattle_profile.md).

## Details

**Natural-scale aliases.** An alias implies `transform = "log"` and, for
`rec_pars`, fills in the column so `slots` needs only the species index.
`transform` is ignored with a warning when an alias is used. A raw slot
instead has `slots` index the full array and `transform` set the scale.

- `"sigmaR"`, `"R_sd"`:

  `R_log_sd`

- `"M1"`:

  `log_M1`

- `"R0"`, `"alpha"`, `"beta"`:

  `rec_pars` columns 1, 2 and 3

- `"q"`:

  `index_log_q`, base catchability, one entry per fleet; `slots` takes a
  fleet index or name. Fleets sharing a `Catchability_index` share one
  q, so naming any of them profiles the group: the rest are added with a
  message and `joint` becomes `"value"`.

**Joint modes.** A joint mode profiles a vector-valued parameter such as
an empirical age-based M, where `"none"` over ten ages and 13 values
would be `13^10` fits.

- `"none"`:

  (default) each slot gets its own grid and they cross.

- `"multiply"`:

  one grid multiplying every cell, so the schedule keeps its shape; `1`
  is the fitted model.

- `"add"`:

  one grid added to every cell, so the schedule keeps the differences
  between cells; `0` is the fitted model.

- `"value"`:

  one grid; every cell takes that value.

A joint mode needs `values` to be one vector, and `"multiply"`/`"add"`
work in natural units, so `transform` must be `"log"` or `"identity"`. A
value that would take a log-scale parameter to zero or below is an error
naming the cell, not a failed fit.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = FALSE, msmMode = 0, avgnMode = 0,
    phase = FALSE, verbose = 0)
#> Warning: Passing ‘phase’, ‘verbose’ directly to fit_mod() is deprecated and will be removed in a future release. Bundle these into fit_control() instead, e.g. fit_control(phase = ..., verbose = ...). Forwarding for now.
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s

# 1-D profile of sigmaR for species 1, natural scale via the alias
p1 <- profile(ss_run, param = "sigmaR", slots = list(1),
    values = list(seq(0.1, 1.5, by = 0.1)))

# 2-D cross-profile of M1 across species 1 and 2 (sex 1, age 1)
p2 <- profile(ss_run, param = "M1",
    slots  = list(c(1, 1, 1), c(2, 1, 1)),
    values = list(seq(0.1, 0.4, length.out = 3),
                  seq(0.1, 0.4, length.out = 3)))

# An empirical age-based M moved as a whole: one slot per age bin on ONE grid
# of multipliers, so the schedule keeps its shape and 1 is the fitted model.
sp <- 1
p3 <- profile(ss_run, param = "M1",
    slots  = lapply(seq_len(ss_run$data_list$nages[sp]),
                    function(a) c(sp, 1, a)),
    values = list(seq(0.6, 1.4, by = 0.05)),
    joint  = "multiply")
# }
# Profiling catchability by fleet name, and reading a profile for where the
# data disagree, are in vignette("model-diagnostics").
```
