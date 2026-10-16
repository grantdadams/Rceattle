# Assessment reporting tables from one or more Rceattle fits

Collects the quantities a stock assessment reports into one set of tidy
tables, so a SAFE chapter or a model comparison comes from a single call
rather than a dozen ad-hoc extractions. Every table holds a `model`
column, so passing several fits gives a like-for-like comparison.

## Usage

``` r
report_tables(
  object,
  model_names = NULL,
  retro = NULL,
  jitter = NULL,
  osa = NULL,
  ci_level = 0.95,
  quantities = c("biomass", "ssb", "R", "biomass_depletion", "ssb_depletion", "F_spp")
)
```

## Arguments

- object:

  An Rceattle fit from
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md),
  or a list of them.

- model_names:

  Names for the models, defaulting to the names of `object` or to
  `Model_1`, `Model_2`, ...

- retro:

  A
  [`retrospective()`](https://afsc-assessments.github.io/Rceattle/reference/retrospective.md)
  result per model, or `NULL` (default) to omit that section.

- jitter:

  A
  [`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md)
  result per model, or `NULL` (default) to omit that section.

- osa:

  An
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
  result per model, summarized here with
  [`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md),
  or `NULL` (default) to omit that section.

- ci_level:

  Confidence level for the `timeseries` interval, default `0.95`.

- quantities:

  Time-series quantities to include, from the set
  [`as.data.frame.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/as.data.frame.Rceattle.md)
  accepts;
  [`quantity_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/quantity_dictionary.md)
  says what each one means.

## Value

A list of data frames with class `"rceattle_report"`, one per section
above, each holding a `model` column.

## Details

The sections follow the AFSC Alaska Groundfish Stock Assessment
Guidelines for what a chapter reports:

- `model`:

  One row per fit: dimensions, switches, the marginal and joint negative
  log-likelihoods, estimated parameter count, AIC, maximum gradient,
  whether the Hessian was positive definite, and run time.

- `parameters`:

  Every estimated parameter with its standard error and its
  natural-scale name and process from
  [`parameter_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/parameter_dictionary.md),
  which is where `sigma_R` and an estimated M are found. Estimates are
  on the parameter's own scale, so a `log_` name needs
  [`exp()`](https://rdrr.io/r/base/Log.html). A **fixed** M is absent,
  never having been estimated: read it off `M_at_age`.

- `likelihood`:

  Negative log-likelihood by component and fleet or species, weighted
  and unweighted.

- `timeseries`:

  Biomass, female spawning-stock biomass, recruitment, depletion and
  fishing mortality by species and year, with standard errors and
  intervals, split into hindcast (`era = "time"`) and projection
  (`era = "fore"`).

- `reference_points`:

  The executive-summary quantities: SPR-based F proxies, unfished and
  target female spawning-stock biomass, the biomass proxies implied by
  `Ptarget` / `Plimit`, and terminal status. A `basis` column says
  whether each was estimated, and if not why.

- `fits`:

  Observed against predicted index and catch, with the standard
  deviation of normalized residuals (SDNR) per fleet.

- `retrospective`, `jitter`, `osa`:

  Present only when the corresponding object is supplied.

Nothing here refits: the three diagnostics are tens to hundreds of
optimizations each, so they are passed in already computed, and a
section whose object is `NULL` is absent from the result.

The standard harvest scenarios of guideline section 4.11.3 are **not**
produced; they need a standard projection module, which Rceattle does
not have. Projected biomass under the model's own control rule is in
`timeseries` at `era = "fore"`.

**Reference points a fit never estimated come back as `NA`, with the
reason in `basis`.** The reported array leaves a plausible *number* in
each, which is how a wrong figure reaches an executive summary:

- `Ftarget` / `Flimit` are estimated only under a control rule that
  defines them, and are off for a species with no projected fishery or
  with fixed numbers-at-age. Unestimated, they sit at `exp(0) = 1`,
  reading as an F of 1.0/yr.

- `SB0` / `B0` under `msmMode > 0` come from the `MSSB0` / `MSB0`
  inputs, which hold a placeholder until
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  derives them from a no-fishing projection. `B_target` and `B_limit`
  are fractions of `SB0` and go with it.

- The per-recruit quantities (`SPR0`, `SPRtarget`, `SPRlimit`) are
  computed only under `msmMode = 0`.

The depletions are **not** blanked alongside `SB0`: under a no-fishing
control rule in multispecies mode the model divides by
last-projection-year biomass, the equilibrated unfished reference, so
the series is meaningful.

## Two negative log-likelihoods

`model` reports both, and **both are minimized**, so smaller is the
better fit. `marginal_nll` is what the optimizer minimized, random
effects integrated out by Laplace, and is what `AIC` is built from.
`joint_nll` is what the template evaluated at the conditional modes, so
it is what `likelihood` sums to, on the same scale as `jnll_comp`. They
are equal when `n_random` is 0 and differ by the Laplace correction
otherwise.

## Supplying diagnostics for several models

A diagnostics list is matched to models **by name**, so
`list(alt = ..., base = ...)` pairs correctly whatever the order. An
unnamed list is paired positionally, with a message. A name that is not
a model name is an error, which catches passing one model's
[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
result stored as a list of parts.

## See also

[`quantity_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/quantity_dictionary.md)
for what each quantity means and its units,
[`as.data.frame.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/as.data.frame.Rceattle.md)
for the time series alone, and
[`standard_output()`](https://afsc-assessments.github.io/Rceattle/reference/standard_output.md)
to relabel the result into the NOAA standardized assessment output that
`stockplotr` and `asar` consume.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_mod(data_list = BS2017SS, estimateMode = 1)
tabs <- report_tables(fit)
tabs$reference_points

# Compare two models, with diagnostics attached
tabs <- report_tables(list(base = fit0, alt = fit1),
                      retro = list(retro0, retro1))
} # }
```
