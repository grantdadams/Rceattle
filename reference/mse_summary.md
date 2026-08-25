# Management strategy evaluation performance metric summary

Management strategy evaluation performance metric summary

## Usage

``` r
mse_summary(mse, om_only = FALSE)
```

## Arguments

- mse:

  MSE runs from
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  or
  [`load_mse`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)

- om_only:

  If TRUE, report only operating-model (true) status and skip the
  estimation-model (EM) perception metrics.

## Value

A named list, one element per entity dimension so each metric lives only
where it applies (no NA padding):

- `species` – a data.frame with one row per species (keyed by `Species`)
  of the conservation / status metrics: per-species average catch, catch
  inter-annual variability (IAV), and P(Closed) = the probability the
  fishery is closed (catch ~ 0); the relative mean-squared error of
  terminal and average SSB; the estimation-model- and operating-model-
  perceived overfishing / overfished probabilities (via
  [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md))
  and the probability each status is misclassified (estimation model
  disagrees with the operating-model truth); and terminal biomass, SSB,
  dynamic SB0, SSB depletion (equilibrium and dynamic), average SSB
  depletion, and SSB-collapse counts.

- `fleet` – a data.frame with one row per fishery fleet (keyed by
  `Fleet_code` / `Fleet_name`) of average catch, catch IAV, and
  P(Closed).

- `total` – a named numeric of the across-fleet total average catch and
  catch IAV.

- `meta` – run provenance (`nsim`, `nspp`, `nflts`, `HCR`,
  projection-year range).

All metrics are averaged across projection years and simulations.

Metric columns carry syntactic names, so they can be typed without
backticks: `avg_catch`, `catch_iav`, `p_closed`, `ssb_rmse_avg`,
`ssb_rmse_terminal`, `em_p_overfishing`, `em_p_overfished`,
`om_p_overfishing`, `om_p_overfished`, `p_overfishing_false_pos`,
`p_overfishing_false_neg`, `p_overfished_false_pos`,
`p_overfished_false_neg`, `om_terminal_biomass`, `om_terminal_ssb`,
`om_terminal_dynamic_sb0`, `om_terminal_depletion`,
`om_terminal_depletion_dynamic`, `om_avg_depletion`,
`om_sims_collapsed`, `om_no_f_sims_collapsed`,
`om_sims_collapsed_from_f`. An `om_` prefix is the operating model's
truth and `em_` the estimation model's perception; the four
`*_false_pos` / `*_false_neg` metrics are the probability the two
disagree. The three `*_sims_collapsed` metrics are **counts of
simulations**, not probabilities.

`om_terminal_depletion` is `NA` for a multispecies run that derived no
unfished reference, which is any run without a harvest control rule
(`HCR = "NoFishing"`): under `msmMode > 0` the model reads spawning
biomass against the `MSSB0` input, and
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
only fills that in by projecting under no fishing when an HCR is
present. Dividing by the placeholder instead reported SSB/999 as a
depletion – on the Pacific hake three-species model, 2.68e3. Use
`om_terminal_depletion_dynamic`, which is computed against the model's
own `DynamicSB0` and is unaffected.

Each frame carries a `"labels"` attribute mapping those names to the
long display strings (e.g. `om_terminal_depletion_dynamic` -\>
`"OM: Terminal SSB Depletion (Dynamic)"`) for plots and tables:
`attr(summ$species, "labels")`.

## Examples

``` r
if (FALSE) { # \dontrun{
mse <- load_mse(dir = "mse_runs")
s <- mse_summary(mse)
s$species   # one row per species
s$fleet     # one row per fleet
s$total     # across the whole system
s$meta      # nsim, nspp, nflts, HCR
} # }
```
