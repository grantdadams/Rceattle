# 0. Model options and functionality

This vignette is an overview of the parameterization choices and
capabilities available in `Rceattle`. It complements:

- vignette 1, [Rceattle: An
  Introduction](https://grantdadams.github.io/Rceattle/articles/introduction.md),
  which walks through fitting a model end-to-end;
- vignette 8, [Model
  parameterizations](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.md),
  the equation-level reference (in progress); and
- the per-topic vignettes 3 (diagnostics), 4 (projections / reference
  points), 5 (single- vs. multi-species), 6 (building a data object in
  R), and 7 (Stock Synthesis migration).

If you’re building a new application, this is the menu of options you
choose between when configuring `data_list` and the arguments to
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).

## 1. Model dimensions and operating mode

> Equations: see vignette 8 [§Population
> dynamics](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.html#population-dynamics).

The top-level switches that decide the observation model live on the
`data_list` and the population model and estimation switches live on the
call to
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md):

| Switch | Where set | Role |
|----|----|----|
| `nspp`, `nages`, `nsex`, `nlengths` | `data_list` | Population structure |
| `styr`, `endyr`, `projyr` | `data_list` | Time horizon (hindcast + projection) |
| `msmMode` | [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md) | Single- vs. multi-species |
| `initMode` | [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md) | Initial age structure |
| `estimateMode` | [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md) | Estimate / project / debug |
| `random_rec`, `random_q`, `random_sel` | [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md) | Random-effects toggles |

### `msmMode` — single-species vs. multi-species

| Value | Meaning | Status |
|:--:|----|----|
| `0` | Single-species (M = M1 only) | Validated |
| `1` | Multi-species, MSVPA Type II predation | Validated |
| `2` | Multi-species, MSVPA Type III predation | Validated |
| `3:9` | Kinzey & Punt (2009) functional responses (Holling I/II/III, predator interference, predator preemption, Hassell-Varley, Ecosim) | **Blocked** — [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md) will refuse these until they are validated against the current parameter set |

### `initMode` — initial age structure

Codes (also accept the string aliases listed):

| Value | Alias | Description |
|:--:|----|----|
| `0` | `"FreeParams"` | Initial N-at-age estimated as free parameters |
| `1` | `"Equilibrium"` | Unfished equilibrium, no init devs, F_(init) = 0 |
| `2` (default) | `"NonEquilibrium"` | Unfished equilibrium with initial deviations |
| `3` | `"FishedNonEquilibrium"` | Fished non-equilibrium with init devs; F_(init) is added to M inside geometric series |
| `4` | `"FishedNonEquilibriumScaled"` | Fished non-equilibrium with init devs and R₀-scaled by F_(init) |
| `5` | `"OffsetEquilibrium"` | Unfished equilibrium, F_(init) = 0, displaced by the year-1 recruitment deviation instead of sitting at R_(init); init devs off, no init-dev penalty (Cole Monnahan / AFSC GOA pollock convention) |

### `estimateMode` — what does `fit_mod()` actually do?

| Value | Behavior |
|:--:|----|
| `0` | Estimate hindcast and projection parameters via `nlminb` / [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) (default). |
| `1` | Estimate only the hindcast parameters via `nlminb` / [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). |
| `2` | Hold the hindcast fixed and run only the projection under the supplied `HCR`. |
| `3` | Evaluate the report once at `inits` (debug / mapping check). |
| `4` | Run through `MakeADFun` and `nlminb` with all parameters mapped out (debug / mapping check). |

### Random-effects toggles

`random_rec`, `random_q`, and `random_sel` move log-recruitment, the
log-catchability deviations, and the time-varying selectivity
deviations, respectively, from penalized fixed effects to random
effects. The deviations are marginalized via the Laplace approximation,
and their distribution parameters are estimated.

## 2. Predation: `suitMode` (per predator)

> Equations: see vignette 8 [§Predation
> sub-model](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.html#predation-sub-model)
> (suitability, predation mortality, consumption / bioenergetics).

`suitMode` chooses how prey preference is modeled. Length 1 (same for
all predators) or length `nspp`.

| Value | String | Description | Available? |
|:--:|----|----|:--:|
| `0` | `"Empirical"` | Empirical proportions from observed stomach data (MSVPA-style; Holsman et al. 2015) | Yes |
| `1` | `"GammaLength"` | Gamma suitability on prey length-at-age | **Not yet** |
| `2` | `"GammaWeight"` | Gamma suitability on prey weight-at-age | Yes |
| `3` | `"LognormalLength"` | Log-normal suitability on prey length-at-age | **Not yet** |
| `4` | `"LognormalWeight"` | Log-normal suitability on prey weight-at-age | Yes |
| `5` | `"NormalLength"` | Normal suitability on prey length-at-age | **Not yet** |
| `6` | `"NormalWeight"` | Normal suitability on prey weight-at-age | Yes |

Suitability is averaged across years `suit_styr:suit_endyr` (defaults to
the full hindcast) if year-specific diet data are not included
(e.g. `Year = 0` in `diet_data`). Per-predator `suitMode` allows mixing
parametric and empirical suitability across species.

## 3. Selectivity (per fleet)

> Equations: see vignette 8 [§Selectivity
> parameterizations](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.html#selectivity-parameterizations).

Set the `Selectivity` column of `fleet_control` to one of the values
below (integers or strings both work):

| Value | String | Notes |
|:--:|----|----|
| `0` | `"Fixed"` | Selectivity supplied via the `emp_sel` data sheet |
| `1` | `"Logistic"` | Two-parameter ascending logistic |
| `2` | `"NonParametric"` | Age- or length-bin-indexed; penalties via `Sel_curve_pen1` / `Sel_curve_pen2` |
| `3` | `"DoubleLogistic"` | Dome-shaped (ascending + descending logistic) |
| `4` | `"DescendingLogistic"` | Descending only |
| `5` | `"Hake"` | Non-parametric à la Taylor et al. (Pacific hake) |
| `6` | `"2DAR1"` | 2D AR(1) age × year (Cheng et al. 2024) |
| `7` | `"3DAR1"` | 3D AR(1) age × year × cohort |

Time-varying behavior is set by `Time_varying_sel` (per fleet):

| Value | String | Behavior |
|:--:|----|----|
| `0` | `"Off"` | Time-invariant |
| `1` | `"IID"` | Independent annual deviations |
| `2` | `"AR1"` | Not yet implemented |
| `3` | `"Block"` | Time block |
| `4` | `"RandomWalk"` | Random walk (not allowed with `NonParametric` or `Hake`) |
| `5` | `"RandomWalkAscending"` | Random walk on ascending portion of double logistic (not allowed with `NonParametric` or `Hake`) |

Bin indexing: with non-parametric and AR(1) forms, selectivity is
indexed over `N_sel_bins`.

Selectivity can be normalized relative to a specific age by setting
`Sel_norm_bin != NA` or to the max (non-differentiable) by setting
`Sel_norm_bin < 0`. Selectivity can also be normalized relative to the
average across an age-range by specifying the lower (`Sel_norm_bin`) and
upper (`Sel_norm_bin_upper`) bound.

`Sel_norm_bin` is a bin index on the fleet’s own
`Selectivity_dimension`: an absolute **age** for an age-based fleet (`6`
means age 6, not the sixth bin) or a 1-based **length-bin ordinal** for
a length-based one.

### Sex structure and relative selectivity

Populations are modelled as one or two sexes via `nsex`. **Every
selectivity form is sex-specific by default** when a species has two
sexes: each sex gets its own free copy of every selectivity parameter,
for all forms.

What differs between configurations is the *relative* scale between the
sexes — whether males and females can be selected at genuinely different
levels, or only with different shapes. Normalization makes two
independent decisions, carried by two columns:

- **`Sel_norm_bin` — where the reference is taken.** `NA` normalizes
  nothing; `>= 0` takes it at that bin (or, with `Sel_norm_bin_upper`,
  the mean over that range); `< 0` takes it at the curve’s maximum over
  bins.
- **`Sel_norm_scope` — whose scale it sets.** `"AcrossSexes"` (default)
  pools one reference over the sexes, so the less-selected sex stays
  below 1 and **relative sex-specific selectivity is retained**.
  `"WithinSex"` divides each sex by its own reference, so both reach 1
  and only the *shape* differs by sex.

| `Sel_norm_bin` | `Sel_norm_scope` | Result |
|----|----|----|
| `NA` | — | nothing normalized; the relative scale is free |
| `< 0` | `"AcrossSexes"` | maximum over bins **and sexes** — one sex peaks at 1, the other keeps its relative level |
| `< 0` | `"WithinSex"` | each sex scaled at **its own plateau**, wherever that falls |
| `>= 0` | `"AcrossSexes"` | anchored at that bin, sexes kept comparable |
| `>= 0` | `"WithinSex"` | both sexes forced to 1 at that bin; only the shape differs |

For strongly dimorphic stocks the relative sex selectivity is usually
the quantity of interest, so the `"AcrossSexes"` default is the one to
keep. If instead the sexes plateau at different ages and you want each
scaled at its own plateau, use `Sel_norm_bin < 0` with
`Sel_norm_scope = "WithinSex"` — that is more robust than naming a bin,
since a named bin stops being the plateau once selectivity is
time-varying while the maximum tracks it.

Note that `Sel_norm_scope` has no effect on a one-sex species, or where
`Sel_norm_bin` is `NA`. Before Rceattle 5.8.0 a named bin always implied
`"WithinSex"` and max-normalization always implied `"AcrossSexes"`; a
two-sex fleet normalizing at a named bin should set
`Sel_norm_scope = "WithinSex"` explicitly to keep its old behaviour.

Two things constrain what is actually estimable:

- **The composition data must be joint.** `comp_data$Sex = 3` stacks
  female and male composition into a single vector that is normalized to
  sum to 1 as a whole (as Stock Synthesis does), so the sex ratio is fit
  and the relative selectivity is informed. `Sex = 1` or `2` put each
  sex in its own row, each summing to 1 independently — those carry
  **no** sex-ratio information, and `Sel_norm_bin` then makes no
  difference however it is set. There is no separate sex-ratio
  likelihood component.
- **`Selectivity = "Hake"` always normalizes within sex**, regardless of
  `Sel_norm_bin`, so relative sex selectivity is unavailable with that
  form.

Related settings elsewhere: the sex ratio at recruitment is **fixed
input**, not estimated — it comes from the `sex_ratio` sheet at the
first modelled age. Sex-specific natural mortality is
`build_M1(M1_model = )` (§7) and sex-specific growth is
[`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
(§8). The per-sex hooks in the linkage grammar are
`linkage_spec(by = ~ sex)` and its `sex =` filter.

Note that `fleet_control` also carries a `Sex` column. It is **inert** —
read nowhere in the model — and is retained only so older workbooks
still read. Sex is set per observation on `comp_data`, not per fleet.

## 4. Catchability (per fleet)

> Equations: see vignette 8 [§Catchability
> parameterizations](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.html#catchability-parameterizations).

Set the `Catchability` column of `fleet_control`:

| Value | String | Description |
|:--:|----|----|
| `0` | `"Fixed"` | q supplied directly |
| `1` | `"Estimated"` | q estimated as a free parameter |
| `2` | `"Estimated-with-prior"` | Estimated with a Normal prior |
| `3` | `"Analytical"` | Closed-form analytical q (concentrated likelihood) |
| `4` | `"PowerEquation"` | Power-function catchability (q · B^(α)) |
| `5` | `"Environmental"` | q linked to an environmental index (`Catchability_index`) |
| `6` | `"AR1"` | Annual AR(1) deviations on log-q fit to a single index (Rogers et al. 2024) |

Time-varying behavior is set by `Time_varying_q` (per fleet):

| Value | String         | Behavior                      |
|:-----:|----------------|-------------------------------|
|  `0`  | `"Off"`        | Time-invariant                |
|  `1`  | `"IID"`        | Independent annual deviations |
|  `2`  | `"AR1"`        | Not yet implemented           |
|  `3`  | `"Block"`      | Time block                    |
|  `4`  | `"RandomWalk"` | Random walk                   |

If `Catchability` is `AR1` or `Environmental`, `Time_varying_q` can be
an integer or character of integers specifying the environmental series
in `env_data` to use (e.g. `Time_varying_q = 1` or
`Time_varying_q = "1,2,3"`).

## 5. Composition likelihoods

Set `Comp_distribution` (or `CAAL_distribution`) per fleet:

| Value | String                                       |
|:-----:|----------------------------------------------|
| `-1`  | `"MultinomialAFSC"` (legacy AFSC convention) |
|  `0`  | `"Multinomial"`                              |
|  `1`  | `"DirichletMultinomial"`                     |

Diet composition uses an analogous switch on the bioenergetics control
sheet: `Diet_distribution = 0` (multinomial) or `1`
(Dirichlet-multinomial).

Conditional age-at-length (CAAL) data go through their own data path and
are weighted via `CAAL_weights` per fleet.

## 6. Recruitment / stock-recruit (`recFun = build_srr()`)

| `srr_fun` | String           | Description                      |
|:---------:|------------------|----------------------------------|
|    `0`    | `"mean"`         | Mean recruitment with deviations |
|    `2`    | `"BevertonHolt"` | Beverton-Holt                    |
|    `4`    | `"Ricker"`       | Ricker                           |

Priors and environmental relationships can be added to stock-functions
via
[`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md).
See
[`?build_srr`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
for options and
[`vignette("environmental-linkages-and-priors")`](https://grantdadams.github.io/Rceattle/articles/environmental-linkages-and-priors.md)
for details.

## 7. Natural mortality (`M1Fun = build_M1()`)

### `M1_model` — fixed-effects structure of M1

| Value | String                | Description                              |
|:-----:|-----------------------|------------------------------------------|
|  `0`  | `"fixed"`             | use the input `M1_base` (no estimation). |
|  `1`  | `"sex_age_invariant"` | estimate one `M1_{spp}`.                 |
|  `2`  | `"sex_specific"`      | estimate `M1_{spp, sex}`.                |
|  `3`  | `"sex_age_specific"`  | estimate `M1_{spp, sex, age}`.           |

Priors and environmental relationships can be added to M1 via
[`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md).
See
[`?build_M1`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
for options and
[`vignette("environmental-linkages-and-priors")`](https://grantdadams.github.io/Rceattle/articles/environmental-linkages-and-priors.md)
for details.

### `M1_re` — random-effects structure on M1

|     Value     | String           | Description                       |
|:-------------:|------------------|-----------------------------------|
| `0` (default) | `"none"`         | No random effects                 |
|      `1`      | `"iid_age"`      | IID by age, constant over years   |
|      `2`      | `"iid_year"`     | IID by year, constant over ages   |
|      `3`      | `"iid_age_year"` | IID across both year and age      |
|      `4`      | `"ar1_age"`      | AR(1) by age, constant over years |
|      `5`      | `"ar1_year"`     | AR(1) by year, constant over ages |
|      `6`      | `"ar1_age_year"` | 2D AR(1) over both year and age   |

Variance and correlation parameters are species-specific but
sex-invariant. `M2_use_prior = TRUE` adds a log-normal prior with mean
`M_prior` and SD `M_prior_sd` on total M (M1 + M2) in multi-species
mode.

## 8. Growth and weight-at-age (`growthFun = build_growth()`)

### `growth_model` — fixed-effects structure of growth

| Value | String | Description |
|:--:|----|----|
| `0` (default) | `"empirical"` | Empirical weight-at-age supplied via the `weight` data sheet. Forecasts re-use the terminal-year weight schedule. |
| `1` | `"vonBertalanffy"` | von Bertalanffy length-at-age + length-weight allometry (`alpha_wt_len`, `beta_wt_len` on the data control sheet) |
| `2` | `"Richards"` | Sex-specific Richards growth |

Priors and environmental relationships can be added to growth via
[`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md).
See
[`?build_growth`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
for options and
[`vignette("environmental-linkages-and-priors")`](https://grantdadams.github.io/Rceattle/articles/environmental-linkages-and-priors.md)
for details.

### `growth_re` — random effects on growth

|     Value     | Description             |
|:-------------:|-------------------------|
| `0` (default) | No random effects       |
|     `≥ 1`     | **Not yet implemented** |

The length-based suitability modes (`suitMode = 1/3/5`) are not yet
enabled; use a weight-based mode (`2/4/6`) for parametric suitability.
When the length-based modes are enabled, prey length-at-age is taken
from the estimated growth model (`growth_model > 0`).

## 9. Harvest control rules (`HCR = build_hcr()`)

| `HCR` | String | Description | Multi-species? |
|:--:|----|----|:--:|
| `0` | `"NoFishing"` | No fishing — estimate hindcast only | Yes |
| `1` | `"CMSY"` | Maximize joint catch across species (option to keep depletion ≥ `Plimit`) | Yes |
| `2` | `"ConstantF"` | Constant F set at `Ftarget` per species | Yes |
| `3` | `"ConstantFSSB"` | F that achieves `Ftarget`% of SB₀ at the end of the projection | Yes |
| `4` | `"ConstantFSPR"` | F_(SPR) at `Ftarget` (NEFSC convention; scale via `Fmult`) | No |
| `5` | `"NPFMC"` | NPFMC Tier 3 SPR-based HCR (`Ftarget`, `Flimit`, `Ptarget`, `Plimit`, `Alpha`) | No |
| `6` | `"PFMC"` | PFMC Category 1 `Ptarget`-`Plimit` ACL with uncertainty buffer around `Flimit` (`Pstar`, `Sigma`) | Yes |
| `7` | `"SESSF"` | SESSF Tier 1 SPR-based HCR | No |

`DynamicHCR = TRUE` switches from static to dynamic SB₀ reference
points. `HCRorder` controls the order in which species F is solved for
iterative multi-species HCRs (e.g. predators before prey). See
[`?build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
for more details.

## 10. Projection and MSE

| Function | Purpose |
|----|----|
| `fit_mod(estimateMode = 2, HCR = build_hcr(...))` | Forward projection under a fixed HCR |
| [`sample_rec()`](https://grantdadams.github.io/Rceattle/reference/sample_rec.md) | Bootstrap historical recruitment deviations (with optional `rec_trend`) into the projection |
| [`remove_F()`](https://grantdadams.github.io/Rceattle/reference/remove_F.md) | Refit a model with F = 0 (used internally for dynamic reference points) |
| [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md) | Closed-loop MSE: refit EM at intervals of `assessment_period`, sample new data at `sampling_period` |
| [`load_mse()`](https://grantdadams.github.io/Rceattle/reference/load_mse.md) | Reload per-simulation `.rds` files written when `dir =` is supplied |
| [`check_mse()`](https://grantdadams.github.io/Rceattle/reference/check_mse.md) | Validate which OM/EM simulations converged |
| [`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md) | Per-fleet performance metrics: mean catch, IAV, P(closed), P(F \> Flimit), P(SSB \< SSBlimit), terminal depletion |

[`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
parallelizes over simulations on a PSOCK cluster. See the [MSE
vignette](https://grantdadams.github.io/Rceattle/articles/hcrs-and-mses.md)
for the seed plumbing and reproducibility guarantees.

## 11. Diagnostics and inference

| Function | Purpose |
|----|----|
| [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md) | Peel terminal years and refit; reports Mohn’s ρ for each quantity / species, optionally with forecast-skill peels (`nyrs_forecast`) |
| [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md) | Refit from random-perturbed starting values; check optimum stability |
| [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md) | Parametric simulation from a fitted model |
| [`compare_sim()`](https://grantdadams.github.io/Rceattle/reference/compare_sim.md) | Estimation-vs-truth comparison across simulated replicates |
| [`model_average()`](https://grantdadams.github.io/Rceattle/reference/model_average.md) | Weighted average of derived quantities across multiple fits, with optional bootstrap uncertainty |
| `calc_mcall_ianelli()` / `calc_mcall_ianelli_diet()` | McAllister-Ianelli composition-data reweighting |

## 12. Plotting

Pop / mortality / fit / predation / biology plots — all named after the
quantity they show: `plot_biomass`, `plot_ssb`, `plot_recruitment`,
`plot_depletion`, `plot_depletionSSB`, `plot_exploitable_biomass`,
`plot_catch`, `plot_f`, `plot_mortality`, `plot_m_at_age`,
`plot_m2_at_age_prop`, `plot_index`, `plot_indexresidual`, `plot_comp`,
`plot_data`, `plot_form`, `plot_b_eaten`, `plot_b_eaten_prop`,
`plot_diet_comp`, `plot_ration`, `plot_selectivity`,
`plot_selectivity_vs_maturity`, `plot_maturity`, `plot_stock_recruit`,
`plot_timeseries`. Most accept a single fitted model or a list of models
with `model_names`.

These are built with **ggplot2** (colorblind-safe Okabe-Ito palette for
discrete series, viridis for continuous ones, `theme_classic`) and
**return the `ggplot` object**, so they print when called interactively
and can be further customized by adding ggplot2 layers
(e.g. `plot_biomass(fit) + ggplot2::ggtitle("My title")`). Pass
`file = "stem"` to also write the figure to `stem_<plot>.png`.

**Display units.** The model works in mt and thousands of fish:
numbers-at-age are in thousands and weight-at-age in kg, so biomass and
SSB come out in mt and recruitment in thousands of fish. The timeseries
plotters rescale to `million mt` and `millions` for the axis, and derive
the recruitment age from `minage` rather than assuming age 1. Supply the
model’s inputs on that convention – catch and index in mt, weight-at-age
in kg – or the axis labels will not describe what is plotted. Pass
`ylab` to override the derived label.

**Confidence intervals.** `add_ci = TRUE` draws a 95% interval. Strictly
positive series take it on the log scale, `exp(log(x) ± 1.92 · sd_log)`,
using the model’s `log_biomass` / `log_ssb` / `log_R` where they exist
and the delta-method identity `sd(log x) = sd(x)/x` where they do not.
The interval is therefore right-skewed and cannot reach zero for a weak
year class or a depleted stock. See
[`?plot_timeseries`](https://grantdadams.github.io/Rceattle/reference/plot_timeseries.md)
for the full derivation.

## References

- Adams, G.D., Holsman, K.K., Barbeaux, S.J., Dorn, M.W., Ianelli, J.N.,
  Spies, I., Stewart, I.J., Punt, A.E. (2022). An ensemble approach to
  understand predation mortality for groundfish in the Gulf of Alaska.
  *Fisheries Research*, 251, 106303.
- Holsman, K.K., Ianelli, J., Aydin, K., Punt, A.E., Moffitt, E.A.
  (2016). A comparison of fisheries biological reference points
  estimated from temperature-specific multi-species and single-species
  climate-enhanced stock assessment models. *Deep Sea Research Part II*,
  134, 360-378.
- Wassermann, S.N., Adams, G.D., Haltuch, M.A., Kaplan, I.C., Marshall,
  K.N., Punt, A.E. (2025). Even low levels of cannibalism can bias
  population estimates for Pacific hake. *ICES Journal of Marine
  Science*, 82(1), fsae064.
- Cheng, M.L.H., et al. (2024). 2D / 3D AR(1) selectivity (citation in
  vignette 8).
- Rogers, L.A., et al. (2024). AR(1) catchability for GOA pollock
  (citation in vignette 8).
