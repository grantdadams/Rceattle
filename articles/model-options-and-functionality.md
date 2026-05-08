# 0. Model options and functionality

This vignette is an overview of the parameterization choices and
capabilities available in `Rceattle`. It complements:

- vignette 1, [Rceattle: An
  Introduction](https://grantdadams.github.io/Rceattle/articles/introduction.md),
  which walks through fitting a model end-to-end;
- vignette 8, [Model
  parameterizations](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.md),
  which is the in progress equation-level reference; and
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

### `estimateMode` — what does `fit_mod()` actually do?

| Value | Behavior |
|:--:|----|
| `0` | Estimate hindcast and projection parameters via `nlminb` / [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) (default). |
| `1` | Estimate only the hindcast parameters via `nlminb` / [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).. |
| `2` | Hold the hindcast fixed and run only the projection under the supplied `HCR`. |
| `3` | Evaluate the report once at `inits` (debug / mapping check). |
| `4` | Run the TMB objects through MakeADFun with out estimation (debug / mapping check). |

### Random-effects toggles

`random_rec`, `random_q`, and `random_sel` move log-recruitment, the
log-catchability deviations, and the time-varying selectivity deviations
(respectively) from penalised fixed effects to random effects
marginalised via the Laplace approximation and estimate their associated
distribution parameters.

## 2. Predation: `suitMode` (per predator)

> Equations: see vignette 8 [§Predation
> sub-model](https://grantdadams.github.io/Rceattle/articles/model-parameterizations.html#predation-sub-model)
> (suitability, predation mortality, consumption / bioenergetics).

`suitMode` chooses how prey preference is modelled. Length 1 (same for
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
`Sel_norm_bin1 != NA` or to the max (non-differentiable) by setting
`Sel_norm_bin1 < 0`. Selectivity can also be normalized relative to the
average across an age-range by specifying the lower (`Sel_norm_bin1`)
and upper (`Sel_norm_bin2`) bound.

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
| `5` | `"Environmental"` | q linked to an environmental index (`Q_index`) |
| `6` | `"AR1"` | Annual AR(1) deviations on log-q fit to a single index (Rogers et al. 2024) |

Time-varying behavior is set by `Time_varying_q` (per fleet):

| Value | String         | Behavior                      |
|:-----:|----------------|-------------------------------|
|  `0`  | `"Off"`        | Time-invariant                |
|  `1`  | `"IID"`        | Independent annual deviations |
|  `2`  | `"AR1"`        | Not yet implemented           |
|  `3`  | `"Block"`      | Time block                    |
|  `4`  | `"RandomWalk"` | Random walk                   |

If `Selectivity` is `AR1` or `Environmental`, `Time_varying_q` can be a
integer or character of integers specifying the environmental series in
`env_data` to use (e.g. `Time_varying_q = 1` or
`Time_varying_q = "1,2,3"`).

## 5. Composition likelihoods

Set `Comp_loglike` (or `CAAL_loglike`) per fleet:

| Value | String                                       |
|:-----:|----------------------------------------------|
| `-1`  | `"MultinomialAFSC"` (legacy AFSC convention) |
|  `0`  | `"Multinomial"`                              |
|  `1`  | `"DirichletMultinomial"`                     |

Diet composition uses an analogous switch on the bioenergetics control
sheet: `Diet_loglike = 0` (multinomial) or `1` (Dirichlet-multinomial).

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
for options and `vignette("environmental-linkages")` for details.

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
for options and `vignette("environmental-linkages")` for details.

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
sex-invariant. `M2_use_prior = TRUE` extends adds a log-normal prior
with mean `M_prior` and SD `M_prior_sd` on total M (M1 + M2) in
multi-species mode.

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
for options and `vignette("environmental-linkages")` for details.

### `growth_re` — random effects on growth

|     Value     | Description             |
|:-------------:|-------------------------|
| `0` (default) | No random effects       |
|     `≥ 1`     | **Not yet implemented** |

Length-based suitability (`suitMode = 1/2/3/4/5/6`) is wired through to
the estimated growth model when `growth_model > 0`.

## 9. Harvest control rules (`HCR = build_hcr()`)

| `HCR` | String | Description | Multi-species? |
|:--:|----|----|:--:|
| `0` | `"NoFishing"` | No fishing — estimate hindcast only | Yes |
| `1` | `"CMSY"` | Maximise joint catch across species (option to keep depletion ≥ `Plimit`) | Yes |
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
parallelises over simulations on a PSOCK cluster. See the [MSE
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
| [`calc_mcall_ianelli()`](https://grantdadams.github.io/Rceattle/reference/calc_mcall_ianelli.md) / [`calc_mcall_ianelli_diet()`](https://grantdadams.github.io/Rceattle/reference/calc_mcall_ianelli_diet.md) | McAllister-Ianelli composition-data reweighting |

## 12. Plotting

Pop / mortality / fit / predation / biology plots — all named after the
quantity they show: `plot_biomass`, `plot_ssb`, `plot_recruitment`,
`plot_depletion`, `plot_depletionSSB`, `plot_exploitable_biomass`,
`plot_catch`, `plot_f`, `plot_mortality`, `plot_m_at_age`,
`plot_m2_at_age_prop`, `plot_index`, `plot_logindex`,
`plot_indexresidual`, `plot_comp`, `plot_data`, `plot_form`,
`plot_b_eaten`, `plot_b_eaten_prop`, `plot_diet_comp`, `plot_ration`,
`plot_selectivity`, `plot_selectivity_vs_maturity`, `plot_maturity`,
`plot_stock_recruit`, `plot_timeseries`. Most accept a single fitted
model or a list of models with `model_names`.

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
