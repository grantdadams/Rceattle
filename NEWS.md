# Rceattle 4.13.0

## New features

* **Save / load a run configuration (`save_config()` / `load_config()`).** A full
  run configuration -- the [model_config()] structure plus the estimation controls
  and `fit_control()` bundle -- round-trips to a documented, git-diffable YAML file
  (each field carries its doc string as a comment; fields at their default are
  omitted so two runs diff to only their real differences). Apply a saved config
  with `fit_mod(data_list, config = load_config("run.yaml"))`; it fills only the
  arguments the caller did not pass (an explicit argument always wins), and every
  fit now records the run configuration it used (`fit$run_config`, also reachable
  via `run_config(fit)`). Follows SAM's `saveConf`/`loadConf` and the
  config-separate-from-weights idiom; adds a `yaml` dependency.

* **Composition tail accumulation (AFSC `ac_yng` / `ac_old`).** Two new
  `fleet_control` columns, `Comp_accum_young` and `Comp_accum_old`, fold the
  young and old tails of a fleet's age/length composition into a boundary bin
  before the composition likelihood (per sex block for joint-sex comps), for
  every composition family including the default `MultinomialAFSC`. They are
  1-based bin ordinals on the fleet's composition dimension; leaving them unset
  (or `young = 1` / `old` at the last bin) applies no accumulation, so every
  existing model is bit-identical. `data_check()` rejects out-of-range or
  inverted (`young > old`) bins, and OSA residuals are not available for a fleet
  with active accumulation (the residuals are built on the un-accumulated bins).

## Bug fixes

* **`self_test()` no longer errors with `object 'getsd' not found`.** Its
  per-simulation refit closure referenced a `getsd` value that the function --
  unlike `retrospective()`, `jitter()`, and `profile()` -- never defined, so
  every simulation died on both the sequential and the parallel-cluster path.
  `self_test()` now takes a `getsd` argument (default `NULL`, inheriting the
  input model's setting) like its sibling functions.

# Rceattle 4.12.0

## New features

* **Intuitive `fleet_control` / data column names.** The workbook and data-list
  columns were renamed to say what they hold, each keeping a back-compat alias
  that is upgraded on read, so existing data lists, `.rda`, workbooks, and
  scripts keep working (with a one-time deprecation message) and fit
  identically: `Q_index`/`Q_init`/`Q_sd_prior` →
  `Catchability_index`/`Catchability_init`/`Catchability_prior_sd`;
  `Comp_loglike`/`CAAL_loglike`/`Index_loglike` →
  `Comp_distribution`/`CAAL_distribution`/`Index_distribution`;
  `Sel_norm_bin1`/`Sel_norm_bin2` → `Sel_norm_bin`/`Sel_norm_bin_upper`;
  `weight1_Numbers2` → `Observation_units`; `proj_F_prop` → `Proj_F_proportion`;
  the control scalar `sigma_rec_prior` → `sigma_rec`; the bioenergetics
  `Diet_loglike` → `Diet_distribution`; and the post-fit output column
  `Est_weights_mcallister` → `Comp_weights_mcallister` (the old name is still
  written for downstream readers).

* **`write_template()`.** A new export that writes a minimal, structurally
  complete single-species starter workbook on the canonical column names; it
  round-trips through `read_data()` and builds under `fit_mod(estimateMode = 3)`.

* **Single-source workbook schema.** The column dictionary now lives once in the
  package (`R/0-column_schema.R`) and drives `switch_check()` defaults,
  `write_data()` object order, the embedded `meta_data` documentation sheet, and
  the roxygen field dictionary, which are kept in sync by guard tests. ~16
  previously used-but-undocumented columns are now documented.

## Bug fixes

* **Robust workbook reads.** `read_data()` now errors clearly when a required
  sheet (`control`, `fleet_control`) is missing and skips optional sheets when
  absent, so a minimal single-species workbook reads cleanly instead of failing
  with a cryptic error. A non-numeric cell in the control or bioenergetics rows
  now errors by name instead of silently becoming `NA`. `rearrange_data()` fails
  clearly on a malformed `fleet_control` rather than via a cryptic `dplyr` error.

* **Removed the dead accumulation-age feature.** `Accumulation_age_lower/upper`
  were only range-validated and never applied to composition data; they are
  dropped (with a soft-deprecation message if an old workbook carries them).

# Rceattle 4.11.0

## New features

* **`initMode = "FishedEquilibrium"` (5).** A new initial-age-structure mode: an
  F = 0 equilibrium seeded by the *first-year* recruitment
  (`exp(rec_pars + rec_dev[year 1])`) decayed by natural mortality, with initial
  deviates turned off and no init-dev penalty. Unlike `initMode = "Equilibrium"`
  (which seeds off the mean-recruitment equilibrium `R0`), the initial numbers
  track the year-1 recruitment deviation, matching the AFSC GOA pollock (Cole
  Monnahan) convention — the first-year cohort and the initial age composition
  share one deviation. Every other mode is unchanged (golden references
  bit-identical).

* **Priors on base selectivity parameters through the linkage grammar.** A
  selectivity linkage with an intercept-only formula and a `priors` entry now
  places a prior on the base logistic parameter that carries the level — the
  ascending/descending slope (`slp_asc` / `slp_desc`, on `log_sel_slp`, log
  scale) or inflection (`inf_asc` / `inf_desc`, on `sel_inf`, natural scale),
  e.g. `build_selectivity(linkages = list(inf_asc = linkage_spec(~ 1, priors =
  list(\`(Intercept)\` = normal(0, 3)))))`. Use `lognormal()` for the log-scale
  slopes and `normal()` for the natural-scale inflections. This mirrors the
  prior-only `build_composition()` path and reproduces the AFSC GOA pollock
  selectivity priors. (Previously such an intercept prior was silently inert for
  selectivity — the `sel` process was missing from the prior re-targeting.)

* **Per-predator suitability reference years** (`suit_styr` / `suit_endyr`).
  These `fit_mod()` arguments (and the underlying `data_list` fields) now accept
  a vector of length `nspp` so each predator can average its suitability over a
  different set of years — e.g. a California Current model with hake 1980–2019,
  arrowtooth 2013–2018, and sablefish 2005–2008. A single scalar is recycled to
  every predator, exactly reproducing the previous global-window behaviour
  (`BS2017MS` and the golden references are unchanged to numerical tolerance).
  Internally `suit_styr` / `suit_endyr` became per-predator `DATA_IVECTOR`s and
  the suitability-averaging and stomach-content prediction loops index them by
  predator.

* **`plot_diet_comp2()` and `plot_diet_comp1()`.** `plot_diet_comp2()` adds
  aggregation-aware diet-composition diagnostics (line plots when prey- or
  predator-age is aggregated, dodged bars when both are, bubble grids when fully
  disaggregated), built on `residuals(source = "diet")`. `plot_diet_comp1()` is
  an alias of `plot_diet_comp()` (the bubble/grid diagnostic).

## Bug fixes

* **`read_data()` tolerates trailing empty age columns in `NByageFixed`.** Older
  writers pad the fixed numbers-at-age sheet to a wider age range than
  `max(nages)`; the all-`NA` trailing columns are now dropped on read instead of
  tripping the `data_check()` column-count validation.

* **Parallel workers now run the in-session package.** `run_mse()`,
  `retrospective()`, and `jitter()` use a FORK cluster on non-Windows platforms,
  so workers inherit the parent's loaded namespace (and the large fitted OM/EM
  objects) via copy-on-write instead of each reloading `library(Rceattle)` and
  receiving a `clusterExport()` of those objects. This fixes silently running a
  stale *installed* package on the workers during `pkgload::load_all()` /
  development sessions, and removes the per-worker startup cost (~3 s for 6
  workers with the Bering Sea multispecies objects). Results are unchanged
  (parallel and serial runs are bit-identical for a given seed). Windows keeps
  the previous PSOCK path as a fallback.

# Rceattle 4.10.0

## New features

* **Random-effect linkages (`~ (1 | group)`).** A linkage formula may now carry
  an IID random effect, so a parameter can vary year to year (or over any
  grouping) as a set of deviations damped by an estimated
  `N(0, sigma)` density, rather than fixed covariate slopes:

    ```r
    fit_mod(...,
      qFun = build_catchability(linkages = list(
               q = linkage_spec(~ (1 | Year), by = ~ fleet))))
    ```

  Each level of the grouping variable gets one deviation (`beta_linkage_re`);
  each distinct group estimates its own log-SD (`log_sigma_linkage`). The
  density is reported in the new `jnll_comp` row *"Linkage random effects"*.
  The deviation SD is routed through the same `linkage_spec()` arguments as
  every other parameter: `init = list(sigma = v)` **fixes** it at an input
  value (reproducing the legacy `Time_varying_*_sd_prior` fixed input), and
  `priors = list(sigma = lognormal(...))` places a **prior** on it and
  estimates it — the first prior on a deviation SD anywhere in the model.

* **`data_requirements()` — see which inputs a model configuration needs.** A new
  exported reader reports, for a given model spec, which top-level data inputs are
  **Required**, **Optional** (used if supplied, otherwise default-filled by
  `clean_data()`), or **Ignored** (not consulted because the feature is switched
  off) — the same conditions `data_check()` enforces at fit time, surfaced up
  front instead of buried in the validator:

    ```r
    data_requirements(msmMode = 1)          # multispecies: diet/ration/bioenergetics Required
    data_requirements(BS2017SS, msmMode = 0) # preview a data object as single-species
    ```

  It accepts either an existing (possibly partial) data list or the convenience
  switch arguments; an explicit switch argument overrides the data list's stored
  value (matching `fit_mod()` precedence). Internally, `data_check()`'s conditional
  presence-requirement gates were refactored to consume one declarative
  requirement table, so the reader and the validator can never drift apart.

* **`build_data()` — assemble a data list in R.** A code-first constructor
  complementing `read_data()`: supply only the blocks a model uses and the
  optional blocks a single-species model does not need are default-filled by
  `clean_data()`. Three combinable entry points cover the workflows real
  assessments use:

    ```r
    build_data(base = BS2017SS, projyr = 2060)          # copy-and-edit a dataset
    build_data(file = "model.xlsx", fleet_control = fc) # read a workbook, override
    build_data(nspp = 1, styr = 1977, endyr = 2023,     # assemble from blocks
               fleet_control = fc, catch_data = catch, index_data = survey)
    ```

  Overrides are checked against the recognised schema, so a typo (`maturty`) is
  caught at construction with a suggestion rather than surfacing later in a fit;
  legacy names (`fsh_biom`, `srv_biom`, `wt`, `pmature`, `Pyrs`) are mapped to
  their canonical equivalents. Validation is deferred to `data_check()` at fit
  time (one source of truth); `build_data()` runs only a light presence
  pre-check so a missing required block is reported early. The result is the
  same bare list `read_data()` returns and round-trips through `write_data()`
  unchanged — a `build_data(base = X)` object fits bit-identically to `X`.

* **`model_config()` — a model configuration that travels with the data.** The
  model-structure arguments of `fit_mod()` (`msmMode`, `initMode`, the HCR and
  the `build_*()` process specifications) can now be bundled into a slot on the
  data list, so a data object records how it is meant to be fit:

    ```r
    dat <- build_data(base = BS2017MS, model_config = model_config(msmMode = 1))
    fit_mod(dat)                     # fits as multispecies without passing msmMode
    ```

  `fit_mod()`'s signature and defaults are unchanged; when a data list carries a
  `model_config`, `fit_mod()` reads each field only for arguments the caller did
  **not** pass (detected with `missing()`), and an explicitly-passed argument
  always overrides the slot. With no slot present the behaviour is byte-identical
  (a `BS2017SS` fit is bit-identical). A call that passes an argument — even at
  its default — overrides the slot, so omit the argument to let the configuration
  take effect. The slot is code-side structure, not a workbook sheet, so it does
  not persist through a `write_data()`/`read_data()` round-trip (a warning fires).

* **Spec-tree `print()` / `summary()` for data objects and fits.** A
  `build_data()` object now carries the class `"Rceattle_data"` and prints as an
  indented specification tree — dimensions → fleets (with their selectivity /
  catchability forms and mirroring) → configured processes → active linkages →
  any attached `model_config()` — instead of dumping the ~40-element list. The
  same tree is shown by `print()` on a fitted model above its fit statistics, so
  "read 600 lines of switch tables" becomes "read the printout". The class is a
  thin tag: every consumer still treats the object as a plain list, so it does
  not change a fit.

## Bug fixes

* **`M1_re = 6` (separable age × year 2D-AR1 on M) now estimates its correlations.**
  A gate bug in `build_map()` (`if(M1_re_model == 5)` inside a block reachable only
  when the mode is 3 or 6) plus a reference to an undefined index left mode 6's age
  and year AR1 correlations unmapped, so the separable AR1 silently collapsed to IID
  (identical to `M1_re = 3`). Mode 6 now frees both `rho` hyperparameters, so it is a
  genuine separable 2D-AR1. This changes results only for models using `M1_re = 6`
  (none of the bundled examples do); golden references are unaffected.

* **The deprecated `est_M1` name is now recognized everywhere `M1_model` is.**
  The natural-mortality estimation switch was renamed `est_M1` → `M1_model`, but
  the data dictionary still documented `est_M1` and no alias existed, so a data
  list carrying `est_M1` was silently ignored — it fell through to the `M1_model`
  default with no warning. `est_M1` is now folded into `M1_model` in `fit_mod()`
  (before the `build_M1()` reconciliation), `switch_check()`, and `combine_data()`,
  with a deprecation message. As with `M1_model`, the value set on the data list
  is a default that an explicit `build_M1(M1_model = ...)` argument overrides — but
  a data-list value that differs from the `build_M1()` setting now *warns* rather
  than being dropped silently. The recommended way to request M1 estimation remains
  `fit_mod(..., M1Fun = build_M1(M1_model = ...))`. The dictionary now documents
  `M1_model`.

* **A catchability linkage on a non-estimated q now errors.** A `q` linkage
  (environmental or random-effect) on a fleet whose `Catchability` is `"Fixed"`
  or an analytical form (`"Analytical"` / `"AnalyticalArith"`) is rejected up
  front, naming the fleet. Previously only the analytical forms were caught, so
  a linkage on a `"Fixed"` fleet silently turned a fixed catchability
  time-varying — contrary to the assessor's setting. Set `Catchability` to
  `"Estimated"` / `"Estimated-with-prior"` to link q.
  This is the same statistical model as the existing `Time_varying_q`/`_sel`
  deviate processes, expressed through the linkage grammar.

* **Random-walk linkages (`rw(1 | group)`).** A deviation that follows a random
  walk — `N(0, sigma)` on successive first differences, the first deviate
  pinned so the walk's level stays with the base parameter — reproducing the
  Dorn-style `Time_varying_* = "RandomWalk"` process through the grammar. The
  grouping variable must be numeric (a real elapsed-time lag).

* **AR1 linkages (`ar1(1 | group)`).** A stationary first-order autoregressive
  deviation, `SCALE(AR1(rho), sigma)` with `sigma` the marginal SD and `rho`
  the correlation (the glmmTMB convention, and the same form as the Rogers et
  al. (2024) QAR1 catchability process). The correlation is routed through the
  same grammar as `sigma`: `init = list(rho = 0.7)` fixes it,
  `priors = list(rho = normal(0, 0.3))` places a prior on it, else it is
  estimated free. Reduces to the IID density at `rho = 0`.

* **State-space environmental covariate (Rogers et al. 2024 QAR1).** An
  `ar1(1 | Year)` term can be turned into a measured latent covariate via
  `observe` / `obs_sd`: the AR1 latent is observed as an `env_data` column with
  a fixed measurement SD, and enters the linked parameter through an estimated
  effect size. This reproduces the Rogers-2024 QAR1 catchability model
  (`Estimate_q = 6`) through the grammar:

    ```r
    build_catchability(linkages = list(
      q = linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 1,
                       observe = "QcovPol", obs_sd = 0.1)))
    ```

  The AR1 process, the covariate observation, and the effect size are all
  carried by the one term; `observe` requires an `ar1` term and a positive
  fixed `obs_sd`.

* **Priors on Dirichlet-multinomial data weights (`build_composition()`).** The
  DM likelihood self-tunes each composition dataset's effective sample size, but
  that weight can be poorly identified when a fleet has few comp years. A new
  `build_composition()` linkage attaches a prior to the DM weight, keeping the
  implied effective sample size in a believable range without reverting to
  Francis / McAllister–Ianelli hand-tuning:

    ```r
    fit_mod(...,
      compFun = build_composition(linkages = list(
                  theta_comp = linkage_spec(~ 1, by = ~ fleet, fleet = c(1, 2),
                                 priors = list(`(Intercept)` = gamma(2, 0.5))))))
    ```

  Keys `theta_comp` / `theta_caal` / `theta_diet` cover age/length composition,
  conditional age-at-length, and predator diet. The process is prior-only (the
  weight is a scalar, not year-varying): a covariate, an `init`, or an
  `est_phase` on a composition spec, or a linkage on a fleet/predator not fit
  with the `DirichletMultinomial` likelihood, errors up front.

* **Formula-linkage effect sizes are now reported.** The estimated linkage
  coefficients (`beta_linkage`), the random-effect deviations (`beta_linkage_re`),
  and the Rogers-2024 QAR1 effect size (`beta_linkage_obs`) appear in
  `fit$quantities`, and the coefficients / effect size are `ADREPORT`'d so they
  carry a standard error in the `sdreport`. Previously the QAR1 effect size — the
  quantity the model exists to estimate — was buried in the raw parameter vector
  with no readable exposure or uncertainty.

* **QAR1 covariates need not span the whole series.** A state-space `observe`
  covariate (e.g. a survey index that began mid-series) is now handled without
  hand-padding: `env_data` is auto-extended to start at `styr` with `NA` for
  missing years, the latent AR1 spans all hindcast years, and the observation is
  applied **only** where the covariate is present (un-observed years are masked
  out — no fabricated observations, unlike the legacy mean-fill). Fixed-effect
  covariates must still be finite over the model range: a missing year in a
  fixed-effect (non-`observe`) covariate now errors clearly rather than silently
  producing an NaN.

* **QAR1 observation SD can be fixed or estimated.** The Rogers-2024 `observe` /
  `obs_sd` state-space covariate holds the measurement SD **fixed** by default;
  pass `linkage_spec(..., obs_sd_est = TRUE)` to **estimate** it (one per observed
  group, started from `obs_sd`), as the reference `Estimate_q = 6` / GOApollock
  model does. Estimation is opt-in because a freely-estimated `obs_sd` collapses
  toward 0 on a smooth covariate (the effect-size / `obs_sd` identifiability
  degeneracy) — keep it fixed unless the covariate is informative (a prior on
  `obs_sd` to regularise this is future work).

* **Selectivity penalties can be given as standard deviations.** The cryptic
  penalty *weights* in `Sel_curve_pen1/2/3` can instead be supplied as SDs via
  `Sel_shape_sd` (+ `Sel_shape_dir` `"Decreasing"`/`"Increasing"` for the
  directional sign), `Sel_curvature_sd`, and `Sel_devmag_sd`. Each penalty is a
  Gaussian SSQ, so `switch_check()` converts `weight = 1/(2*sd^2)`. Legacy
  `Sel_curve_pen` values are never overwritten (existing models are bit-identical);
  a fleet supplying the SD columns fits equivalently. `Sel_shape_sd` /
  `Sel_devmag_sd` apply to `NonParametric` (2/9) and `LogisticPM` (11);
  `Sel_curvature_sd` is `NonParametric`-only (LogisticPM does not use
  `Sel_curve_pen2`). Setting an SD column on a form that doesn't use that slot as a
  weight (e.g. 2D/3D-AR1, which reuse `Sel_curve_pen` for logit-scale correlations),
  or a non-positive SD, errors clearly.

* **Readable string aliases for integer-coded switches.** Following the CIE
  review's "rename options to improve interpretability" (e.g. *constant → not
  estimated*), integer-only switches now also accept self-explanatory strings,
  resolved by `switch_check()` (or `build_srr()`) to the same integer codes (so
  fits are identical). New aliases follow a consistent convention — the
  not-estimated value is `"Fixed"`, estimated is `"Estimated"`:
    - `estDynamics`: `"Estimated"` (0) / `"Fixed"` (1) / `"FixedScaled"` (2) /
      `"FixedScaledByAge"` (3);
    - `Estimate_index_sd` / `Estimate_catch_sd`: `"Fixed"` (0) / `"Estimated"` (1) /
      `"Analytical"` (2);
    - `srr_est_mode` (via `build_srr()`): `"Fixed"` (0) / `"Estimated"` (1) /
      `"LognormalPrior"` (2) / `"BetaPrior"` (3);
    - `suitMode`: the documented string map (`"Empirical"`, `"GammaWeight"`, …) is
      now actually applied (previously defined but never wired).

## Deprecations

* **`Selectivity_block` is now optional; `Q_block` is deprecated and ignored.** The
  per-observation time-block columns in `index_data`/`catch_data` are read only for
  `Time_varying_{sel,q} = "Block"`; every other configuration (including the
  random-effect linkages) ignores them. `clean_data()` now default-fills a missing
  `Selectivity_block` with `1` (a single block), so you need only supply it for
  Block-mode fleets. `Q_block` was never read (q time-blocking reuses
  `Selectivity_block`) — supplying it now warns, and it has been dropped from all
  bundled example datasets (as has the unused `Selectivity_block`). Existing models
  are bit-identical.

* **`Time_varying_q`, `Time_varying_sel`, and `M1_re` are soft-deprecated in favour
  of random-effect linkages.** `fit_mod()` now warns (naming the fleets/species and
  the grammar equivalent) when a model uses these legacy time-variation switches,
  pointing at `build_catchability()` / `build_selectivity()` / `build_M1()` with
  `(1 | Year)` / `rw(1 | Year)` / `ar1(1 | Year)` — which additionally allow a prior
  on, or free estimation of, the deviation SD. **The legacy switches still fit with
  their exact numerics** (they keep their own C++ path); this is only a nudge. The
  warning is not raised where no grammar equivalent exists yet: the environmental /
  Rogers-AR1 catchability modes (which overload `Time_varying_q` to name env
  columns), the non-parametric selectivity forms, and the separable `M1_re = 6`
  (age × year).

* **`fleet_control` columns `Q_prior` → `Q_init`, `Index_sd_prior` → `Index_sd`,
  `Catch_sd_prior` → `Catch_sd`.** These are start/input values, not priors (the
  prior on q lives in `Q_sd_prior`), so the names misled. As with the other
  renames, the old names are accepted and upgraded in place by `switch_check()` /
  `read_data()`, and the bundled datasets were regenerated — existing scripts keep
  fitting identically.

* **`fleet_control` columns `Time_varying_q_sd_prior` / `Time_varying_sel_sd_prior`
  renamed to `Time_varying_q_sd` / `Time_varying_sel_sd`.** They hold the input
  *value* of the time-varying deviate SD, not a prior on it (no density is placed
  on the SD), so the `_prior` suffix was misleading. The old names are still
  accepted — `switch_check()` and `read_data()` upgrade them in place with a
  one-time message — so existing data lists, saved xlsx files, and the bundled
  example datasets keep fitting identically. Update scripts to the new names at
  your convenience.

# Rceattle 4.9.0

## New features

* **Environmental linkages on catchability and selectivity.** The
  formula-driven linkage grammar (`linkage_spec()`) now extends to survey
  catchability and to the parametric selectivity forms, both indexed by a new
  `fleet` stratum (`by = ~ fleet`):

    ```r
    fit_mod(...,
      qFun   = build_catchability(linkages = list(
                 q = linkage_spec(~ temp, by = ~ fleet))),
      selFun = build_selectivity(linkages = list(
                 inf_asc = linkage_spec(~ cold_pool, by = ~ fleet))))
    ```

  Selectivity keys are the shared parameter slots `slp_asc`/`slp_desc`/
  `inf_asc`/`inf_desc`, with DoubleNormal aliases `sigma_asc`/`sigma_desc`/
  `peak`/`right_floor`. Supported forms: `Logistic`, `DoubleLogistic`,
  `DescendingLogistic`, `DoubleNormal`, `LogisticPM`. Every parameter accepts
  `link = "log"` (multiplicative) or `link = "identity"` (additive). A linkage
  on an unsupported form or the non-parametric `coff` errors at fit time
  naming the fleet, rather than being estimated to no effect.

* **`build_catchability()` and `build_selectivity()`** exported, mirroring
  `build_growth()` / `build_M1()` / `build_srr()`, each carrying a `linkages`
  argument for its process.

* Linkage formulas are now parsed with the **`reformulas`** package (shared by
  lme4 and glmmTMB), so `(1 | Year)` / `ar1(Year + 0 | fleet)` syntax is
  recognised and an unknown covariance-structure wrapper errors instead of
  silently degrading to unstructured.

## Behavior changes

* **`Catchability = "Environmental"` (`Estimate_q = 5`) is deprecated.** It
  still fits with its existing numerics but emits a warning pointing at
  `build_catchability(linkages = ...)`, which names covariates by formula and
  carries priors, bounds, and an estimation phase.

# Rceattle 4.8.0

## New features

* **AMAK "avgsel" base-level selectivity penalty** (`fleet_control$Sel_avgsel_pen`).
  Non-parametric (type 9 / `NonParametricPM`) fleets can now carry the AMAK/ebswp
  base-level regulariser `weight * (log(mean(exp(base coefficients))))^2` — ADMB's
  `10 * square(avgsel_*)` — which mildly pins the overall level of the base
  selectivity coefficients that the per-year mean-centering otherwise leaves free.
  The per-fleet weight defaults to `0` (off), so existing models and the `BS2017SS`
  golden reference are unchanged; set `Sel_avgsel_pen = 10` to match AMAK. The
  penalty is accumulated once per shared-selectivity block (on the lead fleet).

## Bug fixes

* **Non-parametric (AMAK "pm", type 9 / `NonParametricPM`) selectivity deviate
  penalties.** Two corrections so the type-9 selectivity and its deviate penalty
  reproduce ADMB/AMAK exactly when the fleet excludes a first bin (e.g. the
  acoustic-survey age-1) and/or starts after the model start year:
    - Excluded bins (below `Bin_first_selected`) are now held at 0 before each
      year's mean-centering instead of being carried through the random walk.
      Previously their log-selectivity accumulated the per-year centering offset
      and drifted, inflating the curvature / random-walk penalty on a bin that is
      zeroed out of the fit anyway.
    - For years up to a fleet's `Sel_start_year`, the curve is now rebuilt from the
      base coefficients each year rather than iterating the running random walk
      over the (data-free) pre-survey years. Iterating instead converged the
      excluded first bin to a different fixed point and perturbed the start-year
      base selectivity, inflating the deviate penalty on the first change year.
  A fleet that starts at `styr` with no excluded bins is unaffected (the reset
  reduces to the original start-year behaviour); the `BS2017SS` golden reference is
  unchanged.
* **`hessian_conditioning` diagnostic now always names the flat direction.** When
  the Hessian's least-identified direction was spread diffusely over many
  coefficients (no single coefficient above the reporting threshold), the message
  read `loads on: .` with nothing after it. The check now aggregates the
  eigenvector's squared loadings by parameter block and reports the block(s) making
  up the direction with their percentage share (e.g. `loads on: rec_dev (69%) +
  ln_srv_sel (31%)`), falling back to `par.fixed` names when `cov.fixed` carries no
  dimnames.

* **Mohn's rho was computed from the wrong row for forecast years beyond the
  terminal year.** `retrospective()` accumulated relative error into
  `mohns[ind, ]` while reading from `mohns[j, ]`; the two indices coincide only
  at forecast year 0, so retrospective bias was wrong for every additional
  forecast year. The observation counter on the adjacent line was correctly
  indexed, which masked the problem.

* **Projection-year conditional age-at-length was missing from every MSE.**
  `run_mse()` built `proj_caal` and then appended a different object
  (`proj_comp`) to a non-existent `proj_caal` field instead of appending the
  CAAL rows to `caal_data`. When the composition branch had not run,
  `proj_comp` did not exist at all.

* **`sim_mod()` reproduced the previous observation's draw instead of
  simulating.** The CAAL branch tested `Comp_loglike` (rather than
  `CAAL_loglike`) against integer codes, and neither the composition nor the
  CAAL branch handled `"MultinomialAFSC"` — the default `Comp_loglike`. In both
  cases no branch fired, leaving the simulation variable holding its value from
  the previous row. Unrecognised likelihood codes now raise an error rather
  than falling through.

* **`sample_rec()` and retrospective forecasts produced `NaN` recruitment
  deviations.** Both took `log()` of the mean of `log(R) - log(R_hat)`, a
  log-scale quantity centred near zero whose mean is routinely negative. The
  trend adjustment is now applied additively on the log scale, matching the
  sibling sampling branch.

* **`combine_data()` discarded merged environmental data** (it assigned the
  merge to its input rather than to the returned object) and appended two
  `NA`-named elements by iterating past the end of its column-name vector.

* **`read_data()` overwrote a user-supplied `fleet_control$Month`** with zero.
  The guard was inverted relative to its own message, so per-fleet survey
  months set in the workbook were silently reset.

* **`build_hcr_map()` disabled SPR reference points too eagerly.** It turned off
  `Ftarget`/`Flimit` when *any* fishery for a species had zero `proj_F_prop`,
  though its comment and message both describe the all-zero case.

* **`set_phases()` declared `log_M1` twice**, misaligning the positional
  phase-to-parameter pairing in `TMBphase()` for every subsequent parameter.

* **`fit_mod()` could fail with `object 'opt' not found`** for
  `estimateMode = 2` combined with `HCR = "NoFishing"`.

* **`Time_varying_q` environmental-index lists were destroyed.** For fleets with
  `Catchability = "Environmental"`, a comma-separated list of `env_data` column
  indices (e.g. `"1,3"`) was coerced with `as.integer()` and became `NA`, both in
  `fleet_control` itself and in the `index_varying_q` vector passed to TMB.
  `process_residuals()` tests that vector against a mode code, so the `NA`
  silently poisoned a branch rather than failing.

* **Natural-mortality environmental linkages used the wrong coefficients.**
  Assigning a scalar to the vector `beta_M1_tmp` broadcast it, so the dot
  product collapsed to `beta_last * sum(env)` rather than `sum_j env_j * beta_j`.
  Dormant in practice because `M1_beta` is mapped out by default.

* **Observation years earlier than `styr` produced a large negative array
  index.** Two consecutive `if` statements (rather than `if`/`else if`) applied
  both year transformations in sequence. Five sites.

* **Selectivity max-normalisation branched on an AD type**, freezing the
  location of the maximum at its initial-parameter position for the whole
  optimisation. If the peak moved during fitting, the curve was normalised by
  the wrong bin and the gradient with respect to the true peak was absent. Now
  folded with `max2()`, matching the neighbouring single-bin normalisation.

* **`plot_form()` errored on every call.** It plots the Kinzey & Punt (2009)
  functional responses, whose parameters (`logH_1`, `logH_1a`, `logH_1b`,
  `logH_2`, `logH_3`, `H_4`) are commented out in the TMB template,
  `build_params()`, `build_map()` and `build_bounds()`, and whose modes
  (`msmMode` 3-9) `data_check()` blocks. `params$logH_1` was therefore always
  `NULL` and the first line failed with "non-numeric argument to mathematical
  function". The function remains exported but now raises a clear message
  naming the supported alternatives; its body is retained, commented, so it can
  be revived with the C++ if these formulations are validated.

* **Predation suitability could become infinite.** The zero-guard tested
  `diet_prop_sum`, but the division is by `suma_suit + other_food_diet_prop`,
  which can reach zero or go negative independently.

## Behavior changes

* **`estimateMode = 3` now returns the real objective and gradient.** It
  previously shared mode 4's placeholder objective (`dummy^2`), so `obj$fn()`
  returned zero and every gradient element was identically zero — a build-only
  object appeared healthy while carrying no information. Mode 3 is the
  build-without-optimizing entry point (the analogue of WHAM's
  `fit_wham(do.fit = FALSE)` and SAM's `sam.fit(run = FALSE)`) and can now be
  used to inspect a model before committing to a fit. **`estimateMode = 4` is
  unchanged** and still returns the placeholder: it maps out every hindcast
  parameter, so it is a plumbing smoke test rather than a likelihood.

## Internal

* **Observed-proportion columns are now found by name rather than by counting
  from a fixed position.** `comp_data` and `caal_data` begin with identifying
  columns (fleet, species, sex, year, sample size) followed by the observed
  proportion at each age or length bin. `sim_mod()` located those proportions
  by a fixed offset (`[, 9:ncol(x)]`, `[, 7:ncol(x)]`), so adding or reordering
  an identifying column would have written simulated values into the wrong
  columns without any error. Both tables now resolve the proportion columns by
  name, as `rearrange_data()` already did.

* **Added an internal parameter dictionary** (`R/0-parameter_dictionary.R`)
  mapping each of the 43 TMB parameters to its natural-scale name, its process
  group, a one-sentence meaning, and its dimensions — so diagnostics and error
  messages can describe `log_sel_slp_dev` or `index_q_rho` in terms a reader
  can act on. A test enforces that the dictionary and the parameter list stay
  in exact correspondence.

* The non-parametric selectivity curvature penalty is computed in one pass
  rather than recomputing the full second difference for each age (O(n) instead
  of O(n^2)). Numerically identical.

# Rceattle 4.7.0

## New features

* **Natural-scale normal survey-index likelihood.** `fleet_control$Index_loglike`
  gains `"Normal"`: the index residual `obs - q*pred` is normal with an *absolute*
  observation SD (the `index_data$Log_sd` column is read as a natural-scale SD, not
  a log-scale CV), i.e. `0.5*(obs - q*pred)^2 / sd^2`. This matches the AMAK/ebswp
  `avo_like` / `cpue_like` survey likelihoods (which use an absolute sigma) so those
  indices can be reproduced exactly, rather than approximated by the default
  lognormal. Fixed alongside a latent crash: the covariance (MVN/MVNORM) block was
  gated on `Index_loglike >= 1`, which now also matched `"Normal"` and applied
  `MVNORM()` to a 1x1 dummy Sigma; it is now restricted to `"MVN"` / `"MVNORM"`.
* **Multivariate-normal (covariance) survey-index likelihood.** A survey/index
  fleet can now use a full variance-covariance matrix for its biomass index
  instead of independent lognormal errors, via the new `fleet_control` column
  `Index_loglike` (`"Lognormal"`, the default, or `"MVN"` / `"MVNORM"`). Supply
  the covariance matrix (e.g. a VAST-derived Sigma) as a named element of the new
  `index_cov` data list, keyed by `Fleet_name`. The likelihood uses TMB's native
  `density::MVNORM(Sigma)` on the natural-scale residual `r = obs - q*pred`:
  `"MVN"` reports the bare quadratic form `0.5 * r' Sigma^-1 r` (the AMAK/ebswp
  `DoCovBTS` bottom-trawl survey value, matching ADMB's reported likelihood),
  while `"MVNORM"` reports the full normalized density
  `0.5 * (r' Sigma^-1 r + logdet(Sigma) + n*log(2*pi))` — the two give an
  identical fit and differ only by a fixed constant. A companion catchability
  option `Catchability = "AnalyticalArith"` gives the arithmetic-mean analytical q
  (`mean(obs)/mean(pred)`) that the AMAK covariance survey uses (the existing
  `"Analytical"` q remains the geometric mean). Defaults are fully back-compatible:
  existing models get `Index_loglike = "Lognormal"` and are numerically unchanged.
* **AMAK-style non-parametric and logistic selectivity forms.** Two new
  `fleet_control$Selectivity` options reproduce the ADMB AMAK ("pm") selectivity:
  `NonParametricPM` (Ianelli coefficient selectivity with the decreasing,
  curvature, and deviation-magnitude penalties, whose weights are set by
  `Sel_curve_pen1` / `Sel_curve_pen2` / `Sel_curve_pen3`) and `LogisticPM`
  (logistic + a free age-1 selectivity). Both take `Time_varying_sel =
  "RandomWalk"`. Four new per-fleet `fleet_control` columns tune the shape
  penalty: `Sel_pen_first_bin` and `Sel_pen_last_bin` (the bin range of the
  adjacent-bin shape penalty, letting it span a narrower range than the first
  selected bin), `Sel_shape_mode` (`"Directional"` one-sided decreasing/increasing,
  the AMAK default, or `"Smooth"` two-sided curvature), and `Sel_cap_bin`
  (hold the realized non-parametric curve flat at/after a bin). Each is given on
  the fleet's own selectivity dimension — an age for age-based fleets, a 1-based
  length-bin ordinal for length-based fleets — and is range-checked accordingly.
  All default to the previous behavior when unset.
* The plotting functions have been overhauled to **ggplot2**. Every exported
  `plot_*()` function now builds its figure with ggplot2 (colourblind-safe
  palettes — the Okabe-Ito qualitative palette for series identity and viridis
  for ordered magnitude such as year — and `theme_classic`) and **returns the
  `ggplot` object**, so
  plots print when called interactively and can be further customised
  (`plot_biomass(fit) + ggplot2::ggtitle(...)`) or saved with `ggplot2::ggsave()`.
  Pass `file = "stem"` to also write the figure to `stem_<plot>.png`. Plot
  conventions follow r4ss / WHAM / SAM (line + 95% ribbon time series; observed
  points + error bars + predicted line for index/catch fits; year-coloured
  at-age curves for selectivity and mortality surfaces).
* `plot_index()` gains a `log` argument for the log-scale survey-index fit.
* The time-series plotters (`plot_biomass()`, `plot_ssb()`, `plot_recruitment()`,
  and the other `plot_timeseries()` wrappers) again honour user-supplied
  `line_col`, `lwd`, and `lty`: pass a colour, line width, and/or line type per
  model to override the defaults
  (e.g. `plot_biomass(list(m1, m2), line_col = c("black", "red"), lty = c(1, 2))`).
  `lwd` keeps the base-graphics convention where the default (3) renders as a
  standard-weight line.
* `plot_stock_recruit()` adds a 95% data ellipse of the SSB–recruitment cloud
  (`add_ci`, default `TRUE`).
* The test suite is reorganised into a flat, navigable `tests/testthat/` (see
  the new `tests/testthat/README.md`), runs fully under continuous integration,
  and has automated code-coverage reporting.

## Breaking changes

* `plot_*()` functions now **return a `ggplot` object** instead of drawing as a
  base-graphics side effect (returning `NULL`). Scripts that only called them
  for the side effect still work (the object prints); scripts that depended on
  the base-graphics device state or on a `NULL` return may need updating.
* `plot_logindex()` has been **removed**; use `plot_index(..., log = TRUE)`.
* The `gplots` and `oce` dependencies have been dropped (no longer used).
* **Time-varying non-parametric selectivity now uses `"RandomWalk"`, not
  `"IID"`.** The Ianelli non-parametric form (`Selectivity = "NonParametric"`)
  previously fired its time-varying coefficient deviations on `Time_varying_sel =
  "IID"`; it now requires `"RandomWalk"` (matching the random-walk structure the
  penalty implements) and rejects `"IID"` with an error at `build_map()`. A model
  using `NonParametric` + `IID` must switch to `RandomWalk`.
* **A selectivity shared across fleets is now penalized once.** When two or more
  fleets share a `Selectivity_index` *and* selectivity type (a mirrored curve),
  the selectivity shape/deviation penalty is accumulated only on the lead fleet
  rather than once per fleet, matching ADMB. This changes the objective only for
  models that both mirror a selectivity **and** put a penalty on it (non-parametric
  or time-varying); models with a unique selectivity index per fleet, and all
  bundled examples, are numerically unaffected.

## Bug fixes

* **`Catchability = "AnalyticalArith"` left an unused free catchability parameter,
  making the Hessian singular.** The arithmetic-mean analytical q solves q from the
  data (like the geometric `"Analytical"`), so its `index_log_q` is never used —
  but `build_map()` excluded only `"Analytical"` from estimation, so the
  `AnalyticalArith` fleet's `index_log_q` was still freed. That parameter never
  entered the objective, leaving a zero-gradient flat direction that prevented
  `sdreport()` from inverting the Hessian (`pdHess = FALSE`). It is now mapped
  out, and such models converge with an invertible Hessian.
* **The catchability prior and deviate penalties were counted once per fleet
  sharing a `Q_index`, not once per estimated parameter.** Fleets sharing a
  `Q_index` estimate one catchability and one deviate vector, but the template
  looped over every fleet, so a mirrored pair applied the `Q_prior` twice to the
  same parameter — tightening an intended prior SD of 0.2 to 0.2/sqrt(2) — and
  penalized the shared `index_q_dev` vector once per fleet for the IID, random
  walk and AR1 forms. A new `flt_q_lead` (the catchability analogue of
  `flt_sel_lead`) accumulates them on one fleet per group. Models whose fleets
  all have distinct `Q_index` are unchanged.
* **`Catchability = "PowerEquation"` is now rejected as not yet implemented.** It
  was accepted as a valid switch, but the power coefficient (`index_q_pow`) is
  not built as a parameter and the template does not apply it, so the fleet
  silently got a plain estimated q. `data_check()` now errors, matching how
  length-based `suitMode` values are handled.
* **`flt_sel_lead` could put the selectivity penalty on an `Off` fleet.** The
  lead was the first fleet in a `Selectivity_index` group by row order. When that
  fleet was `Fleet_type = "Off"` the template's `flt_type > 0` gate then skipped
  the penalty for the whole group, leaving the shared selectivity unpenalized.
  The lead is now the first *estimated* fleet in the group, matching the map
  donor.
* **A mirrored group led by an `Off` fleet stopped estimating selectivity and
  catchability.** `adjust_map_shared_params()` copied the first sharing fleet's
  map slice onto the rest. When that fleet was `Fleet_type = "Off"` its slice is
  all `NA`, so every fleet sharing the index silently had its selectivity /
  catchability parameters fixed at their starting values. The donor is now the
  first *estimated* fleet in the group, and the copy is skipped when the group
  has none. Groups led by an estimated fleet are unchanged.
* **An explicitly set `Sel_start_year` was not shared across mirrored fleets,
  making the fit depend on `fleet_control` row order.** The default derived from
  the data is already the earliest first-observation year across each
  `Selectivity_index` group, but a value the user sets directly was used
  per fleet. Since fleets sharing an index share one deviation block,
  `adjust_map_shared_params()` then overwrote the mirrored fleet's mask with the
  lead fleet's, so whichever fleet appeared first governed the group: when that
  fleet started later, a sharing fleet with earlier data silently lost those
  deviations (12 years in a 1982/1994 pair). `Sel_start_year` now resolves to the
  group minimum for both the map mask and the template's penalty anchor, however
  it was set. `data_check()` warns when a mirrored group has differing
  `Sel_start_year`, and `build_map()` warns when `Bin_first_selected` or
  `N_sel_bins` differ within a group (those are likewise taken from the lead
  fleet). Unmirrored fleets and derived defaults are unchanged.
* **`fleet_control$Fleet_code` is now required to equal the row number.** It is
  used directly as the fleet slot of the per-fleet parameter and map arrays,
  which are built in `fleet_control` row order; a mismatch silently attached
  parameters to the wrong fleet. `data_check()` now rejects it, and the
  remaining places that read a `fleet_control` column by `Fleet_code` instead of
  row index were corrected.
* **Selectivity bin columns were converted to model indices using the species'
  `minage`, which is wrong for length-based fleets.** `Sel_norm_bin1`,
  `Sel_norm_bin2`, and the shape-penalty range/cap columns subtracted
  `minage` to get the 0-based template index. That is correct for an age-based
  fleet (the value is an age), but a length-based fleet's value is a 1-based
  length-bin ordinal and must be offset by 1 — so those columns silently pointed
  at the wrong length bin whenever `minage != 1`. The offset is now chosen per
  fleet from `Selectivity_dimension`. Age-based fits, and length-based fits with
  `minage == 1` (where the two offsets coincide), are unchanged.
* **Non-parametric selectivity penalties now span length bins for length-based
  fleets.** The shape, curvature, and random-walk penalties for the
  `NonParametric` (type 2), `NonParametricPM` (type 9), and `LogisticPM`
  (type 11) forms were hard-coded to the number of ages, so a length-based fleet
  (`Selectivity_dimension = "Length"`) with more length bins than ages left the
  upper length bins unpenalized (and read the wrong array). The penalties now run
  over the fleet's own selectivity dimension (`nlengths` for length-based,
  `nages` for age-based), so every parametric and non-parametric selectivity form
  works on both age and length. Age-based fits are numerically unchanged.
* **Survey-index covariance matrices were not re-aligned when the fitted year
  range changed.** An `index_cov` (MVN/MVNORM) Sigma is positionally keyed to a
  fleet's fitted survey observations, so any workflow that changes that set —
  a `retrospective()` peel, an `endyr` / `styr` subset, or a `run_mse()`
  assessment step that appends survey observations — left the Sigma at its
  original dimension and tripped `rearrange_data()`'s dimension check
  (`"N x N but the fleet has M fitted survey observations"`). `clean_data()` now
  tags each Sigma with its fitted years the first time it is seen and, on every
  subsequent pass, re-keys it to the current fitted set: retained years keep
  their full covariance block, and new (future/simulated) years are added as an
  independent diagonal block with variance `(Observation * Log_sd)^2`. Because
  every re-fit routes through `clean_data()`, retrospective, MSE, and jitter now
  all work with covariance-survey models; fresh fits and non-MVN fleets are
  numerically unchanged.
* **Time-varying selectivity deviations were estimated before a fleet had any
  data.** `build_map()` never consulted `fleet_control$Sel_start_year`, so a fleet
  with time-varying (`"RandomWalk"`) selectivity had deviations estimated across
  *every* hindcast year — including years before its first observation. Those
  deviations are informed by no data and constrained by no penalty (every
  selectivity penalty in the objective is anchored at `Sel_start_year`), leaving
  unidentified flat directions that stall the optimizer. In the EBS pollock model
  this left ~54 free pre-survey deviations on a bottom-trawl survey starting in
  1982 and ~240 on an acoustic survey starting in 1994 — a total parameter count
  of 1483 against 1225 for the equivalent ADMB (AMAK) model, which never declares
  them. Deviations before `Sel_start_year` are now fixed at 0, giving 1249. The
  two selectivity parameterizations differ in where the base curve lives and are
  handled accordingly: `LogisticPM` (and other curve-based forms) estimate a
  separate base, so deviations are fixed *through* the start year; the
  non-parametric random walk maps its mean (`sel_coff`) off and lets the start-year
  deviation carry the base, so only the deviations strictly *before* it are fixed.
* **`Sel_start_year` now defaults to the fleet's first year of data** rather than
  `styr`. It is an optional `fleet_control` column, so the fix above only took
  effect for users who knew to set it — a model with a late-starting survey would
  silently carry unidentified deviations. The default is derived from
  `catch_data` / `index_data` / `comp_data` / `caal_data`, consistent with how
  `switch_check()` already auto-`"Off"`s fleets with no observations. Fleets
  sharing a `Selectivity_index` share one selectivity curve, so the start year is
  the *earliest* first-observation year across the whole group: a fleet whose own
  data starts late but which mirrors an earlier fleet's curve (e.g. an AVO index
  starting in 2006 mirroring an acoustic survey starting in 1994) must not drop
  the deviations the mirrored fleet's data informs. Set the column explicitly to
  override. Only models with time-varying selectivity on a fleet whose data starts
  after `styr` are affected; for those, the previous behaviour is recovered by
  setting `Sel_start_year = styr`.
* **`LogisticPM` selectivity started with an unusable age-1 selectivity.**
  `build_params()` initializes `sel_inf[2]` to `10`, which is correct for its
  usual meaning (the descending-limb inflection *age*), but `LogisticPM`
  (type 11) repurposes that slot as the free first-bin (age-1) **log**-selectivity.
  The shared default therefore started age-1 selectivity at `exp(10)` = 22026
  (for reference, the ADMB AMAK "pm" estimate is `sel_age_one_bts` = -3.19, i.e.
  0.04). When age-1 was selected (`Bin_first_selected = 1`) this swamped the
  survey-index prediction and the optimizer diverged with a gradient blow-up on
  `sel_inf`; when it was not selected the bad value was silently masked by the
  zeroed first bin. `sel_inf[2]` now defaults to `0` (age-1 selectivity = 1) for
  `LogisticPM` fleets.
* `plot_maturity()` read a non-existent `pmature` field and errored on real
  fits; it now reads `data_list$maturity`.
* `plot_ration()` failed for single-sex models (a dropped array dimension); the
  sex dimension is now indexed explicitly.
* `plot_stock_recruit()` drew the mean-recruitment reference line a factor of
  1e6 too high (recruitment points are in millions but the line was not scaled),
  which under ggplot's free y-scale collapsed the SSB–recruitment cloud to the
  axis; the line is now scaled to match the points.
* `plot_index()`, `plot_catch()`, and `plot_indexresidual()` drew
  prediction-only rows (`Year < 0`, excluded from the likelihood) as if they
  were fitted observations; these rows are now omitted.
* `plot_mortality(M2 = TRUE)` labelled the y-axis "M1 + M2" while plotting M2;
  the label now reads "M2".
* `osa_residuals()` now includes fleet in the one-step-ahead conditioning order
  (source, year, **fleet**, then bin), making the sequence fully deterministic
  when several fleets report in the same year. Previously the within-year,
  multi-fleet order fell to the incidental row order of the observation table.
  This does not change results for fixed-effects fits (`random_rec = FALSE`),
  where OSA residuals are order-invariant, and ascending fleet order matches
  WHAM's within-stage sequencing so the WHAM cross-check is unaffected; for
  random-effects fits it makes individual residuals reproducible (the overall
  N(0,1) properties were already order-invariant).

# Rceattle 4.6.0

## New features

* `osa_residuals()` now builds the composition / conditional-age-at-length /
  diet one-step-ahead observation data on demand, so OSA residuals can be
  computed from any fit. Previously the composition OSA data had to be assembled
  during the fit via `fit_control(osa = TRUE)`; that pre-build is no longer
  required (the reconstruction reads only the `*_ctl` / `*_obs` arrays already
  carried by the fitted object, and is identical to an `osa = TRUE` fit).
  Aggregate index/catch OSA residuals are unaffected.
* **Breaking:** the `osa` argument to `fit_control()` has been removed; it is no
  longer needed (see above), and passing it now raises an "unused argument"
  error. The `$osa` element is no longer stored on fitted `Rceattle` objects.
* `fit_control()` gains `bias_adjust_obs` and `bias_adjust_proc` (both default
  `TRUE`) to toggle the lognormal bias correction (`-sigma^2/2`) applied to the
  observation (index / catch) and process (recruitment, initial abundance, and
  the Beverton-Holt steepness prior) likelihoods, respectively. The defaults give
  the standard mean-unbiased behaviour (`E[R] = R0`); set a flag to `FALSE` to
  drop the correction for that likelihood group. The values are passed to the TMB
  template as `DATA_SCALAR`s.

* Added one-step-ahead (OSA) residuals for model validation via the new
  exported `osa_residuals()` (Thygesen et al. 2017; Trijoulet et al. 2023).
  Unlike Pearson residuals, OSA residuals are iid standard normal under a
  correctly specified model even for correlated composition data. All fitted
  observation types are supported: survey indices, fishery catch, age/length
  composition, conditional age-at-length, and predator diet (stomach-content)
  composition. Composition data are residualized with the conditional binomial
  / beta-binomial decomposition (multinomial / Dirichlet-multinomial), and the
  `MultinomialAFSC` likelihood is residualized under the full multinomial. The
  fitted objective is unchanged: a new `obsvec`/`keep`/`osa_mode` machinery
  feeds `TMB::oneStepPredict()` while leaving normal fitting bit-for-bit
  identical. The `oneStepPredict()` call is split by observation type so the
  continuous (lognormal) index/catch series and the composition series can be
  residualized with the correct settings; `discrete = TRUE` treats composition
  as discrete (randomized quantile residuals; Dunn and Smyth 1996) while the
  aggregate series stay continuous.
* `fit_control()` gains `comp_offset`, the small proportion offset added to the
  observed and predicted age/length composition, conditional-age-at-length, and
  predator diet (stomach-content) bins before the multinomial /
  Dirichlet-multinomial likelihood (and to the matching OSA observation vector,
  so fitting and OSA residuals use the same offset). It applies to every
  composition likelihood, including the "Martin's" (`comp_ll_type = -1`) form. It
  defaults to `1e-5` (the historical CEATTLE value, which avoids `log(0)` for
  empty bins); set `comp_offset = 0` for a standard WHAM-style multinomial. The
  value is stored on `data_list$comp_offset` (filled by `switch_check()`), so
  internal re-fits (projections, `retrospective()`, `jitter()`, `run_mse()`)
  inherit it; `fit_control(comp_offset = ...)` overrides the stored value.
* `osa_residuals()` gains a `parallel` argument (default `TRUE`) that computes
  the per-observation one-step-ahead loop with `parallel::mclapply()` -- a
  near-linear speedup for models with random effects, where each observation
  triggers a Laplace re-evaluation (set `options(mc.cores = )` to choose cores;
  serial fallback on Windows). The discrete randomized-quantile path stays serial
  so it remains reproducible given `seed`.
* Added `osa_diagnostics()` for the Stewart and Monnahan (2025) statistical
  diagnostics -- the standard deviation of the normalized residuals (SDNR) and
  the lower/upper tail statistics, each with its standard-normal null interval
  -- and a `plot()` method for `rceattle_osa` objects, styled after the
  NOAA-AFSC `afscOSA` package. The plot draws up to two figures: an *aggregate*
  Q-Q figure for the index/catch series (which have no age/length bin, so no
  bubbles), and a *composition* figure pairing the Q-Q plot with signed
  OSA-residual and Pearson-residual bubble plots, with age-based bins (age
  composition and conditional age-at-length) in the left column and length-based
  bins in the right column, each on its own bin axis. Panel headers use the
  fleet name from `fleet_control`. The `plot()` method takes `source` and
  `species` arguments to subset the data shown (mirroring `residuals()`), and
  `combine = FALSE` to draw the age and length composition as separate figures.
  `osa_residuals()` now also carries a `fleet_name` column and (for composition)
  the matching Pearson residuals.
* Added `process_residuals()` for SAM-style process residuals on the model's
  random-effect deviations (recruitment, initial abundance, and catchability),
  validating the process model as a complement to the observation residuals.
* `residuals()` gains `type = "osa"` and `type = "process"`; `plot_comp()` and
  `plot_indexresidual()` gain a `residual_type = "osa"` option that draws the
  OSA diagnostics through the familiar plotting functions.
* `plot_comp()` was re-implemented in ggplot2 for a consistent look with the OSA
  plots: composition Pearson-residual bubbles plus observed-vs-fitted annual and
  aggregated composition figures. The observed area and fitted line now span only
  the observed bins (they no longer extend past the first / last bin), zero
  observed proportions are retained (only `NA` is dropped), joint-sex (Sex == 3)
  data are drawn on a single age/length axis with females above and males below
  zero, and a `species` argument subsets the species shown.

## Bug fixes

* Fixed the unweighted conditional-age-at-length (CAAL) log-likelihood being
  recorded in the composition slot of `unweighted_jnll_comp` instead of the
  CAAL slot. This affects the reported diagnostic matrix only; the fitted
  objective was unaffected.

## API

* `projection_uncertainty` moved from `fit_mod()` into `fit_control()`,
  consolidating it with the other optimizer / reporting controls. Passing it
  directly to `fit_mod()` still works but emits a deprecation warning and
  forwards the value into `fit_control()` -- the same backward-compatible path
  used by the other former `fit_mod()` control arguments (`phase`, `getsd`,
  ...).
* `residuals()` now follows the `stats::residuals.glm()` convention: `type`
  selects the residual *kind* -- `"response"` (default), `"pearson"`, `"osa"`,
  or `"process"` -- and a `source` argument selects the data source(s)
  (`"index"`, `"catch"`, `"comp"`, `"caal"`, or `"all"`, the default), so by
  default residuals for every applicable source are returned. A `species`
  argument subsets to particular species. Pearson residuals are now available
  for the aggregate index/catch series (standardized by the realized observation
  log-SD) and for predator diet via `source = "diet"` (returned on its own in a
  predator/prey schema); `plot_diet_comp()` now draws its diet residuals from
  this single `residuals()` path. `type` selects the residual kind only; data
  sources are selected with `source`.
* The `source` argument is shared across `residuals()`, `osa_residuals()`, and
  `plot()` (it replaces the earlier `types` argument of `osa_residuals()`),
  accepting `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, and `"all"`, so
  the three entry points select data sources with one consistent vocabulary.

## Documentation

* Reviewed and revised the user vignettes for accuracy, clarity, and
  consistency, targeting assessment scientists rather than developers. Trimmed
  developer-oriented internals and repetition, and corrected the option tables
  so they match the code: `estimateMode`, `suitMode`, selectivity, catchability,
  harvest-control-rule, `M1_model`, and bioenergetics (`Ceq`) codes now agree
  across `model-options-and-functionality`, `single-vs-multispecies`,
  `projections-and-reference-points`, `data-without-excel`, and
  `model-parameterizations`.
* Clarified that the length-based suitability modes (`suitMode = 1/3/5`) are not
  yet enabled; use a weight-based mode (`2/4/6`) for parametric suitability.
* Standardized data-structure terminology across vignettes: weight-at-age is in
  kilograms, and the diet / stomach-content input is `diet_data` (removing the
  legacy `UobsWtAge` / `UobsAgeWt` / `stom_prop_data` names).
* Corrected the `residuals()` examples in the README to the `type` / `source`
  convention (e.g. `residuals(fit, type = "pearson", source = "comp")`).
* Added a website-only "Developer guide" article
  (`vignettes/articles/developer-guide.Rmd`) documenting the fit pipeline, the
  switch system, the TMB module layout, and recipes for adding parameters, data
  inputs, switch options, and likelihood components. Consolidates and updates
  the GitHub wiki developer notes.

# Rceattle 4.5.0

## New features

* Added a general post-fit convergence diagnostics framework via the new
  exported `convergence_diagnostics()` function. It runs a battery of checks
  covering the optimizer gradient, Hessian positive-definiteness and
  conditioning, parameters on bounds, phasing, and parameter estimability, and
  returns a structured `"Rceattle_convergence"` object whose `status` reflects
  the worst severity found (`"OK"`, `"NOTE"`, `"WARN"`, or `"FAIL"`).
* `fit_mod()` now runs the convergence battery automatically and attaches the
  result as `fit$convergence`; `convergence_diagnostics()` can also be called
  directly to re-run it on any fit.
* Added a model-diagnostics vignette section and accompanying unit tests for
  the new framework.

## Bug fixes

* Fixed a parameters-on-bounds misalignment in the convergence diagnostics.
* The joint negative log-likelihood (jnll) is now taken from the optimizer
  objective for consistency.

## Internal / R CMD check

* Moved the suppression of the spurious GCC/Eigen `-Warray-bounds` false
  positive from a `-Wno-array-bounds` compiler flag to a source-level
  `#pragma GCC diagnostic ignored "-Warray-bounds"`, so
  `R CMD check --as-cran` no longer warns about non-portable,
  warning-suppressing compilation flags. Real diagnostics still surface.
* Additional C++ cleanup to remove compiler warnings.

# Rceattle 4.4.2

## Code organisation (no change to fitted results)

The pre-fit pipeline files in `R/` were reorganised so they are easier to
navigate. None of these changes alter model output.

* **File prefixes now follow execution order.** `fit_mod()` runs its stages
  as `clean_data()` -> `data_check()` -> `build_params()` -> `build_map()` ->
  `build_bounds()` -> `rearrange_data()` -> fit -> `rename_output()`, so the
  files were renumbered to match (`data_check` is now `1-`, `build_params`
  `2-`, `build_map` `3-`, `build_parameter_bounds` `4-`). A pipeline map was
  added to the top of `fit_mod()`.
* **Switch lifecycle consolidated** into a single `R/0-switches.R`: the
  string<->integer maps (formerly `0-constants.R`) plus `switch_check()`,
  `revert_switches()`, `validate_switches()`, and `convert_switches()`, with a
  header documenting the order in which they run.
* **HCR helpers co-located.** `build_hcr_map()` moved into `R/0-build_hcr.R`
  alongside `build_hcr()` (the separate `2-build_hcr_map.R` was removed).

## Rename / deprecation

* `rearrange_dat()` is renamed **`rearrange_data()`**. The old name still works
  as a deprecated alias (emits a one-time `.Deprecated()` warning) and will be
  removed in a future release.

## Export hygiene

* `check_composition_data()`, `check_caal_data()`, `calc_mcall_ianelli()`, and
  `calc_mcall_ianelli_diet()` are no longer exported; they are internal
  helpers called only from within the package.

## Internal / R CMD check

* Removed `Rceattle:::` self-references in `build_bounds()` (a package should
  not use `:::` for its own objects).
* `profile.Rceattle()` gained `...` for S3 consistency with the
  `stats::profile` generic.

# Rceattle 4.4.1

## Rename

`profile_param()` was renamed `profile()` and turned into an S3 class.

# Rceattle 4.4.0

## New features

* **Double Normal Selectivity (Type 8)**: Added support for the four-parameter Double Normal selectivity curve (Peak, Ascending SD, Descending SD, and Floor). This includes full support for annual deviates on all four parameters and is compatible with the age- and length-based selectivity engines.
* **Growth SD Control**: The linkage system now supports `sd_L1` and `sd_Linf` parameters. This allows users to specify priors, initial values, and bounds for growth SD endpoints using the same formula-driven interface used for mean growth parameters.
* **Natural-scale linkage API**: Renamed the linkage parameter keys from `log_*` to their natural-scale counterparts (`K` / `L1` / `Linf` / `m`; `M1`; `R0` / `alpha` / `beta`). Internal parameters remain on the log scale.
* **Natural-scale priors**: Standardized prior evaluation so that probability densities are applied to parameters on their natural scale by default, unless a lognormal family is explicitly requested.
* **Natural-scale inits**: Standardized init evaluation so initial values for `"(Intercept)"` parameters are applied on their natural scale.
* **Linkage link functions**: Fully implemented dual-path linkage offsets. `link = "log"` (the new default) applies the offset multiplicatively on the natural scale (additive on the log scale); `link = "identity"` applies it additively on the natural scale.
* **Per-species VB anchor age**: `build_growth()` gains a `growth_age_L1` argument (scalar or length-`nspp` vector) for the age at which mean length equals `L1`. Matches SS3's `Growth_Age_for_L1`. Default `NA` inherits `data_list$growth_age_L1` if set, else falls back to `max(0.5, minage[sp])` so `minage = 0` models pick up an SS3-consistent half-year anchor while `minage >= 1` models stay backwards-compatible.
* **Self-test simulation**: New `self_test()` simulates `nsim` datasets from a fitted model and re-fits the model to each simulated dataset, returning the list of refits. Runs in parallel by default (PSOCK cluster, capped at 2 cores under `R CMD check`) with a per-simulation `seed + i` so results are reproducible under both sequential and parallel execution.
* **Likelihood profile**: New `profile_param()` generalises the legacy `profile_rsigma()` example helper to any parameter slot in `Rceattle$estimated_params`. Supports arbitrary N-D cross-profiles via `expand.grid()` over a list of per-cell value grids (e.g. cross-profiling `log_M1` across sex, or `R_log_sd` across multiple species). Natural-scale aliases `"sigmaR"` / `"R_sd"`, `"M1"`, and `"R0"` / `"alpha"` / `"beta"` apply the implicit `log()` transform and (for the `rec_pars` aliases) auto-fill the column, so `slots` only needs the species index. `slots` defaults to species 1 with a warning. Fits run in parallel on the same PSOCK harness as `jitter()` / `retrospective()` / `self_test()`.
* **Parallel `retrospective()` and `jitter()`**: Both diagnostics now run their independent peels / starts on a PSOCK cluster (same approach as `run_mse()`). New `cores` argument on each (default `parallel::detectCores() - 6`, capped at 2 when `_R_CHECK_LIMIT_CORES_` is set); pass `cores = 1` to force sequential execution.
* **Standard errors in `as.data.frame.Rceattle()`**: The tidy long-format frame now carries a `se` column alongside `value` / `lwr` / `upr`, populated from the TMB `sdreport` for any `ADREPORT`'d quantity. Set to `NA` for non-ADREPORT'd quantities and for fits produced with `getsd = FALSE`.

## Bug fixes

* **SRR logic**: Fixed a bug in `build_srr()` where the `Bmsy_lim` penalty was incorrectly disabled for current Ricker implementations due to an index mismatch.
* **Selectivity RW prior scaling**: Corrected the random walk prior scaling in the TMB template to ensure consistent $4 \times$ SD multipliers for both ascending and descending limb slope/SD parameters.
* **`last_par` returned wrong vector**: Fixed `fit_mod()` so the value stored on the returned fit is the optimizer's last parameter vector rather than a stale prior reference, removing the need for the surrounding `try()` guards in downstream callers.
* **Growth at `minage = 0`**: Fixed a segfault in `growth.hpp` when `minage = 0` by guarding `current_age`, `age_L1`, and the cohort-boundary `age_L1_ceil` against the zero-age anchor. Also corrected length-at-age at `minage = 0` so the closed-form anchor at `L1` is honored on both the within-year and cohort recursion paths.
* **Length-bin midpoint for weight-at-age**: `calculate_weight()` now computes each bin's midpoint as `(lengths[ln] + lengths[ln+1]) / 2` (with the plus-group extended by half the final interior width) rather than `lengths[ln] + (lengths[1] - lengths[0]) / 2`. The previous formula assumed uniform bin widths and silently mis-located the midpoint for non-uniform length bins; for uniform bins (the common case) the two are algebraically identical, so existing fits are unchanged.
* **Plus-group SD-at-age pinned to the upper anchor**: both growth builders in `growth.hpp` now pin the oldest age class's `length_sd` to `exp(growth_log_sd(sp, sex, 1))` (SD at `Linf`), matching WHAM (`SDAA` plus group `= SD_len(1)`), instead of the length-based interpolation between SD(`L1`) and SD(`Linf`) used previously. This makes the plus-group convention consistent across von Bertalanffy and Richards growth and across both builders; expect a small numerical change in `length_sd` at the oldest age relative to prior versions.
* **`fit_mod()` bounds ordering**: `fit_mod()` now indexes parameter bounds by name rather than positional order when assembling `L` / `U` for `nlminb`. Previously, when `map$mapFactor` and `bounds$lower` were not in identical order, parameters could be paired with another parameter's bounds, producing silently wrong constraints. `start_par` is now also subset by name with `drop = FALSE`.
* **`mse_summary()` Tier-3 Flimit**: The internal `flimit_tier3_fun()` returned `Flimit` (its argument) instead of the depletion-adjusted `tier3_flimit` it had just computed, so the Tier-3 (HCR = 5) branch of P(F > Flimit) reduced to the base-Flimit check. Now returns the adjusted vector.
* **`mse_summary()` HCR coercion**: `HCR` is now normalized to its integer code before downstream comparisons (`HCR == 5`, etc.). `build_hcr()` accepts either an integer or a string alias (e.g. `"NPFMC"`); `mse_summary()` previously assumed integer form and silently produced wrong status flags when fits carried the string form.
* **`mse_summary()` OM status at assessment years**: P(F > Flimit) and P(SSB < SSBlimit) are now reported for the OM evaluated at the same assessment years as the EM (previously only the EM's perceived status was returned), and the SSB-limit threshold dispatch is consolidated in one helper so the Tier-3 / Category-1 / dynamic-vs-static cases stay aligned across the EM and OM paths.
* **`clean_data()` inactive-fleet handling**: The auto-Off branch no longer nulls out `proj_F_prop` and `Catchability` on fleets it flips to `"Off"`. The downstream TMB code already ignores those columns for off fleets, and clearing them lost information users needed when re-enabling a fleet.
* **Selectivity block indexing**: Renamed the within-loop `biom_yrs` index vector to `block_yrs` in the selectivity and catchability block branches of `build_map_*()`. The previous name was a leftover from the index-data path and shadowed nothing, but read as if it referred to biomass-observation years.
* **`plot_index()` / `plot_catch()` / `plot_logindex()` warnings**: Wrapped the internal `gplots::plotCI()` calls in `suppressWarnings()` so plotting a fit no longer prints the recurring `"In arrows(...): zero-length arrow is of indeterminate angle and so skipped"` noise when CI half-widths are zero.

## Data checks

* **Empirical growth + CAAL**: `data_check()` now errors when `growth_model == 0` (empirical weight-at-age) is combined with non-empty `caal_data` for a given species. The C++ growth matrix is not populated from the age-transition matrix in the empirical branch, so `pred_CAAL` collapses to ~0 and the multinomial NLL becomes uninformative.
* **Selectivity identifiability**: Fleets with estimated `Selectivity` and `Fleet_type != "Off"` now require at least one positive-sample `comp_data` or `caal_data` row. Either provide composition / CAAL data, mark the fleet as `Selectivity = "Fixed"` with `emp_sel`, or set `Fleet_type = "Off"`.
* **Auto-Off inactive fleets**: `clean_data()` automatically flips `Fleet_type` to `"Off"` for fleets that carry no catch or index observations, preventing the optimizer from drifting on unconstrained selectivity / catchability blocks.
* **`minage` guard**: `data_check()` errors when any species has `minage < 0`.

## Documentation
* Added `vignette("environmental-linkages-and-priors")` (and updated `_pkgdown.yml`) to cover the new linkage intercept behavior, link-function semantics, growth SD endpoints, and Double Normal selectivity.
* Updated all cross-references in `build_srr()` / `build_M1()` / `build_growth()` (deprecation warnings, soft-deprecated arg docs, and the `vignette()` pointers in the model-options vignette) from the old `environmental-linkages` slug to the renamed `environmental-linkages-and-priors` vignette so the soft-deprecation warnings now resolve.
* Expanded the Double Normal (selectivity type 8) doxygen / roxygen so the four estimable parameters (peak, ascending SD, descending SD, and logit right-floor) and their TV deviates are documented in one place; `sel_inf(1)` is the right-tail floor (analogous to SS3 P6 / `end_logit`), not a fixed-by-map placeholder.

## Deprecations

* The soft-deprecated `srr_indices` / `M1_indices` arguments and the
  legacy env-driven integer codes (`srr_fun %in% c(1, 3, 5)`,
  `M1_model %in% c(4, 5)`) continue to work in 4.4.0 with a one-time
  warning that points users at the linkage table. **Removal has been
  rescheduled from v4.2.0 to v4.5.0** to extend the migration window;
  see the "Scheduled removal" section under 4.1.0 below for the
  unchanged cleanup checklist.

# Rceattle 4.3.1

## Bug fixes

* Fixed a segfault in `MakeADFun` triggered by Beverton-Holt / Ricker
  fits with a recruitment linkage. Section 6.9 (expected recruitment)
  was indexing `R0`, `alpha`, and `Beta` with the function-scope `yr`
  variable, which the preceding forecast loop left equal to `nyrs` --
  one column past the end of the `(nspp, nyrs)` matrices. The year-0
  block now uses an explicit `first_yr = 0` constant. Mean-recruitment
  fits (`srr_fun = 0`) happened to dodge the segfault because
  `calculate_recruitment()` doesn't dereference `alpha`/`Beta` for
  that case.
* Fixed the recruitment-parameter offset formula at the top of section
  5.6: `alpha` and `Beta` now apply the linkage offset on the log
  scale (`exp(rec_pars + offset)`) to match the documented contract
  and the `log_R0` formula on the same line. Previously the offset
  was added on the linear scale, which silently corrupted BH/Ricker
  recruitment whenever a non-zero `log_alpha` or `log_beta` linkage
  offset was active.
* Added a `R0/alpha/Beta` column-count assertion at the entry of
  section 6.9 so future stale-loop / sizing regressions surface as a
  clean R-level error from TMB rather than an opaque segfault.

## Reparameterised intercept handling for the linkage system

Linkage rows whose `design_col == "(Intercept)"` no longer carry the
parameter level themselves. Instead:

* The base parameter (`rec_pars`, `log_M1`, `log_growth_pars`)
  remains estimable and holds the level. Phasing and the per-process
  map machinery operate on the base parameter as they would without
  any linkages.
* The `(Intercept)` row's `beta_linkage` slot is fixed at `0` and
  mapped NA. It exists in the table for bookkeeping plus as a hook
  for `init` and `priors`.
* `init = list(`(Intercept)` = X)` on the spec **flows to the base
  parameter** -- it sets the base parameter's starting value rather
  than the linkage row.
* Priors attached to an `(Intercept)` row are **re-targeted to the
  base parameter** in the slot 19 contribution; the prior density
  evaluates against `rec_pars(sp, 0)` / `log_M1(sp, sx, ag)` /
  `log_growth_pars(sp, sx, par)` rather than the (zero) linkage
  value.

For slope-only formulas (`~ 0 + temp`) the behaviour is unchanged:
the base parameter is still mapped NA at its `build_params()`
default, and the linkage row carries the year-by-year offset.

### Recruitment offset semantics

Year 0 of the recruitment block no longer bakes the year-0 covariate
contribution into `R0`. `R0` is computed from `rec_pars(sp, 0)`
alone, and the linkage offset multiplies against `R0` each year
(including year 0):

```
R(yr) = R0 * exp(rec_dev(yr) + linkage_offset(yr))
```

This makes the legacy `srr_fun = 1 / 3 / 5` quirk (which double-counted
`Temp[0]`) obsolete. Users migrating from the legacy paths should now
get clean log-linear behaviour without surprise offsets.

### Schema additions

* New `init_supplied` (logical) column on `Rceattle_linkage_table`
  tracks whether the user explicitly supplied an `init` for that row.
  Used by `build_params()` to decide whether to push a base-parameter
  init.
* New `linkage_is_intercept` IVECTOR in the TMB encoding (set from
  `design_col == "(Intercept)"`) used by the slot-19 prior dispatch
  to evaluate intercept priors against the base parameter.

### Tests

`tests-Dynamics/test-linkage-auto-map.R`,
`tests-Dynamics/test-recruitment-linkage.R`, and
`tests-Dynamics/test-growth-linkage-species.R` were updated to assert
the new contract (base parameters estimable, intercept rows fixed at
0, slope-only offsets in the year-by-year tensor).

# Rceattle 4.3.0

## Tidy long-format extraction: `as.data.frame.Rceattle()`

A new S3 method on `as.data.frame()` flattens derived population
quantities into a long data.frame with columns
`year, species, sex, age, quantity, value, lwr, upr` so that custom
plotting and post-processing don't have to walk the nested
`quantities` list or rely on the dimnames decisions in
`rename_output()`. Two shapes are supported and combined into one
frame:

* **Species-by-year** (default `which`): `biomass`, `ssb`, `R`,
  `biomass_depletion`, `ssb_depletion`, `F_spp`. Other species/year
  series (`B0`, `SB0`, `DynamicB0`, `DynamicSB0`, `DynamicSBF`,
  `exploitable_biomass`, `proj_F`, `fT`) are available by name; pass
  `which = "all"` to get every known quantity present on the fit.

* **Species-by-sex-by-age-by-year**: `N_at_age`, `biomass_at_age`,
  `Z_at_age`, `M_at_age`, `M1_at_age`, `M2_at_age`, `F_at_age`,
  `consumption_at_age`, `B_eaten_as_prey`, `NByage0`, `NByageF`. The
  `age` column is biological age (offset by `data_list$minage`), and
  cells padded out to `max(nsex)` / `max(nages)` for species with
  fewer sexes or ages are dropped rather than returned as `NA`.

`lwr` / `upr` are populated from the TMB `sdreport` for any quantity
that was `ADREPORT`'d (currently `biomass`, `ssb`, `R`); other
quantities and fits produced with `getsd = FALSE` get `NA` for the
band. The `ci_level` argument (default `0.95`) controls width.

## Optional data fields, continued (Phases B, C, D)

Continuing the Phase A work from 4.2.0, three more classes of inputs
that were previously required as non-NULL can now be omitted, with
`data_check()` enforcing them only when the model actually needs them:

* **Phase B: bioenergetics scalars.** `Ceq`, `Cindex`, `Pvalue`,
  `fday`, `CA`, `CB`, `Qc`, `Tco`, `Tcm`, `Tcl`, `CK1`, `CK4` may be
  `NULL` in single-species mode. `switch_check()` fills them with
  safe sentinels so TMB's length-`nspp` `DATA_VECTOR` requirements
  are satisfied. When `msmMode > 0` the scalars are required;
  `data_check()` reports which ones are missing or wrong-length in
  a single grouped error.

* **Phase C: `env_data`.** May be `NULL`. `clean_data()` defaults it
  to a Year-only `data.frame(Year = styr:projyr)` with zero indices.
  Existing checks still error when a feature actually needs an
  index (env-dependent catchability, temperature-dependent
  consumption, env linkages, `srr_indices`, `M1_indices`).

* **Phase D: `emp_sel`.** New requirement check: when any fleet has
  `Selectivity = "Fixed"`, `emp_sel` must be supplied. Other fleets
  do not need it.

* **Tests.** A new `tests-Data-processing/test-optional-fields.R`
  file exercises 25 NULL / requirement scenarios across the four
  phases.

# Rceattle 4.2.0

## Optional data fields & data_check cleanup

Several fields in `data_list` that were previously required as non-NULL
data.frames are now truly optional. Users who do not need composition
data, conditional age-at-length, empirical selectivity, fixed
numbers-at-age, ration data, or diet data can omit them entirely;
`clean_data()` default-fills the missing fields with empty data.frames
that carry the metadata columns the downstream code expects, and
`data_check()` enforces the field only under the conditions where the
model actually needs it.

* **Phase A (this release): `comp_data`, `caal_data`, `emp_sel`,
  `NByageFixed`, `ration_data`, `diet_data` may be `NULL`.** Conditional
  requirements are still enforced (`caal_data` when
  `any(growth_model > 0)`; `NByageFixed` when `any(estDynamics > 0)`;
  `diet_data` when `msmMode > 0`).

* **`data_check()` reorganisation.** The validation function has been
  reorganised into eight topical sections (top-level scalars; per-species
  dimensions; biology; fleet control; observation tables; diet &
  predation; environmental data; switches), with shared `has_data()` /
  `fc_num()` helpers and consolidated duplicate guards. New checks were
  added for year-scalar ordering, lognormal SDs, sample sizes,
  probability ranges, observation values, fleet referential integrity,
  selectivity bin bounds, predation cross-checks, duplicate observations,
  and probability-matrix row sums. Several pre-existing dead branches and
  matrix-`$` access bugs were fixed at the same time.

* **`transpose_fleet_control()` removed.** The deprecated long-format
  fleet_control transposer has been removed from `clean_data()`,
  `read_data()`, and the package namespace.

# Rceattle 4.1.0

## Environmental linkages: a unified, formula-driven API

A new long-format **linkage table** lets users express how process
parameters depend on environmental covariates and on stratifying
factors (species, sex, age) through a single formula-driven helper,
`linkage_spec()`. Each row of the table corresponds to exactly one
estimated coefficient. `fit_mod()` pools every spec into a shared
design matrix `X` and a per-row parameter vector
`beta_linkage`; the TMB template iterates the table once and
accumulates per-process offsets on the linear predictor of the
underlying parameter.

* **New constructor: `linkage_spec()`.** Captures
  `(formula, by, species, link, init, bounds, priors, est_phase)`
  for one process parameter. Anything `model.matrix()` understands
  works: `~ 1`, `~ temp + PDO`, `~ poly(temp, 4)`, `~ I(temp^2)`,
  `~ splines::ns(temp, df = 4)`, `~ temp * PDO`, etc.

* **Per-species formulas.** Register multiple specs against the
  same parameter via `linkages = list(log_K = list(spec_a, spec_b))`
  with each spec's optional `species = ...` argument scoping it
  to a subset of stocks. The pooler unions the design columns
  across specs so there's no duplication when species share
  covariates.

* **Priors.** First-class via `prior_normal()`,
  `prior_lognormal()`, `prior_gamma()`, `prior_beta()`. The same
  constructors are available unprefixed (`normal()` /
  `lognormal()` / ...) **only inside** the `priors = ...` argument
  via a private NSE data mask, so user code stays close to
  mathematical notation without masking
  `base::gamma()` / `base::beta()` at the package level. Priors
  can be a single value applied to every species, or a named list
  keyed by species id (and shortly, by `(species, sex)`).

* **Bounds.** Per-row `lower` / `upper` flow into
  `build_bounds()$lower$beta_linkage` /
  `build_bounds()$upper$beta_linkage`.

* **Growth** (von Bertalanffy / Richards) is the first process
  fully wired to the new pipeline. `build_growth()` gains a
  `linkages` argument and a string-named `fun`
  (`"empirical"` / `"vonBertalanffy"` / `"Richards"`); integer
  codes still work (`fun = 1` is shorthand for
  `fun = "vonBertalanffy"`) so existing scripts don't need to be
  rewritten apart from substituting `fun =` for `growth_model =`.

  ```r
  build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = linkage_spec(
        formula = ~ temp,
        by      = ~ species + sex,
        priors  = list(temp = normal(0, 1))
      )
    )
  )
  ```

* **TMB plumbing.** New `src/TMB/linkage.hpp` accumulator;
  `ceattle_v01_11.cpp` reads parallel `DATA_IVECTOR(linkage_*)`
  inputs plus a `DATA_MATRIX(linkage_X)` and writes a
  `growth_linkage_offset` tensor that is added (additively, on the
  log scale) to `growth_parameters`. Per-row prior densities
  contribute to slot 19 of the joint NLL ("Linkage-table priors"
  in `fit$quantities$jnll_comp`).

* **Documentation.** New vignette
  `vignette("environmental-linkages", package = "Rceattle")` walks
  through the API, prior families, species-keyed priors,
  per-species formulas, basis-expansion formulas, and the
  underlying pipeline.

* **Natural mortality** is the second process wired to the
  pipeline. `build_M1()` gains a `linkages` argument keyed by
  `log_M1`; the offset is added on the log scale to `log_M1` inside
  the `M1_at_age` compute. A row's `age_bin == NA` broadcasts the
  offset across ages; specific values pin it to that age slice.
  `build_M1()` also gains string-form acceptance for `M1_model`
  and `M1_re` (parity with `build_growth(fun)`):

  ```r
  build_M1(M1_model = "sex_age_invariant",   # or 1
           M1_re    = "ar1_age",             # or 4
           linkages = list(
             log_M1 = linkage_spec(formula = ~ temp, by = ~ species)
           ))
  ```

  Growth and M can be linked in the same fit; their rows share the
  same global linkage table and the same `beta_linkage`
  parameter vector.

* **Per-(species, sex) priors.** In addition to scalar and
  species-keyed priors, each `priors[[col]]` value may be a
  two-level nested list keyed first by species id then by sex id:

  ```r
  priors = list(
    temp = list(
      `1` = list(`1` = normal(0, 0.1),
                 `2` = normal(0, 0.2)),    # sp 1 by sex
      `2` = normal(0, 0.5)                  # sp 2, both sexes
    )
  )
  ```

  Missing keys at either level resolve to "no prior" for that
  stratum. The validator checks the structure recursively and
  emits actionable error messages keyed by
  `priors$<col>$<species>[$<sex>]` paths.

* **Default `by = ~ species`.** `linkage_spec()` now defaults the
  `by` argument to `~ species`, so each linkage produces one
  coefficient per species without the user having to spell it out.
  Pass `~ species + sex` for per-(species, sex) coefficients, or
  `by = NULL` to share a single coefficient across every
  species/sex (the prior default). This matches the typical
  multispecies assessment use case where each stock has its own
  environmental sensitivity.

* **Recruitment** is the third process wired to the pipeline.
  `build_srr()` gains a `linkages` argument keyed by `log_R0`,
  `log_alpha`, or `log_beta`; the offset is added on the log scale
  to the corresponding parameter at every recruitment compute
  call site (hindcast, BRPs, projections, expected R). `log_R0`
  is meaningful for any `srr_fun`; `log_alpha` and `log_beta`
  only do work for SRRs that consume them
  (`srr_fun in c(2, 3, 4, 5)` -- Beverton-Holt and Ricker).

  ```r
  build_srr(
    srr_fun  = 2,
    linkages = list(
      log_alpha = linkage_spec(formula = ~ temp,
                               priors  = list(temp = normal(0, 0.5)))
    )
  )
  ```

  Growth, M, and recruitment can be linked in the same fit; their
  rows share the same global linkage table and the same
  `beta_linkage` parameter vector. End-to-end tests in
  `tests/testthat/tests-Dynamics/test-linkage-auto-map.R` verify that
  the base parameter (e.g., `log_growth_pars`, `log_M1`, `rec_pars`)
  is automatically mapped out (set to `NA`) when a linkage is active
  for that parameter, allowing the linkage intercept to define the base value.
  
  `tests/testthat/tests-Dynamics/test-recruitment-linkage.R`
  cover the analytical relations
  `R = R0 * exp(beta * temp[yr])` (mean R) and the
  `growth + M + recruitment` composition.

* **Soft deprecation in `build_M1()`.** The legacy column-index
  argument `M1_indices` and the env-driven structural integer
  codes `M1_model %in% c(4, 5)` are subsumed by the new
  `linkages = list(log_M1 = ...)` argument. Both still work for one
  release cycle, but emit a one-time warning that points users at
  the equivalent linkage-table call. The string aliases
  `"env_sex_invariant"` / `"env_sex_specific"` (added briefly on
  the dev branch) are removed; they were never released. The cpp's
  `M1_beta` / `M1_mult` code path stays in place for now -- both
  paths add additively to `log_M1` on the log scale -- but do not
  use both for the same coefficient or you will double-count.

* **Roadmap.** Recruitment is next on the same pipeline; then
  random-effects pooling on `re_group` for hierarchical
  shrinkage. The legacy `M1_indices` / `M1_model = 4|5` paths
  retire when recruitment migrates.

## Scheduled removal (v4.5.0)

> **Schedule update (v4.4.0):** The removal originally targeted for
> v4.2.0 has been pushed to **v4.5.0** to give downstream users a
> longer migration window for the natural-scale linkage API rolled
> out in v4.4.0. The soft-deprecation warnings continue to point
> users at the equivalent linkage-table call. The cleanup checklist
> below is unchanged.

The soft-deprecated API surfaces below remain functional and emit
one-time warnings pointing users at the linkage table. They will be
**removed entirely in 4.5.0**. To migrate, replace:

| Legacy                                | New                                                      |
|---------------------------------------|----------------------------------------------------------|
| `build_srr(srr_indices = ...)`        | `build_srr(linkages = list(log_R0 = linkage_spec(...)))` |
| `build_srr(srr_fun = 1)`              | `build_srr(srr_fun = 0)` + linkage on `log_R0`           |
| `build_srr(srr_fun = 3)`              | `build_srr(srr_fun = 2)` + linkage on `log_alpha`        |
| `build_srr(srr_fun = 5)`              | `build_srr(srr_fun = 4)` + linkage on `log_alpha`        |
| `build_M1(M1_indices = ...)`          | `build_M1(linkages = list(log_M1 = linkage_spec(...)))`  |
| `build_M1(M1_model = 4)`              | `build_M1(M1_model = 1)` + linkage on `log_M1`           |
| `build_M1(M1_model = 5)`              | `build_M1(M1_model = 2)` + linkage on `log_M1`           |

**Cpp cleanup checklist** (search for `LEGACY` in
`src/TMB/ceattle_v01_11.cpp`):

* Drop `PARAMETER_MATRIX(beta_rec_pars)`, `PARAMETER_ARRAY(M1_beta)`,
  and the scratch vectors `srr_mult`, `beta_rec_tmp`, `env_rec_tmp`,
  `M1_mult`, `beta_M1_tmp`, `env_M1_tmp`.
* Delete the five `srr_env_mult` blocks (hindcast, BRPs, dynamic
  BRPs, projection, R_hat) and the `M1_mult.sum()` term inside the
  M1_at_age compute.
* Pass `Type(0.0)` for `srr_env_mult` at each
  `calculate_recruitment()` call site (or drop the parameter from
  the function signature in `recruitment.hpp` if no caller still
  needs it).

**R-side cleanup**:

* Remove `srr_indices` and `M1_indices` arguments from
  `build_srr()` and `build_M1()`.
* Reject `srr_fun %in% c(1, 3, 5)` and
  `M1_model %in% c(4, 5)` as unknown integer values in
  `.coerce_srr_fun()` and `.coerce_M1_arg()` (drop the
  `.SRR_DEPRECATED_FUNS` / `.M1_DEPRECATED_MODELS` constants).
* Remove the `suppressWarnings()` wrappers around internal
  `build_srr()` / `build_M1()` re-callers in `sim_mod()`,
  `retrospective()`, `jitter()`, `run_mse()`, `project_no_F()`.

# Rceattle 4.0.3

## API

* New `fit_control()` constructor bundles the optimizer / sdreport /
  phasing knobs that previously cluttered `fit_mod()`'s signature
  (`phase`, `getsd`, `bias.correct`, `use_gradient`, `rel_tol`,
  `loopnum`, `newtonsteps`, `getJointPrecision`, `getReportCovariance`,
  `verbose`, `TMBfilename`, `nlminb_control`). Pass the result via the
  new `fit_control` argument:

  ```r
  fit <- fit_mod(
    data_list   = BS2017SS,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, loopnum = 1)
  )
  ```

  `fit_mod()`'s visible argument list shrinks from ~33 to ~22 args, so
  calls now read as model spec rather than a pile of optimizer flags.

* `fit_mod()` emits a deprecation warning if any of the legacy
  control args are passed directly and forwards them into
  `fit_control` for the duration of the deprecation window. Truly
  unknown arguments still error with `Unused arguments to fit_mod():
  ...` (no silent drops).

* Internal callers (`run_mse()`, `retrospective()`, `jitter()`,
  `sim_mod()`, `project_no_F()`) now wrap their control args in
  `fit_control(...)` rather than passing them positionally.

## New methods

* S3 methods on the `"Rceattle"` class so a fit behaves like an R
  model object: `plot()`, `coef()`, `vcov()`, `logLik()`,
  `residuals()`. With `df` set on `logLik`, [stats::AIC()] also
  works without a dedicated method. `nobs()` is intentionally
  *not* defined: counting "observations" in a stock-assessment
  likelihood (composition cells, indices, catches, priors) is not
  well-defined, so [stats::BIC()] does not work — use AIC or
  domain-specific information criteria.

* `plot.Rceattle()` is a thin dispatcher: `plot(fit, what = "biomass")`
  / `"ssb"` / `"recruitment"` / `"depletion"` / `"index"` / `"catch"` /
  `"selectivity"` / `"mortality"` / `"data"`. `...` is forwarded to
  the underlying `plot_*()` function.

* `residuals.Rceattle(type = ...)` returns a long-format data frame
  with rows from one or more of the four fitted data sources:
  `"index"` and `"catch"` (log-scale by default; switch with
  `scale = "natural"`), `"comp"` (Pearson on fitted proportions, with
  the `Age0_Length1` flag preserved), and `"caal"` (Pearson on
  fitted proportions, with both the conditioning `Length` and the
  age `Bin`). `type = "all"` returns all four stacked.

## Documentation

* README now has a self-contained *Getting started* block that fits a
  bundled model and exercises every new S3 method, so first-time
  users on the pkgdown site / CRAN no longer have to bounce to the
  GitHub wiki to see a working example.

* Vignette 8 ("Model parameterizations") is being expanded to fill in
  coverage gaps surfaced during the 4.0.2 release audit:
  * M1 random-effects modes (`M1_re = 0..6`).
  * Full stock-recruit section (Mean / BevertonHolt / Ricker, env-linked
    forms, Beta prior on steepness via `srr_est_mode = 3`).
  * Composition likelihoods (Multinomial, Dirichlet-multinomial, CAAL).
  * Catchability equations for `Catchability = 4` (Power),
    `5` (Environmental), and `6` (AR1).
  * Selectivity equations for `Selectivity = 6` (2DAR1) and
    `7` (3DAR1).
  * Internal growth model (`growth_model = 1` von Bertalanffy /
    `2` Richards).
  * Initial-age-structure modes (`initMode = 0..4`).
  * Harvest control rules (HCR = 0..7) — possibly cross-linking
    vignette 0 §9 rather than duplicating.

# Rceattle 4.0.2

## Bug fixes

* `data_check()` now stops with a clear message if a user requests
  `msmMode = 3:9` (Kinzey & Punt 2009 functional responses --
  Holling I/II/III, predator interference, predator preemption,
  Hassell-Varley, Ecosim). Those code paths are unvalidated against
  the current parameter set and should not be used for advice.
* Plotting functions (`plot_timeseries`, `plot_selectivity`,
  `plot_mortality`, `plot_maturity`, `plot_b_eaten`, `plot_b_eaten_prop`,
  `plot_m_at_age`, `plot_m2_at_age_prop`, `plot_f`, `plot_ration`,
  `plot_index`, `plot_catch`, `plot_logindex`, `plot_indexresidual`,
  `plot_comp`, `plot_selectivity_vs_maturity`, `plot_stock_recruit`)
  now restore graphics `par()` on exit instead of mutating the
  user's device permanently.
* Replaced `T` / `F` shortcuts with `TRUE` / `FALSE` throughout the
  package source (~60 occurrences).
* Replaced `.data$Bin` / `.data$Length` inside `tidyr::pivot_wider`
  arguments with quoted strings (tidyselect 1.2.0 deprecation).
* `examples/Georges_bank_example.R` now calls `plot_mortality()`
  instead of the long-removed `plot_mort()`.
* `_pkgdown.yml` "Get started" / overview navbar links now point at
  the actual generated `articles/Rceattle-overview.html` (was
  `Rceattle-overview-4_17_2025.html`, which 404'd).
* `README.md` example links now reference the correct
  underscore-separated filenames (`Fit_2018_GOA_multi-species_model.R`
  etc.) -- the previous space-encoded URLs 404'd on GitHub.
* Removed duplicate `\seealso{}` block in `?Rceattle-package` by
  consolidating `R/Rceattle.R` into `R/Rceattle-package.R`. Both
  files declared `_PACKAGE`, so roxygen2 emitted the auto-generated
  links twice.
* Added `graphics::box` to package imports (cleared the lone
  `R CMD check` NOTE for `plot_data`).
* TMB: terminal-age geometric series now includes `Finit` in the
  denominator for fished initial-equilibrium modes
  (`initMode = 3` or `4`), correcting a bias in the plus-group
  N-at-age when `Finit > 0`.

## Documentation

* Added Wassermann et al. (2025) cannibalism / Pacific hake
  reference to `inst/CITATION` and `?Rceattle-package`.
* `initMode` accepts integer codes or string aliases.
* `HCR` accepts integer codes or string aliases.

## Tests

* `tests/testthat/test-subdir-folders.R` no longer calls
  `testthat::test_dir()` for each subdirectory. Nested `test_dir()`
  inside `test_check()` triggered an "evaluation nested too deeply:
  infinite recursion" abort inside `rlang`'s trace deparser whenever
  any test failed, masking the real failure. Subdirectory test files
  are now discovered with `list.files()` and pulled in via
  `source()` so they register against the outer reporter directly.
* `tests/testthat.R` now wraps `library()` calls in
  `suppressPackageStartupMessages()` so transient build-version
  notices (e.g. "package 'dplyr' was built under R version 4.5.2")
  do not get captured as test warnings whose backtraces then crash
  rlang's expr_deparse at end-of-run.
* `tests/testthat/test-Dynamics/test-initial-dynamics` evaluates the different starting conditions.

## Parallelism

* `run_mse()` and `check_mse()` now use `parallel::parLapply` on a
  PSOCK cluster instead of `foreach::foreach(...) %dopar%`.
  - The `%dopar%` path triggered `rlang::expr_deparse` infinite
    recursion under nested `test_that` backtraces because each
    `foreach` invocation captures call frames that recurse during
    error formatting. PSOCK workers are clean R processes with no
    captured promise chains, so the issue does not occur.
  - PSOCK clusters work identically on Windows and macOS/Linux.
  - `run_mse()` gains a `cores` argument (default `NULL` picks
    `parallel::detectCores() - 6`); both functions cap at 2 cores
    when `_R_CHECK_LIMIT_CORES_` is set so they comply with CRAN's
    R CMD check limit.
  - `foreach` and `doParallel` removed from `Imports:`.

## Installation / dependencies

* `TMBhelper` moved from `Imports:` to `Suggests:`. Rceattle now uses an
  internal `.fit_tmb()` wrapper that delegates to `TMBhelper::fit_tmb()`
  when the (optional, GitHub-only) package is installed and otherwise
  falls back to a `stats::nlminb()` + `TMB::sdreport()` path. This
  removes the largest install-friction barrier for new users.
* `dplyr` moved from `Depends:` to `Imports:`. The package no longer
  attaches `dplyr` to the user's search path on load (so it no longer
  masks `stats::filter` / `stats::lag`).
* `kableExtra` dropped from `Suggests:`; vignettes now use
  `knitr::kable()` for table rendering.
* `quarto` dropped as a vignette engine; all vignettes are now `.Rmd`.

## API

* `run_mse()` no longer has `om = ms_run`, `em = ss_run` defaults.
  Both arguments are now required and validated as objects of class
  `"Rceattle"` before the MSE loop runs. Calling `run_mse()` with no
  arguments previously produced a confusing "object 'ms_run' not found"
  error; it now stops with a clear message.

## New methods

* `print.Rceattle()` and `summary.Rceattle()`. Auto-printing a fit
  inside knitr / RStudio / R Markdown previously dumped tens of MB of
  nested data and could trigger deep recursion errors during vignette
  rendering.

## Build / packaging

* Tightened `.Rbuildignore` (excludes `examples/`, `R/dev/`,
  `src/TMB/Dev/`, `.Rhistory`, `.claude/`, `.DS_Store`, build tarballs,
  `.Rcheck` directories).
* Tightened `.gitignore` to catch all `*.o` / `*.so` / `*.dll` and
  `*.Rcheck/` directories.
* Suppressed a benign clang `-Wfixed-enum-extension` warning by
  scoping the diagnostic pragma around `#include <TMB.hpp>` rather
  than via global Makevars flags.

# Rceattle 4.0.1

The 4.0.1 development cycle reorganized several `data_list` columns
and `fit_mod` / `build_*` arguments. Models or data files saved
against earlier 4.x revisions may need updating; see the renames
below. Compiled from
`inst/Running_list_of_updates.qmd` plus the `dev` branch commit log.

## Data renames

* `Pyrs` -> `ration_data` (the old name is still accepted on read,
  but is silently renamed).
* `UobsWtAge` -> `stom_prop_data`.
* `fsh_biom` -> `catch_data`.
* `srv_biom` -> `index_data`.
* `Nselages` -> `N_sel_bins` (in `fleet_control`).
* `Sel_norm_bin1` / `Sel_norm_bin2` <- `Age_max_selected` /
  `Age_max_selected_upper` (selectivity normalization bins).
* `Age_first_selected` -> `Bin_first_selected` (in `fleet_control`).
* `sel` -> `sel_at_age` (model report).
* `fleet_control` now carries a `Month` column (month of observation
  for indices / fisheries).

## API renames

* `build_M1`: `M1_prior_mean` -> `M_prior`,
  `M1_prior_sd` -> `M_prior_sd`.
* `build_srr`: `srr_prior_mean` -> `srr_prior`;
  `R_hat_endyr` replaced by `srr_hat_styr` / `srr_hat_endyr`.
* `fit_mod`: `suit_meanyr` replaced by `suit_styr` / `suit_endyr`.
* `initMode` semantics revised: 0 = free-parameter N-at-age,
  1 = unfished equilibrium with no devs, 2 (default) = unfished
  equilibrium with initial devs, 3 = fished equilibrium with initial
  devs. Type 4 ("non-equilibrium scaled") added later.

## New features -- composition and diet likelihoods

* Dirichlet-multinomial composition likelihood. Selected per fleet via
  `fleet_control$Comp_loglike = 1` (or `"DirichletMultinomial"`).
* Conditional age-at-length (CAAL) data path, with `CAAL_loglike` /
  `CAAL_weights` controls in `fleet_control`. CAAL data also flow
  through `sim_mod()` for simulation testing.
* `Diet_loglike` switch on the bioenergetics control sheet selects
  between multinomial (0) and Dirichlet-multinomial (1) for diet
  composition.
* Other-food diet proportion estimates added to the model report.
* Weighted-mean diet data path (annual proportion of prey-at-age in
  predator-at-age averaged across years).

## New features -- selectivity, catchability, growth

* Hake non-parametric selectivity (`Selectivity = "Hake"` or `5`),
  after Taylor et al.
* `2DAR1` (`= 6`) and `3DAR1` (`= 7`) selectivity
  parameterizations, after Cheng et al. (2024).
* `Catchability = 6` ("AR1"): annual AR1 catchability deviates fit
  to an environmental index, after Rogers et al. (2024) for the GOA
  pollock model. Environmental q-link (`Catchability = 5`) also
  exposed.
* Internal growth model. See `build_growth()` and the `growthFun`
  argument to `fit_mod()`. `alpha_wt_len` / `beta_wt_len` added to
  the data control sheet. Length-based suitability (`suitMode = 1` /
  `2` / `3` / `4` / `5` / `6`) wired through to use the estimated
  growth model. Comparison with WHAM growth implemented under
  `tests/comparison/`.
* Predator-specific suitability mode (different `suitMode` per
  predator).
* Suitability calculation now uses configurable year ranges
  (`suit_styr` / `suit_endyr`) instead of "mean year".

## New features -- recruitment and reference points

* Beta-distributed prior on Beverton-Holt steepness, available via
  `srr_est_mode = 3`.
* M1 random effects with optional environmental linkage; `M_prior` /
  `M_prior_sd` priors carried through `build_M1()`.
* `remove_F()` function returns a fitted model with F set to 0 --
  used internally for dynamic reference point calculation.
* `DynamicHCR = TRUE` in `build_hcr()` to switch from static to
  dynamic SB0 reference points.
* CMSY harvest control rule (`HCR = 1`): maximize joint catch across
  species, optionally constrained to keep depletion above `Plimit`.
* PFMC Category 1 40-10 ABC HCR (`HCR = 6`) using `Pstar` /
  `Sigma` uncertainty buffer.
* SESSF Tier 1 HCR (`HCR = 7`).
* Iterative multi-species HCRs: `HCRorder` controls the order in
  which species F is solved (e.g. predators before prey) inside
  `build_hcr()`.

## New features -- MSE and projection

* `run_mse()` now writes per-simulation `.rds` files when `dir` is
  specified, for streaming-friendly long runs. `load_mse()` reads
  those back.
* `check_mse()` validates which OM/EM simulations converged.
* `mse_summary()` produces a per-fleet performance-metric table
  (mean catch, IAV, P(closed), MSE on SSB, P(F > Flimit),
  P(SSB < SSBlimit), terminal depletion, ...).
* MSE function now supports `cap` (catch cap), `catch_mult` (catch
  multiplier), `rec_trend` (linear projected recruitment trend),
  `fut_sample` (future sampling effort), per-fleet
  `assessment_period` / `sampling_period`, `regenerate_past` (refit
  EM to OM-simulated past data), and `timeout`/`try`-error handling
  per simulation.
* `Recruitment_and_fixed_F_projections.R` and `Simulation_testing.R`
  examples added.

## New features -- diagnostics and tooling

* `jitter()` function to perturb starting values and re-fit, for
  global-vs-local-minimum diagnostics.
* `retrospective()` peels with optional `nyrs_forecast`.
* `model_average()` for averaging derived quantities across multiple
  fitted models, with optional bootstrap uncertainty.
* `compare_sim()` and `sim_mod()` for parametric simulation testing.
* `McAllister-Ianelli-reweighting.R` example for composition
  reweighting.
* TMB log-likelihood pieces (unweighted) added to the report for
  composition diagnostics.
* `Selectivity = "Fixed"` (`= 0`) for empirically supplied selectivity
  blocks via the `emp_sel` data sheet.
* `TMBfilename` argument to `fit_mod()` to point at an alternate
  `.cpp` during development.

## Behavior changes

* Removed accumulation-age switches in `fleet_control`. Selectivity
  normalization is now controlled via `Age_max_selected` (i.e.
  `Sel_norm_bin1`) on a per-fleet basis instead of always
  normalizing by the maximum-selectivity age.
* `NA` values inside the valid age/length range of composition data
  are now coerced to 0 with a warning (previously silently dropped
  or errored).
* Selectivity dimensioning switched from age- to bin-indexed for the
  non-parametric and 2D/3D AR1 forms (driven by
  `Selectivity_dimension` and `N_sel_bins`).
* Age-error and age-transition matrices are now dimension-checked
  against `nages` at `data_check()` time.

# Rceattle 4.0.0

* CEATTLE TMB version 4.0.0. See Adams et al. (2022),
  *Fisheries Research*, 251, 106303.
