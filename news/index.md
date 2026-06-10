# Changelog

## Rceattle 4.5.0

### New features

- Added a general post-fit convergence diagnostics framework via the new
  exported
  [`convergence_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/convergence_diagnostics.md)
  function. It runs a battery of checks covering the optimizer gradient,
  Hessian positive-definiteness and conditioning, parameters on bounds,
  phasing, and parameter estimability, and returns a structured
  `"Rceattle_convergence"` object whose `status` reflects the worst
  severity found (`"OK"`, `"NOTE"`, `"WARN"`, or `"FAIL"`).
- [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  now runs the convergence battery automatically and attaches the result
  as `fit$convergence`;
  [`convergence_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/convergence_diagnostics.md)
  can also be called directly to re-run it on any fit.
- Added a model-diagnostics vignette section and accompanying unit tests
  for the new framework.

### Bug fixes

- Fixed a parameters-on-bounds misalignment in the convergence
  diagnostics.
- The joint negative log-likelihood (jnll) is now taken from the
  optimizer objective for consistency.

### Internal / R CMD check

- Moved the suppression of the spurious GCC/Eigen `-Warray-bounds` false
  positive from a `-Wno-array-bounds` compiler flag to a source-level
  `#pragma GCC diagnostic ignored "-Warray-bounds"`, so
  `R CMD check --as-cran` no longer warns about non-portable,
  warning-suppressing compilation flags. Real diagnostics still surface.
- Additional C++ cleanup to remove compiler warnings.

## Rceattle 4.4.2

### Code organisation (no change to fitted results)

The pre-fit pipeline files in `R/` were reorganised so they are easier
to navigate. None of these changes alter model output.

- **File prefixes now follow execution order.**
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  runs its stages as
  [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
  -\>
  [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  -\>
  [`build_params()`](https://grantdadams.github.io/Rceattle/reference/build_params.md)
  -\>
  [`build_map()`](https://grantdadams.github.io/Rceattle/reference/build_map.md)
  -\>
  [`build_bounds()`](https://grantdadams.github.io/Rceattle/reference/build_bounds.md)
  -\>
  [`rearrange_data()`](https://grantdadams.github.io/Rceattle/reference/rearrange_data.md)
  -\> fit -\>
  [`rename_output()`](https://grantdadams.github.io/Rceattle/reference/rename_output.md),
  so the files were renumbered to match (`data_check` is now `1-`,
  `build_params` `2-`, `build_map` `3-`, `build_parameter_bounds` `4-`).
  A pipeline map was added to the top of
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).
- **Switch lifecycle consolidated** into a single `R/0-switches.R`: the
  string\<-\>integer maps (formerly `0-constants.R`) plus
  [`switch_check()`](https://grantdadams.github.io/Rceattle/reference/switch_check.md),
  `revert_switches()`, `validate_switches()`, and `convert_switches()`,
  with a header documenting the order in which they run.
- **HCR helpers co-located.**
  [`build_hcr_map()`](https://grantdadams.github.io/Rceattle/reference/build_hcr_map.md)
  moved into `R/0-build_hcr.R` alongside
  [`build_hcr()`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
  (the separate `2-build_hcr_map.R` was removed).

### Rename / deprecation

- [`rearrange_dat()`](https://grantdadams.github.io/Rceattle/reference/rearrange_data.md)
  is renamed
  **[`rearrange_data()`](https://grantdadams.github.io/Rceattle/reference/rearrange_data.md)**.
  The old name still works as a deprecated alias (emits a one-time
  [`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) warning) and
  will be removed in a future release.

### Export hygiene

- `check_composition_data()`, `check_caal_data()`,
  `calc_mcall_ianelli()`, and `calc_mcall_ianelli_diet()` are no longer
  exported; they are internal helpers called only from within the
  package.

### Internal / R CMD check

- Removed `Rceattle:::` self-references in
  [`build_bounds()`](https://grantdadams.github.io/Rceattle/reference/build_bounds.md)
  (a package should not use `:::` for its own objects).
- [`profile.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/profile.Rceattle.md)
  gained `...` for S3 consistency with the
  [`stats::profile`](https://rdrr.io/r/stats/profile.html) generic.

## Rceattle 4.4.1

### Rename

`profile_param()` was renamed
[`profile()`](https://rdrr.io/r/stats/profile.html) and turned into an
S3 class.

## Rceattle 4.4.0

### New features

- **Double Normal Selectivity (Type 8)**: Added support for the
  four-parameter Double Normal selectivity curve (Peak, Ascending SD,
  Descending SD, and Floor). This includes full support for annual
  deviates on all four parameters and is compatible with the age- and
  length-based selectivity engines.
- **Growth SD Control**: The linkage system now supports `sd_L1` and
  `sd_Linf` parameters. This allows users to specify priors, initial
  values, and bounds for growth SD endpoints using the same
  formula-driven interface used for mean growth parameters.
- **Natural-scale linkage API**: Renamed the linkage parameter keys from
  `log_*` to their natural-scale counterparts (`K` / `L1` / `Linf` /
  `m`; `M1`; `R0` / `alpha` / `beta`). Internal parameters remain on the
  log scale.
- **Natural-scale priors**: Standardized prior evaluation so that
  probability densities are applied to parameters on their natural scale
  by default, unless a lognormal family is explicitly requested.
- **Natural-scale inits**: Standardized init evaluation so initial
  values for `"(Intercept)"` parameters are applied on their natural
  scale.
- **Linkage link functions**: Fully implemented dual-path linkage
  offsets. `link = "log"` (the new default) applies the offset
  multiplicatively on the natural scale (additive on the log scale);
  `link = "identity"` applies it additively on the natural scale.
- **Per-species VB anchor age**:
  [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  gains a `growth_age_L1` argument (scalar or length-`nspp` vector) for
  the age at which mean length equals `L1`. Matches SS3’s
  `Growth_Age_for_L1`. Default `NA` inherits `data_list$growth_age_L1`
  if set, else falls back to `max(0.5, minage[sp])` so `minage = 0`
  models pick up an SS3-consistent half-year anchor while `minage >= 1`
  models stay backwards-compatible.
- **Self-test simulation**: New
  [`self_test()`](https://grantdadams.github.io/Rceattle/reference/self_test.md)
  simulates `nsim` datasets from a fitted model and re-fits the model to
  each simulated dataset, returning the list of refits. Runs in parallel
  by default (PSOCK cluster, capped at 2 cores under `R CMD check`) with
  a per-simulation `seed + i` so results are reproducible under both
  sequential and parallel execution.
- **Likelihood profile**: New `profile_param()` generalises the legacy
  `profile_rsigma()` example helper to any parameter slot in
  `Rceattle$estimated_params`. Supports arbitrary N-D cross-profiles via
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) over a list
  of per-cell value grids (e.g. cross-profiling `log_M1` across sex, or
  `R_log_sd` across multiple species). Natural-scale aliases `"sigmaR"`
  / `"R_sd"`, `"M1"`, and `"R0"` / `"alpha"` / `"beta"` apply the
  implicit [`log()`](https://rdrr.io/r/base/Log.html) transform and (for
  the `rec_pars` aliases) auto-fill the column, so `slots` only needs
  the species index. `slots` defaults to species 1 with a warning. Fits
  run in parallel on the same PSOCK harness as
  [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md)
  /
  [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md)
  /
  [`self_test()`](https://grantdadams.github.io/Rceattle/reference/self_test.md).
- **Parallel
  [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md)
  and
  [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md)**:
  Both diagnostics now run their independent peels / starts on a PSOCK
  cluster (same approach as
  [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)).
  New `cores` argument on each (default `parallel::detectCores() - 6`,
  capped at 2 when `_R_CHECK_LIMIT_CORES_` is set); pass `cores = 1` to
  force sequential execution.
- **Standard errors in
  [`as.data.frame.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/as.data.frame.Rceattle.md)**:
  The tidy long-format frame now carries a `se` column alongside `value`
  / `lwr` / `upr`, populated from the TMB `sdreport` for any
  `ADREPORT`’d quantity. Set to `NA` for non-ADREPORT’d quantities and
  for fits produced with `getsd = FALSE`.

### Bug fixes

- **SRR logic**: Fixed a bug in
  [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  where the `Bmsy_lim` penalty was incorrectly disabled for current
  Ricker implementations due to an index mismatch.
- **Selectivity RW prior scaling**: Corrected the random walk prior
  scaling in the TMB template to ensure consistent $`4 \times`$ SD
  multipliers for both ascending and descending limb slope/SD
  parameters.
- **`last_par` returned wrong vector**: Fixed
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  so the value stored on the returned fit is the optimizer’s last
  parameter vector rather than a stale prior reference, removing the
  need for the surrounding [`try()`](https://rdrr.io/r/base/try.html)
  guards in downstream callers.
- **Growth at `minage = 0`**: Fixed a segfault in `growth.hpp` when
  `minage = 0` by guarding `current_age`, `age_L1`, and the
  cohort-boundary `age_L1_ceil` against the zero-age anchor. Also
  corrected length-at-age at `minage = 0` so the closed-form anchor at
  `L1` is honored on both the within-year and cohort recursion paths.
- **Length-bin midpoint for weight-at-age**: `calculate_weight()` now
  computes each bin’s midpoint as `(lengths[ln] + lengths[ln+1]) / 2`
  (with the plus-group extended by half the final interior width) rather
  than `lengths[ln] + (lengths[1] - lengths[0]) / 2`. The previous
  formula assumed uniform bin widths and silently mis-located the
  midpoint for non-uniform length bins; for uniform bins (the common
  case) the two are algebraically identical, so existing fits are
  unchanged.
- **Plus-group SD-at-age continuity**: `growth.hpp` no longer
  special-cases the oldest age class when computing `length_sd`.
  Previously the plus group used `exp(growth_log_sd(sp, sex, 1))`
  directly (SD at `Linf`); it now uses the same length-based linear
  interpolation between SD(`L1`) and SD(`Linf`) as all other ages above
  `age_L1`. This restores continuity across the plus-group boundary;
  users with archived fits whose plus group did not reach `Linf` should
  expect small numerical drift in `length_sd` at the oldest age.
- **[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  bounds ordering**:
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  now indexes parameter bounds by name rather than positional order when
  assembling `L` / `U` for `nlminb`. Previously, when `map$mapFactor`
  and `bounds$lower` were not in identical order, parameters could be
  paired with another parameter’s bounds, producing silently wrong
  constraints. `start_par` is now also subset by name with
  `drop = FALSE`.
- **[`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  Tier-3 Flimit**: The internal `flimit_tier3_fun()` returned `Flimit`
  (its argument) instead of the depletion-adjusted `tier3_flimit` it had
  just computed, so the Tier-3 (HCR = 5) branch of P(F \> Flimit)
  reduced to the base-Flimit check. Now returns the adjusted vector.
- **[`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  HCR coercion**: `HCR` is now normalized to its integer code before
  downstream comparisons (`HCR == 5`, etc.).
  [`build_hcr()`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
  accepts either an integer or a string alias (e.g. `"NPFMC"`);
  [`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  previously assumed integer form and silently produced wrong status
  flags when fits carried the string form.
- **[`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  OM status at assessment years**: P(F \> Flimit) and P(SSB \< SSBlimit)
  are now reported for the OM evaluated at the same assessment years as
  the EM (previously only the EM’s perceived status was returned), and
  the SSB-limit threshold dispatch is consolidated in one helper so the
  Tier-3 / Category-1 / dynamic-vs-static cases stay aligned across the
  EM and OM paths.
- **[`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
  inactive-fleet handling**: The auto-Off branch no longer nulls out
  `proj_F_prop` and `Catchability` on fleets it flips to `"Off"`. The
  downstream TMB code already ignores those columns for off fleets, and
  clearing them lost information users needed when re-enabling a fleet.
- **Selectivity block indexing**: Renamed the within-loop `biom_yrs`
  index vector to `block_yrs` in the selectivity and catchability block
  branches of `build_map_*()`. The previous name was a leftover from the
  index-data path and shadowed nothing, but read as if it referred to
  biomass-observation years.
- **[`plot_index()`](https://grantdadams.github.io/Rceattle/reference/plot_index.md)
  /
  [`plot_catch()`](https://grantdadams.github.io/Rceattle/reference/plot_catch.md)
  /
  [`plot_logindex()`](https://grantdadams.github.io/Rceattle/reference/plot_logindex.md)
  warnings**: Wrapped the internal
  [`gplots::plotCI()`](https://talgalili.github.io/gplots/reference/plotCI.html)
  calls in [`suppressWarnings()`](https://rdrr.io/r/base/warning.html)
  so plotting a fit no longer prints the recurring
  `"In arrows(...): zero-length arrow is of indeterminate angle and so skipped"`
  noise when CI half-widths are zero.

### Data checks

- **Empirical growth + CAAL**:
  [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  now errors when `growth_model == 0` (empirical weight-at-age) is
  combined with non-empty `caal_data` for a given species. The C++
  growth matrix is not populated from the age-transition matrix in the
  empirical branch, so `pred_CAAL` collapses to ~0 and the multinomial
  NLL becomes uninformative.
- **Selectivity identifiability**: Fleets with estimated `Selectivity`
  and `Fleet_type != "Off"` now require at least one positive-sample
  `comp_data` or `caal_data` row. Either provide composition / CAAL
  data, mark the fleet as `Selectivity = "Fixed"` with `emp_sel`, or set
  `Fleet_type = "Off"`.
- **Auto-Off inactive fleets**:
  [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
  automatically flips `Fleet_type` to `"Off"` for fleets that carry no
  catch or index observations, preventing the optimizer from drifting on
  unconstrained selectivity / catchability blocks.
- **`minage` guard**:
  [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  errors when any species has `minage < 0`.

### Documentation

- Added
  [`vignette("environmental-linkages-and-priors")`](https://grantdadams.github.io/Rceattle/articles/environmental-linkages-and-priors.md)
  (and updated `_pkgdown.yml`) to cover the new linkage intercept
  behavior, link-function semantics, growth SD endpoints, and Double
  Normal selectivity.
- Updated all cross-references in
  [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  /
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
  /
  [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  (deprecation warnings, soft-deprecated arg docs, and the
  [`vignette()`](https://rdrr.io/r/utils/vignette.html) pointers in the
  model-options vignette) from the old `environmental-linkages` slug to
  the renamed `environmental-linkages-and-priors` vignette so the
  soft-deprecation warnings now resolve.
- Expanded the Double Normal (selectivity type 8) doxygen / roxygen so
  the four estimable parameters (peak, ascending SD, descending SD, and
  logit right-floor) and their TV deviates are documented in one place;
  `sel_inf(1)` is the right-tail floor (analogous to SS3 P6 /
  `end_logit`), not a fixed-by-map placeholder.

### Deprecations

- The soft-deprecated `srr_indices` / `M1_indices` arguments and the
  legacy env-driven integer codes (`srr_fun %in% c(1, 3, 5)`,
  `M1_model %in% c(4, 5)`) continue to work in 4.4.0 with a one-time
  warning that points users at the linkage table. **Removal has been
  rescheduled from v4.2.0 to v4.5.0** to extend the migration window;
  see the “Scheduled removal” section under 4.1.0 below for the
  unchanged cleanup checklist.

## Rceattle 4.3.1

### Bug fixes

- Fixed a segfault in `MakeADFun` triggered by Beverton-Holt / Ricker
  fits with a recruitment linkage. Section 6.9 (expected recruitment)
  was indexing `R0`, `alpha`, and `Beta` with the function-scope `yr`
  variable, which the preceding forecast loop left equal to `nyrs` – one
  column past the end of the `(nspp, nyrs)` matrices. The year-0 block
  now uses an explicit `first_yr = 0` constant. Mean-recruitment fits
  (`srr_fun = 0`) happened to dodge the segfault because
  `calculate_recruitment()` doesn’t dereference `alpha`/`Beta` for that
  case.
- Fixed the recruitment-parameter offset formula at the top of section
  5.6: `alpha` and `Beta` now apply the linkage offset on the log scale
  (`exp(rec_pars + offset)`) to match the documented contract and the
  `log_R0` formula on the same line. Previously the offset was added on
  the linear scale, which silently corrupted BH/Ricker recruitment
  whenever a non-zero `log_alpha` or `log_beta` linkage offset was
  active.
- Added a `R0/alpha/Beta` column-count assertion at the entry of section
  6.9 so future stale-loop / sizing regressions surface as a clean
  R-level error from TMB rather than an opaque segfault.

### Reparameterised intercept handling for the linkage system

Linkage rows whose `design_col == "(Intercept)"` no longer carry the
parameter level themselves. Instead:

- The base parameter (`rec_pars`, `log_M1`, `log_growth_pars`) remains
  estimable and holds the level. Phasing and the per-process map
  machinery operate on the base parameter as they would without any
  linkages.
- The `(Intercept)` row’s `beta_linkage` slot is fixed at `0` and mapped
  NA. It exists in the table for bookkeeping plus as a hook for `init`
  and `priors`.
- `init = list(`(Intercept)`= X)` on the spec **flows to the base
  parameter** – it sets the base parameter’s starting value rather than
  the linkage row.
- Priors attached to an `(Intercept)` row are **re-targeted to the base
  parameter** in the slot 19 contribution; the prior density evaluates
  against `rec_pars(sp, 0)` / `log_M1(sp, sx, ag)` /
  `log_growth_pars(sp, sx, par)` rather than the (zero) linkage value.

For slope-only formulas (`~ 0 + temp`) the behaviour is unchanged: the
base parameter is still mapped NA at its
[`build_params()`](https://grantdadams.github.io/Rceattle/reference/build_params.md)
default, and the linkage row carries the year-by-year offset.

#### Recruitment offset semantics

Year 0 of the recruitment block no longer bakes the year-0 covariate
contribution into `R0`. `R0` is computed from `rec_pars(sp, 0)` alone,
and the linkage offset multiplies against `R0` each year (including year
0):

    R(yr) = R0 * exp(rec_dev(yr) + linkage_offset(yr))

This makes the legacy `srr_fun = 1 / 3 / 5` quirk (which double-counted
`Temp[0]`) obsolete. Users migrating from the legacy paths should now
get clean log-linear behaviour without surprise offsets.

#### Schema additions

- New `init_supplied` (logical) column on `Rceattle_linkage_table`
  tracks whether the user explicitly supplied an `init` for that row.
  Used by
  [`build_params()`](https://grantdadams.github.io/Rceattle/reference/build_params.md)
  to decide whether to push a base-parameter init.
- New `linkage_is_intercept` IVECTOR in the TMB encoding (set from
  `design_col == "(Intercept)"`) used by the slot-19 prior dispatch to
  evaluate intercept priors against the base parameter.

#### Tests

`tests-Dynamics/test-linkage-auto-map.R`,
`tests-Dynamics/test-recruitment-linkage.R`, and
`tests-Dynamics/test-growth-linkage-species.R` were updated to assert
the new contract (base parameters estimable, intercept rows fixed at 0,
slope-only offsets in the year-by-year tensor).

## Rceattle 4.3.0

### Tidy long-format extraction: `as.data.frame.Rceattle()`

A new S3 method on
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) flattens
derived population quantities into a long data.frame with columns
`year, species, sex, age, quantity, value, lwr, upr` so that custom
plotting and post-processing don’t have to walk the nested `quantities`
list or rely on the dimnames decisions in
[`rename_output()`](https://grantdadams.github.io/Rceattle/reference/rename_output.md).
Two shapes are supported and combined into one frame:

- **Species-by-year** (default `which`): `biomass`, `ssb`, `R`,
  `biomass_depletion`, `ssb_depletion`, `F_spp`. Other species/year
  series (`B0`, `SB0`, `DynamicB0`, `DynamicSB0`, `DynamicSBF`,
  `exploitable_biomass`, `proj_F`, `fT`) are available by name; pass
  `which = "all"` to get every known quantity present on the fit.

- **Species-by-sex-by-age-by-year**: `N_at_age`, `biomass_at_age`,
  `Z_at_age`, `M_at_age`, `M1_at_age`, `M2_at_age`, `F_at_age`,
  `consumption_at_age`, `B_eaten_as_prey`, `NByage0`, `NByageF`. The
  `age` column is biological age (offset by `data_list$minage`), and
  cells padded out to `max(nsex)` / `max(nages)` for species with fewer
  sexes or ages are dropped rather than returned as `NA`.

`lwr` / `upr` are populated from the TMB `sdreport` for any quantity
that was `ADREPORT`’d (currently `biomass`, `ssb`, `R`); other
quantities and fits produced with `getsd = FALSE` get `NA` for the band.
The `ci_level` argument (default `0.95`) controls width.

### Optional data fields, continued (Phases B, C, D)

Continuing the Phase A work from 4.2.0, three more classes of inputs
that were previously required as non-NULL can now be omitted, with
[`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
enforcing them only when the model actually needs them:

- **Phase B: bioenergetics scalars.** `Ceq`, `Cindex`, `Pvalue`, `fday`,
  `CA`, `CB`, `Qc`, `Tco`, `Tcm`, `Tcl`, `CK1`, `CK4` may be `NULL` in
  single-species mode.
  [`switch_check()`](https://grantdadams.github.io/Rceattle/reference/switch_check.md)
  fills them with safe sentinels so TMB’s length-`nspp` `DATA_VECTOR`
  requirements are satisfied. When `msmMode > 0` the scalars are
  required;
  [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  reports which ones are missing or wrong-length in a single grouped
  error.

- **Phase C: `env_data`.** May be `NULL`.
  [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
  defaults it to a Year-only `data.frame(Year = styr:projyr)` with zero
  indices. Existing checks still error when a feature actually needs an
  index (env-dependent catchability, temperature-dependent consumption,
  env linkages, `srr_indices`, `M1_indices`).

- **Phase D: `emp_sel`.** New requirement check: when any fleet has
  `Selectivity = "Fixed"`, `emp_sel` must be supplied. Other fleets do
  not need it.

- **Tests.** A new `tests-Data-processing/test-optional-fields.R` file
  exercises 25 NULL / requirement scenarios across the four phases.

## Rceattle 4.2.0

### Optional data fields & data_check cleanup

Several fields in `data_list` that were previously required as non-NULL
data.frames are now truly optional. Users who do not need composition
data, conditional age-at-length, empirical selectivity, fixed
numbers-at-age, ration data, or diet data can omit them entirely;
[`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
default-fills the missing fields with empty data.frames that carry the
metadata columns the downstream code expects, and
[`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
enforces the field only under the conditions where the model actually
needs it.

- **Phase A (this release): `comp_data`, `caal_data`, `emp_sel`,
  `NByageFixed`, `ration_data`, `diet_data` may be `NULL`.** Conditional
  requirements are still enforced (`caal_data` when
  `any(growth_model > 0)`; `NByageFixed` when `any(estDynamics > 0)`;
  `diet_data` when `msmMode > 0`).

- **[`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  reorganisation.** The validation function has been reorganised into
  eight topical sections (top-level scalars; per-species dimensions;
  biology; fleet control; observation tables; diet & predation;
  environmental data; switches), with shared `has_data()` / `fc_num()`
  helpers and consolidated duplicate guards. New checks were added for
  year-scalar ordering, lognormal SDs, sample sizes, probability ranges,
  observation values, fleet referential integrity, selectivity bin
  bounds, predation cross-checks, duplicate observations, and
  probability-matrix row sums. Several pre-existing dead branches and
  matrix-`$` access bugs were fixed at the same time.

- **`transpose_fleet_control()` removed.** The deprecated long-format
  fleet_control transposer has been removed from
  [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md),
  [`read_data()`](https://grantdadams.github.io/Rceattle/reference/read_data.md),
  and the package namespace.

## Rceattle 4.1.0

### Environmental linkages: a unified, formula-driven API

A new long-format **linkage table** lets users express how process
parameters depend on environmental covariates and on stratifying factors
(species, sex, age) through a single formula-driven helper,
[`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md).
Each row of the table corresponds to exactly one estimated coefficient.
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
pools every spec into a shared design matrix `X` and a per-row parameter
vector `beta_linkage`; the TMB template iterates the table once and
accumulates per-process offsets on the linear predictor of the
underlying parameter.

- **New constructor:
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md).**
  Captures
  `(formula, by, species, link, init, bounds, priors, est_phase)` for
  one process parameter. Anything
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  understands works: `~ 1`, `~ temp + PDO`, `~ poly(temp, 4)`,
  `~ I(temp^2)`, `~ splines::ns(temp, df = 4)`, `~ temp * PDO`, etc.

- **Per-species formulas.** Register multiple specs against the same
  parameter via `linkages = list(log_K = list(spec_a, spec_b))` with
  each spec’s optional `species = ...` argument scoping it to a subset
  of stocks. The pooler unions the design columns across specs so
  there’s no duplication when species share covariates.

- **Priors.** First-class via
  [`prior_normal()`](https://grantdadams.github.io/Rceattle/reference/prior_normal.md),
  [`prior_lognormal()`](https://grantdadams.github.io/Rceattle/reference/prior_lognormal.md),
  [`prior_gamma()`](https://grantdadams.github.io/Rceattle/reference/prior_gamma.md),
  [`prior_beta()`](https://grantdadams.github.io/Rceattle/reference/prior_beta.md).
  The same constructors are available unprefixed (`normal()` /
  `lognormal()` / …) **only inside** the `priors = ...` argument via a
  private NSE data mask, so user code stays close to mathematical
  notation without masking
  [`base::gamma()`](https://rdrr.io/r/base/Special.html) /
  [`base::beta()`](https://rdrr.io/r/base/Special.html) at the package
  level. Priors can be a single value applied to every species, or a
  named list keyed by species id (and shortly, by `(species, sex)`).

- **Bounds.** Per-row `lower` / `upper` flow into
  `build_bounds()$lower$beta_linkage` /
  `build_bounds()$upper$beta_linkage`.

- **Growth** (von Bertalanffy / Richards) is the first process fully
  wired to the new pipeline.
  [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  gains a `linkages` argument and a string-named `fun` (`"empirical"` /
  `"vonBertalanffy"` / `"Richards"`); integer codes still work
  (`fun = 1` is shorthand for `fun = "vonBertalanffy"`) so existing
  scripts don’t need to be rewritten apart from substituting `fun =` for
  `growth_model =`.

  ``` r

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

- **TMB plumbing.** New `src/TMB/linkage.hpp` accumulator;
  `ceattle_v01_11.cpp` reads parallel `DATA_IVECTOR(linkage_*)` inputs
  plus a `DATA_MATRIX(linkage_X)` and writes a `growth_linkage_offset`
  tensor that is added (additively, on the log scale) to
  `growth_parameters`. Per-row prior densities contribute to slot 19 of
  the joint NLL (“Linkage-table priors” in `fit$quantities$jnll_comp`).

- **Documentation.** New vignette
  `vignette("environmental-linkages", package = "Rceattle")` walks
  through the API, prior families, species-keyed priors, per-species
  formulas, basis-expansion formulas, and the underlying pipeline.

- **Natural mortality** is the second process wired to the pipeline.
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
  gains a `linkages` argument keyed by `log_M1`; the offset is added on
  the log scale to `log_M1` inside the `M1_at_age` compute. A row’s
  `age_bin == NA` broadcasts the offset across ages; specific values pin
  it to that age slice.
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
  also gains string-form acceptance for `M1_model` and `M1_re` (parity
  with `build_growth(fun)`):

  ``` r

  build_M1(M1_model = "sex_age_invariant",   # or 1
           M1_re    = "ar1_age",             # or 4
           linkages = list(
             log_M1 = linkage_spec(formula = ~ temp, by = ~ species)
           ))
  ```

  Growth and M can be linked in the same fit; their rows share the same
  global linkage table and the same `beta_linkage` parameter vector.

- **Per-(species, sex) priors.** In addition to scalar and species-keyed
  priors, each `priors[[col]]` value may be a two-level nested list
  keyed first by species id then by sex id:

  ``` r

  priors = list(
    temp = list(
      `1` = list(`1` = normal(0, 0.1),
                 `2` = normal(0, 0.2)),    # sp 1 by sex
      `2` = normal(0, 0.5)                  # sp 2, both sexes
    )
  )
  ```

  Missing keys at either level resolve to “no prior” for that stratum.
  The validator checks the structure recursively and emits actionable
  error messages keyed by `priors$<col>$<species>[$<sex>]` paths.

- **Default `by = ~ species`.**
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  now defaults the `by` argument to `~ species`, so each linkage
  produces one coefficient per species without the user having to spell
  it out. Pass `~ species + sex` for per-(species, sex) coefficients, or
  `by = NULL` to share a single coefficient across every species/sex
  (the prior default). This matches the typical multispecies assessment
  use case where each stock has its own environmental sensitivity.

- **Recruitment** is the third process wired to the pipeline.
  [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  gains a `linkages` argument keyed by `log_R0`, `log_alpha`, or
  `log_beta`; the offset is added on the log scale to the corresponding
  parameter at every recruitment compute call site (hindcast, BRPs,
  projections, expected R). `log_R0` is meaningful for any `srr_fun`;
  `log_alpha` and `log_beta` only do work for SRRs that consume them
  (`srr_fun in c(2, 3, 4, 5)` – Beverton-Holt and Ricker).

  ``` r

  build_srr(
    srr_fun  = 2,
    linkages = list(
      log_alpha = linkage_spec(formula = ~ temp,
                               priors  = list(temp = normal(0, 0.5)))
    )
  )
  ```

  Growth, M, and recruitment can be linked in the same fit; their rows
  share the same global linkage table and the same `beta_linkage`
  parameter vector. End-to-end tests in
  `tests/testthat/tests-Dynamics/test-linkage-auto-map.R` verify that
  the base parameter (e.g., `log_growth_pars`, `log_M1`, `rec_pars`) is
  automatically mapped out (set to `NA`) when a linkage is active for
  that parameter, allowing the linkage intercept to define the base
  value.

  `tests/testthat/tests-Dynamics/test-recruitment-linkage.R` cover the
  analytical relations `R = R0 * exp(beta * temp[yr])` (mean R) and the
  `growth + M + recruitment` composition.

- **Soft deprecation in
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md).**
  The legacy column-index argument `M1_indices` and the env-driven
  structural integer codes `M1_model %in% c(4, 5)` are subsumed by the
  new `linkages = list(log_M1 = ...)` argument. Both still work for one
  release cycle, but emit a one-time warning that points users at the
  equivalent linkage-table call. The string aliases
  `"env_sex_invariant"` / `"env_sex_specific"` (added briefly on the dev
  branch) are removed; they were never released. The cpp’s `M1_beta` /
  `M1_mult` code path stays in place for now – both paths add additively
  to `log_M1` on the log scale – but do not use both for the same
  coefficient or you will double-count.

- **Roadmap.** Recruitment is next on the same pipeline; then
  random-effects pooling on `re_group` for hierarchical shrinkage. The
  legacy `M1_indices` / `M1_model = 4|5` paths retire when recruitment
  migrates.

### Scheduled removal (v4.5.0)

> **Schedule update (v4.4.0):** The removal originally targeted for
> v4.2.0 has been pushed to **v4.5.0** to give downstream users a longer
> migration window for the natural-scale linkage API rolled out in
> v4.4.0. The soft-deprecation warnings continue to point users at the
> equivalent linkage-table call. The cleanup checklist below is
> unchanged.

The soft-deprecated API surfaces below remain functional and emit
one-time warnings pointing users at the linkage table. They will be
**removed entirely in 4.5.0**. To migrate, replace:

| Legacy | New |
|----|----|
| `build_srr(srr_indices = ...)` | `build_srr(linkages = list(log_R0 = linkage_spec(...)))` |
| `build_srr(srr_fun = 1)` | `build_srr(srr_fun = 0)` + linkage on `log_R0` |
| `build_srr(srr_fun = 3)` | `build_srr(srr_fun = 2)` + linkage on `log_alpha` |
| `build_srr(srr_fun = 5)` | `build_srr(srr_fun = 4)` + linkage on `log_alpha` |
| `build_M1(M1_indices = ...)` | `build_M1(linkages = list(log_M1 = linkage_spec(...)))` |
| `build_M1(M1_model = 4)` | `build_M1(M1_model = 1)` + linkage on `log_M1` |
| `build_M1(M1_model = 5)` | `build_M1(M1_model = 2)` + linkage on `log_M1` |

**Cpp cleanup checklist** (search for `LEGACY` in
`src/TMB/ceattle_v01_11.cpp`):

- Drop `PARAMETER_MATRIX(beta_rec_pars)`, `PARAMETER_ARRAY(M1_beta)`,
  and the scratch vectors `srr_mult`, `beta_rec_tmp`, `env_rec_tmp`,
  `M1_mult`, `beta_M1_tmp`, `env_M1_tmp`.
- Delete the five `srr_env_mult` blocks (hindcast, BRPs, dynamic BRPs,
  projection, R_hat) and the `M1_mult.sum()` term inside the M1_at_age
  compute.
- Pass `Type(0.0)` for `srr_env_mult` at each `calculate_recruitment()`
  call site (or drop the parameter from the function signature in
  `recruitment.hpp` if no caller still needs it).

**R-side cleanup**:

- Remove `srr_indices` and `M1_indices` arguments from
  [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  and
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md).
- Reject `srr_fun %in% c(1, 3, 5)` and `M1_model %in% c(4, 5)` as
  unknown integer values in `.coerce_srr_fun()` and `.coerce_M1_arg()`
  (drop the `.SRR_DEPRECATED_FUNS` / `.M1_DEPRECATED_MODELS` constants).
- Remove the [`suppressWarnings()`](https://rdrr.io/r/base/warning.html)
  wrappers around internal
  [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  /
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
  re-callers in
  [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md),
  [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md),
  [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md),
  [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md),
  `project_no_F()`.

## Rceattle 4.0.3

### API

- New
  [`fit_control()`](https://grantdadams.github.io/Rceattle/reference/fit_control.md)
  constructor bundles the optimizer / sdreport / phasing knobs that
  previously cluttered
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)’s
  signature (`phase`, `getsd`, `bias.correct`, `use_gradient`,
  `rel_tol`, `loopnum`, `newtonsteps`, `getJointPrecision`,
  `getReportCovariance`, `verbose`, `TMBfilename`, `nlminb_control`).
  Pass the result via the new `fit_control` argument:

  ``` r

  fit <- fit_mod(
    data_list   = BS2017SS,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, loopnum = 1)
  )
  ```

  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)’s
  visible argument list shrinks from ~33 to ~22 args, so calls now read
  as model spec rather than a pile of optimizer flags.

- [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  emits a deprecation warning if any of the legacy control args are
  passed directly and forwards them into `fit_control` for the duration
  of the deprecation window. Truly unknown arguments still error with
  `Unused arguments to fit_mod(): ...` (no silent drops).

- Internal callers
  ([`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md),
  [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md),
  [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md),
  [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md),
  `project_no_F()`) now wrap their control args in `fit_control(...)`
  rather than passing them positionally.

### New methods

- S3 methods on the `"Rceattle"` class so a fit behaves like an R model
  object: [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html). With `df` set
  on `logLik`, \[stats::AIC()\] also works without a dedicated method.
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) is intentionally *not*
  defined: counting “observations” in a stock-assessment likelihood
  (composition cells, indices, catches, priors) is not well-defined, so
  \[stats::BIC()\] does not work — use AIC or domain-specific
  information criteria.

- [`plot.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/plot.Rceattle.md)
  is a thin dispatcher: `plot(fit, what = "biomass")` / `"ssb"` /
  `"recruitment"` / `"depletion"` / `"index"` / `"catch"` /
  `"selectivity"` / `"mortality"` / `"data"`. `...` is forwarded to the
  underlying `plot_*()` function.

- `residuals.Rceattle(type = ...)` returns a long-format data frame with
  rows from one or more of the four fitted data sources: `"index"` and
  `"catch"` (log-scale by default; switch with `scale = "natural"`),
  `"comp"` (Pearson on fitted proportions, with the `Age0_Length1` flag
  preserved), and `"caal"` (Pearson on fitted proportions, with both the
  conditioning `Length` and the age `Bin`). `type = "all"` returns all
  four stacked.

### Documentation

- README now has a self-contained *Getting started* block that fits a
  bundled model and exercises every new S3 method, so first-time users
  on the pkgdown site / CRAN no longer have to bounce to the GitHub wiki
  to see a working example.

- Vignette 8 (“Model parameterizations”) is being expanded to fill in
  coverage gaps surfaced during the 4.0.2 release audit:

  - M1 random-effects modes (`M1_re = 0..6`).
  - Full stock-recruit section (Mean / BevertonHolt / Ricker, env-linked
    forms, Beta prior on steepness via `srr_est_mode = 3`).
  - Composition likelihoods (Multinomial, Dirichlet-multinomial, CAAL).
  - Catchability equations for `Catchability = 4` (Power), `5`
    (Environmental), and `6` (AR1).
  - Selectivity equations for `Selectivity = 6` (2DAR1) and `7` (3DAR1).
  - Internal growth model (`growth_model = 1` von Bertalanffy / `2`
    Richards).
  - Initial-age-structure modes (`initMode = 0..4`).
  - Harvest control rules (HCR = 0..7) — possibly cross-linking vignette
    0 §9 rather than duplicating.

## Rceattle 4.0.2

### Bug fixes

- [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  now stops with a clear message if a user requests `msmMode = 3:9`
  (Kinzey & Punt 2009 functional responses – Holling I/II/III, predator
  interference, predator preemption, Hassell-Varley, Ecosim). Those code
  paths are unvalidated against the current parameter set and should not
  be used for advice.
- Plotting functions (`plot_timeseries`, `plot_selectivity`,
  `plot_mortality`, `plot_maturity`, `plot_b_eaten`,
  `plot_b_eaten_prop`, `plot_m_at_age`, `plot_m2_at_age_prop`, `plot_f`,
  `plot_ration`, `plot_index`, `plot_catch`, `plot_logindex`,
  `plot_indexresidual`, `plot_comp`, `plot_selectivity_vs_maturity`,
  `plot_stock_recruit`) now restore graphics
  [`par()`](https://rdrr.io/r/graphics/par.html) on exit instead of
  mutating the user’s device permanently.
- Replaced `T` / `F` shortcuts with `TRUE` / `FALSE` throughout the
  package source (~60 occurrences).
- Replaced `.data$Bin` / `.data$Length` inside
  [`tidyr::pivot_wider`](https://tidyr.tidyverse.org/reference/pivot_wider.html)
  arguments with quoted strings (tidyselect 1.2.0 deprecation).
- `examples/Georges_bank_example.R` now calls
  [`plot_mortality()`](https://grantdadams.github.io/Rceattle/reference/plot_mortality.md)
  instead of the long-removed `plot_mort()`.
- `_pkgdown.yml` “Get started” / overview navbar links now point at the
  actual generated `articles/Rceattle-overview.html` (was
  `Rceattle-overview-4_17_2025.html`, which 404’d).
- `README.md` example links now reference the correct
  underscore-separated filenames (`Fit_2018_GOA_multi-species_model.R`
  etc.) – the previous space-encoded URLs 404’d on GitHub.
- Removed duplicate `\seealso{}` block in `?Rceattle-package` by
  consolidating `R/Rceattle.R` into `R/Rceattle-package.R`. Both files
  declared `_PACKAGE`, so roxygen2 emitted the auto-generated links
  twice.
- Added [`graphics::box`](https://rdrr.io/r/graphics/box.html) to
  package imports (cleared the lone `R CMD check` NOTE for `plot_data`).
- TMB: terminal-age geometric series now includes `Finit` in the
  denominator for fished initial-equilibrium modes (`initMode = 3` or
  `4`), correcting a bias in the plus-group N-at-age when `Finit > 0`.

### Documentation

- Added Wassermann et al. (2025) cannibalism / Pacific hake reference to
  `inst/CITATION` and `?Rceattle-package`.
- `initMode` accepts integer codes or string aliases.
- `HCR` accepts integer codes or string aliases.

### Tests

- `tests/testthat/test-subdir-folders.R` no longer calls
  [`testthat::test_dir()`](https://testthat.r-lib.org/reference/test_dir.html)
  for each subdirectory. Nested `test_dir()` inside `test_check()`
  triggered an “evaluation nested too deeply: infinite recursion” abort
  inside `rlang`’s trace deparser whenever any test failed, masking the
  real failure. Subdirectory test files are now discovered with
  [`list.files()`](https://rdrr.io/r/base/list.files.html) and pulled in
  via [`source()`](https://rdrr.io/r/base/source.html) so they register
  against the outer reporter directly.
- `tests/testthat.R` now wraps
  [`library()`](https://rdrr.io/r/base/library.html) calls in
  [`suppressPackageStartupMessages()`](https://rdrr.io/r/base/message.html)
  so transient build-version notices (e.g. “package ‘dplyr’ was built
  under R version 4.5.2”) do not get captured as test warnings whose
  backtraces then crash rlang’s expr_deparse at end-of-run.
- `tests/testthat/test-Dynamics/test-initial-dynamics` evaluates the
  different starting conditions.

### Parallelism

- [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  and
  [`check_mse()`](https://grantdadams.github.io/Rceattle/reference/check_mse.md)
  now use
  [`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html)
  on a PSOCK cluster instead of `foreach::foreach(...) %dopar%`.
  - The `%dopar%` path triggered
    [`rlang::expr_deparse`](https://rlang.r-lib.org/reference/expr_print.html)
    infinite recursion under nested `test_that` backtraces because each
    `foreach` invocation captures call frames that recurse during error
    formatting. PSOCK workers are clean R processes with no captured
    promise chains, so the issue does not occur.
  - PSOCK clusters work identically on Windows and macOS/Linux.
  - [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
    gains a `cores` argument (default `NULL` picks
    `parallel::detectCores() - 6`); both functions cap at 2 cores when
    `_R_CHECK_LIMIT_CORES_` is set so they comply with CRAN’s R CMD
    check limit.
  - `foreach` and `doParallel` removed from `Imports:`.

### Installation / dependencies

- `TMBhelper` moved from `Imports:` to `Suggests:`. Rceattle now uses an
  internal `.fit_tmb()` wrapper that delegates to
  [`TMBhelper::fit_tmb()`](https://rdrr.io/pkg/TMBhelper/man/fit_tmb.html)
  when the (optional, GitHub-only) package is installed and otherwise
  falls back to a
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) +
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html) path.
  This removes the largest install-friction barrier for new users.
- `dplyr` moved from `Depends:` to `Imports:`. The package no longer
  attaches `dplyr` to the user’s search path on load (so it no longer
  masks [`stats::filter`](https://rdrr.io/r/stats/filter.html) /
  [`stats::lag`](https://rdrr.io/r/stats/lag.html)).
- `kableExtra` dropped from `Suggests:`; vignettes now use
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) for table
  rendering.
- `quarto` dropped as a vignette engine; all vignettes are now `.Rmd`.

### API

- [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  no longer has `om = ms_run`, `em = ss_run` defaults. Both arguments
  are now required and validated as objects of class `"Rceattle"` before
  the MSE loop runs. Calling
  [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  with no arguments previously produced a confusing “object ‘ms_run’ not
  found” error; it now stops with a clear message.

### New methods

- [`print.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle.md)
  and
  [`summary.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/summary.Rceattle.md).
  Auto-printing a fit inside knitr / RStudio / R Markdown previously
  dumped tens of MB of nested data and could trigger deep recursion
  errors during vignette rendering.

### Build / packaging

- Tightened `.Rbuildignore` (excludes `examples/`, `R/dev/`,
  `src/TMB/Dev/`, `.Rhistory`, `.claude/`, `.DS_Store`, build tarballs,
  `.Rcheck` directories).
- Tightened `.gitignore` to catch all `*.o` / `*.so` / `*.dll` and
  `*.Rcheck/` directories.
- Suppressed a benign clang `-Wfixed-enum-extension` warning by scoping
  the diagnostic pragma around `#include <TMB.hpp>` rather than via
  global Makevars flags.

## Rceattle 4.0.1

The 4.0.1 development cycle reorganized several `data_list` columns and
`fit_mod` / `build_*` arguments. Models or data files saved against
earlier 4.x revisions may need updating; see the renames below. Compiled
from `inst/Running_list_of_updates.qmd` plus the `dev` branch commit
log.

### Data renames

- `Pyrs` -\> `ration_data` (the old name is still accepted on read, but
  is silently renamed).
- `UobsWtAge` -\> `stom_prop_data`.
- `fsh_biom` -\> `catch_data`.
- `srv_biom` -\> `index_data`.
- `Nselages` -\> `N_sel_bins` (in `fleet_control`).
- `Sel_norm_bin1` / `Sel_norm_bin2` \<- `Age_max_selected` /
  `Age_max_selected_upper` (selectivity normalization bins).
- `Age_first_selected` -\> `Bin_first_selected` (in `fleet_control`).
- `sel` -\> `sel_at_age` (model report).
- `fleet_control` now carries a `Month` column (month of observation for
  indices / fisheries).

### API renames

- `build_M1`: `M1_prior_mean` -\> `M_prior`, `M1_prior_sd` -\>
  `M_prior_sd`.
- `build_srr`: `srr_prior_mean` -\> `srr_prior`; `R_hat_endyr` replaced
  by `srr_hat_styr` / `srr_hat_endyr`.
- `fit_mod`: `suit_meanyr` replaced by `suit_styr` / `suit_endyr`.
- `initMode` semantics revised: 0 = free-parameter N-at-age, 1 =
  unfished equilibrium with no devs, 2 (default) = unfished equilibrium
  with initial devs, 3 = fished equilibrium with initial devs. Type 4
  (“non-equilibrium scaled”) added later.

### New features – composition and diet likelihoods

- Dirichlet-multinomial composition likelihood. Selected per fleet via
  `fleet_control$Comp_loglike = 1` (or `"DirichletMultinomial"`).
- Conditional age-at-length (CAAL) data path, with `CAAL_loglike` /
  `CAAL_weights` controls in `fleet_control`. CAAL data also flow
  through
  [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md)
  for simulation testing.
- `Diet_loglike` switch on the bioenergetics control sheet selects
  between multinomial (0) and Dirichlet-multinomial (1) for diet
  composition.
- Other-food diet proportion estimates added to the model report.
- Weighted-mean diet data path (annual proportion of prey-at-age in
  predator-at-age averaged across years).

### New features – selectivity, catchability, growth

- Hake non-parametric selectivity (`Selectivity = "Hake"` or `5`), after
  Taylor et al.
- `2DAR1` (`= 6`) and `3DAR1` (`= 7`) selectivity parameterizations,
  after Cheng et al. (2024).
- `Catchability = 6` (“AR1”): annual AR1 catchability deviates fit to an
  environmental index, after Rogers et al. (2024) for the GOA pollock
  model. Environmental q-link (`Catchability = 5`) also exposed.
- Internal growth model. See
  [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  and the `growthFun` argument to
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).
  `alpha_wt_len` / `beta_wt_len` added to the data control sheet.
  Length-based suitability (`suitMode = 1` / `2` / `3` / `4` / `5` /
  `6`) wired through to use the estimated growth model. Comparison with
  WHAM growth implemented under `tests/comparison/`.
- Predator-specific suitability mode (different `suitMode` per
  predator).
- Suitability calculation now uses configurable year ranges (`suit_styr`
  / `suit_endyr`) instead of “mean year”.

### New features – recruitment and reference points

- Beta-distributed prior on Beverton-Holt steepness, available via
  `srr_est_mode = 3`.
- M1 random effects with optional environmental linkage; `M_prior` /
  `M_prior_sd` priors carried through
  [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md).
- [`remove_F()`](https://grantdadams.github.io/Rceattle/reference/remove_F.md)
  function returns a fitted model with F set to 0 – used internally for
  dynamic reference point calculation.
- `DynamicHCR = TRUE` in
  [`build_hcr()`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
  to switch from static to dynamic SB0 reference points.
- CMSY harvest control rule (`HCR = 1`): maximize joint catch across
  species, optionally constrained to keep depletion above `Plimit`.
- PFMC Category 1 40-10 ABC HCR (`HCR = 6`) using `Pstar` / `Sigma`
  uncertainty buffer.
- SESSF Tier 1 HCR (`HCR = 7`).
- Iterative multi-species HCRs: `HCRorder` controls the order in which
  species F is solved (e.g. predators before prey) inside
  [`build_hcr()`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md).

### New features – MSE and projection

- [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  now writes per-simulation `.rds` files when `dir` is specified, for
  streaming-friendly long runs.
  [`load_mse()`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  reads those back.
- [`check_mse()`](https://grantdadams.github.io/Rceattle/reference/check_mse.md)
  validates which OM/EM simulations converged.
- [`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  produces a per-fleet performance-metric table (mean catch, IAV,
  P(closed), MSE on SSB, P(F \> Flimit), P(SSB \< SSBlimit), terminal
  depletion, …).
- MSE function now supports `cap` (catch cap), `catch_mult` (catch
  multiplier), `rec_trend` (linear projected recruitment trend),
  `fut_sample` (future sampling effort), per-fleet `assessment_period` /
  `sampling_period`, `regenerate_past` (refit EM to OM-simulated past
  data), and `timeout`/`try`-error handling per simulation.
- `Recruitment_and_fixed_F_projections.R` and `Simulation_testing.R`
  examples added.

### New features – diagnostics and tooling

- [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md)
  function to perturb starting values and re-fit, for
  global-vs-local-minimum diagnostics.
- [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md)
  peels with optional `nyrs_forecast`.
- [`model_average()`](https://grantdadams.github.io/Rceattle/reference/model_average.md)
  for averaging derived quantities across multiple fitted models, with
  optional bootstrap uncertainty.
- [`compare_sim()`](https://grantdadams.github.io/Rceattle/reference/compare_sim.md)
  and
  [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md)
  for parametric simulation testing.
- `McAllister-Ianelli-reweighting.R` example for composition
  reweighting.
- TMB log-likelihood pieces (unweighted) added to the report for
  composition diagnostics.
- `Selectivity = "Fixed"` (`= 0`) for empirically supplied selectivity
  blocks via the `emp_sel` data sheet.
- `TMBfilename` argument to
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  to point at an alternate `.cpp` during development.

### Behavior changes

- Removed accumulation-age switches in `fleet_control`. Selectivity
  normalization is now controlled via `Age_max_selected` (i.e.
  `Sel_norm_bin1`) on a per-fleet basis instead of always normalizing by
  the maximum-selectivity age.
- `NA` values inside the valid age/length range of composition data are
  now coerced to 0 with a warning (previously silently dropped or
  errored).
- Selectivity dimensioning switched from age- to bin-indexed for the
  non-parametric and 2D/3D AR1 forms (driven by `Selectivity_dimension`
  and `N_sel_bins`).
- Age-error and age-transition matrices are now dimension-checked
  against `nages` at
  [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  time.

## Rceattle 4.0.0

- CEATTLE TMB version 4.0.0. See Adams et al. (2022), *Fisheries
  Research*, 251, 106303.
