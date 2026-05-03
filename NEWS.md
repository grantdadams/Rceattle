# Rceattle 4.0.3 (in development)

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
