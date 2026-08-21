# Rceattle — Claude Code guide

Rceattle (v4.7.0, GPL) fits the **CEATTLE** single- and multi-species, climate-linked,
age-structured stock assessment model. The likelihood is a **TMB / C++** model
(`src/TMB/`); everything around it (data prep, fitting, projection, MSE, diagnostics,
plotting) is R.

## Dev workflow

This is a TMB package, so the C++ must be compiled before the R can run.

```r
pkgload::load_all(".", quiet = TRUE)   # recompiles the TMB DLL + loads R/ (after any .cpp/.hpp edit)
devtools::document(quiet = TRUE)       # regenerate man/*.Rd + NAMESPACE after roxygen changes
NOT_CRAN=true Rscript -e 'devtools::test()'   # full suite (NOT_CRAN runs the skip_on_cran blocks)
rcmdcheck::rcmdcheck()                 # what CI runs (slow; usually backgrounded)
```

- **Toolchain:** prefix R compile/check commands with `export PATH=/usr/bin:$PATH` (system
  toolchain first — avoids a Homebrew clang/gfortran shadowing the TMB build). The repo's
  sessions do this consistently.
- **Editing `src/TMB/*.cpp` or `*.hpp` has no effect until you reload** — `load_all()`
  recompiles via `src/TMB/compile.R`; add `compile = FALSE` for R-only changes. Compiled
  artifacts (`*.o` ~77 MB, `*.so`) are gitignored — never commit them.
- **Dev builds are compiled at `-O2` (fast), not pkgbuild's default `-O0`.** The repo
  `.Rprofile` sets `options(pkg.build_extra_flags = FALSE)`, so `load_all()` compiles the TMB
  model with the same optimization as a production `R CMD INSTALL` — `fit_mod()` runs ~10x
  faster than an unoptimized debug build, bit-identically (measured).
  A normal install was always `-O2`; only `load_all` was `-O0`. To debug the C++ line-by-line
  (gdb/lldb), start R with `RCEATTLE_DEBUG_CPP=1` to restore the `-O0` build. **Any absolute fit
  timing must state its build** — an `-O0` number overstates real cost ~10x.
- **After a `.cpp`/`.hpp` edit, recompile before testing and run the suite serially.**
  `DESCRIPTION` sets `Config/testthat/parallel: true`, and the parallel workers cannot load
  a freshly rebuilt DLL — `devtools::test()` aborts before running anything with
  `testthat subprocess failed to start … getDLLRegisteredRoutines.DLLInfo`. That is a
  toolchain failure, not a test failure. Run `pkgload::load_all(".")` first, then test with
  `TESTTHAT_PARALLEL=false`. (`test-coverage.yaml` forces serial for a related reason.)
- **Tests** run with `NOT_CRAN=true`. To run one file ad-hoc, make the env's parent the
  package namespace so internal (non-exported) helpers resolve, then source the shared
  helpers into it:
  `e <- new.env(parent = asNamespace("Rceattle")); for (f in list.files("tests/testthat","^helper",full.names=TRUE)) sys.source(f, e); testthat::test_file(f, env = e)`.
  A plain `new.env()` makes tests fail spuriously with `could not find function "data_check"`.
  `options(testthat.max_fails = Inf)` shows all failures.
- **Tests that run a real `fit_mod()` optimization must use `testthat::skip_on_cran()`** so
  plain `R CMD check` / `devtools::test()` stay fast (only `NOT_CRAN=true` runs them); leave
  fast unit tests unguarded.
- CI = `.github/workflows/R-CMD-check.yaml` (r-lib actions, multi-OS matrix) +
  `pkgdown.yaml` + `test-coverage.yaml`. There is **no lint config** and no coverage gate.
  `pkgdown.yaml` runs on push to `main` and on PRs targeting `main` — **not** on `dev`, so a
  pkgdown break only surfaces once the PR is open.
- **Slash commands** wrap these: `/recompile`, `/test [file]`, `/document`, `/check`,
  `/golden-check` (local, not tracked in the repo).

## Layout

- **`R/`** — numbered by pipeline/collation order: `0-*` build/prep helpers,
  `1-*` data checks, `2..5-*` params/map/bounds/rearrange, `6-*` fit + rename output,
  `7-*` plotting, `8-*` sim, `9-*` retro/jitter, `10-*` MSE, `11-*` model averaging.
  **The numeric prefixes are meaningful — don't renumber or rename wholesale.**
- **`src/TMB/`** — `ceattle.cpp` is the main model (~3,810 lines, numbered section
  index); process logic lives in headers (`recruitment.hpp`, `selectivity.hpp`,
  `predation.hpp`, `growth.hpp`, `linkage.hpp`, `comp_osa.hpp`, `helper_functions.hpp`,
  `bioenergetics.hpp`, `diet_data.hpp`).
- **`tests/testthat/`** — **flat layout**: every test is a top-level `test-<area>-<topic>.R`
  (e.g. `test-selectivity-logisticpm.R`, `test-likelihood-index-covariance.R`); the old
  `tests-Dynamics/` / `tests-Selectivity/` … subdirectories were removed. Shared
  `helpers-*.R` / `fixtures/` sit alongside. Fast fixtures: `make_test_data()`
  (single-species) or `make_msm_test_data()` (multispecies, incl. diet) with
  `estimateMode = 3` build a non-optimized object. `tests/comparison/` holds WHAM
  cross-checks (not part of `test_check`).
- **`vignettes/`** are `eval = FALSE` (they only need to render). `data/` has bundled
  example datasets (`BS2017SS`, `BS2017MS`, `GOA2018SS`, …).

## Conventions & traps

- **Commits: plain messages, no AI-attribution / `Co-Authored-By` trailer.**
- Roxygen uses markdown; run `devtools::document()` after touching `@`-docs. **Trap:** the
  repo is documented with roxygen2 **8.1.0** (`DESCRIPTION: Config/roxygen2/version`).
  A *different* local roxygen2 rewrites *every* `man/*.Rd` — an older one also swaps the
  key back to the legacy `RoxygenNote:`, leaving both present and contradicting each other.
  After `document()`, check `git diff DESCRIPTION` first: if the version key moved, the
  `man/` churn is the version, not your change. Either regenerate everything under one
  version deliberately, or `git checkout` the unrelated churn and keep only the `.Rd` you
  meant to change.
  Give internal helpers `@noRd` (not just `@keywords internal`) so they generate no `.Rd`.
  **Second trap:** never insert a helper between a function's roxygen block and its
  definition. Contiguous `#'` lines are ONE block and bind to whichever object follows, so
  the helper silently steals the `@export` / `@importFrom` tags and the original function
  loses them. Tests won't catch it (NAMESPACE isn't regenerated) — only the next
  `document()` reveals it. Put helpers above the block or after the function.
- **Released package — preserve the public API.** Rceattle is shipped (v4.7.0, has users;
  see `cran-comments.md` / `inst/RELEASE-CHECKLIST.md`). Exported `build_*()` args carry
  deprecation paths (e.g. SRR codes `1/3/5`) — **deprecate / keep back-compat rather than
  deleting public surface.** Internal refactors are free as long as golden-reference
  equivalence holds.
- **Match the surrounding style.** Canonical references: `src/TMB/recruitment.hpp` (the
  Doxygen-documented header to emulate) and any `R/*.R` + its `tests/testthat/*` pair. The
  codebase favors explanatory section headers and Doxygen on the C++ — match local comment
  density; don't strip comments.
- **`jnll_comp` likelihood rows are addressed by the `JnllRow` enum** in
  `ceattle.cpp` (`JNLL_INDEX`=0 … `JNLL_LINKAGE_RE`=20, `JNLL_N_ROWS` dimensions
  the matrix) — refer to a row by its constant, never a bare integer. Their **display
  names** live separately in `R/6-rename_output.R` (~L130–151) and are kept in sync by
  hand. If you add/reorder a likelihood component, update both the enum and the name vector.
- **`Index_distribution` has a second hand-synced registry.** A family added to
  `index_distribution_map` (`R/0-switches.R`) must also be classified in
  `.index_rows_natural_scale()` in the same file, which is what `residuals(type =
  "pearson")`, `plot_index()`'s observation interval and `plot_indexresidual()` read to
  pick a scale. A natural-scale family that misses it does not error — it silently gets
  the log-scale formula, whose `sigma^2/2` is then a number the size of the index
  squared. `test-index-natural-scale-downstream.R` enumerates the map and is the net.
- **Fleets sharing a `Selectivity_index` / `Q_index` share ONE parameter block.**
  `adjust_map_shared_params()` copies the donor fleet's map slice over the rest of the
  group, so any per-fleet setting that differs within a group is silently taken from the
  donor (`Sel_start_year`, `Bin_first_selected`, `N_sel_bins`, …). The donor is the first
  *estimated* fleet — an `Off` fleet's slice is all `NA` and must never lead. Penalties and
  priors on a shared block are accumulated once, on the lead fleet (`flt_sel_lead` /
  `flt_q_lead`); without that gate they are counted once per sharing fleet. GOA2018SS is
  the worked example (fleets 1 & 7 share selectivity; 9 & 10 share selectivity *and* q).
- **`fleet_control$Fleet_code` must equal the row number.** Per-fleet parameter and map
  arrays are dimensioned by `nrow(fleet_control)` but indexed by `Fleet_code`, so a
  mismatch attaches parameters to the wrong fleet. `data_check()` enforces it; read
  `fleet_control` columns by row index `i`, never by `flt`.
- **Selectivity bin columns are dimension-aware.** `Bin_first_selected`, `N_sel_bins`,
  `Sel_norm_bin1/2`, `Sel_pen_first_bin`, `Sel_pen_last_bin`, `Sel_cap_bin` are bin indices
  on the fleet's own `Selectivity_dimension`: an absolute **age** for age-based fleets
  (offset by `minage`) or a 1-based **length-bin ordinal** for length-based (offset by 1).
  The cpp penalties loop over `nbins`, not `nages`.
- **A new `data_list` element needs `write_data()` / `read_data()` support too**, or it
  round-trips to nothing and the feature is silently lossy through the standard xlsx
  format — this is how `index_cov` was lost.
- **A new `fleet_control` column is defined once, in `R/0-column_schema.R`.**
  `.rce_column_schema()` is the single source of truth: its rows drive `switch_check()`
  defaults, `write_data()`/`read_data()` ordering, the meta sheet, the field dictionary,
  and (via `aliases`) the auto-upgrade of older column spellings on every entry into the
  pipeline. Add the column there and consume it by its **canonical** name — don't hardcode
  a default or a column order anywhere else. `test-schema-canonical.R` pins the aliases.
- **A slow fit is the model, not a regression.** A full `BS2017SS`
  `estimateMode = 0, phase = TRUE` fit is ~55 s: ~14 s for one optimization, ~3x for
  phasing, +~14 s for `sdreport`. That single fit has needed ~500–700 `nlminb` iterations
  since at least 2023 (bisected back through `main`), so don't hunt a recent cause. For dev
  loops `phase = FALSE, getsd = FALSE` gives ~14 s.
- **Numeric changes need golden-reference equivalence:** any edit that can move the fit
  must keep an example fit (e.g. `BS2017SS`) within tolerance and the suite green.
- **`fit_mod(estimateMode=)`:** 0 = hindcast + HCR projection, 1 = hindcast only,
  2 = projection-only from `inits`, 3 = build (`MakeADFun`) without optimizing, 4 = optimize
  with all params mapped out. **Mode 3 returns the real objective**, so `obj$fn()` /
  `obj$gr()` are usable for diagnosing a model before fitting it — the analogue of WHAM's
  `fit_wham(do.fit = FALSE)` and SAM's `sam.fit(run = FALSE)`. (Before v4.7.1 mode 3 shared
  mode 4's placeholder and returned an identically-zero gradient.) **Mode 4 still returns the
  placeholder `jnll = dummy*dummy`** — `build_map()` maps out every hindcast parameter there,
  so `dummy` is the only free parameter and the objective is a plumbing smoke test, not a
  likelihood. Don't read anything into a mode-4 objective, gradient, or Hessian.
- **`fit_control()` bundles the optimizer / uncertainty knobs** (`getsd`, `bias.correct`,
  `loopnum`, `newtonsteps`, `getJointPrecision`). Pass `getsd = FALSE` for fast dev/test fits
  (skips `sdreport`) — but then `sdrep` is NULL, so `vcov()` returns NULL and uncertainty
  bands are NA.
- **The diagnostic *refit* paths go through `.refit_like()` (`R/6-refit_like.R`) and are
  NOT covered by `/golden-check`.** Eight entry points re-invoke `fit_mod()` via this one
  helper — `retrospective()`, `jitter()`, `self_test()`, `profile()`, `run_mse()`,
  `remove_F()`, `sample_rec()`, and `reweight_comps()` — which rebuilds the HCR / SR / M1 /
  growth specs from a source `data_list` and exposes each per-caller divergence as a named
  override. The four golden models exercise none of these paths, so a change that could move
  a refit needs `tools/verify/verify-refit-like.R` (before/after bit-identity, including a
  multispecies MSE) — golden-check alone proves nothing here. **That harness covers six of
  the eight**: `sample_rec(update_model = TRUE)` and `reweight_comps()` are not in it and
  must be checked by hand.
- **`run_mse()` pins the OM's stock-recruit and suitability reference windows to the
  pristine `om$`** (not the advancing `om_use$`) so the hindcast does not drift through the
  projection — essential for multispecies, whose predation suitability must stay fixed. In
  `.refit_like()` these are the `srr_mse_switchyr` / `srr_hat_styr` / `srr_hat_endyr` /
  `suit_styr` / `suit_endyr` overrides; the EM instead advances `srr_mse_switchyr` to its
  current assessment `endyr` each iteration. The invariant (MSE must not perturb the
  hindcast, under any `simulate_data` / `sample_rec`) is checked by
  `tools/verify/verify-mse-hindcast-invariant.R`.
- **Every observation and process error is drawn in a `SIMULATE{}` block beside the
  density that scores it** (`ceattle.cpp` sections 5.12b, 5.13, and one per likelihood
  slot; the multivariate draws live in `comp_sim.hpp`). `sim_mod()` no longer implements
  any observation model in R — it calls `obj$simulate()` once and writes the result back.
  A new likelihood family therefore owes a draw. Three rules the existing blocks follow:
  draw what the density assumes (bias-correction convention and scale included); REPORT
  under a `*_sim` name, because TMB never clears the report environment and a draw under
  the observed object's own name stays readable as the data; and don't draw what the
  model does not define (two densities on one latent — the AMAK/Ianelli stock-recruit
  penalty — have no distribution to draw from, so leave it and warn). A process draw also
  reports a `*_drawn_sim` mask, since the deviation arrays are REPORTed whole.
  `tools/verify/verify-sim-*.R` are the regression net here, as `verify-refit-like.R` is
  for the refit paths — `/golden-check` reaches none of it.
- Scratch outputs (`Rplots.pdf`, `*_osa.png`, `*.RDS` under `tests/comparison/`) are
  gitignored — don't commit them.

## Data assembly, configuration & the linkage grammar

Three subsystems sit in front of `fit_mod()`. The **developer guide**
(`vignettes/articles/developer-guide.Rmd`) is the deep-dive; the load-bearing facts:

- **Assembling / editing a `data_list` in R:** `build_data()` + `model_config()`
  (`R/0-build_data.R`) supply only what a configuration uses and print the model as a
  readable outline; `data_requirements()` (`R/1-data_requirements.R`, table in
  `R/1-data_requirements_table.R`) reports required/optional/ignored inputs *before*
  fitting; `write_template()` writes a blank correctly-shaped workbook. All of this reads
  the column schema (above) — the schema is authoritative, not these callers.
- **Run configuration:** `save_config()` / `load_config()` / `run_config()`
  (`R/0-save_config.R`) round-trip the full model + estimation + `fit_control()` setup to a
  commented, default-omitting **YAML** file. Apply with
  `fit_mod(data_list, config = load_config("run.yaml"))` — it fills only the args the caller
  didn't pass (**explicit args always win**), and every fit records its own config in
  `fit$run_config`. Formulas, `Rceattle_prior` objects, and nested `build_*()` specs all
  serialize, so a saved config reproduces a non-default fit bit-identically.
- **Linkage & priors grammar:** one formula system (`R/0-build_linkage.R`,
  `linkage_encode.R`, `linkage_table.R`, `R/0-priors.R`) covers time-varying / covariate /
  random-effect / prior structure on **every** process — recruitment, M, growth,
  catchability, selectivity, and (prior-only) Dirichlet-multinomial composition weights —
  attached through the `build_*()` constructors via `linkage_spec(formula, by, init,
  priors, est_phase)`. RHS forms: covariate (`~ temp`), time block (`~ cut(Year, ...)`;
  the fixed part is passed straight to `model.matrix()`, so there is no bespoke
  `block()` helper),
  random effect (`~ (1|Year)` IID, `rw()`, `ar1()`), or `~ 1` + `priors` for an
  intercept-only prior (keyed `` `(Intercept)` ``). A **prior on a selectivity /
  catchability / DM parameter** is the `~ 1` case. Inside `priors =`, the bare `normal()` /
  `lognormal()` / `gamma()` / `beta()` resolve to the `prior_*()` constructors via a data
  mask — intentional, don't "fix" them.
- **C++ side of linkages:** `linkage.hpp` holds five per-process accumulators
  (`rceattle_apply_{growth,M,recruitment,q,sel}_linkages`) whose process/param codes are in
  **lockstep with `R/0-linkage_encode.R`** — change one, change both. Each REPORTs a
  `*_linkage_offset` tensor; the constructive linkage tests
  (e.g. `test-dynamics-recruitment-linkage.R`) pin those offsets, so they — not
  `/golden-check` (whose 4 models carry no linkage rows) — are the regression net for any
  linkage-path edit.

## Domain vocabulary (use these exact terms in plots/docs/messages)

Rceattle implements fisheries stock assessments — match this vocabulary in axis labels,
documentation, and console messages; don't substitute lay phrasing.

- **Reference points:** Amendment-56 SPR proxies — F40% = max FABC, F35% = FOFL, B40% = BMSY
  proxy (Tier 3); Tier 1 uses estimated FMSY/BMSY. Don't write "MSY" where an SPR proxy is meant.
- **SSB** = female spawning-stock biomass; "recruitment deviations" (log-scale), not "recruitment error".
- **Selectivity:** name the form (logistic / double-normal / gamma / nonparametric / semi-parametric) —
  don't call every dome shape "double-normal".
- **Composition:** age comps, length comps, conditional age-at-length (CAAL). An ageing-error matrix
  applies only where age/CAAL data are fit — length-only stocks have none.
- **Data weighting:** Francis (2011) or McAllister–Ianelli or Dirichlet-Multinomial. **Diagnostics:** Mohn's rho (retrospective),
  OSA residuals, likelihood profiles.

## After any change — keep docs & version in sync

When a change affects behavior, the public API, or docs (not just local repo/tooling
files), update all three before considering it done:

- **`NEWS.md`** — add a bullet under the top version section (`# Rceattle <version>`,
  grouped under the right `## New features` / `## Bug fixes` subsection) describing the
  user-facing change.
- **`DESCRIPTION` `Version:`** — bump per semver: **patch** for bug fixes / docs, **minor**
  for new features, **major** for breaking API changes. Changes accumulate under the
  current dev version until release. **"Breaking" here means no back-compat path**: a
  removal that ships with a deprecation message and keeps old fits working is a *minor*
  bump even when `NEWS.md` files it under `## Breaking changes` (5.7.0 and 5.9.0 are the
  worked examples — `growth_re` was removed with a `switch_check()` deprecation message
  and a `fit_mod()` guard that drops retired parameter blocks from `inits`).
- **Don't open a NEWS section for a version `DESCRIPTION` never carries.** Intermediate
  headings written mid-branch are folded into the shipped section before merge, and the
  folding is noted rather than the entries renumbered — see the three folded gaps in
  "Active context". A section for an unshipped version also breaks any `(x.y.z)`
  cross-reference pointing at it.
- **Vignettes** — update any `vignettes/*.Rmd` whose documented behavior or API changed
  (they're `eval = FALSE`, so at minimum they must still render).

Full release / tag process lives in `inst/RELEASE-CHECKLIST.md`.

## Sibling model repo

`../Rceattle-models` holds the real assessments (EBS/GOA pollock, sablefish, arrowtooth,
plaice, POP, …) and is a live consumer of this package's API — a breaking change here
breaks scripts there. After one, sweep it:
`grep -rn "<removed_fn>" --include=*.R "../Rceattle-models"`. The v4.7.0 ggplot migration
is the cautionary case: `plot_*()` now return ggplot objects, so every
`plot_x(...); mtext(...)` chain in those scripts failed with "plot.new has not been called
yet". Caveats: some models there are partially implemented and don't run regardless, so not
every hit needs chasing; and a static sweep catches removed/renamed API but not behavioural
drift. Fitted `*.rds` are ~50 MB each — keep them out of git.

## Converting ADMB models

- **`dev_vector` sum-to-zero can't be replicated in TMB.** ADMB deviation vectors
  (`dev_vector`) carry a built-in *sum-to-zero* constraint with no direct TMB equivalent.
  When porting such a model, **turn off estimation of the first element** of the dev
  vector to recover identifiability — *unless* that dev vector already has an additional
  likelihood penalty (e.g. a normal / random-effects penalty), in which case the penalty
  pins it and all elements can stay estimated. (Turning a parameter off = mapping it to
  `NA` in `R/3-build_map.R`.)

## Active context

- **The 5.6.1 release (`dev` → `main`, PR #106).** `main` released through 4.9.1; `dev`
  carries everything from 4.10.0 onward — the linkage grammar, column schema,
  `build_data()`/`model_config()`, `save_config()`/`load_config()` + `fit_mod(config=)`,
  the `JnllRow` enum, `build_growth(sd_plus_group=)`, the `mse_summary()` per-entity
  reshape, the `.refit_like()` collapse, `reweight_comps()`, and the recruitment /
  stock-recruit work. Roadmap and historical record: the commit log and `NEWS.md`.
  Planning documents are kept locally under the untracked `dev/`.
  - **Three version gaps are deliberate, not lost entries.** `4.14.0` was a real
    `DESCRIPTION` version whose NEWS was folded into the 5.0.0 section; `5.2.0`–`5.2.4`
    likewise folded into 5.3.0; and `main`'s 4.9.0 / 4.9.1 are the same recruitment
    changes this line carries as 5.5.0 / 5.5.1, applied to both lines separately. Nothing
    was dropped, and no tag existed above 4.8.0 while these were in flight, so nobody
    could have installed an intermediate. Note the folding rather than renumbering.
  - **Result-changing changes that are not labelled breaking** — the mode-5 selectivity
    penalty fix (GOA Pacific cod SSB 2050 −14.1%), parameter bounds previously applied to
    the wrong parameters, composition weights warm-starting from `inits`, failed
    `run_mse()` simulations returning only a marker, the `mse_summary()` reshape, the
    recruitment fixes (`initMode = 0` random effects, the α-seeding fix, the Ianelli
    steepness prior), and `sim_mod()` drawing the index under the fleet's own
    `Index_distribution`. A model carrying GOA numbers forward needs a refit.
- **Older paused work:** a multi-PR accessibility / code-review refactor (branch
  `accessibility-and-code-review`), plan in
  `~/Downloads/HANDOFF-accessibility-refactor-implementation.md`. Read it before resuming;
  do not start from scratch.
