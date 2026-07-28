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
  faster than an unoptimized debug build, bit-identically (measured; see `dev/PERF-findings.md`).
  A normal install was always `-O2`; only `load_all` was `-O0`. To debug the C++ line-by-line
  (gdb/lldb), start R with `RCEATTLE_DEBUG_CPP=1` to restore the `-O0` build. **Any absolute fit
  timing must state its build** — an `-O0` number overstates real cost ~10x.
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
  `pkgdown.yaml`. There is **no lint config** and no coverage gate.
- **Slash commands** wrap these: `/recompile`, `/test [file]`, `/document`, `/check`,
  `/golden-check` (defined in `.claude/commands/`).

## Layout

- **`R/`** — numbered by pipeline/collation order: `0-*` build/prep helpers,
  `1-*` data checks, `2..5-*` params/map/bounds/rearrange, `6-*` fit + rename output,
  `7-*` plotting, `8-*` sim, `9-*` retro/jitter, `10-*` MSE, `11-*` model averaging.
  **The numeric prefixes are meaningful — don't renumber or rename wholesale.**
- **`src/TMB/`** — `ceattle_v01_11.cpp` is the main model (~3,810 lines, numbered section
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
  repo was documented with roxygen2 7.3.3 (`DESCRIPTION: RoxygenNote: 7.3.3`), but a newer
  local roxygen2 (e.g. 8.0.0) rewrites *every* `man/*.Rd` and swaps `RoxygenNote:` for
  `Config/roxygen2/version:` — a huge spurious diff. After `document()`, `git checkout`
  the unrelated `man/` + `DESCRIPTION` churn and keep only the `.Rd` you meant to change.
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
  `ceattle_v01_11.cpp` (`JNLL_INDEX`=0 … `JNLL_LINKAGE_RE`=20, `JNLL_N_ROWS` dimensions
  the matrix) — refer to a row by its constant, never a bare integer. Their **display
  names** live separately in `R/6-rename_output.R` (~L130–151) and are kept in sync by
  hand. If you add/reorder a likelihood component, update both the enum and the name vector.
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
  NOT covered by `/golden-check`.** `retrospective()`, `jitter()`, `self_test()`,
  `profile()`, `run_mse()`, and `remove_F()` all re-invoke `fit_mod()` via this one helper,
  which rebuilds the HCR / SR / M1 / growth specs from a source `data_list` and exposes each
  per-caller divergence as a named override. The four golden models exercise none of these
  paths, so a change that could move a refit needs `dev/verify-refit-like.R` (before/after
  bit-identity across all six entry functions, including a multispecies MSE) — golden-check
  alone proves nothing here.
- **`run_mse()` pins the OM's stock-recruit and suitability reference windows to the
  pristine `om$`** (not the advancing `om_use$`) so the hindcast does not drift through the
  projection — essential for multispecies, whose predation suitability must stay fixed. In
  `.refit_like()` these are the `srr_mse_switchyr` / `srr_hat_styr` / `srr_hat_endyr` /
  `suit_styr` / `suit_endyr` overrides; the EM instead advances `srr_mse_switchyr` to its
  current assessment `endyr` each iteration. The invariant (MSE must not perturb the
  hindcast, under any `simulate_data` / `sample_rec`) is checked by
  `dev/verify-mse-hindcast-invariant.R`.
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
  priors, est_phase)`. RHS forms: covariate (`~ temp`), time block (`~ block(Year)`),
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
  current dev version until release.
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

- **Data-workflow + linkage-grammar effort** (PRs 0–7): **complete and merged onto
  `dev-data-workflow`** (bumped to 5.0.0; not yet released to `main`). The linkage grammar,
  column schema, `build_data()`/`model_config()`, `save_config()`/`load_config()` +
  `fit_mod(config=)`, the C++ legibility pass + `JnllRow` enum, the
  `build_growth(sd_plus_group=)` WHAM/SS3 feature, the `mse_summary()` per-entity reshape,
  the `.refit_like()` collapse, and the developer-guide expansion all shipped. Tier D2 (the
  `linkage.hpp` accumulator merge) was **declined** (net-negative). Roadmap + historical
  record: `dev/PLAN-data-workflow-and-linkage-grammar.md`; forward backlog + archive index:
  `dev/README.md`.
- **Documentation-quality + roxygen-accuracy pass** (on `dev-data-workflow`): the PR 1–7
  user-facing doc surface was reviewed against three criteria (not AI-verbose,
  current-capabilities-not-changelog, scientist-legible), and a technical-accuracy audit
  cross-checked the `fit_mod` / `build_srr` / `build_M1` / `build_growth` /
  `build_catchability` / `build_selectivity` / `build_composition` / `linkage_spec` /
  `build_hcr` roxygen math against the C++ — fixing accuracy issues incl. the `build_hcr`
  SESSF/NPFMC/PFMC reference-point formula errors, the Ricker β/1e6 reparametrization note,
  and the inert `avgnMode` switch. A companion **code fix** removed the HCR-4 `Fmult`
  double-application (`ceattle_v01_11.cpp` case 4). Doc-only otherwise.
- **Older paused work:** a multi-PR accessibility / code-review refactor (branch
  `accessibility-and-code-review`), plan in
  `~/Downloads/HANDOFF-accessibility-refactor-implementation.md`. Read it before resuming;
  do not start from scratch.
