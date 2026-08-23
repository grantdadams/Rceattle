# CLAUDE.md — Rceattle

Rceattle fits **CEATTLE**, a single- and multi-species, climate-linked, age-structured stock
assessment model. Its output sets **US federal catch limits under Magnuson-Stevens**. A
silently-wrong number here becomes a wrong quota, and it will not announce itself as a crash.

**If a switch code, a unit, or a default is not stated in the schema or a vignette, ask. Do not
infer it.** Accuracy beats speed, every time.

The likelihood is a **TMB / C++** model (`src/TMB/`); everything around it — data prep, fitting,
projection, MSE, diagnostics, plotting — is R.

## Read first

1. `inst/dev/SESSION_HANDOFF.md` — what is in flight, what is verified, where to resume.
2. `R/0-column_schema.R` — the source of truth for every workbook column and switch.
3. `vignettes/articles/developer-guide.Rmd` — the pipeline, the switch system, the schema, the
   linkage grammar.
4. `inst/dev/TRAPS.md` — verified traps with the measured numbers behind them.
5. `inst/RELEASE-CHECKLIST.md` — the release and tag process.

---

## Hard rules

1. **Preserve the public API. Deprecate, never delete.** Rceattle ships and has users
   (`cran-comments.md`). Exported `build_*()` arguments carry deprecation paths. Internal
   refactors are free as long as golden-reference equivalence holds.
2. **A numeric change needs `/golden-check`.** Any edit that can move a fit must keep the four
   reference models within tolerance and the suite green.
3. **`/golden-check` does not cover the refit paths, the simulation draws, or any figure.**
   Use `/verify` to pick the right `tools/verify/*.R` harness. For a plotting change the net is
   `test-plot-*.R` plus a before/after `ggplot_build()` diff. See `inst/dev/TRAPS.md`.
4. **The column schema is the source of truth** for switch values, defaults, and column order.
   Add a column with `/new-column`; consume it by its canonical name. Don't hardcode a default,
   an allowed-value set, or a column order anywhere else.
5. **Behaviour, API, or doc change ⇒ `NEWS.md` + `DESCRIPTION` `Version:` + the affected
   vignette, in the same commit** — plus `_pkgdown.yml` when a documented topic appears or
   disappears. `/doc-sync` checks this. Repo tooling (`.claude/`, `tools/`, `.github/`) and
   developer notes (`inst/dev/`) are exempt: they change no behaviour a user can observe.
   **"Breaking" means no back-compat path**: a removal that ships with a deprecation message and
   keeps old fits working is a *minor* bump, even when NEWS files it under `## Breaking changes`.
   `growth_re` is the worked example — removed with a `switch_check()` deprecation message and a
   `fit_mod()` guard dropping retired parameter blocks from `inits`, and shipped as a minor.
6. **Never hand-edit `man/*.Rd` or `NAMESPACE`.** Run `/document`. Check `git diff DESCRIPTION`
   **first** — if the roxygen version key moved, the `man/` churn is the version, not your change.
7. **TMB source is inert until `pkgload::load_all(".")`.** Then test with
   `TESTTHAT_PARALLEL=false`; the parallel workers cannot load a freshly rebuilt DLL.
8. **Comments state current behaviour, units, and the assessment reason — never bug history.**
   See "Comments" below.
9. **Don't invent a switch code, a default, or a unit. Ask.** A plausible placeholder that
   survives into a fit is the worst failure mode this repo has.
10. **Fleet invariants.** `fleet_control$Fleet_code` must equal the row number — arrays are
    dimensioned by `nrow()` but indexed by `Fleet_code`, so read columns by row index `i`, never
    by `flt`. Fleets sharing a `Selectivity_index` / `Catchability_index` share ONE parameter block.
    Selectivity bin columns are indices on the fleet's own `Selectivity_dimension`.
11. **`nages` is a count of age bins, not the oldest age.** Ages run
    `minage .. minage + nages - 1`; age `a` sits at index `a - minage + 1`. `minage = 1` hides
    every confusion, and that is every bundled dataset and all three live assessments.
12. **`linkage.hpp` and `R/0-linkage_encode.R` are in lockstep.** Their process and param codes
    must match — change one, change both. This is the seam DSEM extends. In a linkage formula the
    fixed part goes straight to `model.matrix()`, so a time block is `~ cut(Year, ...)` — there
    is **no bespoke `block()` helper**, and adding one would be a second grammar.
13. **Commits: plain messages, no `Co-Authored-By` trailer.** Imperative subject, ≤72 chars; the
    body says *why*, and gives the numbers that changed.
14. **`../Rceattle-models` and `../GOA-ATF-ESP` consume this API.** Sweep after a breaking change
    (`/ecosystem-sweep`), and refit a real assessment when a sweep is not enough.
    See `inst/dev/SIBLING-REPOS.md`.

---

## Dev workflow

This is a TMB package, so the C++ must be compiled before the R can run.

```r
pkgload::load_all(".", quiet = TRUE)   # recompiles the TMB DLL + loads R/ (after any .cpp/.hpp edit)
devtools::document(quiet = TRUE)       # regenerate man/*.Rd + NAMESPACE after roxygen changes
NOT_CRAN=true Rscript -e 'devtools::test()'   # full suite (NOT_CRAN runs the skip_on_cran blocks)
rcmdcheck::rcmdcheck()                 # what CI runs (slow; usually backgrounded)
```

- **Toolchain:** prefix R compile/check commands with `export PATH=/usr/bin:$PATH` — system
  toolchain first, so a Homebrew clang/gfortran does not shadow the TMB build.
- `load_all()` recompiles via `src/TMB/compile.R`; add `compile = FALSE` for R-only changes.
  Compiled artifacts (`*.o` ~77 MB, `*.so`) are gitignored — never commit them.
- **To run one test file**, make the env's parent the package namespace so internal helpers
  resolve: `e <- new.env(parent = asNamespace("Rceattle"))`, then source the shared helpers into
  it. A plain `new.env()` fails with `could not find function "data_check"`.
- **A test that runs a real `fit_mod()` optimization needs `testthat::skip_on_cran()`** so plain
  `R CMD check` stays fast. Leave fast unit tests unguarded.
- **CI:** `.github/workflows/R-CMD-check.yaml` (multi-OS) + `pkgdown.yaml` + `test-coverage.yaml`.
  No lint config, no coverage gate. **`pkgdown.yaml` triggers on `main` only**, so a PR to `dev`
  gets no pkgdown CI — run `/pkgdown-check` yourself.
- **Slash commands:** `/recompile`, `/test [file]`, `/document`, `/check`, `/golden-check`,
  `/verify`, `/new-column`, `/doc-sync`, `/pkgdown-check`, `/ecosystem-sweep`, `/handoff`.

## Layout

- **`R/`** — numbered by pipeline order: `0-*` build/prep helpers, `1-*` data checks,
  `2..5-*` params/map/bounds/rearrange, `6-*` fit + rename output, `7-*` plotting, `8-*` sim,
  `9-*` retro/jitter, `10-*` MSE, `11-*` model averaging. **The numeric prefixes are meaningful
  — don't renumber or rename wholesale.**
- **`src/TMB/`** — `ceattle.cpp` is the main model (numbered section index); process logic lives
  in headers (`recruitment.hpp`, `selectivity.hpp`, `predation.hpp`, `growth.hpp`, `linkage.hpp`,
  `comp_osa.hpp`, `comp_sim.hpp`, `helper_functions.hpp`, `bioenergetics.hpp`, `diet_data.hpp`).
  `jnll_comp` rows are addressed by the **`JnllRow` enum** — refer to a row by its constant,
  never a bare integer. Display names live separately in `R/6-rename_output.R` and are kept in
  sync by hand; adding or reordering a component means updating both.
- **`tests/testthat/`** — **flat**: every test is a top-level `test-<area>-<topic>.R`. Shared
  `helpers-*.R` / `fixtures/` sit alongside. Fast fixtures: `make_test_data()` (single-species)
  or `make_msm_test_data()` (multispecies, incl. diet) with `estimateMode = 3` build a
  non-optimized object. `tests/comparison/` holds WHAM cross-checks (not part of `test_check`).
- **`vignettes/`** are `eval = FALSE` — they only need to render, so **nothing executes their
  code and no CI job checks it**. An API break in a vignette is invisible until a user hits it.
  `data/` has the bundled example datasets.
- **`inst/dev/`** — committed developer notes (handoff, traps, sibling repos, ADMB conversion).
  The untracked `dev/` is scratch and does not survive a clone.

## Plotting

- **The shared argument vocabulary lives in `R/7-plot_helpers.R`** and is documented once, in
  `?"rceattle-plot-args"` (`@inheritParams rceattle-plot-args`). `line_col`, `lwd`, `lty`,
  `alpha`, `species`/`spnames`, `minyr`/`maxyr`, `incl_proj`, `incl_mean`, `add_ci`,
  `model_names` each go through one resolver — `.as_colour()`, `.rce_line_params()`,
  `.rce_check_alpha()`, `.resolve_species()`, `.rce_year_filter()`, `.rce_proj_divider()`,
  `.rce_mean_line()`. Adding or converting a plotter means calling those, not writing a second
  reading of the same argument; that divergence is what `plot_f()` and `plot_selectivity()` were.
- **`line_col` and `lty` supply values for whatever the figure separates** — predators in
  `plot_b_eaten_prop()`, sex in `plot_ration()`, the year fan in `plot_selectivity()` — not
  always the model. Say which in the function's own `@details`.
- **Base-graphics arguments from before the ggplot migration are accepted and ignored**
  (`right_adj`, `top_adj`, `mod_cex`, `legend.pos`, `single.plots`, `theta`, `ymax`, `cex`), so
  the assessment scripts keep running. Keep them; document them as ignored, with the ggplot
  equivalent.

## Comments: write for a fisheries scientist, not a programmer

The reader knows fish and knows management. They may not know R or C++ idiom, and they were not
in the room when the decision was made. So:

- **Explain the assessment reason, not the code.** Why this bin edge, this year window, this
  constant.
- **Always give units and the convention.** "metric tons", "log scale", "female SSB".
- **State current behaviour, never bug history.** A comment is read by someone deciding how the
  code behaves *now*.
- **One or two lines.** If it needs a paragraph, it belongs in a vignette or `inst/dev/`.
- **Exceptions that stay:** a comment explaining why an old input path still executes is
  *behaviour*, not history — phrase it as behaviour. Literature citations (AMAK, Ianelli, ADMB,
  Punt, Holsman, Francis, Methot, Kinzey & Punt) are the **specification**; never strip them.
  In a regression test, provenance is what stops the next person deleting the test — put it in
  one header block above `test_that()`, never repeated inline.

```r
# Bad -- narrates the code
comp_weights[flt] <- 0        # set weight to 0

# Good -- states the assessment reason
# A Dirichlet-multinomial estimates its own weight inside the likelihood, so
# tuning it externally would compete with that estimate. 0 on the log scale
# is a weight of 1, i.e. no external multiplier.
comp_weights[flt] <- 0
```

**Match the surrounding style.** Canonical references: `src/TMB/recruitment.hpp` (the
Doxygen-documented header to emulate) and any `R/*.R` + its `tests/testthat/*` pair. The codebase
favours explanatory section headers and Doxygen on the C++ — match local comment density.

**Roxygen:** markdown, regenerated with `/document`. A `@param` is one sentence: what the
argument means, its allowed values, its default. Anything longer belongs in `@details` or a
vignette. Give internal helpers `@noRd` (not just `@keywords internal`) so they generate no
`.Rd`. **Never insert a helper between a function's roxygen block and its definition** —
contiguous `#'` lines are ONE block and bind to whichever object follows, so the helper silently
steals the `@export` and the `@importFrom` tags, and the original loses them. Tests won't catch
it (NAMESPACE isn't regenerated); only the next `document()` will. **Put helpers above the block
or after the function.**

## Domain vocabulary (use these exact terms in plots, docs and messages)

Match this vocabulary in axis labels, documentation, and console messages; don't substitute lay
phrasing.

- **Reference points:** Amendment-56 SPR proxies — F40% = max FABC, F35% = FOFL, B40% = BMSY
  proxy (Tier 3); Tier 1 uses estimated FMSY/BMSY. Don't write "MSY" where an SPR proxy is meant.
- **SSB** = female spawning-stock biomass. **"Recruitment deviations"** (log-scale), not
  "recruitment error".
- **Selectivity:** name the form — logistic / double-normal / gamma / nonparametric /
  semi-parametric. Don't call every dome shape "double-normal".
- **Composition:** age comps, length comps, conditional age-at-length (CAAL). An ageing-error
  matrix applies only where age/CAAL data are fit; length-only stocks have none.
- **Data weighting:** Francis (2011), McAllister–Ianelli, or Dirichlet-multinomial.
- **Diagnostics:** Mohn's rho (retrospective), OSA residuals, likelihood profiles.

## Reference implementations

What to consult when documenting a switch, shaping a workflow, or naming a process argument:

- [**WHAM**](https://github.com/timjmiller/wham) — `basic_info`/`input` argument documentation;
  how a large option surface is presented to an assessment author.
- [**SAM**](https://github.com/fishfollower/SAM) — the `conf` table: a compact, complete
  configuration reference.
- [**Stock Synthesis**](https://github.com/nmfs-ost/ss3-source-code) — control-file reference
  style; exhaustive per-switch documentation with allowed values.
- [**dsem**](https://github.com/James-Thorson-NOAA/dsem) — the DSEM-linked models: formula/path
  grammar, and how linkage structure is specified and reported. DSEM lives on the
  `dsem-v5-integration` branch, not here.

## Known traps

One line each; the evidence and the measured numbers are in `inst/dev/TRAPS.md`.

- **`Index_distribution` has a second hand-synced registry** — a family added to
  `index_distribution_map` must also be classified in `.index_rows_natural_scale()`
  (`R/0-switches.R`), or it silently gets the log-scale residual formula.
- **A `data_list` element with no `write_data()`/`read_data()` support round-trips to nothing** —
  this is how `index_cov` was lost.
- **A `Comp_weights` of 1 under a Dirichlet-multinomial is a starting weight of e** — that
  likelihood reads the column as a log.
- **`fit_mod(estimateMode=)`** takes a string or the integer behind it: `"Estimate"` (0) =
  hindcast + HCR projection, `"Hindcast"` (1) = hindcast only, `"Projection"` (2) =
  projection-only from `inits`, `"DebugBuild"` (3) = build without optimizing,
  `"DebugOptimize"` (4). Prefer the strings.
- **`fit_mod(estimateMode = 4)` returns a placeholder objective** (`dummy*dummy`), because
  `build_map()` maps out every hindcast parameter. Don't read anything into a mode-4 objective,
  gradient, or Hessian. **Mode 3 returns the real objective**, so `obj$fn()` / `obj$gr()` are
  usable for diagnosing a model before fitting it — the analogue of WHAM's
  `fit_wham(do.fit = FALSE)` and SAM's `sam.fit(run = FALSE)`.
- **`fit_control()` bundles the optimizer and uncertainty knobs** — `getsd`, `bias.correct`,
  `loopnum`, `newtonsteps`, `getJointPrecision`, `nlminb_control`, and the bias-adjustment
  flags. `getsd = FALSE` leaves `sdrep` NULL, so `vcov()` returns NULL and uncertainty bands
  are NA. The refit diagnostics forward only `phase` and `getsd`.
- **`run_mse()` pins the OM's stock-recruit and suitability windows to the pristine `om$`**, not
  the advancing `om_use$`, so the hindcast does not drift through the projection — essential for
  multispecies, whose predation suitability must stay fixed.
  `tools/verify/verify-mse-hindcast-invariant.R` checks it.
- **Every observation and process error is drawn in a `SIMULATE{}` block beside the density that
  scores it.** `sim_mod()` implements no observation model in R — it calls `obj$simulate()` once
  and writes the result back, so a new likelihood family owes a draw. Draw what the density
  assumes (bias-correction convention and scale included), REPORT under a `*_sim` name, and
  don't draw what the model does not define. `tools/verify/verify-sim-*.R` is the net.
- **A slow fit is the model, not a regression** — `BS2017SS` has needed ~500–700 `nlminb`
  iterations since at least 2023.
- Scratch outputs (`Rplots.pdf`, `*_osa.png`, `*.RDS` under `tests/comparison/`) are gitignored.
