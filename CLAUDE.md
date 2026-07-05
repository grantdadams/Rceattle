# Rceattle — Claude Code guide

Rceattle (v4.6.0, GPL) fits the **CEATTLE** single- and multi-species, climate-linked,
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
- **`src/TMB/`** — `ceattle_v01_11.cpp` is the main model (~3,574 lines, numbered section
  index); process logic lives in headers (`recruitment.hpp`, `selectivity.hpp`,
  `predation.hpp`, `growth.hpp`, `linkage.hpp`, `comp_osa.hpp`, `helper_functions.hpp`,
  `bioenergetics.hpp`, `diet_data.hpp`).
- **`tests/testthat/`** — organized by process (`tests-Dynamics/`, `tests-Selectivity/`,
  `tests-Likelihoods/`, `tests-Data-processing/`, …) with shared `helpers-*.R` /
  `fixtures/`. Fast fixtures: `make_test_data()` (single-species) or `make_msm_test_data()`
  (multispecies, incl. diet) with `estimateMode = 3` build a non-optimized object.
  `tests/comparison/` holds WHAM cross-checks (not part of `test_check`).
- **`vignettes/`** are `eval = FALSE` (they only need to render). `data/` has bundled
  example datasets (`BS2017SS`, `BS2017MS`, `GOA2018SS`, …).

## Conventions & traps

- **Commits: plain messages, no AI-attribution / `Co-Authored-By` trailer.**
- Roxygen uses markdown; run `devtools::document()` after touching `@`-docs.
- **Released package — preserve the public API.** Rceattle is shipped (v4.6.0, has users;
  see `cran-comments.md` / `inst/RELEASE-CHECKLIST.md`). Exported `build_*()` args carry
  deprecation paths (e.g. SRR codes `1/3/5`) — **deprecate / keep back-compat rather than
  deleting public surface.** Internal refactors are free as long as golden-reference
  equivalence holds.
- **Match the surrounding style.** Canonical references: `src/TMB/recruitment.hpp` (the
  Doxygen-documented header to emulate) and any `R/*.R` + its `tests/testthat/*` pair. The
  codebase favors explanatory section headers and Doxygen on the C++ — match local comment
  density; don't strip comments.
- **`jnll_comp` likelihood rows are magic integers** in `ceattle_v01_11.cpp`; their names
  live separately in `R/6-rename_output.R` (~L130–151) and are kept in sync by hand.
  If you add/reorder a likelihood component, update both.
- **Numeric changes need golden-reference equivalence:** any edit that can move the fit
  must keep an example fit (e.g. `BS2017SS`) within tolerance and the suite green.
- **`fit_mod(estimateMode=)`:** 0 = hindcast + HCR projection, 1 = hindcast only,
  2 = projection-only from `inits`, 3 = build (`MakeADFun`) without optimizing, 4 = optimize
  with all params mapped out. **Trap: for `estimateMode >= 3` the template returns a
  placeholder objective (`jnll = dummy*dummy`) independent of the real parameters** — so
  `obj$fn()` / gradients are meaningless, and a random-effects model built in mode 3 gives a
  *spurious* NaN / singular Hessian. Read the REPORTed `jnll_comp` for the real likelihood;
  use `estimateMode = 1` to diagnose the actual objective / gradient / random effects.
- **`fit_control()` bundles the optimizer / uncertainty knobs** (`getsd`, `bias.correct`,
  `loopnum`, `newtonsteps`, `getJointPrecision`). Pass `getsd = FALSE` for fast dev/test fits
  (skips `sdreport`) — but then `sdrep` is NULL, so `vcov()` returns NULL and uncertainty
  bands are NA.
- Scratch outputs (`Rplots.pdf`, `*_osa.png`, `*.RDS` under `tests/comparison/`) are
  gitignored — don't commit them.

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

## Converting ADMB models

- **`dev_vector` sum-to-zero can't be replicated in TMB.** ADMB deviation vectors
  (`dev_vector`) carry a built-in *sum-to-zero* constraint with no direct TMB equivalent.
  When porting such a model, **turn off estimation of the first element** of the dev
  vector to recover identifiability — *unless* that dev vector already has an additional
  likelihood penalty (e.g. a normal / random-effects penalty), in which case the penalty
  pins it and all elements can stay estimated. (Turning a parameter off = mapping it to
  `NA` in `R/3-build_map.R`.)

## Active context

A multi-PR accessibility / code-review refactor is **planned but paused** (branch
`accessibility-and-code-review`). The self-contained plan and locked decisions are in
`~/Downloads/HANDOFF-accessibility-refactor-implementation.md`. Read it before resuming
that work; do not start editing from scratch.
