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
- **Tests** run with `NOT_CRAN=true` and the shared helpers sourced into an env `e`
  (`for (f in list.files("tests/testthat","^helper",full.names=TRUE)) sys.source(f, e)`),
  then `test_file(f, env = e)`. `options(testthat.max_fails = Inf)` shows all failures.
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
  `fixtures/`. `tests/comparison/` holds WHAM cross-checks (not part of `test_check`).
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
