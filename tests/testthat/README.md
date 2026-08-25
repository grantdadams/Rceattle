# Rceattle test suite

This directory holds the automated [`testthat`](https://testthat.r-lib.org) tests for
Rceattle. It is meant to be **navigable for fisheries scientists**: every test file is named
`test-<theme>-<topic>.R`, so the theme prefix tells you which part of the model it exercises.
Files are flat (no subfolders) because `testthat` only auto-discovers tests directly under
`tests/testthat/`.

## Themes at a glance

| Prefix | What it covers | Main source under test |
|--------|----------------|------------------------|
| `test-data-*`        | Input assembly & validation: `clean_data` → `switch_check` → `data_check`, optional fields, combining data sets | `R/0-clean_data.R`, `R/1-data_check.R`, `R/0-combine_data_sets.R`, `R/0-switches.R` |
| `test-dynamics-*`    | Population dynamics: recruitment, initial age structure, biological reference points (BRPs), spawning biomass per recruit, multi-species fits, recruitment linkage | `R/6-fit_mod.R`, `R/5-rearrange_data.R`, `src/TMB/ceattle.cpp`, `src/TMB/recruitment.hpp`, `src/TMB/spr.hpp` |
| `test-growth-*`      | Growth curves (von Bertalanffy / Richards), CAAL, plus-group continuity, minage = 0, size-based selectivity | `src/TMB/growth.hpp`, `R/0-build_growth.R` |
| `test-likelihood-*`  | Observation likelihoods: index, composition (multinomial / Dirichlet-multinomial), diet, expected composition, OSA residuals | `src/TMB/comp_osa.hpp`, `src/TMB/diet_data.hpp`, `R/6-osa_residuals.R`, `jnll_comp` rows |
| `test-linkage-*`     | Covariate/environmental linkage: formula parsing, encoding, linkage table, priors & bounds, pooling, per-sex / per-species formulas | `R/0-build_linkage.R`, `R/0-linkage_encode.R`, `R/0-linkage_table.R`, `R/0-priors.R`, `src/TMB/linkage.hpp` |
| `test-mortality-*`   | Natural mortality: `M1` input/output, M linkage, time-varying M | `R/0-build_m1.R`, `src/TMB/predation.hpp` |
| `test-selectivity-*` | Selectivity curves: logistic, descending/double logistic, double-normal, non-parametric, catchability, normalization | `src/TMB/selectivity.hpp`, `R/3-build_map.R` |
| `test-functions-*`   | User-facing wrappers: `as.data.frame`, MSE, retrospective, parameter profiling, optional data | `R/9-retro_and_jitter.R`, `R/10-run_mse.R`, `R/0-rceattle_class.R` |

Unprefixed top-level files are whole-package system tests: `test-convergence.R`
(`R/0-convergence.R` diagnostics), `test-parameter-recovery.R` (end-to-end estimation), and
`test-tmb-makeadfun_smoke.R` (TMB object builds).

List the files in a theme with, e.g., `ls tests/testthat/test-growth-*`.

## Running the tests

```r
# Full suite (run the previously-skipped integration tests too -- see policy below)
NOT_CRAN=true Rscript -e 'devtools::test()'

# One theme (filter matches the part of the file name after "test-")
devtools::test(filter = "growth")        # all test-growth-*.R
devtools::test(filter = "selectivity")   # all test-selectivity-*.R

# A single file
testthat::test_file("tests/testthat/test-growth-vonbert-year-invariance.R")
```

Editing `src/TMB/*.cpp` / `*.hpp`? Recompile first: `pkgload::load_all(".")` (or `/recompile`).
On macOS prefix compile/check commands with `export PATH=/usr/bin:$PATH` so the system
toolchain is used (avoids a Homebrew clang/gfortran shadowing the TMB build).

## Skip policy (unit vs. integration)

Tests fall into two classes, and the class decides whether they run under a plain
(CRAN-mode) `R CMD check`:

- **Unit tests** — pure R, no model fit. They run **everywhere**, including CRAN and plain
  `R CMD check`. Most `test-linkage-*`, `test-convergence`, and `test-functions-as-data-frame`
  are unit tests.
- **Integration tests** — call `Rceattle::fit_mod()`, which compiles/builds the TMB model.
  These carry a file-top `testthat::skip_on_cran()` so a plain `R CMD check` stays fast and
  CRAN-safe. They **do** run when `NOT_CRAN=true` (the local command above and the coverage
  workflow), so nothing is silently excluded from coverage.

When you add a test that calls `fit_mod()`, put `testthat::skip_on_cran()` at the top of the
file (after the header comment). When you add a pure-R unit test, do not skip it.

## Test-file conventions

Each file follows a lightweight, FIMS-style template so tests are predictable to read:

```r
# <One-line description of the module/behaviour under test.>
# <Optional: the math or invariant being checked.>
testthat::skip_on_cran()          # only if the file fits a model

# --- input / output -----------------------------------------------------------
testthat::test_that("<function> returns the expected shape/values", { ... })

# --- edge cases ---------------------------------------------------------------
testthat::test_that("<function> handles <boundary>", { ... })

# --- errors & validation ------------------------------------------------------
testthat::test_that("<function> errors on <bad input>", { ... })
```

Inside each `test_that()`, guard optional dependencies with
`testthat::skip_if_not_installed("TMB")` / `("Rceattle")`. Apply the section banners to new
files; retrofit existing files opportunistically rather than churning them all at once.

## Helpers (auto-loaded)

`testthat` automatically sources files matching `^helper`, so tests call these unqualified:

- `make_test_data()` — single-species data fixture (`helpers.R`). Fast, deterministic;
  supports `minage`, parametric vs. empirical `growth`, etc.
- `make_msm_test_data()` — multi-species **operating-model simulator** (`helpers-make-msm-data.R`).
  Returns true `model_quantities`, an Rceattle `data_list`, and `base_data`, so tests can check
  the fitted model against known truth.
- `calc_multinom_nll()`, `calc_dirmultinom_nll()`, `calc_nll_ar1_1d()`, `calc_nll_ar1_2d()` —
  independent R reference implementations of likelihood components, used to cross-check the
  TMB `jnll_comp` math.
- `expect_all_true(x)` — custom expectation: passes only if every element of `x` is `TRUE`
  (NA counts as failure). Used by the TMB map / quantity checks.

## Snapshot tests

Snapshots live in `tests/testthat/_snaps/<testfile>.md` (or `.svg`). After an intentional
change, review and accept with:

```r
testthat::snapshot_review()   # interactive diff
testthat::snapshot_accept()   # accept all
```

Numeric/output snapshots should be captured from a **known-good (golden) fit** so they encode
correct behaviour, not whatever the current run happens to produce. See the `/golden-check`
workflow for full-fit equivalence.

## Coverage

```r
NOT_CRAN=true Rscript -e 'covr::report(covr::package_coverage())'
```

Coverage is reported automatically on pull requests via the GitHub Actions workflow (see
`.github/workflows/`). Because the integration tests run under `NOT_CRAN=true` there, the
coverage number reflects the full suite.

## Not part of `R CMD check`

`tests/comparison/` holds WHAM and CEATTLE-classic cross-checks (`WHAM-*-comparison.R`,
`test-match-orginal.R`). These are run by hand for validation against external models and are
**not** discovered by `test_check()`.
