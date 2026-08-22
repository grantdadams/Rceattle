---
description: Pick and run the tools/verify harnesses that actually cover the current diff
argument-hint: "[optional: a harness name, to run just that one]"
---

`/golden-check` fits four models and diffs the result. It does **not** reach the diagnostic refit
paths, the simulation draws, or any figure. `tools/verify/` holds the harnesses that do.

**They are not all the same shape.** Four are before/after digests you run twice; nine are
self-contained checks you run once. Running a self-contained check on two commits and comparing
the two "PASS" lines verifies nothing — that is the failure mode this section exists to prevent.

## Before running

Work out which harnesses the diff needs, and say so before running anything.

### Before/after digests — run twice, compare

These take an output path and an optional `compare` argument.

| Harness | Covers |
|---|---|
| `verify-refit-like.R` | the `.refit_like()` refit path, via 6 of its callers |
| `verify-mse-om-horizon.R` | equivalence + timing for a shortened OM projection horizon |
| `verify-mse-repro.R` | fixed-seed multi-year MSE trajectory fingerprint |

```
git worktree add /tmp/rce-base <merge-base>
export PATH=/usr/bin:$PATH
cd /tmp/rce-base && NOT_CRAN=true Rscript tools/verify/<harness>.R /tmp/before.rds
cd -            && NOT_CRAN=true Rscript tools/verify/<harness>.R /tmp/after.rds compare /tmp/before.rds
git worktree remove /tmp/rce-base --force
```

### Self-contained checks — run once, on your branch

These take **no path argument** and assert an invariant or a recovery property directly. Do not
pass them a `.rds` path.

`verify-mse-hindcast-invariant.R` (the MSE must not perturb the hindcast) ·
`verify-sim-centering.R` · `verify-sim-diet.R` · `verify-sim-index-families.R` ·
`verify-sim-linkage-re.R` · `verify-sim-local-optimum.R` · `verify-sim-process-error.R` ·
`verify-sim-recovery.R` · `verify-sim-recovery-process.R`

```
export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript tools/verify/<harness>.R
```

`verify-sim-recovery-M.R` also takes one argument, but it is **`nrep`** (default 8), not a path.
Passing it a filename silently yields `NA_integer_`.

### Which to run

| If the diff touches | Run |
|---|---|
| `R/6-refit_like.R` or a caller | `verify-refit-like.R` (digest) |
| `run_mse()`, the OM/EM loop, projection horizons | `verify-mse-hindcast-invariant.R` (once), `verify-mse-om-horizon.R` + `verify-mse-repro.R` (digests) |
| a `SIMULATE{}` block, `sim_mod()`, `simulate.Rceattle()` | the `verify-sim-*.R` set (once each) |
| a `plot_*()` function | **none of these** — see below |

## What nothing here covers

- **`sample_rec(update_model = TRUE)` and `reweight_comps()`.** Both refit through
  `.refit_like()`; neither is in any harness. A change to that helper needs them checked by hand.
  Say so rather than implying the harnesses were exhaustive.
- **Figures.** `/golden-check` and `test-golden-regression.R` diff the fit, so a plotter can
  change every number it draws and stay green. The net is `tests/testthat/test-plot-*.R` plus a
  before/after `ggplot_build()` diff across the affected calls at the merge base.

## After running

Report per harness: which shape it was, PASS/FAIL, and for a FAIL the section and the largest
deviation. State which harnesses you judged not applicable and why, and name the gaps above.
**Do not report a change as verified on the strength of harnesses that do not reach it.**
