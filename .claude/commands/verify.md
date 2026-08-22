---
description: Pick and run the tools/verify harnesses that actually cover the current diff
argument-hint: "[optional: a harness name, to run just that one]"
---

`/golden-check` fits four models and diffs the result. It does **not** reach the diagnostic
refit paths, the simulation draws, or any figure. `tools/verify/` holds the harnesses that do.
Each is a before/after harness: capture on the merge base, then require bit-identity.

## Before running

Work out which harnesses the diff actually needs, and say so before running anything:

| If the diff touches | Run |
|---|---|
| `R/6-refit_like.R`, or any of its callers | `verify-refit-like.R` |
| `run_mse()`, the OM/EM loop, projection horizons | `verify-mse-hindcast-invariant.R`, `verify-mse-om-horizon.R`, `verify-mse-repro.R` |
| a `SIMULATE{}` block, `sim_mod()`, `simulate.Rceattle()` | the `verify-sim-*.R` set |
| a `plot_*()` function | **none of these** — see below |

**Not covered by any harness:** `sample_rec(update_model = TRUE)` and `reweight_comps()`. They
also refit through `.refit_like()`, so a change to that helper needs them checked by hand. Say
so rather than implying the harnesses were exhaustive.

**Figures have no harness and no golden net.** `/golden-check` and `test-golden-regression.R`
diff the fit, so a plotter can change every number it draws and stay green. For a plotting
change the net is `tests/testthat/test-plot-*.R` plus a before/after `ggplot_build()` diff
across the affected calls at the merge base.

## Invocation

Capture the baseline from a worktree on the merge base, so the comparison means something:

```
git worktree add /tmp/rce-base <merge-base>
export PATH=/usr/bin:$PATH
cd /tmp/rce-base && NOT_CRAN=true Rscript tools/verify/<harness>.R /tmp/before.rds
cd -            && NOT_CRAN=true Rscript tools/verify/<harness>.R /tmp/after.rds compare /tmp/before.rds
```

Remove the worktree afterwards (`git worktree remove /tmp/rce-base --force`).

## After running

Report per harness: PASS/FAIL, and for a FAIL the section and the largest deviation. State
plainly which harnesses you ran, which you judged not applicable and why, and that
`sample_rec(update_model = TRUE)` / `reweight_comps()` are outside the set. Do not report a
change as verified on the strength of harnesses that do not reach it.
