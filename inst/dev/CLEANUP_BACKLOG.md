# Cleanup backlog

The 79 `TODO` / `FIXME` markers in the source, triaged. **Add to this file; don't fix these
unasked.** Fix the one in the file you were already asked to touch, in the same commit.

**Cite the marker text, not a line number.** These references have gone stale three times --
twice when roxygen was added above them, once from an unrelated `ceattle.cpp` edit. Grep for the
quoted FIXME text instead; it moves with the code.

Three tiers: a **known defect** is a wrong answer waiting for the right input and should become a
GitHub issue; a **design note** is a wish, not a bug; `TODO(review)` is a deliberate convention
marking a judgement call for Grant, and is never resolved by an agent.

Counts by area: `src/TMB/ceattle.cpp` 29 · `src/TMB/predation.hpp` 5 · `src/TMB/Dev/caal.hpp` 5 ·
`R/3-build_map.R` 5 · `R/10-run_mse.R` 5 · `R/9-retro_and_jitter.R` 4 · `R/5-rearrange_data.R` 4 ·
rest singletons. Re-derive with `grep -rn 'TODO|FIXME' R/ src/TMB/`, excluding the `todo <-`
variable in `R/6-process_residuals.R`.

---

## Tier 0 — known defects, with the input that triggers them

These say, in the source, that the code is wrong under a stated condition. Each needs an issue.

| Where | Condition | Consequence |
|---|---|---|
| `src/TMB/ceattle.cpp` (`will bomb if minage > 1`) | `minage > 1` | "will bomb". Same class as the `nages`-is-a-count trap, but in the **template** rather than the plotters — the R-side plotting instances are fixed and fixture-covered; this one is not. Every bundled dataset and all three live assessments have `minage = 1`, so nothing exercises it. |
| `src/TMB/ceattle.cpp` (`will blow up if nlengths is less than nages`) | `nlengths < nages` | "will blow up" in the composition block. |
| `R/10-run_mse.R` (`does not work for assessments that don't occur annually`) | assessment interval ≠ 1 year | the catch fill-in "does not work for assessments that don't occur annually" — i.e. any `assessment_period > 1` MSE. |
| ~~`R/10-mse_summary.R` (`will bug if not survey`)~~ | ~~a species whose fleets are not all surveys~~ | **Resolved in 5.13.0**: not a defect. `spp_rows <- which(flt_spp == sp)` was assigned once and read by nothing, in this file or anywhere else in the package. The FIXME speculated about code that did nothing; both lines are gone. No behaviour change. |
| ~~`R/3-build_map.R` (`QAR1 is inert`)~~ | ~~**`Catchability = "AR1"`** (the QAR1 form, Rogers et al. 2024)~~ | **Resolved in 5.12.0**: `data_check()` now errors on it, so the branch is unreachable. It was inert — the deviate map is gated on `Time_varying_q %in% c("IID","AR1","RandomWalk")`, but under `Catchability = "AR1"` that column holds an `env_data` **column index**, not a mode, so `index_q_dev` stayed mapped out and q was constant. Not repaired: the Rogers form is implemented correctly by a q linkage (`ar1(1 \| Year)` with `observe`), which GOA pollock 2025 runs. Remaining cleanup is to delete the dead branch and code 6 from `q_map`. These are two different switches sharing a string; an earlier draft of this file named the wrong one. |
| `src/TMB/ceattle.cpp` (`caal_ll_type`) | `CAAL_distribution = "MultinomialAFSC"` | Same shape, one severity down: passes `validate_switches()` and errors in the template. The catch-sigma case of this shape was fixed in 5.12.0; this one remains. |
| ~~`R/3-build_map.R` (`will fail if random_sel = TRUE?`)~~ | ~~`random_sel = TRUE` + `Time_varying_sel = "Block"`~~ | **Resolved in 5.13.0**: confirmed real, and worse than "will fail". The block parameters live in `log_sel_slp_dev`/`sel_inf_dev`, and `fit_mod()` declared those arrays random unconditionally — but the template scores selectivity deviates only for `IID`/`AR1`/`RandomWalk`/`RandomWalkAscending`, so blocks were Laplace-integrated against **no density**, `sel_dev_log_sd` mapped out so there was no variance either. Measured: 8 parameters random, `JNLL_SEL_DEV` identically 0, objective `NaN`, real fit dead with TMB's `NA/NaN gradient evaluation`. Now refused with a message naming the fleets and the way out. `test-selectivity-random-sel-block.R` carries the reproduction and a drift guard pinning `Block` as the only mode the template leaves unscored. |
| `R/6-fit_mod.R` (`swallows EVERY warning build_map() raises`) | any | `suppressWarnings()` swallows *every* warning `build_map()` raises, not just the intended one. |
| `R/10-mse_summary.R` (`EM uses fixed-depletion proxy for HCR 2`) | HCR 2 | the EM uses a fixed-depletion proxy (`0.5 * 0.35`) rather than the configured target. |
| `src/TMB/ceattle.cpp` (`dmultinom_osa` renormalization note) | `Comp_distribution` case 0 | fitting routes through `dmultinom_osa()`, which renormalizes `p`. **Reported value only** -- the FIXME states the gradient and MLE are unchanged (an additive constant), so this is a reporting discrepancy, not a wrong fit. |

## Tier 1 — stated limitations, currently by design

Not bugs, but they bound what the model can be asked. Worth documenting in a vignette rather
than fixing.

- **Forecast growth is ignored** by the retrospective and MSE projection paths
  (`R/9-retro_and_jitter.R:235,249`, `R/10-run_mse.R:436`) — the terminal-year growth is carried
  forward.
- **Projection quantities are held at the terminal hindcast year** (`R/10-run_mse.R:460`).
- **`ration_data` is sized for the hindcast only** (`R/5-rearrange_data.R:683`).
- **SPR reference points**: `sex_ratio` is fixed rather than estimated for two-sex models
  (`src/TMB/ceattle.cpp:1542`), and the M used is the terminal-year value
  (`src/TMB/ceattle.cpp:1514`, `:1522`).

## Tier 2 — design notes and refactor wishes

- **Split `R/0-build_srr_and_M.R`** (1,497 lines, 29 functions). It carries the stock-recruit
  constructors, the M1 constructors and the growth constructors in one file, and is the single
  outstanding item from the superseded `accessibility-and-code-review` plan. The `R/` numeric
  prefixes are meaningful, so a split keeps the `0-` prefix (e.g. `0-build_srr.R`,
  `0-build_m1.R`, `0-build_growth.R`).

No user-visible consequence; do them opportunistically.

- `R/10-run_mse.R:510` — extract `run_one_sim()` as a top-level internal helper.
- `R/3-build_map.R:1163` — express the mapping as a formula.
- `R/3-build_map.R:1300`, `:1325` — add checks for survey selectivity / q sigma.
- `R/3-build_map.R:338`, `:388` — sex-varying `M1_rho` (noted as hard to estimate).
- `R/5-rearrange_data.R:404`, `:428`, `:560` — tidier handling of empty comp/CAAL frames and the
  `age_error` data frame.
- `R/2-build_params.R:269` — variance and AR1 parameters.
- `R/0-clean_data.R:174` — possibly redundant now.
- `src/TMB/ceattle.cpp:1216` — hoist out of the iteration loop.
- `src/TMB/ceattle.cpp:2350` — fit the window form directly when a fleet needs it.
- The `logH_*` / `H_4` / `log_gam_*` markers belong to the stubbed Kinzey-Punt predation forms
  (`msmMode` 3–9) and the gamma predator selectivity. They are pinned as stubbed in
  `tests/testthat/test-schema-registries.R`; leave them until that work is picked up.

## `TODO(review)` — Grant's calls, not an agent's

Six, each a judgement about what the right behaviour *is*:

- `R/0-rceattle_class.R:268` — whether `residuals(source = "all")` should include diet, given
  `osa_residuals("all")` does.
- `R/0-rceattle_class.R:420`, `:452` — how held-out rows (`Year <= 0`) carrying a positive
  observation should be treated.
- `R/6-fit_mod.R:780` — what a user-supplied `NA` bias-adjustment should mean.
- `R/7-plot_osa.R:58` — how process-residual objects should be plotted.
- `src/TMB/growth.hpp` — carries one; see the file.

## Deliberately not changed

- **Non-parametric growth** is declared and calls `error("not yet implemented")`.
- **The `msmMode` 3–9 Kinzey-Punt branches are not declared at all** -- the whole block in
  `predation.hpp` is inside a `/* ... */`, so there is no dispatch, live or erroring. The live
  modes are handled by `if (msmMode == ...)` in `ceattle.cpp`.
  `test-schema-cpp-dispatch.R` pins both, and pins the absence of the switch.
- ~~`flt_sel_ind`~~ — removed in 5.12.0. It was computed from `Fleet_code` on every fit and read
  by nothing.
