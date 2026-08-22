# Cleanup backlog

The 74 `TODO` / `FIXME` markers in the source, triaged. **Add to this file; don't fix these
unasked.** Fix the one in the file you were already asked to touch, in the same commit.

Three tiers: a **known defect** is a wrong answer waiting for the right input and should become a
GitHub issue; a **design note** is a wish, not a bug; `TODO(review)` is a deliberate convention
marking a judgement call for Grant, and is never resolved by an agent.

Counts by area: `src/TMB/ceattle.cpp` 29 · `R/3-build_map.R` 7 · `src/TMB/predation.hpp` 5 ·
`R/10-run_mse.R` 5 · `R/9-retro_and_jitter.R` 4 · `R/5-rearrange_data.R` 4 · rest singletons.

---

## Tier 0 — known defects, with the input that triggers them

These say, in the source, that the code is wrong under a stated condition. Each needs an issue.

| Where | Condition | Consequence |
|---|---|---|
| `src/TMB/ceattle.cpp:1943` | `minage > 1` | "will bomb". Same class as the `nages`-is-a-count trap, but in the **template** rather than the plotters — the R-side plotting instances are fixed and fixture-covered; this one is not. Every bundled dataset and all three live assessments have `minage = 1`, so nothing exercises it. |
| `src/TMB/ceattle.cpp:641` | `nlengths < nages` | "will blow up" in the composition block. |
| `R/10-run_mse.R:556` | assessment interval ≠ 1 year | the catch fill-in "does not work for assessments that don't occur annually" — i.e. any `assessment_period > 1` MSE. |
| `R/10-mse_summary.R:458` | a species whose fleets are not all surveys | "will bug if not survey" when selecting rows per species. |
| `R/3-build_map.R:1181` | **`Catchability = "AR1"`** (the QAR1 form, Rogers et al. 2024) | **QAR1 is inert** — DEPRECATED in 5.12.0 with a warning at fit time; still unfixed. The deviate map is gated on `Time_varying_q %in% c("IID","AR1","RandomWalk")`, but under `Catchability = "AR1"` that column holds an `env_data` **column index**, not a mode. It never matches, `index_q_dev` stays mapped out, and q is constant. These are two different switches sharing a string; the source comment says so explicitly, and an earlier draft of this file named the wrong one. |
| `src/TMB/ceattle.cpp` (`caal_ll_type`) | `CAAL_distribution = "MultinomialAFSC"` | Same shape, one severity down: passes `validate_switches()` and errors in the template. | The catch-sigma case of this shape was fixed in 5.12.0; this one remains. |
| `R/3-build_map.R:692` | `random_sel = TRUE` | "will fail" — marked with a question mark, so unconfirmed. |
| `R/6-fit_mod.R:658` | any | `suppressWarnings()` swallows *every* warning `build_map()` raises, not just the intended one. |
| `R/10-mse_summary.R:494` | HCR 2 | the EM uses a fixed-depletion proxy (`0.5 * 0.35`) rather than the configured target. |
| `src/TMB/ceattle.cpp:3401` | `Comp_distribution` case 0 | fitting routes through `dmultinom_osa()`, which renormalizes `p`. **Reported value only** -- the FIXME states the gradient and MLE are unchanged (an additive constant), so this is a reporting discrepancy, not a wrong fit. |

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
