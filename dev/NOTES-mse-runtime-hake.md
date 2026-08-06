# MSE runtime — measured attribution on the Pacific hake DM model

Companion to `PROMPT-mse-runtime-optimization.md` (the brief) and `PERF-findings.md`
(the BS/GOA baseline). This file records what was measured on a *real* multispecies
assessment, so the next session does not have to re-derive it.

Model: `../Rceattle-models/Pacific hake/Analysis_FinalJUN_DM_dev.R`.
`styr` 1980, `endyr` 2019, workbook `projyr` 2100 (the script overrides to 2025),
`nspp` 3 (hake + arrowtooth + sablefish), `nages` 20/35/35, `niter` 3,
Dirichlet-multinomial comps with lognormal priors on the DM weights.
All timings `-O2`, single-threaded, one machine — treat ratios, not absolutes.

## 1. `fit_mod()` cost is the model *build*, not the optimization

`estimateMode = 3` (data prep + one `MakeADFun`, nothing optimized) versus
`obj$fn` / `obj$gr` microbenched with an alternating par point (TMB caches on the
last `par`; re-calling with the same vector reads ~0 ms and is meaningless):

| model | projyr | nfree | build (s) | fn (ms) | gr (ms) |
|---|---:|---:|---:|---:|---:|
| SingleSpecies | 2025 | 311 |  1.2 | ~1 | ~1 |
| SingleSpecies | 2100 | 311 |  2.9 | ~1 | ~1 |
| TypeIIMSVPA   | 2020 | 314 | 12.4 | ~1 | ~1 |
| TypeIIMSVPA   | 2025 | 314 | 13.7 | ~1 | ~1 |
| TypeIIMSVPA   | 2050 | 314 | 20.3 | ~1 | ~1 |
| TypeIIMSVPA   | 2100 | 314 | 32.8 | ~1 | ~1 |

Multispecies build ≈ **1.3 s fixed + 0.26 s per model year**; it also scales with
`niter` (6.8 / 13.9 / 22.4 s at `niter` 1 / 3 / 5).

`Rprof` on one build (projyr 2100, 31.2 s total) — 98.2 % in `.Call`:

| bucket | s |
|---|---:|
| `MakeADFun` total | 24.5 |
| ├─ `getParameterOrder` (one double-mode template eval) | 6.0 |
| └─ `MakeADFunObject` (the AD taping sweep) | 12.9 |
| `obj$report()` (a second double-mode template eval) | 6.1 |
| R-side (`rearrange_data`, `data_check`, `clean_data`, `rename_output`) | <0.5 |

**The R plumbing is not the problem — do not spend effort caching it.** One
*untaped* evaluation of the template costs ~6 s while a *taped* gradient costs
~1 ms, because the predation block allocates and fills six
`(nspp·max_sex)² × max_age² × nyrs` arrays (`ceattle.cpp:586-611`) and
`calculate_msvpa_predation` runs its 7-deep loop over all `nyrs`. TMBad prunes
those nodes out of the final tape; the build still pays to create and discard them.

## 2. Where `run_mse()` spends its time

`nsim = 1, cores = 1, assessment_period = 1`, `projyr = 2025` → 6 assessment years.
Total `run_mse()` = **70.9 s**:

| phase | per assessment year (s) | share |
|---|---:|---:|
| OM refit (`.refit_like`, `estimateMode = 1`) | 8.3-8.5 | 71 % |
| EM refit (`.refit_like`, `estimateMode = 0`) | 1.7-2.0 | 15 % |
| `remove_F()` (once per sim) | ~8.5 | 12 % |
| `sim_mod()` | 0.00-0.02 | 0.1 % |

Two conclusions that close out items from the brief:

- **`sim_mod()` is not worth optimizing.** The brief flagged its per-row
  `[<-.data.frame` loop and the fact that it simulates the whole table to keep one
  year. Measured at 0.01 s/year here — 0.1 % of the loop. Dropped.
- **The OM step is the target.** It builds a full multispecies tape (~8.4 s) to run
  an `nlminb` over only the new years' `log_F` (~10 ms of actual optimization). It
  is ~100 % overhead by construction.

## 3. Verdict on shrinking the OM's `projyr` each assessment year

Directionally right, with two corrections:

1. **`projyr = endyr` is one step too aggressive.** `10-run_mse.R:410` reads
   `om_use$quantities$max_catch_hat` for `assess_yrs[k]` — a year that was still a
   projection year in the *previous* OM fit — to cap the TAC at exploitable
   biomass. Minimum safe horizon: `projyr = endyr + assessment_period`.
2. **OM only, never the EM.** The EM's HCR reads `ref_SB0 = SB0(sp, nyrs - 1)`
   (`ceattle.cpp:1628`), the terminal *projection* year's unfished spawning
   biomass. Shortening the EM's horizon moves the `Ptarget`/`Plimit` denominator
   and therefore the catch advice — an experimental change, not a speedup.

**This was built** — `run_mse()` now refits the operating model on a horizon of
`assess_yrs[k+1]` and restores the full projection afterwards. What it is
actually worth, measured back-to-back on an idle machine rather than
extrapolated:

| model | horizon | measured |
|---|---|---|
| `BS2017MS` (multispecies) | projyr 2040 | **−9 %** |
| `BS2017SS` (single-species) | projyr 2040 | ~0 % (+0.5 %, i.e. noise) |

Single species sees nothing because its tape cost barely varies with year count;
the saving is a multispecies effect and grows with `projyr − endyr`. **The hake
model has never been measured at its own 2100 horizon** — every figure here and
above is a Bering Sea measurement. An earlier draft of this note projected ~33 %
for hake by extrapolation; treat that as unverified, and measure before relying
on it.

Equivalence is gated by `dev/verify-mse-om-horizon.R`: with `simulate_data =
FALSE` no random numbers are drawn inside the assessment loop, so the shortened
run must be bit-identical, and it is for single- and multi-species across three
scenarios (including a two-year assessment step, where the look-ahead has zero
margin). With `simulate_data = TRUE` the operating model carries fewer rows
during the refit, so `sim_mod()` advances the random stream differently and
results move by ~5 % of peak SSB; runs stay reproducible from a given seed.

Implementation notes (all inside `run_mse()`, around
`10-run_mse.R:427-521`): only `rec_dev`, `log_M1_dev`, `log_growth_par_devs` (plus
DSEM `x_tj`) are `projyr`-dimensioned, so those are the blocks to slice and
restore; `rec_dev` carries the per-simulation recruitment realizations drawn once
by `sample_rec()` at `10-run_mse.R:375`, so slicing must be lossless. `clean_data()`
pads and filters `catch_data` to `projyr` (`0-clean_data.R:177`), so a shortened
horizon drops future rows and silently misaligns `dat_fill_ind`
(`10-run_mse.R:387`) against `max_catch_hat` — restore the full-length table on the
returned object. Only the final OM is retained (`sim_list$OM`), and at the last
assessment year `endyr` already equals `projyr`, which is what makes intermediate
truncation safe.

## 4. Ranked remaining levers

1. **Replicate parallelism.** `nsim = 1` disables it outright
   (`10-run_mse.R:353`: `use_parallel <- nsim > 1L && cores > 1L`). At ~50 min per
   replicate for a `projyr = 2100` run this is the difference between hours and
   days. No package change needed — but verify per-replicate RNG streams so an
   N-worker run equals a 1-worker run for a fixed seed.
2. **OM `projyr` truncation** — §3. Built; measured −9 % on a multispecies MSE to
   2040, ~0 % single-species. Grows with the projection length, unmeasured on hake.
3. **Stop building projection-year predation nodes at all.** `estimateMode` gates
   *nothing* in the C++ (it appears only at `ceattle.cpp:4228`), so a hindcast-only
   fit still tapes every projection year of predation and dynamic-B0.
   `calculate_msvpa_suitability` / `calculate_msvpa_predation` are handed `nyrs`,
   not `nyrs_hind` (`ceattle.cpp:1878`, `:1888`). Gating those when no species is
   forecasting removes the projection cost from *every* fit — the EM's hindcast
   tape and `remove_F` included, which §3 cannot touch. It changes REPORTed
   projection quantities for `estimateMode = 1` fits, so it must be opt-in via
   `fit_control()`.

Already spent, per `PERF-findings.md` and confirmed here: `sdreport` is off
throughout the MSE, phasing is off, `loopnum` is 1, and R-side data prep is <1 s.

## 5. Reproducing

The script needs one edit before it runs: it passes `msmMode = "MSVPA"` (lines 116
and 171), which 5.1.0 rejects — the accepted names are `TypeIIMSVPA` /
`TypeIIIMSVPA` (`0-switches.R:273`). Objective values reproduce exactly:
`ss_run_DM` 2139.032, `ms_run_DM` 2142.280, OM 2262.318.
