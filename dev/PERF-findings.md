# Rceattle runtime — profile, scaling, and ranked recommendations

**Branch:** `perf-profiling` (off `dev-data-workflow`). **Scope (as directed):** *measure &
recommend only* — no package behavior was changed. This document is the deliverable: the
baseline profile at several model sizes, the scaling analysis, and a ranked lever list for
sign-off. Implementation of any lever is a separate, gated follow-up.

**Priority axis (as directed):** general fit + scaling (the C++ per-objective-eval cost that
underlies every fit and every MSE refit).

## How to reproduce

```
export PATH=/usr/bin:$PATH
NOT_CRAN=true Rscript dev/bench-perf.R dev/bench-perf-baseline.rds     # timings + scaling table
NOT_CRAN=true Rscript dev/bench-perf.R dev/bench-perf-after.rds compare dev/bench-perf-baseline.rds
```

Harness = `dev/bench-perf.R` (read-only; wraps public/internal entry points with
`proc.time()`/`Rprof()` and microbenchmarks `obj$fn`/`obj$gr` via an `estimateMode = 3`
build, which keeps the AD object alive and returns the real objective). Raw timings are in
`dev/bench-perf-baseline.rds`. Model recipes mirror `.claude/commands/golden-check.md`
(multispecies warm-started from the SS MLEs).

## Method and honest limitations

- **Rprof cannot look inside the compiled DLL.** All of the objective, gradient, Hessian, and
  AD-tape work funnels through one `.Call`, so Rprof attributes ~85% of every fit to a single
  "objective" bucket and cannot subdivide optimize-eval vs sdreport-eval vs taping, nor C++
  section vs section. Localizing a hot C++ loop needs a *native* sampling profiler
  (macOS Instruments / `perf` / gperftools) — see the recommended next diagnostic.
- **The per-fit `obj$fn`/`obj$gr` microbench is the clean C++-cost instrument.** TMB caches
  fn/gr on the last `par`; the harness alternates between two nearby par points to force a
  fresh sweep each call (otherwise fn reads ~0 ms).
- **Confounded size axes.** The ladder has only two distinct parameter counts (BS = 310, GOA =
  621) and the SS→MS step changes `niter` *and* adds the predation likelihood. So the
  cross-ladder wall-clock exponent is indicative, not a clean single-variable fit; the two
  clean exponents are gr-vs-niter (5 points, shared data) and per-eval-vs-nfree (BS-SS vs
  GOA-SS, both SS/niter=1).

## Baseline profile (single-threaded, this machine)

### Single fit — `estimateMode = 0` (hindcast + projection), `getsd = TRUE`, phased

| model | wall (s) | nlminb iters | nfree | niter | fn (ms) | gr (ms) |
|---|---:|---:|---:|---:|---:|---:|
| BS-SS  |  55.1 |  27 | 310 | 1 |  7 |  16 |
| BS-MS  | 130.3 | 531 | 310 | 5 | 15 |  37 |
| GOA-SS | 376.1 |  77 | 621 | 1 | 29 |  69 |
| GOA-MS | 904.4 |  39 | 621 | 3 | 49 | 121 |

Rprof buckets put **~85% of each fit's wall-clock in the compiled objective/gradient/Hessian
`.Call`**; R-side data-prep, `MakeADFun` assembly, and post-processing are each <1 s and do
not scale with problem size. **The fit cost is the C++, not the R plumbing.**

### Objective evaluation — `obj$fn` / `obj$gr` (median ms per call)

Per-eval cost and its two scaling axes:

- **vs model size (fleets/params), SS, niter=1:** BS-SS gr = 16 ms (nfree 310) → GOA-SS gr =
  69 ms (nfree 621). 2.0× the parameters → **4.3× the gradient time → exponent ≈ 2.1**. fn
  behaves the same (7 → 29 ms, ≈2.05). **The objective evaluation is roughly quadratic in
  model size.** This is the headline scalability finding.
- **vs `niter` (BS-MS, shared data):** fn = 7/9/11/15/21 ms at niter = 1/2/3/5/8 → **k ≈ 0.53**
  (r² 0.98); gr **k ≈ 0.58** (r² 0.99). Section-6 population dynamics runs `niter`× per eval,
  but it is only a fraction of the whole eval, so the multispecies `niter` premium is
  **sublinear**.

### MSE — serial (`cores = 1`), 3 assessment years, simulate + resample on

| model | wall (s) | per assessment-yr (s) |
|---|---:|---:|
| BS-SS  |  43.7 | 14.6 |
| BS-MS  | 250.3 | 83.4 |

Each assessment year ≈ one OM refit + one EM refit. Rprof shows **no R-side hot spot** (the
per-iteration data reshaping / `sim_mod` / alias-upgrade work is <0.2 s combined); the cost is
entirely the two refits' compiled evaluations. MSE speedups therefore come from **cheaper
fits (the C++ objective) or fewer/parallel fits**, not from optimizing the R loop — and the R
loop is already well-optimized (no phasing, no sdreport, `loopnum = 1`, warm-started from the
previous MLE, replicate parallelism).

### OSA residuals — on a converged fit

| model | build_osa_data (s) | retape .osa_build_obj (s) | comp OSP net (s) | index OSP net (s) | all (s) | n resid |
|---|---:|---:|---:|---:|---:|---:|
| BS-SS | 0.3 |  2.4 | ~21 | ~2 | 25.2 | 4783 |
| BS-MS | 0.1 | 29.1 | ~42 | ~23 | 72.8 | 4783 |

("net" = the per-source `osa_residuals()` wall minus the retape it repeats; the realistic
single `osa_residuals()` call pays the retape once.) Two facts:
1. **Composition `oneStepPredict` is the dominant term** (≈21 s SS, ≈42 s MS) — as users
   report. `build_osa_data` is already cheap and linear.
2. **The `.osa_build_obj` retape is a large fixed cost that grows steeply with the model**
   (2.4 s SS → 29 s MS, because the MS retape carries the predation machinery).

## Ranked levers (re-ranked by the evidence)

Each: mechanism → estimated win → risk → correctness gate. Honest headline: **the MSE loop and
OSA are already optimized on the obvious axes; the remaining wins are in the compiled model.**

### #1 — C++ objective evaluation (the priority; ~85% of all fit time, ~O(n²) in size)

The evaluation dominates every fit and every MSE refit, and it scales ≈quadratically between
BS (7 fleets) and GOA (16 fleets). That super-linearity — not the absolute cost — is the
highest-value target: it is why GOA-SS (376 s) is ~7× BS-SS (55 s) for ~2× the parameters, and
it compounds through the ~500-iteration optimizations and the MSE refit loop.

- **1a. Localize the super-linear loop with a native profiler** *(do this first — it's a
  measurement, zero code risk).* Rprof can't see inside the DLL; a sampling profiler
  (Instruments/`perf`/gperftools) on a BS vs GOA fit will show whether the O(n²) lives in the
  Section-13 likelihood (per-fleet selectivity-penalty loops `ceattle_v01_11.cpp:2754-3004`,
  RE-linkage loops `:3496-3598`), the predation/suitability section 7, or a fleet×obs /
  bin×bin accumulation. **Win:** identifies the real target; **risk:** none (profiling only).
- **1b. Flatten the identified O(n²) hot loop.** Once localized, replace the quadratic
  accumulation with a linear pass (precompute/group once). **Win:** potentially large and
  *scaling* (helps big models most); **risk:** high — C++, federal template, AD graph must be
  identical. **Gate:** `/recompile` → `/golden-check` (4 models ~1e-8) → OSA suites →
  `dev/verify-mse-repro.R` bit-identical, one hot loop per sub-commit.
- **1c. Hoist Section-6 loop-invariants out of the `for(iter<niter)` wrapper.** Real but
  **modest** — the measured `niter` exponent is only ~0.55, so Section-6 is a minority of the
  eval; hoisting helps MS/MSE but will not move SS (niter=1). Same gate as 1b.

### #2 — OSA residuals

- **2a. Composition `oneStepPredict` (dominant term).** Confirm the cheapest valid method is
  used per type (default `oneStepGaussianOffMode` for continuous; generic only for discrete —
  already correct) and that the tail-accumulation fold reaches large folded stocks so fewer
  bins are residualized. **Win:** proportional to comp bins; **risk:** low if method unchanged.
- **2b. Parallel engagement.** `parallel = TRUE` silently falls back to serial whenever
  multiple TMB DLLs are loaded (dev `load_all`, Rceattle + WHAM, full test suite). In a clean
  single-DLL production session it already parallelizes. **Win:** session-dependent, up to
  ~Ncore× on the continuous group; **zero** in a clean install. **Honest:** low priority.
  **Gate:** OSA suite determinism (1e-8) with and without a second DLL loaded.
- **2c. Retape cost (29 s for MS).** Inherent to the current design (a fresh `MakeADFun` at
  `last.par.best` with the obsvec/keep structure). Folding the OSA machinery into the main
  tape would remove it but slow *every* fit — unfavorable. **Recommend: leave as-is.**

### #3 — MSE loop

No R-side hot spot; already optimized (phasing off, sdreport off, `loopnum = 1`, warm-start,
replicate parallelism). Remaining wins are **downstream of #1** (cheaper fits) plus:
- **3a. AD-object / tape reuse across MSE iterations** — largest *potential* single MSE win
  (avoids re-taping ≥2 `MakeADFun` builds per refit), but **highest correctness risk**: each
  iteration appends projection years so array dims change, which generally forces a retape; a
  partially-reused stale tape is a quota-setting catastrophe. **Recommend deferring** to its
  own reviewed effort with a dimension-change adversarial suite. Given how optimized the loop
  already is, the risk/reward is unfavorable unless #1 is exhausted.

### #4 — `getJointPrecision = TRUE` default (interactive `getsd = TRUE` fits only)

Verified off in both projection fits already, so it never touches MSE; it costs only on an
interactive `fit_mod(getsd = TRUE)` hindcast, where it scales with the parameter count and may
be a large share of GOA-SS's 376 s. Flipping the *default* is user-visible (changes the
returned object; needs NEWS + DESCRIPTION bump). **Recommend:** document that
`getJointPrecision = FALSE` roughly halves large-model `sdreport` cost, rather than a silent
default flip. Size it first with the targeted `getsd`/`getJointPrecision` on-off timing
(recommended diagnostic #2 below) — the harness has hooks for it.

### Non-issues (verified — do not spend effort here)

- Reduce MSE `loopnum` — already `loopnum = 1`.
- Warm-start EM instead of phasing — already done (`inits = em_use$estimated_params`, `phase = FALSE`).
- Skip `sdreport` in MSE — already done (`getsd = FALSE` in `.refit_like`).
- `rearrange_data` `sum()` fleet×obs scan — lives inside the MVN `index_ll_type %in% c(1,2)`
  branch; does not execute for BS/GOA (or any non-covariance-index model).
- R-side data-prep / alias-upgrade in the fit path — Rprof shows it is <1 s and does not scale.

## Recommended next diagnostics (measurement, before any implementation)

1. **Native C++ profile of a BS vs GOA fit** (Instruments/`perf`/gperftools) to localize the
   ~O(n²) hot loop — this is the single most valuable next step and unblocks lever #1.
2. **Targeted `getsd` on/off + `getJointPrecision` on/off timing** on BS-SS and GOA-SS
   (`estimateMode = 1`, run isolated so the timing is clean) to size lever #4 — how much of
   GOA-SS's 376 s is `sdreport` vs the joint-precision matrix. Cheap (~4 fits); not run here
   because #4 is low priority, but a quick, decisive input if the default flip is considered.

## Correctness gates captured at the branch point

No package code was changed in this effort (`git diff` touches only `dev/`), so the gates are
green by construction; the runs confirm the environment reproduces the pinned references and
capture the "before" digests that any future optimization must reproduce bit-identically.

- **OSA suites** `test-likelihood-osa-residuals.R` (PASS 107 / FAIL 0) +
  `test-composition-accumulation.R` (PASS 73 / FAIL 0) — **green** (warnings are deprecation
  notices only). The pinned OSA golden jnll (BS2017SS = 1537036.2876…) reproduces, an
  independent confirmation the fit is bit-faithful.
- **`dev/verify-refit-like.R`** → `dev/refit-before.rds` — captured; all 9 refit paths
  (retro/jitter/self_test/profile/remove_F + SS/regen/MS MSE) fingerprinted with no errors.
- **`dev/verify-mse-repro.R`** → `dev/mse-repro-before.rds` (net-new full-MSE digest) —
  captured; bs_ss / bs_ms / goa_ss = 10 fit-nodes each, no errors.
- **`dev/verify-mse-hindcast-invariant.R`** — **green**; OM hindcast SSB max|dev| = 0.0e+00
  for SS and MS under both (simulate,sample) = (F,F) and (T,T).
- **`/golden-check`** (4 models vs `dev/golden-ref.rds`) — **green**; max|dev| = 0.0e+00 on
  every field (par/obj/jnll/ssb/R) for all four models (ss / ms / goa_ss / goa_ms), objectives
  exact to the pinned reference.

These `*-before.rds` digests plus `dev/bench-perf-baseline.rds` are the reference set for the
`--compare` mode of every future optimization commit.
