# Handoff — `integrate = FALSE` for linkage random effects (penalized fixed-effect deviations)

## Goal

Let a `rw()` / `ar1()` / `(1|Year)` linkage term be estimated as a **penalized
fixed effect** (deviations in the objective as a plain penalty, *not* integrated
out by the Laplace approximation) instead of a proper random effect — **but only
when its SD is fixed**. This makes the grammar reproduce ADMB/AMAK-style
penalized time-varying deviations (`Time_varying_sel`, legacy `Time_varying_q`,
goa_pk's convention) exactly, which the current always-integrated behaviour
cannot.

Proposed API:

```r
linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05), integrate = FALSE)
```

`integrate = TRUE` (default) keeps today's behaviour.

## Why this is needed (the concrete case)

During the GOA pollock 2025 work we converted the survey-3 (adfg) time-varying q
to `rw(1 | Year)` and it reproduced the legacy fit exactly (free fit: jnll diff
0, SSB 4.6e-6; forward pass: q3 to 5.6e-17). That worked because survey-3 q was
**already** integrated in the legacy config (`random_q = TRUE` puts `index_q_dev`
in the random set).

The **fishery ascending-selectivity random walk** (`Time_varying_sel =
"RandomWalkAscending"`) did **not** convert: the legacy config treats it as a
**penalized fixed effect** (`random_sel = FALSE`, matching goa_pk), whereas the
grammar's `rw()` is a Laplace-**integrated** random effect. Integrating it moved
the free fit to jnll 1934 (marginal, includes the Laplace log-det) vs the legacy
1815, and shifted SSB ~14% — a genuine modelling change, not a reparametrization,
and it breaks exact goa_pk reproduction.

That conversion is deferred with a `TODO` in
`../Rceattle-models/GOA pollock/2025/03-model-comparison.R` and
`04-fit-and-diagnostics.R` pointing at this feature.

## The statistical principle (why the fixed-SD guard is required, not optional)

A variance parameter is only consistently estimated by **integrating** the random
effect (marginal likelihood). Estimating the deviations *and* their SD jointly as
fixed effects is degenerate (drive σ→0 and devs→0). This is the same reason
goa_pk fixes `sigmaR = 1.3` when it treats recruitment deviations as penalized
fixed effects — and we confirmed it empirically this session: with `sigmaR`
**estimated** (i.e. rec devs integrated), goa_pk lands on 1.0158, matching
Rceattle's estimated `R_sd` 1.0157 to four significant figures.

So `integrate = FALSE` **must error unless the group's σ is fixed**
(`init = list(sigma = v)` with no σ prior → `sigma_fixed = TRUE`, decided today at
[R/3-build_map.R:1503](../R/3-build_map.R#L1503)). A penalized fixed effect with
a free SD is not a well-posed model.

## Design — the one real constraint

TMB integrates by **parameter name**, not per-element, and `beta_linkage_re`
([src/TMB/ceattle.cpp:413](../src/TMB/ceattle.cpp#L413)) is a single vector
spanning every RE group. We need a *mix* — e.g. the Shelikof QAR1 integrated and
the fishery selectivity penalized in the same model — so the penalized groups
need their **own parameter vector**:

- `beta_linkage_re`      — integrated groups (goes in `random`).
- `beta_linkage_re_pen`  — penalized groups (a fixed effect; **not** in `random`).

Both carry the same row-20 density; only the first is integrated. This mirrors
how the model already keeps `index_q_dev` (integrated via `random_q`) separate
from the penalized `*_dev` selectivity arrays.

## Implementation

### R side
1. **`linkage_spec()`** ([R/0-build_linkage.R:141](../R/0-build_linkage.R#L141)) —
   add `integrate = TRUE`; store on the spec; document it near the σ/ρ reserved-key
   note ([L119-121](../R/0-build_linkage.R#L119)).
2. **Validation** (spec build / `pool_linkages`) — if `integrate = FALSE` and the
   group's σ is not fixed, stop with a clear message ("penalized fixed-effect
   linkages (`integrate = FALSE`) require a fixed sigma: `init = list(sigma = …)`
   with no sigma prior"). Reuse the `sigma_fixed` logic
   ([R/3-build_map.R:1503](../R/3-build_map.R#L1503)).
3. **Linkage-table schema** ([R/0-linkage_table.R] `LINKAGE_COLS`) — add a per-row
   `re_integrate` flag; set it in `pool_linkages`.
4. **Encoding** ([R/0-linkage_encode.R] `encode_linkage_for_tmb()` + the empty
   case) — emit a per-group (or per-slot) `linkage_re_integrate` DATA vector, and a
   second slot index space so each RE slot points into `beta_linkage_re` **or**
   `beta_linkage_re_pen`.
5. **Params** ([R/2-build_params.R]) — split the RE init vector into
   `beta_linkage_re` (integrated slots) and `beta_linkage_re_pen` (penalized
   slots), each init 0.
6. **Map** ([R/3-build_map.R]) — both estimable; keep the existing rw first-element
   pin and the `log_sigma_linkage` NA-fixing for `sigma_fixed` groups
   ([L1502-1519](../R/3-build_map.R#L1502)) for **both** vectors.
7. **fit_mod** ([R/6-fit_mod.R:643](../R/6-fit_mod.R#L643)) — add **only**
   `beta_linkage_re` to `random_vars` (never `beta_linkage_re_pen`). Ensure
   `build_bounds()` emits bounds for `beta_linkage_re_pen` (it is now a mapped
   non-random estimated parameter — the bounds-assembly loop misaligns L/U if an
   entry is missing).

### C++ side ([src/TMB/ceattle.cpp])
1. `PARAMETER_VECTOR(beta_linkage_re_pen);` beside `beta_linkage_re`
   ([L413](../src/TMB/ceattle.cpp#L413)).
2. `DATA_IVECTOR(linkage_re_integrate);` (per group, 1/0) — or reuse a per-slot
   index that already says which vector a slot lives in.
3. **Effective-beta swap** ([L657-660](../src/TMB/ceattle.cpp#L657)) — route each
   RE row to `beta_linkage_re(idx)` or `beta_linkage_re_pen(idx)` by the flag.
4. **Density block** (row `JNLL_LINKAGE_RE`,
   [L4111-4143](../src/TMB/ceattle.cpp#L4111)) — apply the same rw first-difference
   / AR1 `SCALE` density to whichever vector the group uses. For penalized groups
   the density is a pure penalty; for integrated groups it is the RE prior. Keep
   the `+= SCALE(AR1(...))` sign and the `-= dnorm(re(t)-re(t-1), …)` rw form
   exactly.
5. **Ecov OSA** ([L4160](../src/TMB/ceattle.cpp#L4160)) — the `observe=` path reads
   `re(t)`; make it read from the correct vector. (Likely moot: `observe=` requires
   `ar1`, and an observed AR1 latent should stay integrated — consider simply
   disallowing `integrate = FALSE` with `observe=`.)
6. `random = beta_linkage_re` only (set on the R side).

## Verification (the acceptance bar — this is federal advice)

There is **no golden reference** for a new penalized-density path, so verify
constructively:

- **Golden bit-identical** — no RE ⇒ both vectors empty, row 20 = 0. Run
  `/golden-check`.
- **The pollock fishery-selectivity equivalence** — the decisive test. Convert the
  fishery ascending RW to `slp_asc`/`inf_asc` `rw(1|Year)` with `integrate = FALSE`
  and **fixed** σ (see the exact SDs below); the free fit must reproduce the legacy
  `Time_varying_sel = "RandomWalkAscending"` fit (jnll, SSB, selectivity), the way
  survey-3 q reproduced its legacy fit to ~1e-6. This is the case that *failed*
  (1934 vs 1815) with `integrate = TRUE`.
- **Forward pass** — with the fishery sel on `integrate = FALSE`, the 02-bridge
  forward pass must still reproduce goa_pk (fold each RW's first deviate into its
  base level; the grammar pins the rw first element — see traps).
- **Guard test** — `integrate = FALSE` with an *estimated* σ must error.
- **Density-reviewer refute pass** — obs/process split, the two-vector routing, the
  map/random consistency, and the AD graph. (Use the `density-reviewer` agent.)

## Traps discovered this session (do not rediscover)

- **4× multiplier on the ascending inflection.** The legacy `RandomWalkAscending`
  penalizes the ascending *slope* deviate at `sel_dev_sd` and the ascending
  *inflection* deviate at **`4 * sel_dev_sd`**
  ([src/TMB/ceattle.cpp:3189-3190](../src/TMB/ceattle.cpp#L3189)). For the pollock
  fishery (`sel_dev_sd = 0.05`) that is **`slp_asc` σ = 0.05, `inf_asc` σ = 0.20**.
  The descending limb is not varied in mode 5. Survey-3 q RW has **no** multiplier
  (σ = 0.05, [L3355](../src/TMB/ceattle.cpp#L3355)).
- **rw pins the first deviate to 0** ([R/3-build_map.R:1502-1509](../R/3-build_map.R#L1502)),
  but the legacy leaves it free. In a free fit the base intercept absorbs it (why
  survey-3 q matched to 1e-6). In a **forward pass** you must fold goa_pk's first
  deviate into the base level explicitly — e.g. survey-3:
  `index_log_q[3] <- log_q3_mean + log_q3_dev[1]` and
  `beta_linkage_re[adfg] <- log_q3_dev - log_q3_dev[1]` (see the worked example in
  `../Rceattle-models/GOA pollock/2025/02-bridge.R`). ar1 does **not** pin its first
  element.
- **Multi-spec base composition.** The fishery `slp_asc`/`inf_asc` already carry an
  intercept prior shared across the ascending-limb fleets. Register the RW as a
  *second* spec on the same parameter (a `list()` of `linkage_spec`s), fishery-only;
  its `(Intercept)` row is mapped NA (inert), so the base stays with the `~1` prior
  spec and the RW adds deviations around it. `init = list(sigma = v)` with no σ
  prior is what fixes the SD ([R/0-build_linkage.R:1481](../R/0-build_linkage.R#L1481)).
- **`sigma_index` vs whole-vector overwrite.** When a model has multiple linkage σ
  groups, set only the target group's slot of `log_sigma_linkage` in any inits
  mapping — overwriting the whole vector length-mismatches (hit this in 02-bridge).

## Coordination

This edits the linkage-grammar files (`0-build_linkage.R`, `0-linkage_encode.R`,
`0-linkage_table.R`, `2-build_params.R`, `3-build_map.R`, `6-fit_mod.R`,
`src/TMB/ceattle.cpp`) that another session has been actively committing to on
`mse-debug`. Coordinate / confirm that session is quiet before starting.

## On landing

- `NEWS.md` (new feature), `DESCRIPTION` version bump.
- Convert the pollock fishery selectivity to `integrate = FALSE` and remove the
  `TODO`s in `03-model-comparison.R` / `04-fit-and-diagnostics.R`; re-verify the
  forward pass and free fit.
- Note in the environmental-linkages vignette that `integrate = FALSE` reproduces
  penalized (ADMB/AMAK-style) deviations, permitted only with a fixed σ.
