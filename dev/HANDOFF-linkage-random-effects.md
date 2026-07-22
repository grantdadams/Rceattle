# Handoff: random-effect linkage density (resume point)

Branch **`dev-data-workflow`**, checkpoint commit **`97775cb3`**. Everything below the
"Remaining" line is not yet built. The RE density is the most delicate piece of the linkage
work — Laplace-approximation AD taping — so it was deliberately deferred to a fresh session.
**Read this whole file before writing any C++.**

> **This code sets US federal fishing quotas.** A silently-wrong likelihood term can misset a
> catch limit — it will not announce itself as a crash. Treat statistical accuracy, not
> "it converged," as the acceptance bar. This density has no golden reference (a new RE fit is
> genuinely new), so its correctness rests entirely on the constructive verification below —
> do not skip it, and do not report the feature done until simulation self-consistency and the
> limiting cases pass. Note and test every edge case; an unsupported structure errors loudly.

## Why this piece is different from everything before it

Every earlier linkage increment (catchability, all of selectivity, the RE foundation) was
verifiable by **bit-identity** against a golden reference — the offsets were inert when no
linkage was supplied, so `BS2017SS = 10241.0304274978` and `BS2017MS = 10304.2212737654`
proved correctness. **The RE density has no golden reference**: a new `~ (1|Year)` model is a
genuinely new fit. A subtle taping bug produces a plausible-but-wrong answer that no
bit-identity check catches. So verification here must be *constructive*, not comparative:

1. **`estimateMode = 3`** (build only) then `obj$fn()` / `obj$gr()` — confirm finite, no NaN.
2. **Simulation self-consistency:** simulate deviates from a known sigma, refit, check the
   estimated sigma recovers it within sampling error (this is the real test).
3. **IID-with-sigma-fixed-tiny ≈ fixed-effect** limit: a `(1|Year)` with sigma pinned very
   small should collapse toward the no-deviate fit; sigma pinned large ≈ free deviates.
4. **`jnll_comp` row 20 ("Linkage random effects")** is non-zero and finite when RE present,
   exactly zero when absent (the latter IS golden-checkable — keep it bit-identical).
5. The M1_re machinery is a *working, tested* analogue — cross-check behaviour against a
   `M1_re`-based model with the same structure.

## What is already built and verified (do not redo)

- **`97775cb3`** — IID materialization. `(1 | group)` expands to one indicator column per
  level of `group`, appended after the fixed columns in the design matrix `X`. Each RE row
  carries `re_group` (the grouping variable name) and `re_struct = "us"`. See
  `.materialize_re_design()` in [R/0-build_linkage.R](../R/0-build_linkage.R). The fixed
  intercept of `~ (1|Year)` is a *separate* row that re-targets the base parameter (carries
  the mean); the RE rows are the deviations around it. Verified: `~ (1|Year)` on M1 gives
  1 intercept + N year rows; `rw()`/`ar1()` rejected "not yet wired".
- **`54615b3e`** — the three parameter vectors exist and are length-0/inert:
  `PARAMETER_VECTOR(beta_linkage_re)`, `PARAMETER_VECTOR(log_sigma_linkage)`,
  `PARAMETER_VECTOR(trans_rho_linkage)` at [cpp:383-387](../src/TMB/ceattle_v01_11.cpp#L383).
  `jnll_comp` and `unweighted_jnll_comp` are sized **21** (row 20 = RE density). Row 20 is
  named "Linkage random effects" in [R/6-rename_output.R:151](../R/6-rename_output.R#L151).
  `build_params` creates all three as `numeric(0)`
  ([R/2-build_params.R:205-215](../R/2-build_params.R#L205)).
- **Fit-time guard** in [R/6-fit_mod.R](../R/6-fit_mod.R) (just after
  `.check_q_linkage_support`) errors on any materialized RE row ("not yet consumed by the TMB
  template"). **Removing this guard is the last step**, once the density fits.

## The template: copy M1_re verbatim

[cpp:3388-3435](../src/TMB/ceattle_v01_11.cpp#L3388) already implements exactly this pattern:

```cpp
Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_y, 2)), 0.5);
// ... fill M_re_yr vector ...
jnll_comp(15, sp) += SCALE(AR1(rho_M_y), Sigma_M)(M_re_yr);   // 2D: SCALE(SEPARABLE(AR1,AR1), .)
```

For the linkage version, indexed by **re_group** instead of species:
- `sigma = exp(log_sigma_linkage(group))`
- `rho   = rho_trans(trans_rho_linkage(group))`  (only for ar1; IID this commit ignores it)
- **IID:** `jnll_comp(20, col) -= sum(dnorm(re_slice, Type(0), sigma, true))` — the simplest
  form; do NOT use AR1() for IID (avoid the `1/(1-rho^2)` scaling). Get IID working and
  simulation-verified before touching rw()/ar1().

## Remaining — implement in this order, golden/constructive-check each

1. **Group registry + encoding** (R, in `pool_linkages()` / `encode_linkage_for_tmb`,
   [R/0-linkage_encode.R](../R/0-linkage_encode.R)). After pooling, assign each RE row an
   `re_index` (0-based position in `beta_linkage_re`) and each distinct `re_group` value a
   `sigma_index` (position in `log_sigma_linkage`). Emit DATA vectors:
   `linkage_re_index` (per table row: -1 for fixed rows, else its beta_linkage_re slot),
   `linkage_re_sigma_index` (per RE row: its group's sigma slot),
   `linkage_re_struct` (per group: 0=us/iid). Make groups unique **per
   process|param|species|sex|age_bin|fleet** so different fleets/params get their own sigma
   (mirror the `group_key` in `map_linkage_adjuster`, which already includes `fleet`).
2. **Accumulator split** (cpp `linkage.hpp`). Every accumulator reads `beta(i)` from
   `beta_linkage`. Pass `beta_linkage_re` and `linkage_re_index` in; for a row with
   `re_index >= 0` read `beta_linkage_re(re_index)`, else `beta_linkage(i)`. The offset math
   is otherwise identical (indicator design → adds the deviation in that level's year only).
3. **Density** (cpp, new block before the `jnll = jnll_comp.sum()` at
   [cpp:3870](../src/TMB/ceattle_v01_11.cpp#L3870)). Loop over RE groups; IID `dnorm` sum into
   `jnll_comp(20, col)`. Choose `col` consistently (e.g. group's species, or column 0).
4. **Populate the vectors** (R `build_params`): `beta_linkage_re` length = number of RE rows
   (init 0); `log_sigma_linkage` length = number of RE groups (init from `linkage_spec(init=)`
   or a sensible default like `log(0.3)`); `trans_rho_linkage` length = number of ar1 groups.
   `build_map` maps them (all estimable by default; `init`-fixed → NA).
5. **Variance via the grammar** (Grant's ask): `log_sigma_linkage` init comes from the spec's
   `init` (fixed input value) and its prior from the spec's `priors` — evaluated in the same
   re-targeting prior loop that handles intercepts
   ([cpp:3399-3438](../src/TMB/ceattle_v01_11.cpp#L3399)), adding a
   `RCEATTLE_PROC`-agnostic sigma branch. This is what makes "input a variance" and "prior on
   a variance" fall out of one mechanism. **No prior on the deviate SD exists anywhere in the
   model today** — this is genuinely new.
6. **`random_vars`** ([R/6-fit_mod.R:429-451](../R/6-fit_mod.R#L429)): append
   `"beta_linkage_re"` when any pooled spec is random.
7. **Remove the fit-time guard.** Now RE fits are real.
8. **Then** rw()/ar1() structures, then the legacy translators
   (`Time_varying_sel`/`Time_varying_q`/`M1_re` → grammar, bit-identical) and the
   `Time_varying_*_sd_prior` → `_sd` rename (it is an input value, not a prior — see the plan).

## Local-dev gotcha (cost 30 min this session — do not rediscover it)

The parallel functions (`retrospective`, `jitter`, `run_mse`, ...) have their workers run
`library(Rceattle)`, which loads the **installed** package, not the `load_all` dev version.
So after any change to `jnll_comp`'s row count or the linkage-table schema, the dev parent
produces a fitted object the stale installed workers can't process → **"subscript out of
bounds" in `checkForRemoteErrors`**, only in the parallel path (serial `cores = 1` is fine).
It is a version-mismatch artifact, NOT a code bug — CI installs first, so CI is green.
**Fix: `R CMD INSTALL .` (with `export PATH=/usr/bin:$PATH`) before running the parallel
tests locally.** The RE work changes `jnll_comp` again (already did: 20→21) and the schema,
so reinstall before every local `devtools::test()` that reaches `test-functions-retrospective`
/ `-jitter` / `-mse`.

## Traps specific to this work

- **The numeric-Year / true-lag rule** (from the plan): when rw()/ar1() land, the RE vector's
  index order must be *real elapsed time*, not factor-column order. IID (this commit) is
  order-invariant, so it is safe — but bake the correct ordering in when you build the
  registry so ar1() inherits it.
- **`SEPARABLE`/`AR1` are for correlated structures only.** IID is a plain `dnorm` sum; using
  the AR1 density with rho=0 introduces a `1/(1-0)` = 1 no-op but obscures intent — keep them
  separate.
- **Row 20 must stay exactly 0 when no RE is present** (bit-identity gate still applies to
  every existing model). Guard the density loop on `beta_linkage_re.size() > 0`.
- The plan file (`dev/PLAN-data-workflow-and-linkage-grammar.md`, canonical copy in
  `~/.claude/plans/shiny-hugging-yeti.md`) has the full PR 1 "Random effects" section and the
  variance-prior design under "Variances become linkage parameters (Grant's ask)".

## Reference models in ../Rceattle-models (real validation targets)

These are configured, runnable assessments that use the *existing* time-varying q/sel
mechanisms. They are the constructive-verification cases the density work needs — the new
grammar path must reproduce each **bit-identically** through the translator (step 8), and each
exercises a specific RE structure. Codes: `Time_varying_* = 0` Off, `1` IID, `4` RandomWalk;
`Estimate_q = 6` QAR1 (AR1 catchability).

| Model / file | Config | Structure it validates |
|---|---|---|
| `EBS pollock/Dev/2024 EBS pollock re models.R` | `Time_varying_sel[1] <- 1` + `random_sel = TRUE` (+ `random_rec = TRUE`) on the Fishery | **IID selectivity**, estimated sigma — the cleanest IID target |
| `EBS pollock/scratchpad/random_sel.R` | `random_sel = TRUE` on non-parametric `sel_coff` fishery, `getsd = TRUE`, reports pdHess/maxgrad | IID sel on `coff` (existing mechanism; note the grammar excludes `coff` linkages — this stays on the legacy path) |
| `GOA pollock/2021 pollock build data.R` | `Time_varying_sel = c(rep(0,nindex), 1)`, `Time_varying_sel_sd_prior[FISH] = 0.1`; `Time_varying_q = c(4,0,4,0,0,0,NA)`, `Time_varying_q_sd_prior = c(0.038,NA,0.05,...)` | **IID sel + RandomWalk q** with **input SDs** — the richest single case, and proof the `_sd_prior` column is an input value, not a prior |
| `GOA pollock/2024 pollock bridging.R` | `Time_varying_q[1] <- 1` + `random_q = TRUE` | **IID catchability**, estimated sigma |
| `GOA pollock/2021 pollock build data.R` (`pollock23`) | `Estimate_q[1] <- 6` | **AR1 catchability (QAR1)**, Rogers et al. 2024 — the `est_index_q == 6` path in cpp |

Verification recipe per case: run the legacy config, capture `obj$fn()` / sigma estimate /
`jnll_comp`; run the equivalent new-grammar spec; assert bit-identical. The `0.1` / `0.038` /
`0.05` input SDs become `linkage_spec(init = list(sigma = log(<value>)))` — which is exactly
Grant's "input a variance value" ask, demonstrated against a real model.

## State to resume from
- Branch `dev-data-workflow` @ `97775cb3`, working tree clean, full suite green.
- Version is **4.9.0** (q + selectivity linkages shipped). The RE work bumps to **4.10.0**
  when it lands, with NEWS + the environmental-linkages vignette's "On the roadmap" updated.
