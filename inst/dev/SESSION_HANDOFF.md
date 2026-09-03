# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**PR #139 is OPEN into `main`** — https://github.com/grantdadams/Rceattle/pull/139
(`dev` -> `main`, releasing **5.26.0 + 5.27.0**; `main` is at 5.25.1). Reviewed 2026-09-02
against CLAUDE.md for accuracy, ease of use, concise language and sibling breakage. Two P0
defects were found, both measured, both fixed on `dev`; the fixes are described below and are
NOT yet committed.

**One root cause behind both.** Under the default `estimateMode = "Estimate"` with an
estimating HCR, `fit$obj` and `fit$sdrep` belong to the **projection**, not the hindcast —
`build_hcr_map()` maps every hindcast parameter off. Measured on `Atka2022` with `HCR = 5`:
`fit$obj$par` is length **2** (`log_Flimit`, `log_Ftarget`) against `.conv_hindcast$par`'s
**584**. Every test the PR added fits at `estimateMode = "Hindcast"` or `3`, so neither new
feature was exercised on the mode every real assessment uses. Full account in `TRAPS.md`,
section "`fit$obj` and `fit$sdrep` are the projection's under an estimating HCR".

1. **`parameter_index()` labelled the wrong parameters.** `fit$identified` and
   `fit$.conv_hindcast` are captured before the projection, so the convergence battery matched
   hindcast positions against a projection index: a flagged `rec_pars` printed as
   `log_Flimit  element 1`. Fixed by storing the hindcast index on `.conv_hindcast$index` (built
   in `.capture_opt_convergence()`, the last point holding the hindcast `obj`) and by
   `.conv_index_ok()` / `.conv_index_for()`, which verify an index names the vector being
   labelled before any coordinate is printed. A mismatch now prints block names alone.
2. **`plot_selectivity(add_ci = TRUE)` drew a zero-width band.** The projection `sdreport`
   holds the selectivity parameters fixed, so all 1,012 `log_sel_at_age` errors come back
   exactly 0 — a hairline reading as certainty, worse than the bare line it replaced, and the
   existing decline never fired because the errors were *present*. `fit_mod()` now warns at the
   point the option is set, and the plotter declines an all-zero error naming
   `estimateMode = "Hindcast"` / `projection_uncertainty = TRUE`.

**Also fixed in the same pass:**

- **`Bin_first_selected` is a 1-based bin ordinal, not an age** (`Sel_norm_bin` beside it *is*
  an age — opposite conventions, and the `Age_first_selected` alias argues for the wrong one).
  The PR's prose read it as an age in four places. Schema description rewritten; docs now say
  "first selected bin". `minage = 1` hides this everywhere it is currently written.
- `.rce_sel_lead()` re-derived the selectivity grouping instead of reading `data_list`'s own
  `flt_sel_lead` / `flt_sel_type`. Now reads them; its test builds the vector with the real
  `.group_lead()` so it proves parity rather than restating it.
- `test-selectivity-standard-errors.R`'s back-compat test asserted through `fit_mod()`, which
  sets `adreport_sel` itself, so it passed either way. Now asserts on `build_osa_data()`.
- Prose: the log-vs-logit argument was written out five times. The 28-line `ceattle.cpp` block
  is now 12 and points at `?fit_control`; `plot_selectivity()`'s roxygen essay moved to
  `model-options-and-functionality.Rmd`, which gained a runnable `estimateMode` example.
- `add_ci` moved after `colour_by` in `plot_selectivity()`'s signature (it had been inserted
  before it, breaking positional `colour_by` callers — no sibling does this today).
- `.rce_sel_se()` now attaches a `(Fleet, Sex, Year)` row lookup, so a 16-fleet model does not
  rescan a 30,000-row frame once per panel and year.

**An adversarial review of these fixes found a HIGH-severity regression I had introduced,
and it is fixed.** `.rce_sel_lead()` had been "simplified" to read `data_list$flt_sel_lead` /
`flt_sel_type` — fields that do NOT exist on a fitted object, because `fit$data_list` is the
pre-`rearrange_data()` list. The helper early-returned on every real fit, so mirrored fleets
drew no confidence band beside a banded lead. Measured: `GOApollock`'s
`Pollock_survey_1_shelikof_acoustic_BS` resolved to lead `7` instead of `1, 7`, and a mirrored
`Atka2022` went from a finite band at ages 1-3 to `NULL`, leaving 506 of 1012 ribbon rows
finite. `GOApollock`, `GOA2018SS`, `GOAatf2023` and `GeorgesBank3spp` all have real selectivity
groups. **The rewritten unit test passed throughout**, because it asserted against a hand-built
list carrying the very field under test — the "guards are not themselves guarded" shape, now in
`TRAPS.md` as its own section. Now resolved with `rearrange_data()`'s own `.group_lead()` on
`fleet_control`, and the test asserts on a fitted object: lead `1 2`, band restored,
**1012/1012** ribbon rows finite.

Two smaller findings from the same review, also fixed:

- The all-zero-sd decline could not tell "all exactly 0" (HCR projection) from "all NaN"
  (non-PD Hessian), and pointed an assessor with a genuine identifiability problem at two
  remedies that cannot help. Split into two messages; the NaN one names
  `fit$convergence`'s `sdreport_failed` / `hessian_conditioning`.
- The docs scoped the `fit$obj` caveat to "an estimating HCR". Measured: `ConstantF` also
  leaves `obj$par` at length 1 against a 584-row `cov.fixed`. The accurate scope is "any HCR
  but `NoFishing`" for `obj`, and "an estimating HCR" only for `sdrep`. CLAUDE.md rule 10 was
  also left contradicting the new `Bin_first_selected` trap, and now names both conventions.

**A non-finding worth recording:** the memory cost of storing `.conv_hindcast$index` on every
fit was measured, not argued — 0.07%-1.5% of a fit (31-109 KB against 3.7-156 MB), ~1 MB over
ten retro peels. No action.

**Verified:** `FAIL 0 | PASS 439` across `plot-selectivity-ci`, `selectivity-standard-errors`,
`schema-parameter-index`, `convergence`, `schema-canonical`, `vignette-api`; `FAIL 0 |
PASS 4147` on the wider `plot|selectivity|schema|convergence|vignette-api|fit-verbose` net;
**golden regression `FAIL 0 | PASS 20`**. `src/TMB/ceattle.cpp` is comment-only (0 non-comment
lines changed). `devtools::document()` touches only `fit_control.Rd`, `parameter_index.Rd`,
`plot_selectivity.Rd`. `inst/extdata/meta_data_names.xlsx` regenerated via
`data-raw/regenerate_meta_xlsx.R` after the schema edit, as its drift guard requires.

**Still to run before this PR merges:** `/pkgdown-check`, and the release checklist in the PR
body (`build_vignettes`, `urlchecker`, `spell_check`, tag the merge commit `5.27.0` bare,
publish a GitHub Release, `deep-checks` on `main`, ecosystem sweep).

**Sibling repos: no breakage found.** `../Rceattle-models/Model comparison/run_rceattle.R`'s
`get_condition_number()` already reads `covariance_condition_number` with a pre-5.26.0
fallback. All 20 `plot_selectivity()` call sites across the siblings are named or
single-argument. `build_osa_data()` defaults `adreport_sel`, so an old saved `data_list` still
builds.

## Release 5.21.0 — two things that did not work

- **The `release: published` event fired no pkgdown run.** `main`'s workflow has
  `release: types: [published]` and 5.20.0 has such a run; 5.21.0 got none. Dispatched manually
  (`gh workflow run pkgdown.yaml --ref main`), which started immediately while push/PR runs sat
  queued 45-57 minutes. Actions was degraded that afternoon: a 0-second `startup_failure`, and a
  run marked `failure` whose job never left `queued`. **Check that a release actually rebuilt
  the site rather than assuming the event fired.**
- **`deep-checks` reports success while its `safebounds` job fails** — it is
  `continue-on-error: true` ("read it, do not gate on it"). See Blocked.

## Done & verified on 5.22.0

- **Suite green**: `FAIL 0 | WARN 229 | SKIP 3 | PASS 7934` (`NOT_CRAN=true`,
  `TESTTHAT_PARALLEL=false`, 2026-08-26), `test-golden-regression.R` running.
- **Before/after `ggplot_build()` diff against the released 5.21.0 source** (rule 3's net for a
  plotting change): layers, labels and series order identical on `relative = "own"`,
  `"minimum"` and `"none"`. Existing figures do not move.
- **Ecosystem sweep for `profile(`**: every caller in `Rceattle-models`, `GOA-ATF-ESP`, the
  vignettes, `examples/` and `tools/verify/verify-refit-like.R` uses NAMED arguments, so
  appending `joint` breaks nothing. It is appended after `getsd` regardless.
- **Four defects caught in review, all in this session's own work.** The one that mattered:
  `profile(param = "q")` on a fleet with `Catchability = "Analytical"` (3) or `"AnalyticalArith"`
  (7) would have come back FLAT. Those forms solve q from the fleet's own index and overwrite
  `index_q` (`ceattle.cpp:2526`), so `index_log_q` never reaches the likelihood — every grid
  point returns the same fit, reading as "the data do not inform q". Now an error naming the
  fleet, and such fleets are excluded from group expansion (the schema says a group containing
  one does not share a catchability). `est_index_q` is **not carried on a fitted object**, so the
  raw column is resolved with `.map_switch(fc$Catchability, q_map, ...)` rather than a second
  reading of `q_map`.

## Done & verified on 5.21.0

- **Suite green on the recompiled DLL**: `FAIL 0 | WARN 229 | SKIP 3 | PASS 7893`
  (`NOT_CRAN=true`, `TESTTHAT_PARALLEL=false`, 2026-08-26). `test-golden-regression.R` ran
  (real `fit_mod()` at `:47` and `:50`) and passed, so the four reference objectives are
  unmoved.
- **No `/golden-check` run separately, and it is not needed.** The only template change moves a
  term between `jnll_comp` rows; the objective is `jnll_comp.sum()` (`ceattle.cpp:5308`), so it
  is bit-identical by construction. The golden regression in the suite is the check.
- **`/pkgdown-check` clean** (exit 0), `_pkgdown.yml` carries the three new topics.
- **New guard**: `test-schema-jnll-rows.R` parses every `jnll_comp(JNLL_*, col)` write in the
  template — 124 writes, all 21 rows covered — and asserts each row's declared axis matches the
  column it is actually indexed by. Confirmed non-vacuous by printing the parsed table.

## Done & verified on #119 (merged)

- **Suite green**: 0 failures, 3 skips (`NOT_CRAN=true`, `TESTTHAT_PARALLEL=false`).
- **Golden bit-identical** on all four reference models — but see the caveat below; it does not
  cover the 5.15.0 changes.
- **Pacific hake MSE run end to end on this branch AND on `main`**, and the hindcast is identical
  across the whole release: 2133.821 / 2134.471 / 2137.443 / 2260.706, vulnerability ATF→hake
  0.8172141 and SBF→hake 0.7685885. `run_mse()` completed both sims; `mse_summary()` returned all
  four blocks. **The `-1.612` on the estimated-suitability stage is not a regression** — it is
  against a 5.6.1 inline comment in `../Rceattle-models/Pacific hake/04-mse.R`, and `main` returns
  the same 2260.706. That comment is stale; the folder `README.md` only pins the first three.
- **Ecosystem sweep clean.** No exported symbol removed or added between `main` and `dev`;
  `NAMESPACE` gains only three `Rceattle_fit_control` replacement methods and `simulate.Rceattle`.
  Every deprecated spelling still live downstream resolves through a shim: `right_adj` ×154,
  `est_M1` ×40, `srr_pred_fun` ×24, `qFun` ×7, `single.plots` ×7, `rearrange_dat` ×6.
- `urlchecker::url_check()` clean (20 URLs); `_pkgdown.yml` covers every export.

## The 5.15.0 work, and what actually covers it

`MSSB0` is the multispecies unfished SSB. The template discards its own equilibrium `SB0` under
`msmMode > 0` and reads this `DATA_VECTOR` instead, so `ssb_depletion` is `ssb / MSSB0`. No
workbook can supply it: `clean_data()` seeds `.RCE_MSSB0_PLACEHOLDER` (999 mt) and `fit_mod()`
section 10.2 derives the real value by projecting under no fishing.

- **It was written only into `data_list_reorganized`**, so the returned `data_list` kept the 999
  and every refit off a fitted object re-entered the template with it — `.refit_like()`,
  `remove_F()`, every `run_mse()` projection. `SB0` is also the HCR threshold and the `posfun`
  floor on `ssb_depletion - Plimit`, so a multispecies refit compared SSB against 999 mt. Fixed.
- **`mse_summary()` divided by the placeholder regardless.** With `HCR = "NoFishing"` no unfished
  reference is ever derived (section 10 is gated on `HCR != "NoFishing"`), and the hake model
  reported a terminal depletion of 2.68e3 where the dynamic reference gave 0.96. `NA` now.
- **`DynamicHCR` was unusable under multispecies**: the no-HCR arm of the depletion block ran
  unconditionally and overwrote `ssb/DynamicSB0`, pinning terminal depletion at exactly 1.
  Gated on `DynamicHCR == 0` now. **That arm is correct as written** — with `HCR = 0` the
  projection is unfished, so terminal-year SSB *is* multispecies SB0, given enough projection
  years to equilibrate.

**`/golden-check` is silent on all of it.** The golden recipe passes no HCR, so `HCR = 0`: the
reference-point penalties never enter the objective and `SB0` is overwritten by `MSSB0` anyway.
The net is `test-dynamics-unfished-reference.R` (9 assertions, mutation-tested — removing the
persistence line leaves `MSSB0` at 999 where the fit derived 27796.7 and 189.7 mt), a before/after
probe against the pre-change build on BS2017MS, and the hake MSE.

## For Grant's review

1. **`SB0`/`SBF` are still carried on `M_at_age`** under predation, whose `M2` comes from the
   *fished* predator field — an unfished reference built on a fished mortality. Correcting it needs
   an equilibrium `M2` that `calculate_msvpa_predation()` does not compute: it carries three
   variants (fished, dB0, dBF), and an equilibrium pair means 6 new arrays and 6 new solver
   parameters inside the `iter` fixed-point loop. Measured on BS2017MS at `estimateMode = 4`:
   switching `NByageF` to `M_at_age_dBF` alone moved `SBF` by 5.27; `NByage0` to `M_at_age_dB0`
   moved summed `NByage0` by 16808.7. Marked `TODO(review)` at the site, deferred as too invasive
   for a release PR. **Settle the definition first**: `NByage0` is not a strict equilibrium — mean
   recruitment and terminal-year weights, but year-specific M1 and R0 — so "equilibrium M2" here
   means dB0 without recruitment deviations.
2. `DynamicSB0`/`DynamicSBF` now decay to spawning on the M they were propagated on. Only bites
   when `spawn_month != 0` **and** `msmMode > 0`; no bundled dataset carries both, so nothing
   measurable moved. Worth confirming against a real assessment that has a spawning month.
3. **PR #119 grew past its original scope.** It began as a release of already-reviewed commits and
   now also carries a template change and a behaviour change to multispecies refits. Those are
   verified, but they were not part of what #109–#118 were reviewed as.

## Known flags

- **`jnll_comp` columns count fleets on rows 1–8 and species on rows 9–20**, so `rowSums()` pools
  two different axes. `.JNLL_ROW_AXIS` (`R/9-profile.R`) is now a third hand-synced partner to the
  `JnllRow` enum and `R/6-rename_output.R`; all three must move together.
- **`unweighted_jnll_comp` is written for 5 of its 21 rows** — composition, CAAL, stomach and the
  two linkage rows, the ones carrying a `Comp_weights` multiplier. Everything else is
  structurally zero there, not small, so `profile_components(weighted = FALSE)` returns a much
  smaller set of series.
- **The QAR1 block in `ceattle.cpp:4269` is dead code.** `data_check()` refuses
  `Catchability = 6` (`R/1-data_check.R:92`) and `est_index_q` is only ever set from that column;
  the live QAR1 form is a q linkage. Its AR1 density was scoring into the "Catchability prior"
  row and now scores into "Catchability deviates" — a reporting fix that no fit can reach.
  `R/6-process_residuals.R:204` still branches on `est_index_q == 6`, which is why the block was
  fixed rather than removed. **Removing the dead QAR1 path — the C++ block and that R branch
  together — is Grant's call and wants its own change.**
- **A model carrying GOA numbers forward needs a refit.** Result-changing changes on this line not
  labelled breaking: the mode-5 selectivity penalty fix (GOA Pacific cod SSB 2050 −14.1%),
  parameter bounds previously applied to the wrong parameters, `remove_F()` zeroing fitted hindcast
  F when `suit_endyr < endyr`, composition weights warm-starting from `inits`, failed `run_mse()`
  simulations returning only a marker, the `mse_summary()` reshape, the recruitment fixes, and
  `sim_mod()` drawing the index under the fleet's own `Index_distribution`.
- **5.10.0 moved three predation figures' numbers** — `plot_m2_at_age_prop()` (a share, not a
  contribution), `plot_ration()` (× average numbers-at-age) and `plot_b_eaten()` (million mt).
  `plot_selectivity()` also renames `p$data$Age` to `Bin`.
- **MSE projection statistics are not comparable across 5.13.0.** `sim_mod()` draws the index
  under the fleet's own `Index_distribution` now, which shifts the RNG stream. At `nsim = 2` the
  hake summary swung `catch_iav` 0.25 vs 0.74 between branches on identical OM and EM fits.
- The golden reference on `dev` is `ss = 10241.0304272585` (`ms = 10267.2478324443`,
  `goa_ss = 12868.0052289274`, `goa_ms = 12932.7931701136`), pinned in
  `test-golden-regression.R`.

## Blocked

Nothing.

**Cleared 2026-08-27: the bounds check runs.** `verify-safebounds.R` now builds the package
SHLIB as well as the TMB model, asserts both DLLs loaded before any case runs, and reports
"could not run" apart from "violation". All **six** configurations pass locally from a clean tree
(all three `.so` deleted first), which is the CI condition — the sixth is
`test-selectivity-catchability.R`, the file that crashed on Windows, which used to be a second CI
step and was never bounds-checked at all. 206 s end to end.

**`continue-on-error` is off now**, so a bounds failure fails the run. It gates nothing —
`deep-checks` is nightly and dispatch-only, never a PR check, and `main` is not protected — but
the failure is visible, which it was not before.

The history below is kept because it is the reason not to trust a green `deep-checks` run
without reading the job.

<details><summary>What was wrong</summary>

**The bounds check had measured nothing since at least 2026-08-25, and reported green.**
`deep-checks`' `safebounds` job is `continue-on-error: true`, so the workflow's overall
conclusion is `success` while the job fails. It fails at the DLL, not on a violation:

```
confirmed: -DTMB_SAFEBOUNDS in the compile line
Failed to load at least one DLL.
  ERROR: 'ceattle' was not found in the list of loaded DLLs.   x 5 cases
BOUNDS VIOLATIONS or errors in 5 case(s):
```

All five configurations error before running, so **zero cases are bounds-checked**. The runs on
2026-08-25 (5.20.0 dispatch) and the 2026-08-26 nightly failed identically, so it is not new.

Cause: `NAMESPACE` loads two DLLs (`useDynLib(Rceattle, .registration=TRUE); useDynLib(ceattle)`).
`verify-safebounds.R` rebuilds only the TMB lib, then calls `pkgload::load_all(compile = FALSE)`
— deliberately, since a recompile would overwrite its bounds-checked model — so nothing ever
builds `src/Rceattle.so`, and `src/Makevars` has `$(SHLIB): tmblib`. On a fresh CI checkout that
file has never existed. **It passes locally only because a developer's tree already has one**,
which is why "five configurations clean on macOS" in the 5.20.0 notes was not the clearance it
looked like.

Both are fixed. It is the only net for the Windows `0xC0000005`, which is still uncured — see the
5.21.0 release notes' known limitation.

</details>

## Resume here

**TODO, in order.**

1. **Read tonight's `deep-checks`** (nightly, 05:00 UTC). It is the first run in which the bounds
   check measures anything, and the first in which a failure shows red — `continue-on-error` is
   off. Green means the Windows `0xC0000005` finally has a net; red is a real finding. This is
   the only thing in PR #128 not verified on a real runner.
2. **Merge PR #128** (`dev` → `main`, 5.22.0) once its CI and the above are green. Then tag
   `5.22.0` — bare, no `v` prefix — and **publish** a GitHub Release;
   `inst/RELEASE-CHECKLIST.md` §3.
3. **Confirm the site actually rebuilt.** The `release: published` event fired no pkgdown run for
   5.21.0; dispatch `pkgdown.yaml` by hand if it happens again. See the 5.21.0 note above.
4. **Dispatch `deep-checks`** (§3b) after the release.

Three loose ends, none blocking:

1. **`inst/dev/SIBLING-REPOS.md` has an uncommitted edit that predates this session** — the
   `../GOA-multispecies-assessment` entry. Coherent and complete; left out of the 5.21.0 and
   5.22.0 commits because it is unrelated to that work, not because it is wrong.
2. **`../Rceattle-models/GOA pollock/2025/04-fit-and-diagnostics.R` gained an M-at-age-6
   component profile** (uncommitted, separate repo, not yet run). `minage = 1`, `nages = 10`, and
   `M1_base` runs 1.39 at age 1 to 0.30 from age 6 on, so age 6 is bin 6 and 0.30 is the value a
   0.18-0.42 grid brackets. It fixes the age-6 cell **alone**, leaving ages 7-10 at base —
   narrower than goa_pk's own M profile, which scales the whole vector. **5.22.0's
   `joint = "multiply"` now expresses the goa_pk version**; switching it is a two-line change and
   probably what the chapter wants.
3. **The dead QAR1 path.** `ceattle.cpp:4269` was fixed in place rather than removed because
   `R/6-process_residuals.R:204` still branches on `est_index_q == 6`. Removing both together is
   Grant's call.

After that, `inst/dev/CLEANUP_BACKLOG.md` has
64 items left (Tier 0 and Tier 2 are cleared); `inst/dev/BACKLOG-PLAN.md` sequences them by who is
exposed. The equilibrium-`M2` item above is the largest and wants its own PR with the hake MSE as
its check, since golden cannot see it.

Queued and not started: `inst/dev/CONTRIBUTOR-EXPERIENCE.md`. **Item 0 comes first and is a
conversation, not code**: ask the sibling-repo authors where they actually stopped. A–E and H are
doc and tooling; G is additive; F is the only item that can break a caller and needs
`/golden-check` + `/ecosystem-sweep`.

## Older paused work

**The `accessibility-and-code-review` refactor is superseded, not outstanding.** The branch is
gone from local and remote, its plan no longer exists, and no commit in any branch mentions it.
Three of its four locked "chosen extras" landed by other routes — the `JnllRow` enum, the repaired
cpp section index, and Doxygen `@file`/`@brief` blocks on the previously bare headers. Its one
concrete leftover, splitting `R/0-build_srr_and_M.R`, is done in 5.14.0. Nothing of that plan is
outstanding. Do not resume from the branch name or the handoff path; both are dead references.
