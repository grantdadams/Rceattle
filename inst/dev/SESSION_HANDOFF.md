# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**In flight: branch `reporting-tables` (worktree `/private/tmp/rce-report`, off `dev` at
93ee38b7), version 5.24.0.** Kalei Shotwell (AFSC) asked for two things: a cheat sheet for
output names, and one call producing every table she reports, so BSAI ATF next round does not
mean re-deriving it. Both are built and the suite is green
(`FAIL 0 | SKIP 3 | WARN 229`, `NOT_CRAN=true`, `TESTTHAT_PARALLEL=false`, 2026-08-29) —
the warning count is unchanged from 5.22.0, so nothing new was introduced.

- **`quantity_dictionary()`** — all 99 `fit$quantities` names with meaning, units, dims,
  `se`, and `standard_label`. `test-schema-quantity-dictionary.R` asserts it against the
  template and a real fit, so a quantity added or renamed in `ceattle.cpp` fails a test.
- **`report_tables()`** — `model`, `likelihood`, `timeseries`, `reference_points`, `fits`,
  plus `retrospective`/`jitter`/`osa` when passed. Never refits.
- **`standard_output()`** — the NOAA 33-column standard for `stockplotr` / `asar`.

**Evaluated against the 2026 GOA three-species assessment** (`../GOA-multispecies-assessment`,
`models/GOA_26_mod_list.RData` + `retrospectives.RData` + `jitter_ms.rds` + `osa_ms.rds`).
Those cached objects are the cheapest real exercise in the ecosystem: 3 species, 16 fleets, a
two-sex stock, a live `sdreport`, 1977-2025 hindcast with projection to 2050, and all three
diagnostic objects already computed. No refit needed. The dictionary covers exactly the 99
quantities both its fits report.

**Five defects the real data found that the fixtures could not.** All fixed, all with
regression tests:

1. `Ftarget` / `Flimit` are unestimated under most HCRs and sit at `exp(0) = 1`, so the
   assessment reported a target F of **1.0/yr**. Gate on `build_hcr_map()` -- NOT on the fit's
   own `map`, which `build_map()` sets to NA in the hindcast whatever the HCR.
2. Under `msmMode > 0`, `SB0`/`B0` are the `MSSB0` 999 mt placeholder, so `B_target` came out
   **399.6 mt** against a true 1e5-1e6 mt. `data_list$MSSB0_derived` is the flag.
3. `standard_output()` dropped every FLEET likelihood row (31 of 38) by filtering `unit` as a
   species when it holds fleet names on rows 1-8.
4. A diagnostics list was paired with models positionally and silently -- and `osa_ms.rds` is
   stored as a list of parts (`index`/`catch`/`comp`), exactly the shape that mispairs. Now
   matched by name; unnamed pairs positionally with a message.
5. A species name containing `", "` could not be selected in `standard_output()`.

**One over-correction caught in review, worth not repeating:** blanking the depletions
alongside `SB0` is WRONG. Section 12.1 of the template does not divide by `SB0` under
`HCR = NoFishing & msmMode > 0`; it divides by biomass in the LAST projection year, the
equilibrated unfished reference. That series is valid on exactly the multispecies fits the
package exists for.

**Model-level findings NOT fixed here** (they are the model's, not the reporting layer's):

- ~~`SPR0` NA for the two-sex species~~ **FIXED in 5.24.1.** SPR is spawning output per TOTAL
  recruit, so the female fraction belongs in it for every species -- once. A one-sex schedule is
  sex-combined, so the age-varying `sex_ratio` applies at every age (`mature_females`, 5.4). A
  two-sex species' sex-0 schedule is already female, so only the recruitment split
  `sex_ratio(sp, 0)` applies. 6.2 applied the age-varying ratio to both, re-applying on a
  two-sex species a split already in its schedule, and going NaN where the table stopped short
  of its own `nages`. Arrowtooth SPR0 NA -> 1.2300187; one-sex species unchanged;
  `ssb`/`biomass`/`R`/`F_spp`/`N_at_age`/`SB0`/`jnll` do not move.
  **Getting this wrong once is easy:** the first attempt dropped the ratio entirely for a
  two-sex species, which reports per FEMALE recruit -- a different quantity from the one the
  other species report. `test-dynamics-brps.R:229` and `test-dynamics-spr.R:168` pin the
  per-total-recruit convention on fixtures where sex_ratio is 0.5 at every age, and caught it.
- `data_check()` relaxed to match (a two-sex species reads `sex_ratio` only at age 1), and
  **`nsex` is now validated** for length and value -- it was read straight into loop bounds, and
  an NA aborted `data_check()` with R's "missing value where TRUE/FALSE needed".
- **Still open, pre-existing:** a NULL or NA `nspp` crashes `data_check()` at
  `if (data_list$nspp != max(data_list$weight$Species))` with the same uninformative abort,
  before any accumulated error is raised. Same family as the `nsex` bug; not fixed.

**VERSION COLLISION -- needs a decision before either branch merges.** `osa-cdf-method` and
`reporting-tables` both branched from `dev` at 5.23.0 and both wrote a **5.24.0** NEWS section:
the former for `osa_residuals(method = "cdf")`, the latter for the reporting tables (plus
5.24.1 for the SPR fix). Whichever merges second must renumber its own section; neither is
tagged, so either order works.
- `ms_mod` has 15 top-level elements and `ss_mod` 17 -- concrete evidence for why stockplotr's
  positional indexing into the fit object cannot work.

**Do not re-derive these** (each cost real time to establish):

- A grep for `REPORT(` over `ceattle.cpp` **over-counts**: `mature_females`,
  `sex_ratio_hat`, `N_at_age_dB0/dBF` and the Kinzey block are commented out. Enumerate from
  a fit. The 99 names are identical single- and multi-species.
- SPR quantities are computed only under `msmMode == 0` and are **exactly zero** on a
  multispecies fit. Reported as NA.
- `data_list$Ftarget` is the SPR percentage (0.40); `quantities$Ftarget` is the F value.
- `residuals()` keys species by integer code; every other section keys by name. That
  mismatch silently emptied the whole `fits` section from `standard_output()`.

**Open, needs Grant:**

- **stockplotr.** `convert_output()` on `nmfs-ost/stockplotr` main already accepts
  `model = "rceattle"` and is broken against current Rceattle: reads `index_ln_q` (ours is
  `index_log_q`) and `catch_h` (ours is `catch_hat`), indexes the fit object positionally,
  and never mentions `nspp`/`spnames` — so a multispecies fit gets `year = NA` on
  biomass/ssb/R. Single-species looks fine, which is what makes it dangerous. An issue was
  offered and NOT sent; nothing has been contacted upstream.
- **SPM is the one required thing missing.** The AFSC guidelines make §4.11.3 Standard
  Harvest Scenarios required for Tiers 1-3. `report_tables()` has no `scenarios` section and
  says so rather than omitting it quietly. Grant already flagged adding SPM for next year.
- Not run: `/golden-check` (the `ceattle.cpp` change is comment-only, so it cannot move a
  fit) and `/ecosystem-sweep` (purely additive — three new exports, no signature changes).

**5.21.0 is released** — tag `5.21.0` on `f9595235`, published 2026-08-26. It carried
`profile_components()`, `plot_profile()`, `plot.Rceattle_profile()` and `profile()`'s `$alias`.

**In flight: `dev` → `main` for 5.22.0, the profile follow-ups Grant asked for.**

- **`relative = "scaled"`** on `profile_components()` / `plot_profile()`. One component routinely
  dwarfs the rest — on a BS2017SS M profile the bottom-trawl composition moves 60.5 objective
  units against the bottom-trawl index's 0.697, a factor of 87 — so the others flatten onto the
  axis and where they prefer the parameter cannot be read. Each series is divided by its own
  change, putting every curve on 0 to 1. It **discards magnitude**, so `minfraction` matters more
  with it; `minfraction` and the legend order are computed on the RAW change, before rescaling.
- **`profile(joint =)`** — `"multiply"`, `"add"`, `"value"` move every cell in `slots` on ONE
  grid. Slots otherwise cross, so a ten-age M schedule over 13 values was `13^10` fits. This is
  what makes an empirical age-based M profilable the way it is normally reported.
- **`profile(param = "q")`** — base catchability, `slots` takes a fleet name. A shared
  `Catchability_index` group is moved together.

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
