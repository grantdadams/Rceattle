# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**In flight (2026-09-23): exact SS3 -> Rceattle bridge for AI and GOA Pacific cod, branch
`cod-bridge`** (off `dev` at `cf82f27e`; WIP commit `f0796442` pushed, 5.42.0 in DESCRIPTION;
paired `Rceattle-models` commits `3c6e9f0`, `99f40e2` and `3b2fc4b` on master). Goal: from a cold start Rceattle reaches the same solution as an SS3 reference
run. Plan: `../Rceattle-models/SS3-bridge/PLAN.md` (phases 0-5). Targets are SS3 runs
adjusted only in *estimation-method* choices (F_Method 2, `max_bias_adj -1`, F_Ballpark off,
InitEQ lambda 0); biology, selectivity and likelihood are built into Rceattle.

- AI target: `Rceattle-models/AI cod - Dev/Data/M24_1_adjusted` (SS3 3.30.22.1, NLL 531.003,
  gradient 6.3e-5). `M24_1_baseline` reproduces the original M24_1 (535.161 vs 535.174).
- GOA target: `Rceattle-models/GOA cod/Data/goa_pcod-no init and ramp` (NLL 2048.07).
- The superseded `origin/dev-cod-bridge` (798 behind `dev`) is a read-only reference; nothing
  is merged from it.

The release notes below still stand; the cod bridge ships in the same `dev` -> `main` release
or after it, Grant's call.

**`dev` is at 5.41.0**, and once this notes consolidation merges nothing is open against it.
`main` is at 5.33.0, so the next step is one `dev` -> `main` release covering 5.34.0 through
5.41.0, per `inst/RELEASE-CHECKLIST.md`. Read that file's pkgdown note before tagging: the
`release: published` event has silently failed to fire once already.

The 2026-09-14 backlog plan is finished. Eight branches, listed below in version order
(#150 merged before #149), each reviewed
adversarially before commit and again by a second session before merge:

| Version | PR | What landed |
|---|---|---|
| 5.34.0 | #144 | MSE dynamic-SB0 and fixed-numbers reference-point masking |
| 5.35.0 | #145 | Silent-wrong-number fixes; `estDynamics = 3` retired |
| 5.36.0 | #146 | Config overlay by field; OSA outliers flagged per panel |
| 5.37.0 | #147, #148 | QAR1 path removed; stored-map guard; `CONTRIBUTING.md`, the Doxygen build and `adding-a-selectivity-form.Rmd` |
| 5.38.0 | #149 | Per-sex apical selectivity offset (`log_sel_apical`) |
| 5.39.0 | #150 | Multispecies stock-recruit bounds and a degenerate-curve check |
| 5.40.0 | #151 | `NonParametricIID` (13) and `NonParametricRW` (14) |
| 5.41.0 | #152 | `osa_residuals(method = "cdf")` |

## Done & verified (cod bridge, 2026-09-23, WIP on `cod-bridge`)

**Phase 1, growth and biology.** `build_growth(pop_lengths, sd_form = "SD"/"CV",
plus_group_length = "M1"/"none"/"SS3.24"/"decay", plus_group_decay)`; control columns
`L50_mat_len` / `slope_mat_len` (added via `/new-column`; round-trip, template, meta xlsx
regenerated); one `spawn_output[nspp, nages, nyrs]` array now feeds SSB, SB0/SBF, dynamic B0
and SPR (`ceattle.cpp` 5.7). Defaults reproduce the old behaviour.
- `test-growth-ss3-options.R`: all pass (checks against an independent R implementation of
  SS3's ALK / weight / fecundity).
- `/golden-check` via `test-golden-regression.R`: **PASS**, 24 expectations, 0 failures, after
  Phase 1.

**Phase 2, selectivity.**
- `Selectivity = "DoubleNormalSS3"` (15): SS3 pattern 24, own array `sel_dn6[6, flt, sex]`
  on SS3's scales, unnormalized, -999 ends -> `sel_dn6_ends` data flag + parameter mapped out.
  Linkage codes 6-11 (`dn_peak`, `top_logit`/`dn_top`, `ascend_se`/`dn_asc`,
  `descend_se`/`dn_desc`, `start_logit`/`dn_init`, `end_logit`/`dn_final`), in lockstep in
  `R/0-linkage_encode.R` and `linkage.hpp`. Blocks = identity-link `~ cut(Year, breaks)`;
  SS3 annual devs = log-link `(1 | Year)` with `integrate = FALSE` and fixed SD (existing
  grammar, no new dev machinery). `Time_varying_sel` must be "Off" for this form.
- Selected body weight (behaviour change, Grant approved): length-selective fleets with
  estimated growth weigh catch and survey biomass by sum P(l|a) s(l) w(l) / sum P(l|a) s(l)
  (`ceattle.cpp`, right after `calculate_selectivity()`).
- `test-selectivity-double-normal-ss3.R`: 12 pass (curve to 1e-10 vs SS3 formula, -999 ends,
  block linkage). Targeted subset (schema, linkage-encode, selectivity, parameter, growth):
  262 tests, the only failures were `test-schema-cpp-dispatch.R` exemptions, now added (not
  re-run since).
- `test-golden-regression.R` after Phase 2 C++: **PASS**, 24 expectations, 0 failures.

**Parity harness** `Rceattle-models/SS3-bridge/parity_check.R` (G1 forward state at SS3 MLE,
G2 gradient/NLL at SS3 MLE, G3 cold start); converter consolidated to
`Rceattle-models/SS3-bridge/ss3_to_rceattle.R`. Driver: `Rscript SS3-bridge/run_parity.R "AI cod - Dev"`
from `Rceattle-models` (sources the stock's forward pass, then `parity_report()`).
Set `RCEATTLE_PKG` to load a worktree of this branch instead of `../../Rceattle`.

**AI cod, run 2026-09-23 with the converter's ageing error and sample-size factor.**
`test-schema-cpp-dispatch.R` re-run with the new exemptions: 111 pass, and the whole
`schema` + `switches-schema` set is green.

G1, 12 rows at tol 1e-5 (Report.sso print precision), **11 PASS, 1 FAIL**: length-at-age
4.5e-6, weight 2.4e-6, fecundity 2.8e-6, ALK 4.3e-7, sel FshComb 7.3e-7, sel Srv 3.7e-7,
N-at-age 7.5e-6, SSB 4.7e-6, R 3.9e-6, length comp FshComb 1.3e-7, length comp Srv 8.7e-8,
**predicted CAAL 1.2e-1**. The two composition rows are new — see "The CAAL gap" below.

G2 still FAILS: **max |gradient| 547**, unchanged from 532 and still on
`log_growth_pars` (547, 472, 233, 141) and `growth_log_sd` (55, 16), then `rec_pars` -26,
`beta_linkage` 10, `index_log_q` -7. Total NLL Rceattle 519.71 vs SS3 531.00. Netting off
the constants SS3 drops from its densities, what is left as a difference in fit is:

| component | Rceattle | SS3 | constant | **residual** |
|---|---|---|---|---|
| Age_comp (CAAL) | 407.4520 | 402.4730 | -- | **+4.9790** |
| Catch | -70.2254 | 0.3080 | 34 x (log 0.05 + 0.5 log 2pi) = -70.6110 | **+0.0776** |
| Length_comp | 140.0588 | 140.0590 | -- | **-0.0002** |
| Recruitment | 28.1199 | -2.9041 | 34 x 0.5 log 2pi = 31.2439 | **-0.2199** |
| Survey | 3.9741 | -8.9361 | 14 x 0.5 log 2pi = 12.8651 | **+0.0451** |
| [Rce only] Initial abundance deviates | 12.8700 | -- | -- | -- |
| [Rce only] Linkage-table priors | -2.5415 | -- | -- | -- |

**The ageing-error fix worked and was not the growth gradient.** CAAL fell from +725 to
+4.98, so `build_ss3_age_error()` closed 99.3% of that component. The growth gradient did
**not** move (532 -> 547), so the previous note's hypothesis is refuted.

**The CAAL gap, measured.** What has been ruled out, each against SS3's own Report.sso:

- Observed CAAL is exact. `caal_data` is row-for-row with SS3's `agecomp` (same order,
  same keys) and the proportions agree to 0 (max |diff| over 1157 x 13 cells).
- Sample sizes are exact: `Nsamp_ss3 / Sample_size_rce` is the constant 6.1289572 for every
  row, which is SS3's variance adjustment 0.163372 divided by the add-to-comp 1 + 13 x 1e-4.
- The age-length key is exact **over the whole matrix**, not just Jan 1: `growth_matrix`
  indices 1-2 match SS3's Sub_Seas 1 ALK to 4.3e-7 and indices 3-4 (the per-fleet keys, at
  SS3 month 7 = Rceattle Month 6) match Sub_Seas 2 to 4.4e-7, across all 143 lengths x 14 ages.
  Note SS3's ALK rows come out of `r4ss` in descending length order; sort before comparing.
- N-at-age is exact (7.5e-6) and the survival to survey time agrees: SS3's own
  `natage` mid/begin ratio is 0.8118 and is flat across ages 0-3, so it cancels in the
  conditional either way.
- No ageing-error matrix can close it. Solving for the matrix that maps Rceattle's reported
  `pred_CAAL` onto SS3's `condbase` (a linear, well-posed fit, since every row sums to 1)
  leaves 0.1148 against 0.1176 for the converter's. The error is in the joint, not the smear.
- It is not the length axis: CAAL rows are single 1-cm bins (`Lbin_hi - Lbin_lo` = 1 for all
  1160), and SS3's population bins are the data bins (143, 1 cm, 0.5-142.5).

What is left: the predicted **length marginal** is exact to 1.3e-7 while the **age split
within a length bin** is out by up to 0.118. The residual is concentrated at 23.5-26.5 cm
(ages 1 v 2) and 36.5-40.5 cm (ages 2 v 3) -- where adjacent ages overlap -- and at 25.5 cm
reads Rceattle 0.5705/0.4225 against SS3 0.4529/0.5377. **SS3's `condbase` expectation is
not reproducible from SS3's own printed ALK and N-at-age**: that reconstruction gives
0.5704, i.e. Rceattle's answer, not SS3's. Two one-parameter fits close most of it and are
confounded on a ridge -- a growth timing of 0.465 yr instead of 0.5 (RMSE 0.0018 against
0.0124) or a +0.5 cm shift of the length axis (0.0026) -- and neither reaches 1e-5, so
neither is the cause. **Do not implement either.** The next step is to read SS3's source
for how the conditional age-at-length expectation is formed for a survey at month 7,
because the printed pieces do not compose into the printed answer.

## Known flags (cod bridge)

- **SS3 months are calendar months** (1 = 1 Jan): Rceattle `Month` = SS3 month - 1; fisheries
  use Month 6 (SS3 mid-season ALK). Fixed in the converter (`ss3_month_to_rce()`); this alone
  took AI survey selectivity from 2.9e-2 to 3.7e-7.
- **SS3 age bins start at 1 for both stocks** (AI 1-13, GOA 1-10) while Rceattle ages start at
  0; SS3's ageing error sends true age 0 into bin 1. The converter's `build_ss3_age_error()`
  does the same (age-0 obs column empty). Run and verified 2026-09-23: it took AI CAAL from
  +725 to +4.98. It is **not** the growth gradient, which did not move.
- **Rceattle keeps one age-length key per fleet, not per data row**: `growth_matrix` is
  indexed `nspp * 2 + flt`, and both the length-comp block and the CAAL block read it at
  `flt_month(flt)`. `caal_data` carries no `Month` column, so a CAAL row's timing comes from
  its fleet. That is right for SS3, whose sub-season ALKs these reproduce to 4.4e-7.
- SS3 multinomial = `MultinomialAFSC` x 1/(1 + n_SS3bins x min_comp); the converter sets
  `comp_offset = addtocomp` and divides `Sample_size` by that factor. A `MultinomialSS3`
  family was written and **reverted at Grant's request** (2026-09-23). Tail compression is
  refused with a commented placeholder in the converter.
- `sd_plus_group = "SS3"` is mislabelled for `Growth_Age_for_L2 = 999` (SS3 then pins the plus
  group, which is `"WHAM"`); documented, behaviour unchanged.
- SS3 fecundity is 0 below `First_Mature_Age`; Rceattle has no cut (AI age 0: 1.2e-9 kg, ~2e-9
  of SSB). The harness compares mature ages only.
- GOA notes in `Rceattle-models/GOA cod/Bridging/*.md` have SS3's F_Method numbering backwards
  and call Pope's the blocker; `dev` already uses Baranov with F as parameters (= F_Method 2).
- Windows: `load_all()`'s debug build overflows the object file ("file too big"). Use
  `pkgbuild::compile_dll(".", debug = FALSE)` then `pkgload::load_all(".", compile = FALSE)`;
  the bridge scripts do. `devtools::document()` on this machine churns unrelated `man/` and
  `NAMESPACE` (R 4.5.1 link targets); keep only the intended `.Rd`.
- Every switch value needs a string alias (Grant, standing rule).

## Blocked (cod bridge)

- GOA cannot run on `dev` until Phase 2 items 3-4 and 3b land: its forward pass uses
  branch-only switches (`BlockDev`, `SS3Robust`, `EnvExp`, `*_addtocomp`,
  `Age_first_selected`).

## Before the release

1. Work `inst/dev/SIMPLIFY-LOG.md`. It is the accumulated list of API, switch and workflow
   simplifications found while doing the above. Every row is logged rather than done, by
   standing rule. Grant picks which become PRs after the release.
2. `/ecosystem-sweep` the four consumer repos in `SIBLING-REPOS.md`. 5.35.0 retired
   `estDynamics = 3` and 5.37.0 refuses unknown names in a stored `map`, so a sweep is not
   optional this cycle.
3. Run the hake `MSE_yr2024.R`. It is the only end-to-end `run_mse()` and the only routine
   exercise of estimated suitability; reference objectives are in `SIBLING-REPOS.md`.

## Open work, by where it is recorded

- `TODO-selectivity.md` — `Hake` ignores `Sel_norm_scope` (and `Sel_norm_bin_upper`), with the
  fix sketch and the double-normalization trap in it; the parked `sel-penalty-form` branch; and
  whether a bias-corrected Laplace or `tmbstan` should be the recommended way to report the
  non-parametric deviation SD.
- `CONTRIBUTOR-EXPERIENCE.md` — items A (two of three recipes), B, E, G and H are open, and so
  is item 0, which is still the one that should reorder the rest.
- `TODO-srr-multispecies.md` — initial ages do not decay with M1 + M2 under predation. Built
  and measured during #150, then reverted: year-1 N and M2 feed back across predation
  iterations, and `BS2017MS`'s default starts reached an infinite objective at `niter = 10`.
- `CLEANUP_BACKLOG.md` — everything found and deliberately not fixed, in tiers. Absorbed
  `TODO-5.34-followups.md`: the PFMC `Ftarget` assignment and the `goa_ss` second minimum
  live in `TRAPS.md`, the unbounded `log_Ftarget` and the dead average-F branch here.
- `TODO-projection-module.md`, `TODO-mse-horizon.md` — unchanged by this batch.
- `TRAPS.md` — verified traps with the measured numbers behind them.

## Parked branches

- `sel-penalty-form` (`Sel_penalty_form`, 5 commits) — parked by decision, not by defect.
- `dsem-v5-integration` — PR #111 closed unmerged 2026-09-09.
- `reporting-tables` — one stray doc commit, `4716968c`, which reached `dev` as `3255fb49`
  via PR #132.

## Resume here

**Cod bridge** (`git checkout cod-bridge` and `git pull`; pull `../Rceattle-models` master too):
1. **The one open question is the predicted CAAL** (G1 fails at 1.2e-1, CAAL NLL +4.98,
   `log_growth_pars` gradient 547). Everything it is built from has been checked against
   Report.sso and is exact — see "The CAAL gap" above for what is ruled out and with what
   numbers. Read SS3's source for how a survey's conditional age-at-length expectation is
   formed; do not fit the timing or the length axis, both of which close most of the gap
   and are confounded.
2. The rest of G2 is within 0.35 nats of SS3 once the densities' constants are netted off
   (`parity_report()` prints the residual column). Two blocks have no SS3 counterpart:
   `init_dev` (+12.87) and the linkage-table prior (-2.54). `rec_pars` (-26) and
   `index_log_q` (-7) are still Phase 4 (InitF / equilibrium catch, analytical q); the
   -0.22 recruitment and +0.045 survey residuals are the same two items.
3. Phase 2 items for GOA: length selectivity on population bins (sel, ALK, comps, selected
   weight), SS3 age pattern 10 with length selectivity; Phase 3b multiple ageing-error
   definitions; then port GOA's forward pass off the branch-only switches.
4. Before committing: NEWS 5.42.0 needs AI before/after numbers for selected body weight;
   `/doc-sync`, `/document` (keep only intended `.Rd`), full `/test`. None of the above has
   touched `R/` or `src/` — the only code change this session was to the harness in
   `Rceattle-models` (`3b2fc4b`).

**Otherwise:** read `inst/RELEASE-CHECKLIST.md` and start the release, or pick from
`SIMPLIFY-LOG.md` first. Both are Grant's call.
