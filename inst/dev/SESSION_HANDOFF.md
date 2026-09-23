# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**In flight (2026-09-23): exact SS3 -> Rceattle bridge for AI and GOA Pacific cod, branch
`cod-bridge`** (off `dev` at `cf82f27e`, nothing committed yet, version bumped to 5.42.0 in the
working tree). Goal: from a cold start Rceattle reaches the same solution as an SS3 reference
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

## Done & verified (cod bridge, 2026-09-23, uncommitted on `cod-bridge`)

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

**AI cod G1: all 9 checks PASS at tol 1e-5 (Report.sso print precision)** — last run before
the Phase 3 converter edits: length-at-age 4.5e-6, weight 2.4e-6, fecundity 2.8e-6, ALK
4.3e-7, sel FshComb 7.3e-7, sel Srv 3.7e-7, N-at-age 7.5e-6, SSB 4.7e-6, R 3.9e-6.
**G2 still FAILS: max |gradient| 532, all on `log_growth_pars`/`growth_log_sd`**; NLL
Rceattle 1248.10 vs SS3 531.00 (CAAL +725, catch -70.5 = lognormal constant
34 x log(0.05 sqrt(2 pi)), survey +12.9, recruitment +43.9).

## Known flags (cod bridge)

- **SS3 months are calendar months** (1 = 1 Jan): Rceattle `Month` = SS3 month - 1; fisheries
  use Month 6 (SS3 mid-season ALK). Fixed in the converter (`ss3_month_to_rce()`); this alone
  took AI survey selectivity from 2.9e-2 to 3.7e-7.
- **SS3 age bins start at 1 for both stocks** (AI 1-13, GOA 1-10) while Rceattle ages start at
  0; SS3's ageing error sends true age 0 into bin 1. The converter's new
  `build_ss3_age_error()` does the same (age-0 obs column empty) — the likely source of the G2
  growth gradient. **Not yet run.**
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

**Cod bridge** (`git checkout cod-bridge`; all work is uncommitted, including in
`../Rceattle-models`):
1. Re-run `test-schema-cpp-dispatch.R` (exemptions added after the last run).
2. Rerun AI parity (G1 + G2) with the new converter ageing error / sample-size factor; expect
   the growth gradient to fall. Then work G2 down: recruitment (+43.9) and survey (+12.9) should
   be constants only; `rec_pars` (-26) and `index_log_q` (-7) point at Phase 4 (InitF /
   equilibrium catch, analytical q).
3. Phase 2 items for GOA: length selectivity on population bins (sel, ALK, comps, selected
   weight), SS3 age pattern 10 with length selectivity; Phase 3b multiple ageing-error
   definitions; then port GOA's forward pass off the branch-only switches.
4. Before committing: NEWS 5.42.0 needs AI before/after numbers for selected body weight;
   `/doc-sync`, `/document` (keep only intended `.Rd`), full `/test`.

**Otherwise:** read `inst/RELEASE-CHECKLIST.md` and start the release, or pick from
`SIMPLIFY-LOG.md` first. Both are Grant's call.
