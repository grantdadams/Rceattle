# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**PR #119 is merged** (`bd9d0fcf` on `main`, 2026-08-26). `main` and `dev` are both 5.20.0.

**In flight: `dev` → `main` for 5.21.0, component likelihood profiles.** `dev` is one prior
commit ahead of `main` (`ccf7e6d2`, the release-checklist edit) plus this work. Small PR, not a
release of accumulated branches.

Prompted by Cole asking for a built-in way to extract the NLL components of a profile and plot
them — the total answers how well the data determine a parameter, not whether the data sources
agree about it.

- `profile_components()` (`R/9-profile.R`) — `quantities$jnll_comp` from every grid fit as one
  long data frame, each cell labelled with the fleet or species it belongs to.
- `plot_profile()` + `plot.Rceattle_profile()` (`R/7-plot_profile.R`) — `r4ss::SSplotProfile()`
  layout: each series re-zeroed at its own minimum with a point marking it, total black and out
  of the colour legend, `minfraction = 0.01`.
- `profile()` gains `$alias`, so `print()` and the axis say `M1` where the grid holds M rather
  than `log_M1`, which is what `param` resolves to.

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

## Resume here

Merge the 5.21.0 PR once CI is green, then tag `v5.21.0` and draft a GitHub Release —
`inst/RELEASE-CHECKLIST.md` §3. The `pkgdown` workflow rebuilds on the `release` event.

Two loose ends from this session, neither blocking:

1. **`inst/dev/SIBLING-REPOS.md` has an uncommitted edit that predates this session** — the
   `../GOA-multispecies-assessment` entry. It is coherent and complete; it was left out of the
   5.21.0 commit because it is unrelated to that work, not because it is wrong. Commit it
   whenever suits.
2. **`../Rceattle-models/GOA pollock/2025/04-fit-and-diagnostics.R` gained an M-at-age-6
   component profile** (uncommitted, separate repo). `minage = 1` and `nages = 10`, so age 6 is
   bin 6; `M1_base` runs 1.39 at age 1 down to 0.30 from age 6 on. It fixes the age-6 cell
   **alone**, leaving ages 7–10 at their base values — narrower than goa_pk's own M profile,
   which scales the whole M-at-age vector together. Not yet run.

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
