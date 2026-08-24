# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**PR #119 releases `dev` onto `main`: 5.8.1 → 5.15.0, 133 commits.** It is the first release on
this line since 5.8.1, and it carries nine merged PRs (#109, #110, #112–#118) plus three commits
added on the PR itself. `main` is 0 commits ahead, so the merge is a fast-forward.

Open as a **draft** until CI reports. Once merged, tag `v5.15.0` and draft a GitHub Release from
it — `inst/RELEASE-CHECKLIST.md` §3. The `pkgdown` workflow rebuilds on the `release` event.

## Done & verified on #119

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

Merge #119 once CI is green, then tag and release. After that, `inst/dev/CLEANUP_BACKLOG.md` has
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
