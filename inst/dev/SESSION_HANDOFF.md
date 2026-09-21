# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

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
- `TODO-5.34-followups.md` — PFMC `Ftarget` assigned inside the F loop; `zero_N_pen`
  over-count; the dead average-F refit in `run_mse()`.
- `TODO-projection-module.md`, `TODO-mse-horizon.md` — unchanged by this batch.
- `CLEANUP_BACKLOG.md` — everything found and deliberately not fixed, in tiers.
- `TRAPS.md` — verified traps with the measured numbers behind them.

## Parked branches

- `sel-penalty-form` (`Sel_penalty_form`, 5 commits) — parked by decision, not by defect.
- `dsem-v5-integration` — PR #111 closed unmerged 2026-09-09.
- `reporting-tables` — one stray doc commit, `4716968c`, which reached `dev` as `3255fb49`
  via PR #132.

## Resume here

Read `inst/RELEASE-CHECKLIST.md` and start the release, or pick from `SIMPLIFY-LOG.md` first.
Both are Grant's call.
