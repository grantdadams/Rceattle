# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**`dev` is at 5.41.0 and `main` at 5.33.0.** The next step is one `dev` -> `main` release
covering 5.34.0 through 5.41.0, per `inst/RELEASE-CHECKLIST.md`. Read that file's pkgdown note
before tagging: the `release: published` event has silently failed to fire once already.

**Everything the checklist asks for before tagging is done** (2026-09-23). Full suite 9,506
assertions / 0 failures / 3 skips; ecosystem sweep clean; hake `MSE_yr2024.R` identical to
5.33.0 on all six fits; a reproducible install into a temporary library driven through them.

**The release sequence, from here:**

1. Merge whatever `inst/dev` PRs are still open (they are documentation only).
2. Open and merge the `dev` -> `main` release PR. Say what forces a refit, what breaks and what
   is new; do not paste `NEWS.md`.
3. Tag the MERGE COMMIT on `main`, then publish a GitHub Release from the tag.
4. Confirm pkgdown actually rebuilt, then `gh workflow run deep-checks.yaml --ref main`.
5. Tell the consumer repos to pin the tag rather than track `main`.

**Two things are red before the release starts, and neither is from this batch.** Do not read
either as evidence against the tag, and do not spend the release chasing them:

- **`deep-checks` `golden` fails on `main` at 5.33.0.** `goa_ss` lands in the second local
  minimum and `goa_ms` inherits it through its warm start. See `TRAPS.md`; the robustness fix
  is the first job after the release, because until it lands this guard cannot gate anything.
- **Windows `R-CMD-check` fails intermittently**, about 2 runs in 30, with an access violation
  that also reproduces on `main`. The file the framework names carries no information. See
  `TRAPS.md`.

macOS was red 2026-09-20 to 09-22 for two unrelated upstream reasons and recovered on its own;
branch `ci/macos-libomp` holds an unmerged remedy if the OpenMP one recurs.

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
| 5.40.0 | #151 | `NonParametricIntegrable` (13) |
| 5.41.0 | #152 | `osa_residuals(method = "cdf")` |

## After the release, in order

1. **Make `golden` robust**, so `deep-checks` can gate a release. Warm-start the reference fits
   from the pinned parameters, or take the lower of two starts. It is a harness change and
   cannot move a fitted number.
2. **Decide on the three inert test guards** (`CLEANUP_BACKLOG.md`): restore or delete. The
   multispecies one is hiding an unexplained disagreement with the old EBS CEATTLE.
3. **Work `SIMPLIFY-LOG.md`.** Fifteen open rows: four change behaviour and need a deprecation
   path or a shim, one moves a golden reference, four are internal, two additive, one doc.
   Every row is logged rather than done, by standing rule; Grant picks which become PRs.

Done for this cycle, so do not repeat them: the ecosystem sweep of the four consumer repos, and
the hake `MSE_yr2024.R` run. Both are recorded above with their results.

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

Read `inst/RELEASE-CHECKLIST.md` and start the release, or pick from `SIMPLIFY-LOG.md` first.
Both are Grant's call.
