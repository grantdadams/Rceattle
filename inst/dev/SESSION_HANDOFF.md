# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**`dev` is at 5.42.1 and `main` at 5.33.0.** The next step is one `dev` -> `main` release
covering 5.34.0 through 5.42.1, per `inst/RELEASE-CHECKLIST.md`. Read that file's pkgdown note
before tagging: the `release: published` event has silently failed to fire once already.
**The tag is the DESCRIPTION version, so it is 5.42.1, not 5.41.0** -- #160 (5.42.0) and the
review of #158 (5.42.1) both landed after the release PR was written.

**The pre-tag measurements are older than this head.** The checklist run recorded 2026-09-23 --
ecosystem sweep clean, hake `MSE_yr2024.R` identical to 5.33.0 on all six fits, a reproducible
install into a temporary library driven through them -- predates 5.42.0 and 5.42.1 entirely,
and the 9,506-assertion suite figure was measured 2026-09-21. **Re-run all four before
tagging.** What HAS been re-done at 5.42.1: the ecosystem sweep, widened to the 5.42.0
refusals (375 workbooks, 183 with a `fleet_control` sheet, across the four consumer repos) --
one negative `Sel_curve_pen1` anywhere, EBS pollock 2024's ATS/AVO fleets, `NonParametricPM`
with no `Sel_shape_mode` column, which the directional exemption keeps legal; no workbook sets
`Sel_shape_dir` or `Sel_devmag_sd`; no selectivity prior or apical linkage sits on an `Off`
fleet or a shared-block follower (the GOA pollock 2025 prior fleets are each their group's
lead); and no group anywhere mixes selectivity forms.

**PR #159 merged into `dev` on 2026-09-24**, #160 on 2026-09-25 (`dev` head `c01ea717`), and
the review of #158 after it. #159 asked four questions of #158: does the language read as
AI-written, is the API frictionless, are the docs concise, can a developer find and change the
model. Full suite green (235 files, 0 failures), golden unchanged to ~1e-11.

What a reviewer should still go at hardest:

- **The selectivity form collapse.** `NonParametricIID` (13) and `NonParametricRW` (14) became
  one `NonParametricIntegrable` (13), with `Time_varying_sel` picking the structure. Code 14 is
  free. **Golden does not cover it**: the four reference models use forms 0-4 (verified: the
  union of their `Selectivity` columns is {0,1,2,3,4}), so golden passing only shows the new
  `sel_case` dispatch is inert for the other forms. What covers the merge is
  `test-selectivity-nonparametric-integrable.R` and its independent `dnorm` oracle.
- **The collapse left stale text in six shipped schema descriptions and in NEWS 5.40.0.** Fixed
  on `fix/release-doc-corrections`; `test-docs-anchors.R` now fails if any schema description
  names a selectivity code `sel_map` does not accept. A schema `doc` string is written verbatim
  into `meta_data_names.xlsx`, so it is user-facing, not a comment.
- **The language sweep touches 69 files** and is isolated in one commit. It recasts clauses
  rather than transliterating dashes; the risk to look for is a `carry` that meant *propagate*
  being flattened to `hold`. Six such were caught in roxygen; assume more exist.

**The release sequence, from here:**

1. Merge the 5.42.1 review branch, then whatever `inst/dev` PRs are still open (they are
   documentation only). `fix/release-doc-corrections` is already in.
2. Re-run the four checklist measurements at the merge head (suite, sweep, hake MSE, install):
   the recorded ones predate 5.42.0. Then merge the `dev` -> `main` release PR #158. Its body
   must say what forces a refit, what breaks and what is new, and must cover 5.42.0 and 5.42.1;
   do not paste `NEWS.md`.
3. Tag the MERGE COMMIT on `main` with the DESCRIPTION version, then publish a GitHub Release
   from the tag.
4. Confirm pkgdown actually rebuilt, then `gh workflow run deep-checks.yaml --ref main`.
5. Tell the consumer repos to pin the tag rather than track `main`.

**The last installable tag is `5.28.0`, not 5.33.0.** `main` carried 5.29.0, 5.30.0, 5.31.0,
5.32.0, 5.32.1 and 5.33.0 without a tag being pushed for any of them. So a consumer who pins
tags, which is what step 5 asks for, moves **5.28.0 -> 5.42.1**, fourteen minor versions, not
eight. Say that in the release body, and treat step 3 as the fragile step it has proven to be:
the pkgdown `release: published` miss at 5.21.0 is the same step failing in a different way.

**Three things are red before the release starts, and none is from this batch.** None is
evidence against the tag. But "known" is not "ignore": read each post-merge run against the
signature below, because that is the only thing separating a known red from a new one:

- **`deep-checks` `golden` fails on `main` at 5.33.0.** `goa_ss` lands in the second local
  minimum and `goa_ms` inherits it through its warm start. Last fully green 2026-09-13.
  **Before dispatching `deep-checks` on `main`, expect exactly this signature: a `goa_ss` delta
  of 52.9 with the other three models bit-identical.** Any other pattern is a real regression
  and stops the release. Diagnose from the gradient at the reference `par`, not the objective
  (`TRAPS.md`). The robustness fix is the first job after the release and ships as 5.42.2; until
  it lands this guard cannot gate anything.
- **`deep-checks` `suite` never finishes.** It is `cancelled` in every recent run at 5h00-5h01
  wall clock, a timeout rather than a pass. So the one job that runs the 140 `skip_on_cran()`
  test files proves nothing about `main` right now. Either raise the ceiling or shard it; until
  then do not cite a `deep-checks` run as suite coverage.
- **`R-CMD-check` fails intermittently, and it is mostly macOS, not Windows.** Over the last 100
  runs: 79 completed, 8 failures, of which 7 were macOS-only and 2 touched Windows. The Windows
  access violation is real and reproduces on `main` (the file the framework names carries no
  information, see `TRAPS.md`), but it is ~2 in 79, not 2 in 30.

macOS was red 2026-09-20 to 09-22 for two unrelated upstream reasons and recovered on its own;
branch `ci/macos-libomp` holds an unmerged remedy if the OpenMP one recurs.

The 2026-09-14 backlog plan is finished. Nine branches across eight versions, listed below in
version order (#150 merged before #149; 5.37.0 took two branches), each reviewed
adversarially before commit and again by a second session before merge:

| Version | PR | What landed |
|---|---|---|
| 5.34.0 | #144 | MSE dynamic-SB0 and fixed-numbers reference-point masking |
| 5.35.0 | #145 | Silent-wrong-number fixes; `estDynamics = 3` retired |
| 5.36.0 | #146 | Config overlay by field; OSA outliers flagged per panel |
| 5.37.0 | #147, #148 | QAR1 path removed; stored-map guard; `CONTRIBUTING.md`, the Doxygen build and `adding-a-selectivity-form.Rmd` |
| 5.38.0 | #149 | Per-sex apical selectivity offset (`log_sel_apical`) |
| 5.39.0 | #150 | Multispecies stock-recruit bounds and a degenerate-curve check |
| 5.40.0 | #151 | Two integrable non-parametric forms (13, 14), collapsed to `NonParametricIntegrable` (13) in #159 |
| 5.41.0 | #152 | `osa_residuals(method = "cdf")` |
| 5.42.0 | #160 | Selectivity-penalty sign refusals; apical/prior `Off` and shared-block gates |
| 5.42.1 | review of #158 | Penalty lead keyed as the template keys it; the code-14 guard widened |

## After the release, in order

1. **Make `golden` robust**, so `deep-checks` can gate a release. Warm-start the reference fits
   from the pinned parameters, or take the lower of two starts. It is a harness change and
   cannot move a fitted number. **This is not just the next cleanup: it gates the NOAA
   transfer** (`PLAN-adoption-and-NOAA-transfer.md` section 0, item 5) and it has a release
   vehicle already chosen, 5.42.2 (`TODO-pre-transfer.md` B3). Do it before anything below.
   While doing it, fix the `deep-checks` `suite` timeout too; a guard that cannot finish is
   the same problem in a different job.
2. **Decide on the three inert test guards** (`CLEANUP_BACKLOG.md`): restore or delete. The
   multispecies one is hiding an unexplained disagreement with the old EBS CEATTLE.
3. **Work `SIMPLIFY-LOG.md`.** Seventeen rows, three struck through, so **fourteen open**: six
   change behaviour (two of them needing a deprecation path or a shim), one moves a golden
   reference, four are internal, two additive, one doc. Every row is logged rather than done,
   by standing rule; Grant picks which become PRs.

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
- `PLAN-adoption-and-NOAA-transfer.md` — moving the package off a personal account to a NOAA
  org, and the two adoption barriers behind it. **Section 0 holds five decisions reserved for
  Grant** (destination org, license, co-maintainer, whether `Rceattle-models` moves, timing);
  agents do not pick these.
- `TODO-pre-transfer.md` — the execution checklist for that plan, stages A-F with owner tags.
  Stage B is this release. **B3 is the `golden` robustness fix**, to ship as 5.42.2 if it lands
  after the tag.

## Parked branches

- `sel-penalty-form` (`Sel_penalty_form`, 5 commits) — parked by decision, not by defect.
- `dsem-v5-integration` — PR #111 closed unmerged 2026-09-09.
- `reporting-tables` — **local only, never pushed**, so it is not on the remote to triage. Its
  one stray doc commit, `4716968c`, reached `dev` as `3255fb49` via PR #132 (which merged from
  `docs/minfraction`, so the content was re-applied rather than merged from this branch).

## Resume here

Read `inst/RELEASE-CHECKLIST.md` and start the release, or pick from `SIMPLIFY-LOG.md` first.
Both are Grant's call.
