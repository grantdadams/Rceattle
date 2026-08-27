# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**PR #111 (`dsem-v5-integration` -> `dev`) is at 5.23.0 and synced with `dev`.** `dev` released
5.20.0, 5.21.0 and 5.22.0 while this branch was open; all are merged in.

**The version keeps colliding.** Both lines have now independently used 5.19.0, 5.21.0 and
5.22.0, and each time the branch had to move up after the fact and renumber a NEWS section that
was already written. `dev`'s number always wins, because `dev` is what releases. Before bumping
this branch again, read `dev`'s DESCRIPTION -- or better, leave the bump until the PR is about
to merge, since anything chosen earlier is a guess about what `dev` will not use.

The last block of work answered a report that **`self_test(process = TRUE)` did not simulate
process error** — on a DSEM or without one. It did not, for five separate reasons; all five are
fixed, plus two features and two crash/consistency fixes that came out of the same conversation.

## Done & verified on #111

Full suite **198 files, 0 failures, 4 skips** (all pre-existing), run after the `dev` merge that
brought component profiles, on the pushed tree.

## Done & verified on #111

Full suite **196 files, 0 failures, 4 skips** (all pre-existing: three `Skipping` guards and
`test-dsem-recruitment.R:23`, whose reference objectives predate v5.0). Verified after the `dev`
merge, on the pushed tree.

- **Golden reference green** — the four fits reproduce their pinned objectives. Every new draw is
  gated on `dsem_on == 1`, so non-DSEM numerics are untouched by design and by measurement.
- **`tools/verify/verify-dsem-equivalence.R` PASSES — for the first time ever on this branch.**
  See "Known flags"; it was dead.
- **End to end on a `family = "normal"` DSEM**: fit, redraw, write back, rebuild, refit, both
  replicates converged.
- Golden, full suite and the equivalence harness were all re-run **after** the merge, not only
  before it.

## What the five defects were

Worth reading before touching `sim_mod()` or `dsem.hpp`; the detail and the measured numbers are
in `inst/dev/TRAPS.md` under "The SIMULATE contract".

1. **`calculate_dsem()` had its own `do_simulate` branch that discarded the R draw.** It assigns
   the whole latent field (`x_tj = draw_tj`), and it ran *after* `sim_mod()` substituted its
   conditional draw — adding 100 to every latent state changed nothing. It consulted neither the
   map nor `dsem_cond_k`, so it also redrew the covariate columns; under `family = "fixed"` those
   columns **are** the `env_data` series, so each replicate was generated under an environment the
   refit never saw (a scaled covariate with sd 1 moved by up to 2.96). The draw is now taken only
   in R and `dsem_do_sim` is pinned to 0.
2. **The mask stayed zero, so no truth came back.** `rec_dev_drawn_sim` is written inside the
   template's IID block, which is gated off under a DSEM, so `attr(x, "process_sim")` was `NULL`
   on exactly the models whose process *had* been redrawn — indistinguishable from nothing
   happening, and enough to stop `compare_sim()` warning it was measuring bias against the wrong
   baseline. The R draw now carries its own mask and hands back `dsem_x_tj` as well as `rec_dev`.
3. **`init_dev` was gated on `dsem_on == 0`**, so a replicate's initial age structure kept its
   fitted values while the hindcast was re-realized — two histories in one data set.
4. **The AMAK/Ianelli gate did not reach the SEM draw.** `srr_fun = 0` with `srr_pred_fun > 0`
   scores the deviations twice; the template has always refused the standard draw for that, the
   SEM draw did not.
5. **`self_test()` swallowed every `sim_mod()` warning.** Replicates run through `parLapply` and a
   worker's warnings never reach the caller, so at the default `cores` a `process` request that
   redrew nothing was completely silent. The unique set is now re-emitted once, and asking for
   process error and getting none back warns in its own right.

## The two features from the same conversation

- **A DSEM covariate with a measurement family is simulated, state and observation.**
  `family = "fixed"` is data and is untouched; any other family makes the covariate a latent state
  observed with error, so the state is drawn with the field and the observation is drawn from its
  own density (all seven families) and written back into `env_data`. Checked as a moment: over 200
  draws, `sd(y_sim - mu)` was 0.9995 against a fitted measurement sd of 1.
- **`build_bounds(dsem = )` bounds the SEM's standard deviations at 0.** A lag-0 two-headed self
  path is a Cholesky diagonal whose sign is unidentified, so the surface is exactly symmetric
  about 0 — harmless at the MLE, fatal for MCMC. Only diagonals that stand alone in their row
  qualify; a cross-covariance or a lagged `<->` keeps unbounded support and `fit_mod()` names those
  variables. A fit at the negative root is flipped rather than rejected (10401.62 either way).
  `vignette("dsem")` shows the `tmbstan` call.

## For Grant's review

- **`NEWS.md` 5.22.0 folds two release blocks into one.** The branch had filed its DSEM work as
  5.19.0 while separate from `dev`, and `dev` released a different 5.19.0. Both existed after the
  merge. `dev`'s keeps its number; the branch's is merged into 5.21.0 under the four normal
  headings (39 bullets in, 39 out). Read it once as a single release — the wording of the two
  halves has not been harmonized.
- **The PR body is out of date in the worst direction.** It still says `sim_mod(process = TRUE)`
  refuses on a DSEM and that a `self_test()` does not re-realize the process. Both are now the
  opposite, and that is the headline of the change.

## Known flags

- **A `tools/verify/` harness runs in NO CI job, so a broken one is invisible.**
  `verify-dsem-equivalence.R` — the only numerical check on the vendored `dsem.hpp` — could not
  run for the whole time it sat on this branch: `cond_k` was added to `calculate_dsem()` and to
  `dsem_standalone.cpp`, but the harness's `dsem_data()` was never taught to supply it, so every
  invocation died in `getParameterOrder()` before the first comparison. **Adding an argument to
  `calculate_dsem()` means updating three places, not two** — the header, `ceattle.cpp`, and
  `tools/verify/dsem_standalone.cpp` *plus* the R harness feeding it. Now in `TRAPS.md`.
- **A file split defeats git's rename detection and merges cleanly wrong.** `dev` moved
  `self_test()` and `profile.Rceattle()` out of `R/9-retro_and_jitter.R` into new files. Those
  arrived *already staged as adds*, never appeared as conflicts, and a `git commit` of the merge
  took `dev`'s copies — dropping every change this branch had made to both functions, silently,
  leaving a package that loads and functions that exist. Fixed in `6ff40df5` by taking the
  branch's version of each and re-applying `dev`'s two small deltas on top. **After any merge that
  moves a function between files, assert each one is defined exactly once and grep for your own
  change inside the new file.**
- `sim_mod(simulate = FALSE)` returns `env_data` as supplied, not as fitted values — the one data
  type that does not return its hat. Deliberate and documented in `?sim_mod`.
- An `env_data` column that is also a `linkage_spec` covariate is an errors-in-variables case
  under a measurement family: the replicate is generated from the series as supplied while the
  refit sees the redrawn one. Faithful to the real assessment, documented, not a defect.
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

**The bounds check has measured nothing since at least 2026-08-25, and reports green.**
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

Two separable fixes: the DLL bootstrap, and a report that says "could not run" rather than
"BOUNDS VIOLATIONS" when no case executed. This matters because it is the only net for the
Windows `0xC0000005`, which is still uncured — see the 5.21.0 release notes' known limitation.

## Resume here

Update the PR body -- it still describes only the original DSEM work, not `bound_sd`, the
zero-variance guard, or `$bounds$par_lower` -- read `NEWS.md` 5.23.0 once as a whole, then merge
#111. `run_mse()` and two `process_residuals()` processes still refuse on a DSEM and are the
remaining gap in that suite.

Loose ends inherited from `dev`, none blocking:

1. **`inst/dev/SIBLING-REPOS.md` has an uncommitted edit that predates these sessions** — the
   `../GOA-multispecies-assessment` entry. Coherent and complete; left out of the 5.21.0 and
   5.22.0 commits as unrelated, not as wrong.
2. **`../Rceattle-models/GOA pollock/2025/04-fit-and-diagnostics.R` gained an M-at-age-6
   component profile** (uncommitted, separate repo, not yet run). `minage = 1`, `nages = 10`, and
   `M1_base` runs 1.39 at age 1 to 0.30 from age 6 on, so age 6 is bin 6 and 0.30 is the value a
   0.18-0.42 grid brackets. It fixes the age-6 cell **alone**, leaving ages 7-10 at base --
   narrower than goa_pk's own M profile, which scales the whole vector. **5.22.0's
   `joint = "multiply"` now expresses the goa_pk version**; switching it is a two-line change.
3. **The dead QAR1 path.** `ceattle.cpp:4269` was fixed in place rather than removed because
   `R/6-process_residuals.R:204` still branches on `est_index_q == 6`. Removing both together is
   Grant's call.

After #111: `inst/dev/CLEANUP_BACKLOG.md` has 64 items left (Tier 0 and Tier 2 cleared);
`inst/dev/BACKLOG-PLAN.md` sequences them by who is exposed. Queued and not started:
`inst/dev/CONTRIBUTOR-EXPERIENCE.md` — **item 0 comes first and is a conversation, not code**: ask
the sibling-repo authors where they actually stopped.

## Older paused work

**The `accessibility-and-code-review` refactor is superseded, not outstanding.** The branch is
gone from local and remote, its plan no longer exists, and no commit in any branch mentions it.
Three of its four locked "chosen extras" landed by other routes — the `JnllRow` enum, the repaired
cpp section index, and Doxygen `@file`/`@brief` blocks on the previously bare headers. Its one
concrete leftover, splitting `R/0-build_srr_and_M.R`, is done in 5.14.0. Nothing of that plan is
outstanding. Do not resume from the branch name or the handoff path; both are dead references.
