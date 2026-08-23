# Working the cleanup backlog

How to turn `inst/dev/CLEANUP_BACKLOG.md` into shipped fixes, in an order that respects what
each item can break. Written 2026-08-22, after the Claude-readiness work landed the drift guards
that make several of these safe to touch.

## The shape of the problem

81 `TODO`/`FIXME` markers, triaged into three tiers. **The tiers are not a priority order.** A
Tier 0 "known defect" that no live assessment can reach is less urgent than a Tier 1 limitation
every MSE run hits. Sequence by *who is exposed*, not by tier label.

Three facts that should drive the order:

1. **Every Tier 0 defect names the input that triggers it.** None is a latent mystery; each says,
   in the source, "this is wrong when X". So each can be reproduced before it is fixed, and the
   reproduction is the test.
2. **Most are unreachable from the bundled data.** `minage > 1`, `nlengths < nages`,
   `assessment_period > 1`, a non-survey species — none occurs in `BS2017SS`, `BS2017MS`,
   `GOA2018SS`, or the three live assessments. That is why they have survived; it also means
   `/golden-check` will stay green through the fix and prove nothing. **Each fix needs a fixture
   that reaches it.**
3. **Two are silently wrong rather than loud.** QAR1 returns a constant q, and the HCR-2 proxy
   uses a fixed depletion. Those produce a number, not an error, which is the failure mode this
   package cannot afford.

## Recommended order

### First — the two that produce a wrong number silently

**1. QAR1 (`Catchability = "AR1"`).** Already deprecated with a warning in 5.12.0, so nobody is
silently affected today. Decide whether to fix or remove: the deviate map is gated on
`Time_varying_q`, which under this form holds an `env_data` column index rather than a mode.
Fixing means gating on `Catchability` instead. Removing means dropping the form and pointing at
the AR1 linkage, which is what the warning already recommends. **Removal is probably right** —
the linkage grammar supersedes it and carries priors, bounds and a phase, which this does not.

**2. `mse_summary()`'s HCR-2 fixed-depletion proxy.** The EM uses `0.5 * 0.35` rather than the
configured target, so a Tier-2 MSE reports performance against the wrong reference point.
Reachable by anyone running HCR 2. Needs a fixture with a non-default `Ptarget`.

### Second — the three that error loudly on a real input

These stop rather than mislead, so they are lower risk but they *block* work.

**3. `minage > 1` in the template.** Same class as the `nages`-is-a-count trap the plotters had,
but in the C++. The R side is already fixture-covered (`test-plot-smoke.R` runs the plot family at
`minage = 3`); this needs the equivalent for the fit. **Do this before any assessment adopts a
`minage != 1` species** — it is the one item with a plausible near-term trigger.

**4. `nlengths < nages`** in the composition block.

**5. `run_mse()` catch fill-in at `assessment_period > 1`.** Every non-annual-assessment MSE.
Given how many strategies assess biennially, this may deserve promoting above (3) — check whether
any live MSE configuration uses it before deciding.

### Third — the two dispatch gaps of the shape already fixed once

**6. `CAAL_distribution = "MultinomialAFSC"`** passes `validate_switches()` and then errors in the
template. Identical in shape to `Estimate_catch_sd = "Analytical"`, fixed in 5.12.0 by
implementing the missing `case`. Decide the same way: implement, or reject in R with a message
that names the alternative. `test-schema-cpp-dispatch.R` pins it as a known gap either way.

**7. `random_sel = TRUE`** in `build_map()` — marked with a question mark, so **confirm it is real
before scheduling it.** It may already work.

### Fourth — the structural one

**8. Split `R/0-build_srr_and_M.R`** (1,497 lines, 29 functions: stock-recruit, M1 and growth
constructors in one file). The last outstanding item from the superseded
`accessibility-and-code-review` plan. Pure refactor, no behaviour, but it touches the
`build_*()` public surface, so it needs `/golden-check` and an `../Rceattle-models` sweep. Keep
the `0-` prefix — the `R/` numbering is meaningful.

### Never — the ones to leave alone

`fit_mod()`'s `suppressWarnings()` swallowing every `build_map()` warning is worth narrowing, but
it is load-bearing for the diagnostic refits (which re-validate a data list the caller already
fitted). The `logH_*` / `H_4` / `log_gam_*` markers belong to stubbed features and are pinned as
such. The six `TODO(review)` markers are Grant's judgement calls and are never resolved by an
agent.

## How to do each one

The same shape every time, because the common failure here is a fix that passes the suite while
the defect survives:

1. **Reproduce first, in a test that fails.** The marker names the triggering input; build the
   fixture that reaches it. A fix whose test passes before the fix is not a test.
2. **Check what actually covers it.** `/golden-check` will be green either way for almost every
   item on this list — none of the four reference models reaches these inputs. Use `/verify` to
   pick the right `tools/verify/*.R` harness, and remember no harness reaches
   `sample_rec(update_model = TRUE)`, `reweight_comps()`, or any figure.
3. **For a C++ change**, recompile before testing (`pkgload::load_all(".")`) and run the suite
   serially (`TESTTHAT_PARALLEL=false`). Golden is required even when you expect no movement —
   the analytical catch sigma was bit-identical and running it was still the right call.
4. **Adversarially review the diff before committing.** Across this work every review found
   something real, including two changes that moved fitted numbers and would otherwise have
   shipped.
5. **NEWS + `DESCRIPTION` + the affected vignette, same commit.** `/doc-sync` checks it.
6. **Move the item to the "Deliberately not changed" section or delete it** — a backlog that only
   grows stops being read.

## Keeping the backlog honest

- **Cite the marker text, not a line number.** These references have gone stale three times: twice
  when roxygen was added above them, once from an unrelated `ceattle.cpp` edit.
- **A Tier 0 claim is a claim about behaviour, so check it against the code, not just the comment.**
  One entry named the wrong switch (`Time_varying_q` instead of `Catchability`) for a whole
  session, and the source comment goes out of its way to warn about exactly that confusion.
- **Re-derive the counts when you touch it.** `Rscript` over `R/` and `src/TMB/` for
  `TODO|FIXME`, excluding the `todo <-` variable in `R/6-process_residuals.R`.
