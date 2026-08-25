# Verified traps

Each entry is something that has actually bitten, with the evidence. `CLAUDE.md` carries a
one-line version of the ones you need before touching the code; this file holds the detail and
the measured numbers so a claim can be checked rather than believed.

Verified against source, 2026-08.

## Build and test

**Dev builds are `-O2`, not pkgbuild's `-O0`.** The repo `.Rprofile` sets
`options(pkg.build_extra_flags = FALSE)`, so `load_all()` compiles the TMB model with the same
optimization as a production `R CMD INSTALL` — `fit_mod()` runs ~10x faster than an
unoptimized debug build, bit-identically (measured). A normal install was always `-O2`; only
`load_all` was `-O0`. To debug the C++ line-by-line (gdb/lldb), start R with
`RCEATTLE_DEBUG_CPP=1`. **Any absolute fit timing must state its build** — an `-O0` number
overstates real cost ~10x.

**A slow fit is the model, not a regression.** A full `BS2017SS` `estimateMode = 0,
phase = TRUE` fit is ~55 s: ~14 s for one optimization, ~3x for phasing, +~14 s for
`sdreport`. That single fit has needed ~500–700 `nlminb` iterations since at least 2023
(bisected back through `main`), so don't hunt a recent cause. For dev loops
`phase = FALSE, getsd = FALSE` gives ~14 s.

**Editing a `.hpp` alone does not rebuild the model.** `TMB::compile()` compares the `.o`
against the `.cpp` and nothing else, so a change confined to a header leaves `load_all()`
reporting success while the loaded DLL is still the old one. Every process in this model
lives in a header — `recruitment.hpp`, `predation.hpp`, `selectivity.hpp`, `dsem.hpp` — so
that is most C++ edits here. It reads as a change that "does not take effect", or as an
error message you already deleted still appearing. Force it:

```sh
rm -f src/TMB/ceattle.o src/TMB/ceattle.so src/ceattle.so
cd src/TMB && Rscript compile.R      # load_all() will not rebuild once the .o is gone
```

Touching `ceattle.cpp` works too. Editing a `.cpp` and a `.hpp` in the same round hides
this, which is why it takes a while to notice.

**The parallel-test DLL failure is a toolchain failure, not a test failure.**
`DESCRIPTION` sets `Config/testthat/parallel: true`, and the workers cannot load a freshly
rebuilt DLL, so `devtools::test()` aborts before running anything with
`testthat subprocess failed to start … getDLLRegisteredRoutines.DLLInfo`. After a `.cpp`/`.hpp`
edit run `pkgload::load_all(".")` first, then test with `TESTTHAT_PARALLEL=false`.
(`test-coverage.yaml` forces serial for a related reason.)

**Roxygen is pinned at 8.1.0** (`DESCRIPTION: Config/roxygen2/version`). A *different* local
roxygen2 rewrites *every* `man/*.Rd`; an older one also swaps the key back to the legacy
`RoxygenNote:`, leaving both present and contradicting each other. After `document()`, check
`git diff DESCRIPTION` **first**: if the version key moved, the `man/` churn is the version, not
your change. Either regenerate everything under one version deliberately, or `git checkout` the
unrelated churn and keep only the `.Rd` you meant to change. A matching version is necessary but
not sufficient — the build environment can rewrite `man/` and `NAMESPACE` too.

## Shared parameter blocks

**Fleets sharing a `Selectivity_index` / `Catchability_index` share ONE parameter block.**
(`Q_index` is a deprecated alias of `Catchability_index`; the schema upgrades it on read.)
`adjust_map_shared_params()` copies the donor fleet's map slice over the rest of the group, so
any per-fleet setting that differs within a group is silently taken from the donor
(`Sel_start_year`, `Bin_first_selected`, `N_sel_bins`, …). The donor is the first *estimated*
fleet — an `Off` fleet's slice is all `NA` and must never lead. Penalties and priors on a shared
block are accumulated once, on the lead fleet (`flt_sel_lead` / `flt_q_lead`); without that gate
they are counted once per sharing fleet.

**Worked example: GOA2018SS.** Fleets 1 and 7 share selectivity; fleets 9 and 10 share
selectivity *and* q.

## Silent-wrong-number traps

**`Index_distribution` has a second hand-synced registry.** A family added to
`index_distribution_map` (`R/0-switches.R`) must also be classified in
`.index_rows_natural_scale()` in the same file, which is what `residuals(type = "pearson")`,
`plot_index()`'s observation interval and `plot_indexresidual()` read to pick a scale. A
natural-scale family that misses it **does not error** — it silently gets the log-scale formula,
whose `sigma^2/2` is then a number the size of the index squared.
`test-index-natural-scale-downstream.R` enumerates the map and is the net.

**`nages` is a COUNT of age bins, not the oldest age.** Ages run
`minage .. minage + nages - 1`, and the array's 3rd dimension is a 1-based bin index, so age `a`
is at index `a - minage + 1`. `minage = 1` makes the two coincide — which is every bundled
dataset and all three live assessments — so a bin/age confusion passes every test and every real
model, and only shows up on the next `minage != 1` species. Write
`seq_len(nages[sp]) - 1 + minage[sp]` and mean it.

**A `data_list` element with no `write_data()`/`read_data()` support round-trips to nothing.**
The feature is then silently lossy through the standard xlsx format. This is how `index_cov` was
lost.

**A `Comp_weights` of 1 under a Dirichlet-multinomial is a starting weight of e.** That
likelihood reads the column as a log (`DM_pars_comp = exp(comp_weights)`); a multinomial reads it
as a natural-scale multiplier. `write_template()` seeds 1. `fit_mod()` reports this since 5.11.0.

**Golden reference: `BS2017SS = 10241.0304272585`,** measured on `dev` 2026-08-21 with the
`/golden-check` recipe (every fit at `newtonsteps = 3`). Other values have been quoted in this
repo's tooling and were wrong. The constants are platform-independent (verified across Windows
and macOS/arm64 to 16 digits), but on a new machine capture the baseline **before** the change
you want to test — a failure you cannot attribute is worse than no check.

## The refit path, and the eight callers

**Eight entry points re-invoke `fit_mod()` through `.refit_like()`** (`R/6-refit_like.R`), which
rebuilds the HCR / SR / M1 / growth specs from a source `data_list` and exposes each per-caller
divergence as a named override:

`retrospective()` · `jitter()` · `self_test()` · `profile()` · `run_mse()` · `remove_F()` ·
`sample_rec()` · `reweight_comps()`

`tools/verify/verify-refit-like.R` covers **six**. `sample_rec(update_model = TRUE)` and
`reweight_comps()` are **not in it** and must be checked by hand.

**`run_mse()` pins the OM's stock-recruit and suitability reference windows to the pristine
`om$`**, not the advancing `om_use$`, so the hindcast does not drift through the projection —
essential for multispecies, whose predation suitability must stay fixed. In `.refit_like()` these
are the `srr_mse_switchyr` / `srr_hat_styr` / `srr_hat_endyr` / `suit_styr` / `suit_endyr`
overrides. The EM instead advances `srr_mse_switchyr` to its current assessment `endyr` each
iteration. `tools/verify/verify-mse-hindcast-invariant.R` checks the invariant.

## The SIMULATE contract

Every observation and process error is drawn in a `SIMULATE{}` block **beside the density that
scores it** (`ceattle.cpp` sections 5.12b and 5.13, plus one per likelihood slot; the
multivariate draws live in `comp_sim.hpp`). `sim_mod()` implements no observation model in R — it
calls `obj$simulate()` once and writes the result back. **A new likelihood family therefore owes
a draw.** Three rules, with their reasons:

1. **Draw what the density assumes** — bias-correction convention and scale included. A draw that
   does not match its own density makes a self-test look biased when the model is fine.
2. **REPORT under a `*_sim` name.** TMB never clears the report environment, so a draw REPORTed
   under the observed object's own name reads back as the data. A process draw also reports a
   `*_drawn_sim` mask, since the deviation arrays are REPORTed whole.
3. **Don't draw what the model does not define.** Two densities on one latent — the AMAK/Ianelli
   stock-recruit penalty — have no distribution to draw from. Leave it, and warn. This gate is
   the model's, not the parameterization's: a DSEM under the same stock-recruit penalty has the
   same two densities on one latent and has to refuse the same way.

**One process, one draw — and the mask has to come from whichever one ran.** A DSEM's field is
drawn in R (`.dsem_draw_process()`), because the draw must be conditional on the cells the map
pins. `calculate_dsem()` also carried a `do_simulate` branch that assigned the whole field
(`x_tj = draw_tj`), and it ran *after* the substitution, so the R draw was inert — adding 100 to
every latent state changed nothing — while the template's own draw redrew the covariate columns
too. Under `family = "fixed"` those columns are the `env_data` series, so every replicate was
generated under an environment the refit never saw (a scaled covariate with sd 1 moved by up to
2.96). Two consequences worth carrying forward:

- **A second draw is invisible to a two-seed test.** "Two calls with different seeds differ" is
  satisfied by the observation draw alone, so it passes with the process draw a complete no-op.
  Test a draw by *substituting a value that cannot be confused with one* and checking it reaches
  the dynamics, at one fixed seed.
- **A draw taken in R still owes a mask.** `rec_dev_drawn_sim` is written inside the template's
  IID block, which is gated off under a DSEM, so it stayed zero and `attr(x, "process_sim")` came
  back `NULL` on exactly the models whose process *had* been redrawn. That is indistinguishable
  from nothing happening, and it silenced `compare_sim()`'s guard against measuring bias from an
  operating model that is no longer the truth.

**A covariate is data only under `family = "fixed"`.** Any other family makes it a latent state
observed with error: the state belongs to the process draw, and the observation belongs to the
observation draw (`ceattle.cpp` 5.5c) and must be written back into `env_data`, or the replicate
carries the original series and the SEM's paths look more recoverable than they are.

## MSE draws and the assessment schedule

`sim_mod()` draws **once per observation row**, so anything that changes how many rows the
operating model carries changes how far the random stream advances — and therefore every draw
after it. That makes the random stream sensitive to things that look like pure bookkeeping.

**The trap it already caused.** 5.14.0 refit the operating model only as far as `assess_yrs[k+1]`,
the *next* assessment year, because that is all the loop reads (`max_catch_hat`, for the
exploitable-biomass cap). `clean_data()` filters every data frame to `styr:projyr`, so the
shortened model carried fewer projection rows and consumed fewer draws — a count set by **when the
next assessment fell**. Two runs differing only in their assessment schedule then drew different
observation error from the *first* assessment onward, before diverging in any other way.

Measured on BS2017SS (endyr 2017, projyr 2020, `nsim = 1`, `seed = 666`), dropping the 2019
assessment moved 2019 catch by **2.1%** — a year whose advice is identical by construction, set by
the same 2018 assessment under both schedules. Per-year seeding took it to 0.53%; refitting over
the whole projection took it to 0.003%, which is optimizer noise rather than the draws. Both
landed in 5.18.0, and the 9% the shortening was worth on a Bering Sea multispecies MSE went with
it — see `inst/dev/TODO-mse-horizon.md` for how to get it back.

**What is still open.** Seeding is keyed to the **assessment** year, but one `sim_mod()` call
draws every year in the interval. A year sitting inside a longer interval is drawn under that
interval's seed, while the same year in a schedule that assessed it directly gets its own. So two
schedules realize different observation error in the years *between* assessments even where the
stock is in the same state, and that enters the next assessment as noise rather than as the
schedule. Closing it needs the draw keyed to the observation's year, not the assessment's.

**The rule.** Before changing anything that alters the operating model's row count, its horizon,
or the order of draws, ask what it does to a comparison of two schedules — not just to a single
run's reproducibility. `tools/verify/verify-mse-schedule.R`'s `gap_crn_clean` check is the gate:
it asserts a year both schedules take from the same assessment comes out the same. Note it
measures **catch**, so it does not see the between-assessment observation gap above.

## The `nages` age-vs-index class, and what closes it

`plot_ration(minage=)`, `plot_m_at_age(age=)` and `plot_m2_at_age_prop(age=)` once indexed the
age array directly while labelling the axis with an age. They now resolve the argument per
species through `.rce_age_index()` / `.rce_age_plus_index()`.

The class is closed by fixtures, not by those three fixes. `test-plot-predation-args.R` shifts a
fitted model to `minage = 2` and requires `age = 2` to draw what `age = 1` drew before;
`test-plot-smoke.R` runs the rest of the exported plotters against a `minage = 3` model. Without
a `minage != 1` fixture the whole suite passes on a bin index read as an age, because every
bundled dataset and all three live assessments have `minage = 1`.

## The guards are not themselves guarded

Two of this repo's most valuable checks were, until 5.16.0, guarded off in a way nothing
verified. Both are the same shape: a speed optimization that silently removes coverage.

**The golden regression ran in NO automated job.** `test-golden-regression.R` carries both
`skip_on_cran()` and `skip_on_covr()`. `R-CMD-check` sets `NOT_CRAN=false` deliberately (to keep
the multi-OS matrix fast), so the first fires; `test-coverage` runs under covr, so the second
does. The committed backstop for "a numeric change needs golden equivalence" therefore only ever
ran when someone typed `/golden-check` locally. Its measured cost of `0` in
`tools/ci/coverage-costs.tsv` is the tell — four phased fits do not cost nothing. The
`deep-checks` workflow now runs it nightly, and **asserts the run produced assertions**, because
a skipped file otherwise reports as a pass.

**`NOT_CRAN=false` was set through `$GITHUB_ENV` and did not hold reliably.** On `7e16fa73` the
skip fired (testthat recorded `test-selectivity-catchability.R:3:1`); on `7611bd1c`, differing
only in a Markdown file, the same file executed and the worker died with exit code
`-1073741819` = `0xC0000005`, a Windows **access violation**. Set it as step-level `env:` on the
`check-r-package` step instead — step env beats both the job env and anything a prior step
exported, and it is printed in the step's own env block, so the log says which mode ran. 1 of 7
recent Windows runs slipped; when it does not slip you learn nothing.

**An access violation is memory corruption, not a bad optimum.** A fit in a different basin
gives a huge gradient and a failed `sdreport` — not a fault. The model is built
`safebounds = FALSE`, so an out-of-range access is a silent write into adjacent memory; that is
precisely what the 5.15.0 `age_hat` overflow was. Build with `RCEATTLE_SAFEBOUNDS=true` to turn
that into an R error naming the object, and run `tools/verify/verify-safebounds.R`. That harness
carries a **positive control** — it deletes the artifacts, rebuilds, and asserts
`-DTMB_SAFEBOUNDS` is in the compile line — because `pkgload` only recompiles when sources
change, so setting the flag alone guarantees nothing and a clean result against a stale
unchecked `.so` is meaningless. As of 5.16.0 it reports no violation over five configurations
(ragged comps both directions, joint sex, predation arrays, the CI crash config, `sim_mod()`
draws) on macOS. **That is not a clearance** — the fault is intermittent and was seen on Windows.

Two further sharp edges found while building that harness, both worth knowing:
`src/TMB/compile.R` resolves `ceattle.cpp` against the **working directory**, so running it from
the repo root skips the compile, prints its messages, and exits 0; and `on.exit()` registered at
a script's top level has no function frame and never fires, which is how the harness first left
a slow bounds-checked `.so` in the tree for every later fit to use.

## Coverage gaps in the golden check

`/golden-check` fits four models and diffs the result. It does not reach:

- **The refit paths.** Eight entry points re-invoke `fit_mod()` through `.refit_like()`. Use
  `tools/verify/verify-refit-like.R`, which covers six of them —
  `sample_rec(update_model = TRUE)` and `reweight_comps()` are **not** in it.
- **The simulation draws.** `tools/verify/verify-sim-*.R`.
- **Any figure.** A plotter can change every number it draws and stay green. The net is
  `tests/testthat/test-plot-*.R` plus a before/after `ggplot_build()` diff at the merge base.
- **Linkage rows.** The four golden models carry none. The constructive linkage tests
  (`test-dynamics-recruitment-linkage.R`, `test-linkage-random-effects.R`) are the net.
- **Selectivity normalization.** No bundled dataset enters the normalization block in
  `selectivity.hpp` — every one has `Sel_norm_bin` all-`NA` except one-sex `Atka2022` — so the
  four reference fits pass trivially through any change to it. `GOA2018SS` *is* reachable
  (**species 2**, arrowtooth, is the two-sex one -- `nsex = c(1, 2, 1)`; 9 is the `Fleet_code` of
  `ATF_bottom_trawl`, not a species index) once you set a named `Sel_norm_bin`; that is exactly how
  `681f3ba0` shipped with `test-selectivity-normalization.R` failing.
- **Growth.** `growth_fun` is unset on BS2017SS and GOA2018SS, so `growth.hpp`'s SDAA path —
  `sd_plus_group` included — is never exercised. Covered constructively by the
  `test-growth-*.R` files.

**`goa_ms` (fixed-M GOA multispecies) sits on a flat likelihood ridge:** the same objective at
different `par`/`ssb` across *different* code, though deterministic on same-code re-runs. Judge
it on `obj`/`jnll`, not `par`/`ssb`.
