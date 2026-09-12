# Cleanup backlog

The 72 `TODO` / `FIXME` markers in the source, triaged. **Add to this file; don't fix these
unasked.** Fix the one in the file you were already asked to touch, in the same commit.

**Cite the marker text, not a line number.** These references have gone stale three times --
twice when roxygen was added above them, once from an unrelated `ceattle.cpp` edit. Grep for the
quoted FIXME text instead; it moves with the code.

Three tiers: a **known defect** is a wrong answer waiting for the right input and should become a
GitHub issue; a **design note** is a wish, not a bug; `TODO(review)` is a deliberate convention
marking a judgement call for Grant, and is never resolved by an agent.

62 remain as of 5.29.0. Counts by area: `src/TMB/ceattle.cpp` 24 ·
`src/TMB/predation.hpp` 5 · `src/TMB/Dev/caal.hpp` 5 · `R/10-run_mse.R` 4 ·
`src/TMB/growth.hpp` 3 · `R/3-build_map.R` 3 · `R/9-retro_and_jitter.R` 3 ·
`R/0-rceattle_class.R` 3 · rest 1–2. Re-derive with
`grep -rnE 'TODO|FIXME' R/ src/TMB/ | grep -v 'todo <-' | grep -v 'TODO-'` -- the `-E` is
needed for the alternation, the first filter drops a variable in `R/6-process_residuals.R`, and
the second drops pointers to `inst/dev/TODO-*.md` notes, which are not markers.

---

## Tier 0 — known defects, with the input that triggers them

These say, in the source, that the code is wrong under a stated condition.

Three further defects of the same class were found reviewing the fixes below, and are resolved
in 5.13.0 alongside them. None carried a marker, which is why none appeared in this file: they
are what the markers pointed *near*, not what they said.

Rows marked **Open** were found reviewing PR #143 (2026-09-12) and have not been reproduced with
a fit. None carries a source marker.

| Where | Condition | Consequence |
|---|---|---|
| `R/3-build_map.R` (`# Age-independent scalar`) | **Open.** `estDynamics = 3` in a multispecies model (`msmMode > 0`) | Gated on `estDynamics[sp] == 2 \| msmMode != 0`, so under `msmMode > 0` columns 2..`nages` of `log_pop_scalar` are mapped out for every species, and the age-specific scalar of `estDynamics = 3` is never estimated (it stays `exp(init)`, 1 without `inits`). Gating on `estDynamics[sp] == 2` alone is safe for codes 0-2. Check the Hessian before freeing them: they are informed only through predation and that species' index. The `# Age-dependent scalar` branch has the same `\|` but maps only padding beyond `nages`, which is harmless. |
| `R/3-build_map.R` (`# Don't estimate the scalar`) | **Open, question for Grant.** `msmMode = 0` | Every column is mapped out, so `estDynamics` 2 and 3 behave as 1 in single-species models. The 5.30.0 NEWS says so; the schema and `?BS2017SS` do not. |
| ~~`JNLL_Q_PRIOR` against a q linkage intercept prior~~ | ~~`Catchability = "Estimated-with-prior"` plus a q linkage `(Intercept)` prior on the same fleet~~ | **Resolved in 5.33.0**: both penalized that fleet's log q, so the prior counted twice. `.check_q_linkage_support()` now refuses the pair, taking fleet 1 for a row with no fleet, as the template does. `test-linkage-double-prior-guards.R`. |
| ~~`JNLL_M_PRIOR` against an M1 linkage intercept prior~~ | ~~`M1_use_prior = TRUE` with `M2_use_prior = FALSE`, plus an M1 linkage `(Intercept)` prior~~ | **Resolved in 5.33.0**: both penalized that species' log M1. `.check_M_linkage_prior()`, run from `fit_mod()`, now refuses the pair, taking species 1 for a row with no species. `test-linkage-double-prior-guards.R`. |
| `srr_terms_on` | **Open, low.** A recruitment linkage intercept prior on an `estDynamics > 0` species | The gate covers the stock-recruit prior and the curve penalty, not the linkage-prior loop. `rec_pars` is mapped out for such a species, so the prior only adds a constant: it moves the objective and `JNLL_LINKAGE_PRIOR` (likelihood tables, AIC), but no estimate. Not yet checked: slope rows on such a species, which `build_map_fixed_natage()` does not map out. |
| `// Input SB0 (if running in multi-species mode)` | **Open.** `msmMode = 0` with `estDynamics > 0` | Equilibrium `SB0`/`SBF` are still built on placeholder recruitment; from 5.30.0 only the dynamic runs keep the input numbers. HCRs 5, 6 and 7 read `SB0` when `DynamicHCR = FALSE` (5 also reads `SBF`), and with `DynamicHCR = FALSE` `ssb_depletion` and `biomass_depletion` divide by `SB0` and `B0` under every HCR. Projected numbers come from `NByageFixed` and are right; the reported depletion, and F and catch advice under those HCRs, are not. |
| ~~`int spawn_yr = yr - minage(sp);`, `int rp_yr = yr - minage(sp);`~~ | ~~`minage = 0` with a stock-recruit curve~~ | **Resolved in 5.33.0** by refusal: each lag read the same year's `ssb`/`SB0`/`DynamicSB0` before it was accumulated, so the curve gave R = 0 in the hindcast (`srr_fun >= 2`), the reference points (`srr_pred_fun >= 2`) and a curve-based projection. `data_check()` now refuses `minage = 0` with any curve; computing SSB before recruitment would allow it, but the GOA and AI cod `minage = 0` models fit no curve. `test-data-check-srr-guards.R`. |
| ~~`src/TMB/ceattle.cpp` (male slot writes)~~ | ~~every species one-sex (`max_sex == 1`)~~ | **Resolved in 5.13.0**: ten lines wrote sex index 1 unconditionally, but arrays are dimensioned `max_sex`, so that index does not exist when no species has two sexes. Value written is 0 (`sex_ratio` is set to 1 for a one-sex species first), but the write is out of range and lands on `(sp, 0, age + 1, yr)` — the next age — surviving only because the age loop overwrites it. Fires on BS2017SS and BS2017MS every evaluation. Reproduced with `TMB::compile(safebounds = TRUE)`, which raises Eigen's range assertion; guarded, the fit is clean at an unchanged 1537036.287629372. `test-dynamics-sex-index-bounds.R`. |
| ~~`R/1-data_check.R` (no `comp_data$Sex` check)~~ | ~~`Sex` 2 or 3 on a one-sex species~~ | **Resolved in 5.13.0**: `M1_base`, `weight` and `ration_data` are all checked against `nsex`; composition was not. Two registries disagree on what "joint" means — `check_composition_data()` uses `nsex == 2 & Sex == 3`, the template uses `flt_sex == 3` alone — so a joint row on a one-sex species was sized at `nages` and written to `nages * 2`, corrupting the NEXT observation's predicted composition and its likelihood. Refused at the boundary rather than reconciled in the template. `test-data-check-comp-sex.R`. |
| ~~`src/TMB/ceattle.cpp` (reference-point recruitment arms)~~ | ~~a stock-recruit curve with `proj_mean_rec = TRUE` (the default)~~ | **Resolved in 5.13.0**: the mean-rec arm required `proj_mean_rec == 1 & srr_pred_fun < 2` and the curve arm required `proj_mean_rec == 0`, so that combination matched neither and reference-point recruitment stayed 0 after year 1. `SB0` became the initial cohort decaying (3.344 → 1.230 over six years), and `SB0` in the terminal year is what HCR 5 and 6 read as the depletion reference — so perceived depletion and the resulting catch advice were both wrong. The curve arm now fires whenever a curve exists; the projection switch is read separately and is unchanged. `build_srr(proj_mean_rec =)` was also documented backwards. `test-dynamics-refpoint-mean-rec.R`. |

---

**All nine are resolved as of 5.13.0.** The table is kept struck through rather than deleted:
each row records what the defect actually was once reproduced, which in six of the nine differed
from what its FIXME claimed. Add new rows above it.

| Where | Condition | Consequence |
|---|---|---|
| ~~`src/TMB/ceattle.cpp` (`will bomb if minage > 1`)~~ | ~~`minage > 1`~~ | **Resolved in 5.13.0**: it did not bomb, it read adjacent memory. `ssb(sp, yr - minage(sp))` went negative for the first `minage - 1` years and Eigen does not bounds-check in a release build. Measured (BevertonHolt, nages 5): R came back `8.6e-314` with `is.finite()` TRUE, objective 14407.38 → 15162.60 → 17532.15 across minage 1/2/3. Those years now take `R_init * exp(rec_dev)` -- equilibrium recruitment at `F = Finit`, already what year 0 uses -- following Stock Synthesis's equilibrium-plus-early-devs treatment of the pre-start period (WHAM fixes the lag at one year instead). NOT `R0`: under a stock-recruit hindcast `build_map()` maps the mean-recruit parameter out and only `R0[, 1]` is overwritten with the derived `(alpha - 1/SPR0)/Beta`, so a guard reading `R0[, yr]` gets `exp(9) = 8103.08`. The first fix did exactly that, giving 14407.38 &rarr; 20421.68 &rarr; 30920.21; corrected it is 14407.38 &rarr; 13959.97 &rarr; 13587.68. Expected recruitment (`R_hat`) shares the anchor, or the stock-recruit penalty reads the mean-versus-curve level gap as signal (0.83 nats/year under the Ianelli configuration). All four sites guarded; `minage = 1` cannot fire, golden unmoved. `test-dynamics-recruitment-minage.R`. |
| ~~`src/TMB/ceattle.cpp` (`will blow up if nlengths is less than nages`)~~ | ~~`nlengths < nages`~~ | **Resolved in 5.13.0**: `age_hat`/`age_obs_hat` are written at AGE indices (to `nages*2` joint-sex) but were sized `= comp_obs`, whose width is the workbook's `Comp_` columns — `nlengths` for a length-only model. Silent out-of-bounds WRITE, not a crash. Both now sized from the widest age index the model's own `comp_ctl` says will be written -- `nages`, or `nages*2` for a joint-sex row -- and never narrower than the observations. Reserving `max_age*2` unconditionally (the first fix) also closes the overrun but widens every model's REPORTed `age_hat`: BS2017SS went 25 &rarr; 42 columns, 17 permanently zero, and assessment scripts `cbind` that array. No end-to-end reproduction is possible from R (an unchecked write is not observable), so `test-composition-age-hat-width.R` pins the invariant and the absence of `= comp_obs`. |
| ~~`R/10-run_mse.R` (`does not work for assessments that don't occur annually`)~~ | ~~assessment interval ≠ 1 year~~ | **Resolved in 5.13.0**: only the scalar-`cap` branch was affected. `dat_fill_ind` spans the whole interval, so `sum(Catch[dat_fill_ind]) > cap` held a multi-year total to a one-year ceiling — at `assessment_period = 2`, roughly halving projected catch. Now applied per projection year. The same line carried a second and larger defect: `ifelse()` returns the shape of its length-1 test, so whichever branch was taken was truncated to its first element and recycled across every row. Two species at 80/20 t came back 40/40 against a 50 t cap (80 t total, over the ceiling) and 80/80 against a non-binding 500 t cap (160 t total, catch invented). **Any stored MSE using a scalar `cap` moves, at every `assessment_period` including 1, binding or not** -- only an exactly equal split was safe. The species-specific vector branch was always per row. `test-mse-cap-and-hcr2-threshold.R`. |
| ~~`R/10-mse_summary.R` (`will bug if not survey`)~~ | ~~a species whose fleets are not all surveys~~ | **Resolved in 5.13.0**: not a defect. `spp_rows <- which(flt_spp == sp)` was assigned once and read by nothing, in this file or anywhere else in the package. The FIXME speculated about code that did nothing; both lines are gone. No behaviour change. |
| ~~`R/3-build_map.R` (`QAR1 is inert`)~~ | ~~**`Catchability = "AR1"`** (the QAR1 form, Rogers et al. 2024)~~ | **Resolved in 5.12.0**: `data_check()` now errors on it, so the branch is unreachable. It was inert — the deviate map is gated on `Time_varying_q %in% c("IID","AR1","RandomWalk")`, but under `Catchability = "AR1"` that column holds an `env_data` **column index**, not a mode, so `index_q_dev` stayed mapped out and q was constant. Not repaired: the Rogers form is implemented correctly by a q linkage (`ar1(1 \| Year)` with `observe`), which GOA pollock 2025 runs. The dead `build_map()` branch is deleted in 5.13.0. Code 6 **stays in `q_map`**: `validate_switches()` runs before `data_check()`, so dropping it would replace that migration message with a generic "invalid value" for exactly the workbooks that need the recipe (GOA pollock 2024/2025 still carry a 6). These are two different switches sharing a string; an earlier draft of this file named the wrong one. |
| ~~`src/TMB/ceattle.cpp` (`caal_ll_type`)~~ | ~~`CAAL_distribution = "MultinomialAFSC"`~~ | **Resolved in 5.13.0**: implemented, as the catch-sigma case of this shape was in 5.12.0. The AFSC multinomial pseudo-likelihood is a published AMAK form already implemented for age comps, so extending it to CAAL was mechanical rather than inventing a likelihood. Verified against the form computed by hand from the reported CAAL proportions (2.8e-14). The `test-schema-cpp-dispatch.R` exemption is removed. `test-likelihood-caal-afsc.R`. |
| ~~`R/3-build_map.R` (`will fail if random_sel = TRUE?`)~~ | ~~`random_sel = TRUE` + `Time_varying_sel = "Block"`~~ | **Resolved in 5.13.0**: confirmed real, and worse than "will fail". The block parameters live in `log_sel_slp_dev`/`sel_inf_dev`, and `fit_mod()` declared those arrays random unconditionally — but the template scores selectivity deviates only for `IID`/`AR1`/`RandomWalk`/`RandomWalkAscending`, so blocks were Laplace-integrated against **no density**, `sel_dev_log_sd` mapped out so there was no variance either. Measured: 8 parameters random, `JNLL_SEL_DEV` identically 0, objective `NaN`, real fit dead with TMB's `NA/NaN gradient evaluation`. Now refused with a message naming the fleets and the way out. `test-selectivity-random-sel-block.R` carries the reproduction and a drift guard pinning `Block` as the only mode the template leaves unscored. |
| ~~`R/6-fit_mod.R` (`swallows EVERY warning build_map() raises`)~~ | ~~any~~ | **Resolved in 5.13.0**: the comment named the wrong warnings — the shared-block ones it cited are raised by `data_check()`, not `build_map()`. What was actually swallowed was `build_map()`'s own set (M1 sex mismatch, selectivity-form incompatibilities), each of which changes what is estimated. Now de-duplicated via `withCallingHandlers()` and passed through, so `.refit_like()`'s per-peel re-entry prints each distinct warning once instead of hundreds of times. |
| ~~`R/10-mse_summary.R` (`EM uses fixed-depletion proxy for HCR 2`)~~ | ~~HCR 2~~ | **Resolved in 5.13.0**: in single-species mode the HCR 2 arm now reads `ssb_limit_thresh()`, the helper the operating model already uses, so both sides of the cross-tab score one criterion (absolute `0.5 * SBF`) and branch on the same scale flag. `Plimit` is NOT the answer -- `build_hcr()` defaults it to 0, so reading it reports a default ConstantF run as never overfished; the first fix did that. Under `msmMode > 0` both sides still fall through to `Plimit`, which is the operating model's own multispecies rule, so they agree and it is left alone. `test-mse-cap-and-hcr2-threshold.R`. |

## Tier 1 — stated limitations, currently by design

Not bugs, but they bound what the model can be asked. Worth documenting in a vignette rather
than fixing.

- **Forecast growth is ignored** by the retrospective and MSE projection paths
  (`ignores forecasted growth`, twice in `R/9-retro_and_jitter.R`, once in `R/10-run_mse.R`) —
  the terminal-year growth is carried forward.
- **Projection quantities are held at the terminal hindcast year** (`R/10-run_mse.R`,
  `assuming same as terminal year of hindcast`).
- **`ration_data` is sized for the hindcast only** (`R/5-rearrange_data.R`,
  `Change for forecast`).
- **SPR reference points**: `sex_ratio` is an input rather than estimated for two-sex models,
  and the M used is the terminal-year value (`src/TMB/ceattle.cpp`, `rates for a reference point
  are the terminal hindcast year's`). The `sex_ratio` marker went in 5.24.1, when SPR was
  corrected to apply only the recruitment split `sex_ratio(sp, 0)` to a two-sex species
  (`female_split`); the ratio itself is still read from data.
- **Linkage random-effect priors are penalties, not proper densities** (`src/TMB/ceattle.cpp`,
  `FIXME(jacobian)`, twice). The sigma and rho priors sit on the natural scale without the
  Jacobian of the `log` / `rho_trans` transform. Fine as a penalty under maximum likelihood;
  under a Bayesian (`tmbstan`) run the stated density is not the prior actually applied. A
  `lognormal` sigma prior is exempt; rho has only normal and beta families.

## Tier 2 — design notes and refactor wishes

Cleared in 5.14.0 except where noted. As in Tier 0, three of these were not what their marker
said, so each struck row records what it actually turned out to be.

- ~~**Split `R/0-build_srr_and_M.R`**~~ — **Done in 5.14.0**, but not as described. The file was
  1,497 lines and **52** top-level objects, not 29, and the three-way srr/M1/growth split named
  here would have stranded 612 lines (41%): catchability, selectivity and composition linkage
  machinery that touches none of those three. Split six ways instead, one file per process, all
  keeping the `0-` prefix. `.coerce_switch_arg` went to `0-switches.R` beside `.canon_switch()`
  (it is generic over `srr_fun`/`M1_model`/`sd_plus_group`); `.default_stratum`,
  `.resolve_auto_by` and `.stamp_param` went to `0-build_linkage.R`, which already called the
  last of them. The roxygen block for `.check_q_linkage_support` was bound to
  `.message_auto_fleet_linkages` and is reattached. Pure relocation: 441 object bodies unchanged,
  and the multiset of non-blank lines across `R/` is identical.
- ~~`R/5-rearrange_data.R` (empty comp/CAAL frames, `age_error`)~~ — **Done in 5.14.0.** The comp
  and CAAL blocks were the same five lines twice; now one `.normalise_rows()`. The zero-row guard
  stays — `t(apply())` does not preserve the shape of a matrix with no rows. `age_error` keeps its
  `as.data.frame()` coercion, which is load-bearing because the loop mixes `$` and positional
  `[i, ]` access, and its `1:nrow()` becomes `seq_len()`.
- ~~`R/2-build_params.R` (`variance and AR1 parameters`)~~ — **Not work.** Both already exist
  directly above the marker (`log_sigma_linkage`, `trans_rho_linkage`), and the `else` branch
  zeroes every name the `if` branch sets. A leftover placeholder at a section boundary; deleted.
- ~~`src/TMB/ceattle.cpp` (`can probably outside iter loop`)~~ — **Not work.** Section 5.12 is
  already ~260 lines above the only `iter` loop. The comment now states what the placement means,
  including that the forecast arm is provisional and recomputed in section 6.7.
- ~~`R/0-clean_data.R` (`may be redundant now?`)~~ — **Half right, and not the half that matters.**
  The template does derive `SB0`/`B0` itself, and reads `MSSB0`/`MSB0` only to overwrite that under
  `msmMode > 0`. But neither has `read_data()`/`write_data()` support, so no workbook can supply
  them and this default is the only thing that creates the required `DATA_VECTOR`s. Kept, with the
  reason and with `999` named as the placeholder `fit_mod()` fills in. **Correction, 5.15.0:**
  section 10.2 filled it into `data_list_reorganized` only, so the *returned* `data_list` kept the
  999 and every refit off a fitted object re-entered the template with it. Fixed by carrying
  `MSSB0`/`MSB0` onto `mod_objects$data_list`.
- ~~`R/3-build_map.R` (`add checks for surveys sel sigma`)~~ — **Done in 5.14.0.** Real: fleets
  sharing a `Selectivity_index` estimate one `sel_dev_log_sd` between them, and a differing
  `Time_varying_sel_sd` was reconciled silently. Not to the first member's value, which is the
  intuition to unlearn: TMB's `updateMap()` collapses a shared parameter with
  `tapply(par, map, mean)`, and this one is held on the log scale, so the group starts at the
  **geometric mean** of its estimated members' values — `sqrt(0.3 * 0.7) = 0.4583` for a
  two-fleet group at 0.3 and 0.7, a value neither row asks for. `.warn_shared_dev_sd()` reports
  it, once per group, and runs at the END of `build_map()`: `build_map_f_and_data_weights()`
  maps the parameter out for `Off` fleets and `build_map_fixed_natage()` for a fixed-dynamics
  species, both after the sharing pass, so a check placed inside
  `adjust_map_shared_params()` counts fleets that end up estimating nothing.
- ~~`R/3-build_map.R` (`add checks for surveys q sigma`)~~ — **Done in 5.14.0, after fixing the
  reason there was nothing to check.** `index_q_dev_log_sd` was mapped out for every fleet and
  never turned back on, so no shared group could discard it. That was itself the defect:
  `random_sel` frees `sel_dev_log_sd` alongside the selectivity deviates it integrates out, but
  `random_q` integrated the catchability deviates out and left their sd fixed at
  `Time_varying_q_sd`. Now symmetric, so `random_q = TRUE` estimates it — **and any fit using that
  flag moves.** With the sd estimable the shared-group copy is meaningful, so a shared
  `Catchability_index` goes through the same `.warn_shared_dev_sd()` check, on the same
  geometric-mean footing described in the row above.

  Caveat, measured rather than assumed: on a 40-year index with the observation sd FIXED and q
  deviations injected at a true sd of 0.4, the marginal MLE pins to its lower bound at observation
  sd 0.1, 0.3 and 1.0, and reaches only 0.06 at 0.05. A short or noisy index can therefore return
  an sd that reads as a constant q. That is a diagnostic to check, not a reason to withhold the
  parameter -- the estimate and its gradient show it, and whether the series informs it is the
  assessor's call. `index_q_log_sd`, the prior sd on q itself, stays fixed: estimating the width of
  one's own prior is not meaningful.

Still open. No user-visible consequence; do them opportunistically.

- `R/10-run_mse.R` (`extract run_one_sim as a top-level internal helper`) — 455 lines, 18 free
  variables. Three hazards, all verified: `%!in%` is defined **inside `run_mse()`** and is not a
  package-level object; the `<<-` in the OM-no-F handler currently resolves to the closure's own
  `sim_list` and would walk to the namespace if the handler moves out; and `estimate_mode_base`
  and `sim_dat` exist as both `run_mse()` locals and closure-local re-assignments, so the closure
  does not depend on the outer copies. Note also that the closure's environment is passed to
  `.parallel_lapply()` for the PSOCK export, so the free-variable set is load-bearing for
  parallelism, and that `test-mse-cap-and-hcr2-threshold.R` greps this file for literal source
  strings and must move in lockstep. Needs `verify-mse-repro.R` and `verify-mse-om-horizon.R` as
  before/after digests.
- `R/3-build_map.R` (`use formula`) — retiring the `Time_varying_q` overload that holds
  comma-separated `env_data` column indices in favour of the q linkage. Deprecates a public switch
  path, so it needs a shim, not an opportunistic edit.
- `R/3-build_map.R` (`may want sex-varying?? Hard to estimate`) — `M1_rho` already has a sex
  dimension, but the mode-6 branch offsets by `nspp` to avoid colliding with the mode-4/5
  `sp`-valued rho, and going sex-varying needs that counter scheme reworked. The markers say the
  obstacle is estimability; decide that first.
- `src/TMB/ceattle.cpp` (`fit the window form directly when a fleet needs it`) — needs `t1` and
  `D` per fleet plumbed through to the template, **plus seasonal dynamics that do not exist**: the
  annual recursion spreads F evenly across the year, so the window predictor would be inconsistent
  with the dynamics fitted against it. The derivation above the marker already quantifies the
  current approximation at ~1.5% trend error over F 0.05–0.8, against −29% to +33% for the
  snapshot it replaced.
- The `logH_*` / `H_4` / `log_gam_*` markers belong to the stubbed Kinzey-Punt predation forms
  (`msmMode` 3–9) and the gamma predator selectivity. They are pinned as stubbed in
  `tests/testthat/test-schema-registries.R`; leave them until that work is picked up.
- `src/TMB/ceattle.cpp` (`penalize every selectivity deviation rather than a sub-range`) —
  would pin the unidentified directions and drop the year/bin indexing of the deviation penalty.
  It would **not** retire the four columns the marker names: `Sel_cap_bin` holds the
  NonParametricRPM curve flat past a bin, `Sel_start_year` builds the curve from the base
  coefficients through that year and pins the random walk's level in `build_map()`, and
  `Sel_pen_first_bin` / `Sel_pen_last_bin` bound the shape penalty, not the deviation penalty.
  Moves every fit with penalized deviations, so it needs `/golden-check`.
- `R/0-osa_data.R` (`switch_check() does not run`) — the comment above it says `comp_offset` is
  filled by `switch_check()`; on the exported `rearrange_data()` path it and the
  `bias_adjust_*` scalars are filled in `build_osa_data()`. Reword only; no behaviour.

- **Single-species hindcast-curve projection double-counts the SSB drop.** `sample_rec()` and
  `retrospective()` set the deviation to `log(mean(...R...)) - log(...R0...)`. Under a curve R0
  is unfished recruitment, so the curve applies the SSB drop again: 8.6% below its mean at
  h = 0.8 and 40% SB0, 27% at h = 0.5. It hits `run_mse(sample_rec = FALSE)` OM years and the
  dynamic reference points. The multispecies branch's `log(mean(R / R_hat))` is the model.
- **Dynamic B0 under the penalty form** ("Dynamic reference points (Includes annual recruitment
  deviation: pass rdev)") applies `exp(rec_dev)` = R/R0, not R/R_hat, to the curve.
- **`.map_switch()` passes a factor through**, so a factor `srr_est_mode` skips `build_srr()`'s
  checks and fits as its level code.
- **The refit warning for retired srr codes is hidden** by the `suppressWarnings()` wrapped
  around `.refit_like()` in `retrospective()`, `jitter()`, `profile()` and `self_test()`.
- **A one-year retrospective peel** averages over that year, though its warning says "after the
  first".

## `TODO(review)` — Grant's calls, not an agent's

Six, each a judgement about what the right behaviour *is*:

- `R/0-rceattle_class.R` (`osa_residuals("all") includes diet`) — whether
  `residuals(source = "all")` should include diet too.
- `R/0-rceattle_class.R` (`held-out rows (Year <= 0) with a positive observation`, twice) — how
  those rows should be treated.
- `R/6-fit_mod.R` (`a user-supplied NA`) — what an `NA` bias-adjustment should mean.
- `R/7-plot_osa.R` (`process-residual objects`) — how `process_residuals()` output should be
  plotted.
- `src/TMB/growth.hpp` (`this branch (and its Richards mirror below) tests`) — see the file.

A seventh, on multispecies SBF, was restated as a known limitation in `5d423172`; it is settled
under "Deliberately not changed".

## Deliberately not changed

- **Multispecies SBF sits on the projection's realized M2** (`src/TMB/ceattle.cpp`, `Multispecies:
  M_at_age carries the projection's realized M2`). Under `msmMode > 0` it is reported but nothing
  live reads it. The one rule that does, NPFMC (5), is refused there along with 4 and 7;
  ConstantFSSB tunes realized SSB against SB0, CMSY reads depletion, and `mse_summary()` reads SBF
  only when `msmMode == 0`. Allowing HCR 5 in multispecies mode would make this a defect;
  `test-switches-hcr-multispecies.R` pins the refusal.
- **Non-parametric growth** is declared and calls `error("not yet implemented")`.
- **The `msmMode` 3–9 Kinzey-Punt branches are not declared at all** -- the whole block in
  `predation.hpp` is inside a `/* ... */`, so there is no dispatch, live or erroring. The live
  modes are handled by `if (msmMode == ...)` in `ceattle.cpp`.
  `test-schema-cpp-dispatch.R` pins both, and pins the absence of the switch.
- ~~`flt_sel_ind`~~ — removed in 5.12.0. It was computed from `Fleet_code` on every fit and read
  by nothing.
- **The `dmultinom_osa()` renormalization under `Comp_distribution` case 0.** Fitting routes
  through `dmultinom_osa()`, which renormalizes `p`, so the *reported* multinomial NLL carries a
  per-row constant the old `dmultinom()` did not. The gradient and the MLE are unchanged, so this
  is a reporting discrepancy, not a wrong fit. Correcting it would move the golden reference
  numbers for a cosmetic gain; reviewed and left in place 2026-08-23.
