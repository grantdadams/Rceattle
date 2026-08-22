# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

Branch `api-names-docs-and-guards` — PR A of three, targeting `dev`. Carries the Stage 6 API
work (done, both numeric gates exact) plus the Claude-tooling, drift-guard and documentation
work. PR B (schema authority, the `nages` figure fix) and PR C (`write_template`, model-switch
table) stack on it.

## Done & verified

- **`dev` merged to `main` in PR #106.** `main` is at 5.8.1; `dev` carries everything since —
  the linkage grammar, column schema,
  `build_data()`/`model_config()`, `save_config()`/`load_config()` + `fit_mod(config=)`, the
  `JnllRow` enum, `build_growth(sd_plus_group=)`, the `mse_summary()` per-entity reshape, the
  `.refit_like()` collapse, `reweight_comps()`, and the recruitment / stock-recruit work.
- **5.11.0 API names** — the fitted model is `object` across ten diagnostics, old names shimmed
  and silent. Verified 2026-08-22: suite 6925 pass / 0 fail on the committed state; golden 4-model capture
  `identical()` TRUE with every field at `0.000e+00`; `verify-refit-like.R` bit-identical across
  all 32 fits; ecosystem sweep found 15 old-name call sites, all shimmed, and zero calls passing
  both spellings.
- **The golden reference on `dev` is `ss = 10241.0304272585`** (with `ms = 10267.2478324443`,
  `goa_ss = 12868.0052289274`, `goa_ms = 12932.7931701136`). Measured, not quoted. A different
  value that had been recorded in the local agent file was wrong and has been corrected.

## Known flags

- **5.10.0 moved three predation figures' numbers** — `plot_m2_at_age_prop()` (a share now, not
  a contribution), `plot_ration()` (× average numbers-at-age, not biomass-at-age) and
  `plot_b_eaten()` (million mt, so `p$data` moves by 1e6). Any figure regenerated from them
  differs from earlier runs of the same model. `plot_selectivity()` also renames `p$data$Age` to
  `Bin`.
- **Result-changing changes on this line that are not labelled breaking** — the mode-5
  selectivity penalty fix (GOA Pacific cod SSB 2050 −14.1%), parameter bounds previously applied
  to the wrong parameters, composition weights warm-starting from `inits`, failed `run_mse()`
  simulations returning only a marker, the `mse_summary()` reshape, the recruitment fixes
  (`initMode = 0` random effects, the α-seeding fix, the Ianelli steepness prior), and
  `sim_mod()` drawing the index under the fleet's own `Index_distribution`. **A model carrying
  GOA numbers forward needs a refit.**
- The open `nages` age-indexing defect in three predation plotters is recorded in
  `inst/dev/TRAPS.md`, not here — a known unfixed defect should not live in the file this
  command rewrites every session. Scheduled for PR B with a `minage != 1` fixture.

## Blocked

Nothing.

## Resume here

PR A: finish the Stage 1 drift guards and the Stage 2 documentation sweep, then open the PR.
The `write_template()` coverage test in Stage 1 is deliberately `skip()`ped until PR C.

## Older paused work

A multi-PR accessibility / code-review refactor on branch `accessibility-and-code-review`. Its
plan is at `~/Downloads/HANDOFF-accessibility-refactor-implementation.md` (outside the repo, so
it does not survive a clone — ask Grant for it). Read it before resuming; do not start fresh.
