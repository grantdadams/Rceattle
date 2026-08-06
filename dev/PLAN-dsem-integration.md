# Plan — integrate DSEM onto the v5.0 (data-workflow) codebase

## Context & goal

`dev-DSEM` (44 commits since the 4.5.0 base) adds a **Dynamic Structural Equation Model on
recruitment**: `build_DSEM()` makes the recruitment deviations the latent states of a
`dsem`-package GMRF, so environmental covariates + lagged/simultaneous effects drive
recruitment. Meanwhile `dev-data-workflow` (294 commits, now merged into `dev`, v5.1.0)
rewrote the surrounding machinery — linkage grammar, column schema, `build_data`/config,
a C++ legibility pass (renamed `ceattle_v01_11.cpp` → `ceattle.cpp`, added the `JnllRow`
enum + recruitment/linkage accumulators), and flattened the test layout.

**Goal:** land DSEM on top of the v5.0 codebase so the model has *both* the linkage grammar
and DSEM, converging bit-identically to `dev-DSEM` on DSEM models and to `dev`/golden on
everything else.

**Branch:** `dsem-v5-integration` (off `dev-DSEM`, merged `dev-data-workflow`).

## Status

- [x] **Tier 0 — mechanical merge (done, commit `b1275313`).** Resolved all 25 conflicts to
  the **v5.0 base**; carried DSEM's standalone files (`R/0-build_DSEM.R`, `src/TMB/dsem.hpp`,
  `examples/DSEM_example.R`, `tests/testthat/test-dsem-recruitment.R`, vignettes,
  `man/*.Rd`); moved `check_dsem_spec()` from `convergence.R` into `R/0-build_DSEM.R`;
  DESCRIPTION `+dsem (== 3.0.0), +utils`; NAMESPACE `+build_DSEM/build_dsem_objects/
  check_dsem_spec`; untracked the `dev/` scratch scripts. **Loads + compiles at v5.0; DSEM
  present but INERT** (`fit_mod()` has no `dsem=` arg; `ceattle.cpp` does not `#include
  dsem.hpp` or call `calculate_dsem`).
  - **Auto-merge leak caught & reset (commit after `b1275313`).** Several v5.0 files
    auto-merged with dev-DSEM edits and moved the fit — `R/3-build_map.R` shifted the
    BS2017SS golden objective ~18 units. Reset to **pure v5.0**: `R/3-build_map.R`,
    `R/4-build_parameter_bounds.R`, `R/0-rceattle_class.R`, `R/0-read_write_excel_data.R`
    (+ the auto-merged `tests/`, `man/`, examples). The branch now equals v5.0 + the inert
    DSEM files (diff vs `dev-data-workflow` is only DESCRIPTION/NAMESPACE/README +
    DSEM-new files).
    - **Correction (2026-08-05):** the note above said dev-DSEM's additions to
      `R/2-build_params.R` / `R/3-build_map.R` / `R/4-build_parameter_bounds.R` "must be
      re-applied in Tier 1". They must **not**. Those diffs are pure *deletions* —
      `param_list$rec_dev` and `param_list$R_log_sd` are commented out, as are their map
      and bounds blocks, and `rec_dev` becomes a local `matrix<Type>`
      (`dev-DSEM:ceattle_v01_11.cpp:420`) rather than a `PARAMETER_MATRIX`. On dev-DSEM,
      DSEM is a **hard replacement** of the recruitment-deviation parameterization, not an
      optional feature: `fit_mod()` defaults to `dsem = build_DSEM()` and *overwrites*
      `random_vars`. There is no "DSEM off" path there, which is why dev-DSEM cannot
      reproduce a v5.0 golden fit. Tier 1 is a redesign into an opt-in path, not a port.
      Keep `rec_dev` and `R_log_sd` as parameters; default `dsem = NULL`.
- [x] **Tier 0b — synced with `dev` (2026-08-05).** Merged `origin/dev` (through
  `3a321ef0`); zero conflicts, as the branch's diff vs `dev` is purely additive. Moved
  `dsem (== 3.0.0)` from `Imports:` to `Suggests:` (the runtime version guard already lives
  in `build_dsem_objects()`, and a hard `==` pin on a released package would tax every user
  and every CI job); rewrote the README note; marked `test-dsem-recruitment.R` skipped with
  a pointer to Tier 3. **BS2017SS golden re-pinned — see the verification gates below; the
  earlier `10241.0304275` predates the `mse-debug` work and the move to `newtonsteps = 3`.**
  The `mse-debug` PR was still open at this point; a second merge is needed once it lands.
- [ ] Tier 1 — R wiring (`fit_mod`, params/map, data I/O)
- [ ] Tier 2 — C++ re-wiring (`ceattle.cpp` recruitment path) — the crux
- [ ] Tier 3 — tests, docs, verification

## The load-bearing design decision (resolve before Tier 2)

DSEM and the v5.0 **linkage grammar both own "structure on recruitment deviations"**:
- DSEM replaces the recdev density with `jnll_dsem` (latent-state GMRF; env covariates + lags).
- The linkage grammar does covariate / `(1|Year)` / `rw()` / `ar1()` on recruitment via
  `build_srr(linkages=)`, the `rceattle_apply_recruitment_linkages` accumulator, and the
  `JNLL_LINKAGE_RE` row.

Decide how they coexist in the rewritten recruitment path:
- **(A) Mutually exclusive** — a model uses DSEM *or* a recruitment linkage, never both; a
  single switch in `ceattle.cpp` routes the recdev density to one or the other. *Less work;
  recommended for the first landing.*
- **(B) Unified** — express DSEM as a linkage-grammar backend so it composes with other
  linkages. *Much larger; defer.*

**Resolved 2026-08-05: (A) now, (B) when expanding to growth and mortality.** The fall
workflows need a working recruitment DSEM on v5.0, so the first landing is bespoke. But (B)
is the committed target, because DSEM is wanted on growth and M too — and repeating (A)'s
glue per process would mean neutering each process's own density and mapping age/sex-
structured deviations onto `x_tj`'s one-column-per-series layout.

**(B) is cheaper than this file previously assumed, and the seam already exists.**
`calculate_dsem()` takes 17 arguments, *none* recruitment-specific — `dsem.hpp` contains zero
occurrences of `rec_` and knows only an `n_t x n_j` stacked state space plus a RAM. All
recruitment semantics are ~12 lines of glue (`ceattle_v01_11.cpp:888-899`) and five hardcoded
`recdevs<sp>` sites in `R/0-build_DSEM.R`. Meanwhile commit `131e6b31` routes every
random-effect read through one slot-space vector (`beta_linkage_re_all` -> `beta_linkage_eff`),
and that single vector feeds all five `rceattle_apply_*_linkages` accumulators. **A joint GMRF
density attached at slot space would therefore drive growth, M, recruitment, q and sel with no
accumulator changes at all.** (B) needs exactly three things the grammar lacks: multi-
`(process, param)` RE groups with a within-group variable axis (`.assign_re_registry()` currently
keys on `process|param|...`, forcing them apart), a `GMRF(Q)` branch in the density dispatch
(`ceattle.cpp` ~L4171, today a 3-way `us`/`rw`/`ar1` scalar switch), and a path-coefficient
parameter vector (nothing analogous to `trans_rho_linkage` exists for cross-variable paths).

**Requirements on the (A) landing so (B) is not a rewrite:**
- Leave `dsem.hpp` generic — no recruitment-specific arguments on `calculate_dsem()`.
- Collapse the five hardcoded `recdevs<sp>` sites into one registry (default-sem generation
  L112-121, column injection L128-141, the `^recdevs[0-9]+$` family regex L164, the SD
  self-loop lookup L253-283, the `rec_dev_col` match L300, plus `check_dsem_spec()` L438).
  Highest-value cheap thing in the landing.
- Isolate the recruitment glue in `ceattle.cpp` to one marked block.
- Give DSEM its own `JnllRow` slot; do not repurpose slot 10 (on v5.0 `JNLL_REC_DEV = 10` is
  the live standard recdev density).
- Keep `dsem = NULL` opt-in and the standard recruitment path textually intact.
- Fix `R/8-sim_mod.R` L206-225, which indexes `x_tj[, sp]` by species number and so assumes
  recdev column == species index — false the moment M/growth columns are added.

## Where DSEM's hooks live (to re-apply) — `dev-DSEM:src/TMB/ceattle_v01_11.cpp`

| what | dev-DSEM location | re-apply onto v5.0 `ceattle.cpp` |
|---|---|---|
| include | `#include "dsem.hpp"` (L28) | add near the other header includes |
| DSEM data | `DATA_IVECTOR(options)` + DSEM DATA block (L71–92) | add to the DATA section |
| DSEM params | DSEM PARAMETERS block (L92+) | add to the PARAMETER section |
| recruitment = latent | `rec_dev.row(sp) = x_tj.col(rec_dev_col(sp))` (L889) | in the recruitment build, replacing the free recdev where DSEM is on |
| density | `Type jnll_dsem=0; calculate_dsem(jnll_dsem, …)` (L862–889) | Section 6 recruitment; **reconcile with the `JNLL_RECRUITMENT` / `JNLL_LINKAGE_RE` rows** — under decision (A), skip the standard recdev density when DSEM owns it |
| report + total | `REPORT(jnll_dsem)` (L3437), `jnll += jnll_dsem` (L3503) | REPORT block + objective accumulation; consider giving DSEM its own `JnllRow` slot |

`src/TMB/dsem.hpp` and `helper_functions.hpp`'s DSEM additions are already in the tree (dsem.hpp
carried; verify the `helper_functions.hpp` DSEM helpers — v5.0's version was taken, so
**re-add the ~34 lines dev-DSEM added to `helper_functions.hpp`** if `calculate_dsem` needs them).

## Tier 1 — R wiring

- **`R/6-fit_mod.R`**: re-add the `dsem = build_DSEM()` arg and the DSEM-object build /
  `data_list$dsem_settings` threading (dev-DSEM `fit_mod` ~L104, L350–369), coexisting with
  the v5.0 `config=` overlay. Mirror how `recFun`/`M1Fun` flow.
- **`R/2-build_params.R` / `R/3-build_map.R` / `R/4-build_parameter_bounds.R`**: re-add the
  DSEM parameters (latent states, SEM path coefficients, obs-SD), their map entries, and
  bounds. All were reset to pure v5.0 in Tier 0 — re-apply dev-DSEM's additions (get them via
  `git diff <base> dev-DSEM -- R/2-build_params.R R/3-build_map.R R/4-build_parameter_bounds.R`).
  **Gate each with the BS2017SS golden check — these are exactly the files whose auto-merge
  moved the objective.**
- **`R/0-rceattle_class.R`**: reset to v5.0; re-apply dev-DSEM's output-class additions
  (the `$dsem` slot, print/summary of DSEM results) — `+126` lines on dev-DSEM.
- **`R/0-read_write_excel_data.R`**: dev-DSEM added `env_data` / dsem-settings I/O (`+42`
  lines). It **auto-merged**, so verify `env_data` round-trips through `read_data`/`write_data`
  on the v5.0 schema (this is the "new element needs read/write support" trap from CLAUDE.md).
- **`R/0-convergence.R`**: v5.0's version was taken; `check_dsem_spec()` (now in
  `0-build_DSEM.R`) calls `.conv_record` / `.conv_overall` — confirm those exist in v5.0's
  convergence.R (same names) or adapt.

## Tier 2 — C++ re-wiring (the crux)

Apply the hooks table above onto `src/TMB/ceattle.cpp`, honoring decision (A). This is the
hardest part (federal-quota model): the recruitment section was rewritten, so DSEM's latent-
state assignment and density must slot into the new structure without disturbing the standard
recruitment path when DSEM is off. Recompile (`/recompile`), iterate.

## Tier 3 — tests, docs, verification

- **`tests/testthat/test-dsem-recruitment.R`**: update to the v5.0 API (flat layout already);
  it is the DSEM regression net.
- `devtools::document()` (mind the roxygen-version churn trap in CLAUDE.md); update NEWS.md
  with a DSEM-on-v5.0 entry; DESCRIPTION version bump.

### Verification gates
1. **`/golden-check` bit-identical** on the 4 non-DSEM models — DSEM being present must not
   move them (with DSEM off, the recruitment path must be byte-for-byte the v5.0 path).
   Read the pinned objectives off `tests/testthat/test-golden-regression.R` rather than any
   value quoted in this file; they were restated when `mse-debug` moved the golden fits to
   `newtonsteps = 3`.
2. **DSEM equivalence — REVISED 2026-08-05.** The original gate ("confirm the same solution
   as `dev-DSEM` today") is **not achievable and must not be chased.** v5.0's 316 commits
   moved the surrounding likelihood substantially (bias adjustment alone is worth several
   hundred jnll units on BS2017SS), and the DSEM GMRF carries normalizing constants the
   standard `dnorm` recdev density does not. Replace with:
   - **2a. DSEM-off = `dev`** (the strong gate): `dsem = NULL` reproduces all four golden
     values bit-identically. This is what proves the Tier-2 surgery is clean.
   - **2b. DSEM-IID ~ standard `random_rec`** (the primary DSEM gate, single-branch): fit the
     same model with `dsem = NULL, random_rec = TRUE` and with `dsem = build_DSEM()` (the
     synthesized IID sem). Recruitment deviations, SSB and R should agree closely; objectives
     differ by the GMRF normalizing constant.
   - **2c. `dev-DSEM` cross-check — sanity only.** Build `dev-DSEM` in a throwaway worktree
     and fit in a *separate* `Rscript` process (only one Rceattle DLL loads at a time, and the
     branches even have different cpp filenames). Compare recdev series shape, `R_sd`, and the
     estimated path coefficients. **Do not compare objectives.**
3. Constructive linkage tests, OSA suites, `dev/verify-refit-like.R`. Note the linkage tests
   are the regression net for the recruitment/linkage path — the 4 golden models carry no
   linkage rows. `tests/testthat/test-linkage-random-effects.R` L103/L288 address the
   linkage-RE row as `jnll_comp[nrow(jnll_comp), ]`, so adding a `JNLL_DSEM` row silently
   re-points them; repair to an explicit index.
4. `dev/verify-mse-hindcast-invariant.R` + `dev/verify-mse-om-horizon.R`. `x_tj` is
   year-dimensioned, so `.mse_proj_param_yrdim` in `R/10-run_mse.R` needs `x_tj = 1L` or a
   shortened OM refit misaligns. Whether `run_mse()` supports DSEM at all is a separate
   decision: `rec_dev` is derived from `x_tj`, so `sample_rec()` must draw into
   `x_tj[, rec_dev_col]`, not `rec_dev`.

## Notes / traps
- The **`JnllRow` enum changed** — the recruitment/linkage rows moved; don't index `jnll_comp`
  by old row numbers. Give DSEM its own row if it carries a distinct density.
- Keep decision (A) simple first; only pursue (B) unification once DSEM lands and is verified.
  Resist starting (B) at Tier 2 just because the slot-space seam is visible.
- `dsem (== 3.0.0)` is now in **`Suggests:`**, not `Imports:` — `build_dsem_objects()` already
  carries the `requireNamespace` + `packageVersion != "3.0.0"` guard with an install hint.
- `src/TMB/dsem.hpp` carries dead code: the `nyrs_hind` and `proj_mean_rec` arguments are never
  referenced in the body, and the `dsem_hindcast_positions` / `dsem_subset_vec` helpers are never
  called. Delete them rather than carrying them into (B).
- `helper_functions.hpp` on this branch is pure v5.0 and is **missing** the `sign()`, `dlnorm()`
  and `devresid_tweedie()` that `dsem.hpp` calls (L543/553/563, L578, L593). Nothing compiles the
  moment `ceattle.cpp` includes `dsem.hpp`. Put them in `dsem.hpp` itself, not in the shared
  header — they are DSEM-only, and a template named `sign` in a header included by every
  translation unit is a needless overload-resolution hazard. Keeping them local also means the
  fix cannot move the fit.
- `env_data` alignment is inconsistent between the two subsystems: the linkage grammar aligns
  **positionally** (`R/0-build_linkage.R`, row r -> year `styr + r - 1`, with
  `.extend_env_data()` gap-filling) while `build_dsem_objects()` aligns **by `Year`** (a
  `full_join`). Benign when `env_data$Year` is complete and sorted (all bundled datasets),
  dangerous otherwise. (B) needs one contract; routing DSEM through `.extend_env_data()` is it.
- `dsem_settings` is a run specification, not data — persist it on `data_list$dsem_settings`
  (read back by `.refit_like()`) and through `save_config()`/`load_config()`. Do **not** add a
  workbook sheet; `env_data` already round-trips on v5.0, so the CLAUDE.md read/write trap does
  not apply here.
- The six diagnostic refit paths collapsed into `.refit_like()` on v5.0, so DSEM needs **one**
  `dsem = data_list$dsem_settings` override there, not dev-DSEM's per-caller edits to
  `9-retro_and_jitter.R` / `OPT-phaser.R` / `10-run_mse.R`. Pass the *settings*, never a prebuilt
  object: `x_tj`/`eps_tj`/`y_tj` are dimensioned over `styr:projyr` and every refit changes the
  horizon.
