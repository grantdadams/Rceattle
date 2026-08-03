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

Everything in Tier 2 depends on this. Start with (A).

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
- **`R/2-build_params.R` / `R/3-build_map.R`**: re-add the DSEM parameters (latent states,
  SEM path coefficients, obs-SD) and their map entries. dev-DSEM diffs: build_params `+4`,
  build_map `+11` — small; graft onto the v5.0 builders.
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
2. **DSEM equivalence** — run `examples/DSEM_example.R` on this branch and confirm the same
   solution as `dev-DSEM` today (objective / SSB / R / recruitment latents), the same
   before/after method used for the GOA-ATF vs main check. Build `dev-DSEM` in a throwaway
   worktree (`git worktree add … dev-DSEM`), fit both at `-O2`, diff.
3. Constructive linkage tests, OSA suites, `dev/verify-refit-like.R`.

## Notes / traps
- The **`JnllRow` enum changed** — the recruitment/linkage rows moved; don't index `jnll_comp`
  by old row numbers. Give DSEM its own row if it carries a distinct density.
- Keep decision (A) simple first; only pursue (B) unification once DSEM lands and is verified.
- `dsem (== 3.0.0)` is a hard pin — the branch needs dsem 3.0.0 installed (it is on this box).
