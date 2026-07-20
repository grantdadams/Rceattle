# Handoff — data-workflow modernization

**Status:** designed, not started. No code written.
**Branch:** work from `dev-data-workflow` (currently at `016d2f69`, same tip `dev-ebs-pk`
had when this was scoped) or a fresh branch off `dev`.

## The problem

As Rceattle's feature set grows, the *input* grows with it, but any one model uses only a
slice. A single-species age-structured model needs nothing from the predation, diet,
bioenergetics, growth, or environmental-linkage blocks — yet the only supported entry
points make the user confront all of it:

- `read_data()` reads a ~22-sheet xlsx workbook.
- The alternative is hand-building a `data_list` with ~56 elements.

There is no code-first way to build a minimal model, and nothing tells a user which blocks
their configuration actually requires. The stated goal: **users should supply only the data
their model needs.**

## Key finding: the logic already exists internally

This is the important part — the required/optional rules are already encoded, just split
across two functions and not exposed:

- **`clean_data()`** (`R/0-clean_data.R`) already default-fills every optional block:
  `comp_data`, `caal_data`, `emp_sel`, `NByageFixed`, `ration_data`, `diet_data`,
  `env_data`, `index_cov`, `MSSB0`/`MSB0`. The engine already tolerates a minimal input.
- **`data_check()`** (`R/1-data_check.R`) already encodes what is *conditionally* required,
  gating on `msmMode`, `growth_model`, `estDynamics`, `suitMode`, per-fleet `Catchability`
  / `Selectivity` / `Time_varying_*`, and environmental-index references.
- **`switch_check()`** (`R/0-switches.R`) fills switch defaults and derives some values from
  the data (e.g. `Sel_start_year` from first year of observations, grouped by
  `Selectivity_index`).

So the gap is **only at the user-facing surface**. A thin layer over existing logic — not a
rewrite. That is what makes this tractable and low-risk on a released package.

## Prior art worth copying

- **WHAM** `prepare_wham_input(basic_info, selectivity = , M = , NAA_re = , ...)` — one
  constructor taking a small named list per process, everything else defaulted.
- **SPoRC** — per-process builders composed into an input object.
- **FIMS** — modular components assembled explicitly.

## Options considered

Four entry points were scoped. They are not mutually exclusive; 2 is a natural precursor to 1.

1. **`build_data()` constructor** *(recommended first deliverable)*
   Pass only the blocks the model uses; everything else defaulted through the existing
   `clean_data()` path.
   ```r
   dat <- build_data(nspp = 1, styr = 1977, endyr = 2023, nages = 10,
                     maturity = mat, sex_ratio = sr,
                     fleet_control = fc, catch_data = catch,
                     index_data = survey, comp_data = ages)
   fit_mod(dat)
   ```
2. **`data_requirements()` introspection** — given a model spec, report which blocks are
   Required / Optional / Ignored, reusing the `data_check()` gates. Lowest effort, useful
   immediately, and it forces the gate logic into a reusable form that (1) can then consume.
3. **Modular per-process builders** — `build_population() |> add_fishery() |> add_survey()
   |> add_composition() |> add_predation()`. Most flexible and self-documenting; largest new
   surface to test and maintain.
4. **Tolerant xlsx + `write_template()`** — emit only the sheets a configuration needs.
   Smallest change; helps existing xlsx users but not code-first users.

## Recommended sequence

**Step 1 — extract the requirement table.** Refactor the `data_check()` gates into a
declarative spec: one row per `data_list` element with the condition under which it is
required. `data_check()` then *consumes* that table instead of hard-coding the conditions.
No behavior change; this is the enabling refactor and should hold golden equivalence.

**Step 2 — `data_requirements()`.** A reader over that table. Pure introspection, no risk.

**Step 3 — `build_data()`.** Constructor validating against the same table, defaulting
everything else via `clean_data()`. Keep `read_data()` / `write_data()` untouched.

**Step 4 (optional) — per-process builders** as sugar over `build_data()`, only if demand
justifies the surface.

## Constraints (from CLAUDE.md and this codebase)

- **Released package — preserve the public API.** `read_data()` / `write_data()` /
  `fit_mod()` signatures must keep working. Add, don't replace; deprecate rather than delete.
- **Golden-reference equivalence.** Any change that can move a fit must hold
  BS2017SS = `10241.030427` and keep the suite green.
- **`Fleet_code` must equal the `fleet_control` row number** (now enforced in `data_check()`).
  A constructor must preserve that invariant.
- Numeric prefixes in `R/` are meaningful; a new file fits the `0-*` prep tier.
- Tests are flat `tests/testthat/test-<area>-<topic>.R`; anything running a real `fit_mod()`
  needs `testthat::skip_on_cran()`.
- Match local comment density; keep comments about *the code*, not about the change.

## Open questions for whoever picks this up

1. Does `build_data()` validate eagerly (error on construction) or defer to `data_check()`
   at fit time? Deferring keeps one source of truth; erroring early gives better messages.
2. Should the requirement table be internal, or exported so users can query it directly?
3. Multispecies: `build_data()` for `msmMode > 0` needs the predation/diet blocks — one
   constructor with optional args, or a separate `add_predation()`?
4. Is `write_template()` (option 4) worth doing for the existing xlsx user base, or does
   the code-first path supersede it?

## Where to start reading

- `R/0-clean_data.R` — the optional-block defaults (and `.align_index_cov`).
- `R/1-data_check.R` — the conditional-requirement gates; the per-fleet loop is the densest part.
- `R/0-switches.R` — `switch_check()` defaults and data-derived values.
- `R/0-read_write_excel_data.R` — the current xlsx surface (`read_data` is already tolerant
  of several missing sheets via `%in% sheetnames`, a useful precedent).
- `R/data.R` — the documented `data_list` element and `fleet_control` column reference.
