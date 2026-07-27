# Handoff — PR 5: optional Excel data + one canonical `fleet_control` schema

Start a **fresh session** for this. Self-contained. Branch **`dev-data-workflow`**.

## Where things stand (PR 4 is done)

PR 4 shipped the user-facing data-construction + introspection layer, all committed and
verified: `data_requirements()` (introspection over a declarative requirement table),
`build_data()` (code-first constructor), `model_config()` (a fit-config slot on the data list),
and a spec-tree `print()`/`summary()` for a new `Rceattle_data` class. The `dev`, MSE, and hake
branches were merged in; DESCRIPTION is at **4.11.0**; the full NOT_CRAN suite is green (447
tests, 0 failures); the 4-model golden (BS + GOA) was re-pinned and cross-checked against `main`
(numerically inert). See the `golden-reference-state` and `intuitive-naming-priority` memories.

**Key enabling asset for PR 5:** `data_requirements()` + the requirement table
(`R/1-data_requirements_table.R`) are the beginnings of a single source of truth for "what a
data element is". PR 5 extends that idea to the **`fleet_control` column schema**, which today
is duplicated across **four** unsynchronized places (this is the drift PR 5 kills):
1. `R/data.R` roxygen `@format` (the field dictionary, ~L56-92) — hand-maintained.
2. `inst/extdata/meta_data_names.xlsx` (the `meta_data` sheet embedded in every workbook).
3. `write_data()`'s ordered `row_labels` / `bioenergetics_control` vectors
   (`R/0-read_write_excel_data.R:31-56, 118-141`).
4. `switch_check()`'s per-column defaults (`R/0-switches.R`), with a literal
   "must be kept in sync by eye" comment at `R/0-switches.R:27-30`.

## PR 5 scope (from `dev/PLAN-data-workflow-and-linkage-grammar.md`, "PR 5")

1. **One canonical `fleet_control` schema.** Define the schema **once** — for each column:
   name, type, default, allowed values, TMB target, doc string — as an internal table (mirror
   the `R/1-data_requirements_table.R` pattern). Then **generate** from it: the `R/data.R`
   roxygen, `switch_check()`'s defaults, `meta_data_names.xlsx`, and the Excel template. This is
   the fix for:
   - the **22 used-but-undocumented** `fleet_control` columns (survey/grep the columns actually
     read in `R/2/3/4/5-*.R` + the cpp vs those documented);
   - **stale documented entries**: `Acuumulation_age_lower/upper` (sic — also the
     `Accumulation` typo), `Catch_units`, `Sex`, and the legacy `Estimate_survey_sd` /
     `Survey_sd_prior` documented *instead of* the current `Estimate_index_sd` / `Index_sd`
     (renamed in PR 4-era work — see `intuitive-naming-priority` memory);
   - the `weight1_Numbers2` case mismatch;
   - the "kept in sync by eye" comment at `R/0-switches.R:27-30`.
2. **Optional Excel sheets.** Guard the **9 unguarded `read_xlsx` sites** in
   `R/0-read_write_excel_data.R` (grep `read_xlsx` / `readxl::`); route a missing sheet through
   the existing `clean_data()` / `switch_check()` defaulting instead of erroring — so a minimal
   workbook (single-species, no diet/predation sheets) reads cleanly, mirroring what
   `build_data()` already tolerates in-memory.
3. **`write_template()` (deferred from PR 4).** Emit a **minimal single-species xlsx template**
   containing only the sheets a config needs, generated from `data_requirements()` + the new
   canonical schema. (This is the PR-4 open question we deliberately punted to here.)
4. **Read robustness.** Fix `read_data()` (~L207) swallowing non-numeric control entries via
   `suppressWarnings(as.numeric(...))` — a typo currently silently becomes `NA` and fails much
   later; make it error at read with a clear message naming the cell.
5. **Add the missing `fleet_control` presence guards in `rearrange_data()`** (it assumes columns
   present that a minimal workbook may lack).

## How to work it (this repo's standard)
- **Plan first, get approval** (EnterPlanMode/ExitPlanMode); then **incremental verified
  sub-commits**. Each keeps `/golden-check` within tolerance (all 4 models) and the fast suite
  green. Schema-generation changes are the risky part — a generated `switch_check()` default
  that differs by a hair moves a fit, so **golden-gate every step**.
- **Constructive verification** (federal-rigor memory): round-trip real datasets through
  `write_data()` → `read_data()` and assert identical; assert the generated `meta_data_names.xlsx`
  / roxygen / defaults match the hand-written ones *before* switching over (prove the generator
  reproduces the status quo, then delete the hand-copies). Spawn an adversarial reviewer after
  each substantive step.
- **Released package** — preserve the public API. `read_data()`/`write_data()`/`fit_mod()`
  signatures unchanged; keep deprecation aliases for renamed columns (don't break the sibling
  `../Rceattle-models` scripts — sweep them after any column rename).
- Commits: **plain messages, no AI-attribution trailer**. Prefix R with `export PATH=/usr/bin:$PATH`.
  After a schema/`.rda` change, regenerate the bundled `.rda` and `R CMD INSTALL .` before
  local parallel tests. After any column rename, sweep `../Rceattle-models` (+ the active repos:
  `Climate_MSE_projections`, `GOA-ATF-ESP`, `GOA_circlulation_study`).
- Per CLAUDE.md, update **NEWS.md** (under `# Rceattle 4.11.0`) + bump **DESCRIPTION** (minor)
  + touch affected vignettes (`data-without-excel.Rmd`, `stock-synthesis-conversion.Rmd`).

## Key files
- `R/0-read_write_excel_data.R` (read/write, the `read_xlsx` guards, `row_labels`), `R/0-switches.R`
  (defaults + the "sync by eye" comment), `R/data.R` (roxygen `@format`), `R/5-rearrange_data.R`
  (presence guards), `inst/extdata/meta_data_names.xlsx` + `*.xlsx` templates, and the new
  canonical-schema table (new `R/0-*` file, mirror `R/1-data_requirements_table.R`).
- Plan: `dev/PLAN-data-workflow-and-linkage-grammar.md` (PR 5 section; PR 6 `save_config`/
  `load_config` reuses this schema; PR 7 contributor docs + C++ legibility follow).

## Deferred / parked (not PR 5)
- **GOA multispecies M-estimation convergence** — see `dev/PROMPT-goa-ms-m-estimation.md`
  (its own session). The `goa_ms` golden uses fixed M until that's solved.
- Per-process builders (`build_population() |> add_fishery() |> ...`) — PR 4 step 4, only if
  demand justifies it.
