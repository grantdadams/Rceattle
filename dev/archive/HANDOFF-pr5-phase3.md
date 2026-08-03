# PR 5 — resume at Phase 3 (renames + docs + read robustness + write_template)

> **STATUS: COMPLETE (v4.12.0).** Shipped in 11 commits on `dev-data-workflow`
> (`f59fee43`..`ec055d13`), each golden bit-identical + adversarially reviewed: 3a accumulation
> deprecation + column docs, 3b.0–3b.5 the schema-driven alias cascade + the full intuitive-naming
> rename pass (all alias-backed) + `.rda`/xlsx regen + back-compat tests, 3c the roxygen
> field-dictionary + drift-guard, 4.7–4.10 read-path robustness + `write_template()`, finalize
> the committed golden-regression test + NEWS + DESCRIPTION 4.12.0. `rcmdcheck`: 0 errors (the
> remaining WARNING/NOTES are pre-existing). See the `pr5-schema-progress` memory. The doc below
> is the original resume brief, kept for the record.

Start a **fresh session** for this. Self-contained. Branch **`dev-data-workflow`**.

## Read first
- Memory `pr5-schema-progress` (index in MEMORY.md) — the one-screen status + resume pointer.
- Plan `~/.claude/plans/we-re-continuing-the-rceattle-twinkly-puffin.md` — **every locked decision**
  (rename table, accumulation-deprecation, doc-pass scope, phase sequencing). Read the "New renames
  introduced by this PR" and "Phase 3a/3b/3c" sections before touching code.

## Where things stand (Phase 1 + 2 shipped)
PR 5 kills the four-way `fleet_control`/workbook schema drift by defining the column dictionary
**once** in `R/0-column_schema.R` (`.rce_column_schema()`) and having every consumer read from it.
Committed on `dev-data-workflow`, all **golden bit-identical**, each adversarially reviewed:
- `e553d32b` Phase 1 — canonical schema + guard tests (`tests/testthat/test-schema-canonical.R`).
- `ae8e7758` Phase 2.1 — `switch_check()` fleet_control defaults via `.rce_apply_default()`.
- `de462d96` Phase 2.2 — `write_data()` control/bio object order from `.rce_schema_names()`.
- `a40394eb` Phase 2.3 — regenerated `meta_data_names.xlsx` (legacy→canonical) + exact drift-guard;
  reproducible recipe `data-raw/regenerate_meta_xlsx.R`.
- `e11b4151` Phase 2 hardening — `has_default` is wired for **fleet_control only** (control/bio
  `default` fields are documentation-only); non-self-referential defaults pin.

The schema is now the single source consumed by (a) `switch_check()` defaults, (b) `write_data()`
control/bioenergetics object order, (c) the generated `meta_data_names.xlsx`. `aliases` is currently
**documentation-only** — the rename cascades in `read_data()`/`switch_check()` are still hand-coded
(Phase 3b wires the schema `aliases` into them).

## The task — Phase 3a → finalize (in order)

**3a — document the gaps + deprecate accumulation ages.**
- **Deprecate accumulation ages ENTIRELY.** `Accumulation_age_lower/upper` is a dead feature: only
  *validated* at `R/1-data_check.R:457-468`, never *applied* (no comp grouping anywhere in R or the
  C++ template). Remove that validation block; drop BOTH spellings (incl. the `Acuumulation` roxygen
  typo) from `R/data.R`; do NOT add it to the schema; emit a one-time soft-deprecation warning if a
  workbook/data list still carries the column. No bundled or sibling data uses it → golden-inert.
- Flip the ~18 used-but-undocumented columns from `meta = FALSE` to `meta = TRUE` (so they appear in
  the xlsx): `Selectivity_dimension, Sel_curve_pen3, Sel_start_year, Sel_shape_sd, Sel_curvature_sd,
  Sel_devmag_sd, Sel_shape_dir, Sel_pen_first_bin, Sel_pen_last_bin, Sel_shape_mode, Sel_avgsel_pen,
  Sel_cap_bin, Index_loglike, Month`, plus bioenergetics `Diet_loglike, Diet_comp_weights`.
- Regenerate the xlsx (`Rscript data-raw/regenerate_meta_xlsx.R`); the drift-guard + alias-pin tests
  update accordingly. Golden bit-identical (doc-only).

**3b — the intuitive-naming renames (the big user-facing step).** Each keeps a back-compat alias via
the double-fire-safe `rename_deprecated_col()` pattern (NOT `dplyr::rename`). Migrate the rename
cascades to consume the schema `aliases`; add data_list-level handlers for the non-fleet_control
ones; update `utils::globalVariables()` (`R/Rceattle-package.R`); regen the 3 bundled `.rda` +
`inst/extdata/*.xlsx`; `R CMD INSTALL .`; sweep `../Rceattle-models` (+ `Climate_MSE_projections`,
`GOA-ATF-ESP`, `GOA_circlulation_study`). **Golden bit-identical** each rename + a back-compat test
(old-name data list fits identically) + old-name-xlsx round-trip test.

Locked rename table (fleet_control input columns):

| Current | New |
|---|---|
| `Sel_norm_bin1` | `Sel_norm_bin` |
| `Sel_norm_bin2` | `Sel_norm_bin_upper` |
| `weight1_Numbers2` / `Weight1_Numbers2` | `Observation_units` (fix the case bug in the alias) |
| `Comp_loglike` | `Comp_distribution` |
| `CAAL_loglike` | `CAAL_distribution` |
| `Index_loglike` | `Index_distribution` |
| `Q_index` | `Catchability_index` |
| `Q_init` | `Catchability_init` |
| `Q_sd_prior` | `Catchability_prior_sd` |
| `proj_F_prop` | `Proj_F_proportion` |

- **KEEP** `Sel_curve_pen1/2/3` (reused containers — shape/curvature/devmag only for NonParametric;
  2DAR1/3DAR1 use them as logit correlations, LogisticPM as RW weights).
- Output column (post-fit, `R/6-rename_output.R`): `Est_weights_mcallister` → `Comp_weights_mcallister`
  (parity with `Diet_weights_mcallister`; not an input, no data alias).
- data_list-level: `Diet_loglike` → `Diet_distribution` (per-species vector; bespoke alias + Excel
  `bioenergetics_control` label + `combine_data_sets`/`build_data` name lists); `sigma_rec_prior` →
  `sigma_rec` (control scalar; misleading `_prior`; alias + row_labels + `.rda` regen).
- Drop orphans `Catch_units`, `Sex`, `R_sexr` (already `status = "orphan"` in the schema).

**3c — intuitive documentation pass (field-dictionary scope only).** Rewrite the ~50 schema `doc`
strings for clarity (plain language, list readable string AND integer code, reference the new names),
which drives BOTH the xlsx and (next) the roxygen. Then ALIGN `R/data.R` `@format` to the schema
`doc` and add the roxygen drift-guard test (this folds in the old Phase 2.4). `document()`; regen
xlsx. Golden-inert. Do NOT touch the data-object `@describe` prose or function docs (PR 7).

**4.7–4.10 — read robustness + template.**
- 4.7: guard the 8 unguarded `read_xlsx` sites in `R/0-read_write_excel_data.R` (hoist
  `excel_sheets()`; `if(sheet %in% sheetnames)`; required-and-absent → clear `stop()`, optional-and-
  absent → NULL/defaulted). A minimal single-species workbook must read cleanly.
- 4.8: fix the silent `suppressWarnings(as.numeric(...))` swallow at `read_data()` ~L245-247 (and the
  bioenergetics loop ~L444-446) — error naming the offending Object + species cell on an unparseable
  non-empty value.
- 4.9: `rearrange_data()` presence guards for the required-identity / hard-`$` fleet_control sites
  (`Species`, `Fleet_code`, `Fleet_type`, `Fleet_name`, `Selectivity`, `Q_init`→`Catchability_init`,
  and the `Sel_start_year`/`Selectivity_index` helpers).
- 4.10: `write_template()` (new export) — minimal single-species xlsx from `data_requirements()` +
  the schema; round-trip → `data_check()` → `fit_mod(estimateMode = 3)` builds.

**Finalize.** Add a committed `skip_on_cran` **golden regression test** keyed on the objective within
~1e-6 of inline reference constants (`ss=10241.0304275402, ms=10267.2478327352,
goa_ss=12807.4375258732, goa_ms=12866.2957391829`) — Grant approved this. Update `NEWS.md` (under
`# Rceattle 4.11.0`), bump `DESCRIPTION` 4.11.0 → **4.12.0** (minor), touch `vignettes/
data-without-excel.Rmd` + `stock-synthesis-conversion.Rmd`, `document()`, background `rcmdcheck`.

## Working standard (non-negotiable — federal quota advice)
- **Constructive verification, not "it runs."** Prove equivalence; enumerate every intended diff.
- **Golden-gate every step**: `/golden-check` (or `dev/golden-ref.rds` recipe) across all 4 models
  (BS + GOA, SS + MS) — bit-identical for renames/docs. Use `getsd = FALSE` (validated inert, faster);
  the SSB field is `quantities$ssb` (NOT `biomassSSB`); MS fits warm-start from SS MLEs; GOA-MS uses
  fixed M. Fast NOT_CRAN suite green.
- **Spawn an adversarial reviewer (Explore agent) after each substantive step**, before committing.
- **Incremental verified sub-commits; Grant confirms before each commit.** Plain messages, **no
  AI-attribution trailer**. Prefix R with `export PATH=/usr/bin:$PATH`.
- **Released package** — preserve `read_data`/`write_data`/`fit_mod` signatures; every renamed column
  keeps a deprecation alias; sweep the sibling repos after a rename.
- Toolchain / dev loop: `pkgload::load_all(".", quiet=TRUE)` recompiles TMB; run one test file via
  `e <- new.env(parent = asNamespace("Rceattle")); source helpers; testthat::test_file(f, env=e)`.

## Key files
- `R/0-column_schema.R` (schema + generators + `.rce_apply_default`), `R/0-switches.R` (defaults +
  rename cascade), `R/0-read_write_excel_data.R` (read/write, read_xlsx guards, `write_template()`),
  `R/data.R` (roxygen `@format`), `R/5-rearrange_data.R`, `R/6-rename_output.R`,
  `R/Rceattle-package.R` (`globalVariables`), `R/0-combine_data_sets.R` + `R/0-build_data.R`,
  `data-raw/regenerate_meta_xlsx.R`, the 3 bundled `data/*.rda` + `inst/extdata/*.xlsx`, and the
  `tests/testthat/test-schema-*.R` / `test-switches-schema-defaults.R` / `test-data-workbook-roundtrip.R`.
