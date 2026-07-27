# Prompt — documentation-quality review of the PR 1–7 doc surface

> **DONE (2026-07-27).** Executed on `pr7-legibility`. All five clusters applied —
> linkage/config roxygen, the data-workflow roxygen + `R/data.R`/`column_schema.R` field
> dictionaries (regenerated `meta_data_names.xlsx`), the `data_check()` console messages,
> the vignettes (substantial rewrite of `environmental-linkages-and-priors.Rmd`; light
> passes elsewhere), and a NEWS trim. Two factual defects fixed along the way: the vignette's
> slot-19-vs-20 contradiction (verified `JNLL_LINKAGE_PRIOR = 19`) and a stale "`obs_sd` …
> is planned" claim that contradicted the shipped `obs_sd_est = TRUE`. `devtools::document()`
> touched only the 16 intended `man/*.Rd`; doc-coupled tests + vignette knit-checks green.
> A separate technical-accuracy audit of the `fit_mod`/`build_*` roxygen math followed.
> Superseded — kept for history only.

Copy the block below into a fresh session. This is a **documentation-only** pass over
everything the data-workflow + linkage-grammar effort (PRs 1–7) added or changed. It does
**not** touch code, the public API, or any fit — the code was adversarially reviewed and
golden-verified throughout; the remaining value is purely in the prose.

---

We're reviewing the **user-facing documentation** produced by PRs 1–7 on `dev-data-workflow`
(the effort from DESCRIPTION `4.6.0` → `5.0.0`; commit set `git log $(git merge-base main
dev-data-workflow)..dev-data-workflow`, ~232 commits). Rceattle is a released package that
drives US federal quota advice and has real assessment-scientist users, so the docs are part
of its operational credibility — worth getting right for the 5.0.0 release.

## The three quality criteria

Judge every piece of documentation against these. Flag violations with `file:line`, the
offending text, and a concrete rewrite.

1. **Not AI-verbose.** Terse, direct, no filler. Cut throat-clearing ("This function
   provides a convenient way to…"), hedging ("it should be noted that", "in general"),
   marketing adjectives ("powerful", "flexible", "seamless", "robust", "comprehensive"),
   and restating the obvious. Lead with the verb.
   - Bad: *"This function provides a convenient and flexible interface for comprehensively
     specifying a stock-recruitment relationship for use in the model."*
   - Good: *"Specify a stock-recruit relationship."*

2. **Current capabilities, not the history of changes.** Documentation describes what the
   software **does now**, not how it got there. Delete "previously / now / new / this PR /
   as of v4.x / used to / we added / replaces the old". A reader should not be able to tell
   which release introduced a feature from its help page. **The sole exception is `NEWS.md`**
   — that is the changelog and *should* describe changes (but still terse + clear).
   - Bad: *"Previously `q` was fixed; now you can make it time-varying via the new linkage
     system introduced in 4.9.0."*
   - Good: *"Make catchability time-varying: `build_catchability(linkages = list(q =
     linkage_spec(~ year)))`."*

3. **Understandable to an average fisheries scientist / ecologist.** Explain in terms of the
   assessment concept, not the implementation. A stock-assessment scientist should recognize
   the task without knowing TMB, AD, or the internal table encoding.
   - Use the domain vocabulary from `CLAUDE.md` ("Domain vocabulary" section): SSB = female
     spawning-stock biomass; F40%/F35%/B40% SPR proxies (not "MSY"); name the selectivity
     form; Francis / McAllister–Ianelli / Dirichlet-multinomial weighting; Mohn's rho; OSA
     residuals; CAAL; ageing-error matrix.
   - Bad (user-facing): *"a prior on the intercept-only case of the linkage table, keyed by
     `(Intercept)`."*
   - Good: *"put a prior on a survey's catchability: `linkage_spec(~ 1, priors =
     list(`(Intercept)` = normal(0.8, 0.2)))`."*
   - Keep TMB/AD/engineering jargon out of user-facing help; it is fine in the
     developer guide (contributor-facing) and `CLAUDE.md` (out of scope here).

## Scope (the doc surface)

Edit the **source of truth**, not generated files — roxygen in `R/*.R` (never `man/*.Rd`
directly) and the `.Rmd`. Then regenerate. Surface:

- **Roxygen** on the PR 1–7 functions (~47 `R/*.R` → 81 `man/*.Rd`): the linkage system
  (`linkage_spec`, `build_srr`/`build_M1`/`build_growth`/`build_catchability`/
  `build_selectivity`/`build_composition`, the `*_LINKAGE_*` constant docs, `Rceattle_priors`);
  the data workflow (`build_data`, `model_config`, `data_requirements`, `write_template`,
  the `R/0-column_schema.R` field docs); run configuration (`save_config`, `load_config`,
  `run_config`, `fit_mod(config=)`); `mse_summary`; and the `?BS2017SS`-style `@format` field
  dictionary. Check `@description`, `@details`, every `@param`, `@return`, and that `@examples`
  are runnable and show a realistic assessment task.
- **Vignettes** (9 touched): `environmental-linkages-and-priors.Rmd`, `data-without-excel.Rmd`,
  `articles/developer-guide.Rmd`, plus the older `hcrs-and-mses`, `model-diagnostics`,
  `model-options-and-functionality`, `model-parameterizations`,
  `projections-and-reference-points`, `stock-synthesis-conversion` (check only where PR 1–7
  changed them). They are `eval = FALSE` — they must still knit.
- **`NEWS.md`** 4.7.0 → 5.0.0: criterion 1 + 3 only (it is allowed to describe changes). Make
  each bullet a terse, scientist-legible statement of the user-facing change.
- **Console / validation messages** in `data_check()`, `data_requirements()`, `switch_check()`,
  `combine_data()`: are they intelligible to a scientist hitting them, and do they name the
  fix? (These are the first docs a user actually reads.)
- **`README.md`** (package landing) and the field dictionary (`meta_data_names.xlsx`).
- Out of scope: `CLAUDE.md` (Claude-facing, curated already), `dev/` notes, code.

## Guardrails (from CLAUDE.md — do not trip these)

- **roxygen2 version churn:** the repo is documented with a specific roxygen2 version; a newer
  local one rewrites *every* `man/*.Rd` + swaps `RoxygenNote`. After `devtools::document()`,
  `git checkout` the unrelated `man/`/`DESCRIPTION` churn and keep only the `.Rd` you meant to
  change.
- Never insert anything between a function's roxygen block and its definition (the block binds
  to the next object — it would steal `@export`/`@importFrom`). Keep internal helpers `@noRd`.
- Preserve the public API and every roxygen tag; this is a **prose-only** edit. Deprecation
  text stays. Domain terms stay exact.
- Keep `@examples` and vignette chunks runnable / knittable (vignettes are `eval = FALSE`).

## Method + deliverable

The surface is large but parallelizable — it is a good fit for a **multi-agent fan-out**
(one reviewer per cluster: linkage roxygen, data-workflow roxygen, config roxygen, each
vignette, NEWS, console messages), each returning a findings list keyed to the three criteria,
then a synthesis + apply pass. Or do it sequentially by cluster if running solo.

For each cluster: (1) report findings (`file:line`, quote, which criterion, suggested rewrite);
(2) after review, apply the accepted rewrites; (3) `devtools::document()` and clean the churn;
(4) confirm the vignettes still knit (`eval = FALSE`, so a render check suffices).

**No golden-check / no version bump** — documentation-only, zero behavior change. Commit in
doc-only batches by cluster with plain messages (no AI-attribution trailer). If a "doc" change
turns out to imply a behavior or API change, stop and surface it — it is out of scope here.

## Is it worth it?

Yes, scoped to prose. The code is verified; this closes the gap between a correct model and a
model an assessment scientist can *use* correctly — the operational-credibility priority. Keep
it strictly documentation; resist re-reviewing the code.
