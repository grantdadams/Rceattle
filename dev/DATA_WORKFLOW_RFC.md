# RFC: Rceattle data-input workflow redesign

Author: Claude (for Grant Adams) · Date: 2026-07-05 · Status: **draft for review**

> Goal: as Rceattle's feature set grows (predation/multispecies, DSEM environmental
> linkages, time-varying selectivity, OSA, growth), the *input* dataset has grown
> too — but any single model configuration needs only a **subset**. This RFC proposes
> how to let users supply only the data their model needs, without breaking the
> released public API, and drawing on SPoRC, WHAM, and FIMS.

---

## 1. The problem, precisely

Investigation of the current codebase surfaced that the pain is **not the monolithic
list itself** — it is three specific things layered on top of it:

1. **An Excel-vs-programmatic asymmetry.** `read_data()` (Excel) reads most sheets with
   **no existence guard** (`R/0-read_write_excel_data.R`), so a single-species user must
   physically supply `diet_data` (117 rows), `ration_data` (2,025 rows),
   `bioenergetics_control`, `emp_sel`, `NByageFixed`, … sheets they never use — the
   `BS2017SS.xlsx` template ships **19 sheets** for a model that reads ~10. The
   *programmatic* path, by contrast, **silently back-fills** missing pieces via
   `clean_data()` / `switch_check()`. Same model, two very different burdens.

2. **A hidden second input surface.** Half the model *configuration* —
   `msmMode`, `initMode`, `HCR`, `srr_*`, `M1_*`, `growth_*`, `comp_offset`,
   `random_rec` — lives in `fit_mod()` args and `build_*()` helpers, **not in the data
   object**. So "what do I need to provide?" has no single answer, and the data you must
   supply depends on switches that live somewhere else.

3. **Late, cryptic validation.** Control-flag ↔ data coupling is ~80% caught up front in
   `data_check()`, but the rest fails deep in `rearrange_data()` with dplyr/`NA` errors
   (e.g. an unguarded `fleet_control` column pull, an ALK/CAAL `stop()`). A user learns
   the requirement only after a confusing traceback.

Secondary: the **schema is split** across four hand-synced sources that already diverge —
`write_data()`/`read_data()`, the `@format` roxygen in `R/data.R`, the switch maps in
`R/0-switches.R`, and the static `inst/extdata/meta_data_names.xlsx` glossary. And
`fleet_control` is **32–48 columns wide**, mostly selectivity plumbing irrelevant to a
simple logistic model.

---

## 2. What the reference packages do

All three ultimately build the same TMB/RTMB `(data, par, map, random)` bundle for
`MakeADFun`. The difference is entirely in the R-side ergonomics of populating it.

| Package | Pattern | "Off" idiom | Validation |
|---|---|---|---|
| **SPoRC** | Incremental **builder chain**: `Setup_Mod_Dim()` seeds an `input_list`, then `Setup_Mod_Rec()`, `Setup_Mod_Biologicals()`, `Setup_Mod_SrvIdx_and_Comps()`, … each take it and return it augmented. One call per active process. | Per-element string sentinels: `fish_idx_type = c("biom","none")`, `M_spec="fix"`. | Each `Setup_Mod_*` validates its arrays against the dims from `Setup_Mod_Dim()`; errors at the offending call. |
| **WHAM** | **Monolithic orchestrator + NULL slots**: `prepare_wham_input(asap3, selectivity=NULL, M=NULL, NAA_re=NULL, ecov=NULL, …)` threads an input list through a fixed sequence of `set_*()` helpers. The **same `set_*()` are public and re-enterable** to overwrite one component of an existing input. | `NULL` = defaults; a *disabled* process is still materialized as a **present-but-inert** dummy (`set_ecov(NULL)` sets `n_Ecov=1` with a dummy obs, `Ecov_use_obs=0`). | `set_selectivity()` is the gold standard: enumerated-option checks, dim-consistency checks naming the exact element, mutual-exclusivity, plus an info log of auto-filled defaults. Errors surface in R. |
| **FIMS** | **Modular object graph + data-driven default table**: independent Rcpp modules linked by integer IDs; on top, `FIMSFrame(data)` → `create_default_configurations(data)` → `create_default_parameters()` produce a **tibble spec whose defaults scale to the data** (3 fleets ⇒ 3 selectivity modules), edited with `dplyr::rows_update()` before `initialize_fims()` materializes anything. | Absence = off (modules simply never instantiated). **Not available to compiled TMB.** | `FIMSFrame()` validates the data object; the spec is sized from it; runtime log surfaces C++ messages. |

**The compiled-TMB constraint that rules the choice:** a compiled TMB template's `DATA_*`
members must *all exist* — "off" cannot be *absent*, it must be a **present, flagged/zeroed
structure**. That kills FIMS's true-absence model for Rceattle and points squarely at
**WHAM's dummy-off**. (This is exactly Rceattle's existing `map → NA` convention.)

---

## 3. Recommendation

Adopt **WHAM's pattern (C), layered on a dimensional backbone (A), with WHAM's dummy-off
(E) and its declarative "how" linkage** — because it is the lowest-friction fit for a
*compiled-TMB, released, public-API* package, and it reuses the defaulting Rceattle
already has (`clean_data`/`switch_check`).

1. **Keep the monolithic list as the internal representation.** It already *is* the TMB
   bundle. Users stop touching it directly; it becomes the *output* of the builders. This
   is what preserves the public API — the big list can stay; only *how users populate it*
   changes. Existing scripts that build the list by hand keep working.

2. **Add a dimensional-backbone constructor** — the one thing all three references have
   and Rceattle lacks (there is no `build_data()` today):

   ```r
   d <- ceattle_data(nspp = 1, styr = 1977, endyr = 2024, projyr = 2100,
                     nages = 15, nsex = 1, fleets = c("Fishery","BTS","ATS"))
   ```
   This fixes species/years/ages/fleets and becomes the single source of truth every
   later component validates against. Especially valuable given Rceattle's multispecies
   dimensions.

3. **Refactor each subsystem into a re-enterable `set_*()` helper**, dispatched from one
   orchestrator via `NULL`-default args (WHAM's dual-use trick):

   ```r
   d <- d |>
     set_catch(catch_df) |>
     set_index(index_df, cov = list(BTS = Sigma)) |>   # <- covariance survey, see §5
     set_composition(comp_df) |>
     set_selectivity(Fishery = "NonParametricPM", BTS = "LogisticPM") |>
     set_weight(waa) |> set_maturity(mat) |> set_M(m1)
   # optional processes — only mentioned if used:
   d <- set_predation(d, ration_df, diet_df)          # pulls in msmMode + bioenergetics
   d <- set_environment(d, env_df, linkages = ...)     # pulls in the DSEM linkage table
   ```
   One call builds sensible defaults; the *same* function re-modifies one slice later
   (`d <- set_selectivity(d, BTS = "Logistic")`). New `NULL` args are non-breaking →
   satisfies the deprecation policy. Each helper *owns its defaults and validation*.

4. **Represent "off" as present-but-inert, never absent** — a dummy structure +
   `use_* = 0` + `map → NA`, so the compiled model always gets a valid slot. This is
   already Rceattle's convention, and it's exactly what the covariance-matrix input in §5
   does (a 1×1 dummy precision matrix for non-MVN fleets).

5. **Move validation into the `set_*` helpers** with `stop()` messages that name the exact
   element and the dimension it violated — killing the late `rearrange_data()` failures.
   `data_check()` remains the aggregate backstop.

6. **For the climate/multispecies niche specifically**, adopt WHAM's declarative "how"
   linkage: a matrix/spec saying *which env covariate drives which process*
   (`recruitment_how`, `M_how`, `q_how`, …) with `"none"` = off, instead of scattering
   env flags across the list. This is essentially your DSEM linkage table, formalized —
   the single most transferable idea for the differentiator.

### What this explicitly does **not** do
- It does **not** delete or rename the monolithic list, the numbered `R/` pipeline, or any
  exported `build_*()` arg. Golden-reference equivalence is preserved throughout.
- It does **not** force the Excel path away; but see Phase 1 for making its sheets
  genuinely optional.

---

## 4. Phased rollout (each phase independently shippable — no big-bang rewrite)

- **Phase 0 — worked exemplar (DONE this session).** The new covariance (MVN) survey
  likelihood input was built following every principle above: an optional, name-keyed
  `index_cov` list (present only if used), a per-fleet `Index_loglike` flag defaulting to
  the back-compatible `"Lognormal"`, a **present-but-inert** 1×1 dummy for non-MVN fleets,
  and **early validation in `data_check()` keyed to the flag** ("fleet X is MVN but no
  covariance matrix / wrong dimension / not symmetric"). It is the template for `set_index(..., cov=)`.
- **Phase 1 — make optional data genuinely optional (low risk, high relief).** Guard every
  sheet read in `read_data()` behind `if (sheet %in% sheetnames)` and route missing
  optional sheets through the existing `clean_data()`/`switch_check()` defaulting (the
  programmatic path already does this). Ship a **minimal single-species template** (~8
  sheets) beside the full one. This alone removes most of the "I had to supply diet data
  for a single-species model" pain.
- **Phase 2 — backbone + `set_*` helpers.** Add `ceattle_data()` and wrap the existing
  `clean_data`/`switch_check`/`data_check` defaulting in `set_*()` helpers (they populate
  slices of the list that already exists). This is additive; the hand-built-list path
  still works.
- **Phase 3 — fold the hidden switches into the data object (the open fork — see §6).**
- **Phase 4 — declarative env "how" linkage table** for the climate niche.

Each phase updates the four schema sources together (or, better, generates the roxygen /
template / glossary from one canonical schema definition — a stretch goal that would end
the drift permanently).

---

## 5. The covariance-matrix input as the pattern in miniature

The Phase-0 work is a concrete, already-merged illustration of the target ergonomics:

- **Optional, keyed by name:** `data_list$index_cov <- list(BTS = Sigma)` — a fleet that
  isn't MVN never mentions it.
- **Flag-driven, back-compatible default:** `fleet_control$Index_loglike` defaults to
  `"Lognormal"`; existing models are numerically unchanged.
- **Present-but-inert for the compiled model:** TMB always receives a length-`n_flt` list
  of precision matrices; non-MVN fleets get a 1×1 dummy that is never read.
- **Validated at the flag, early:** `data_check()` errors with a message naming the fleet
  and the exact dimension mismatch — not a cryptic failure inside `rearrange_data()`.

Under the redesign this becomes `set_index(d, index_df, cov = list(BTS = Sigma))`, and the
`Index_loglike = "MVN"` flag is inferred from the presence of `cov` for that fleet.

---

## 6. The one open design fork (needs Grant's call)

**Is Phase 3 in scope — should the `fit_mod()` run-control switches (`msmMode`, `initMode`,
`HCR`, `srr_*`, `M1_*`, `growth_*`) be absorbed into the data object?**

- **Yes (unify):** one object fully specifies a model; "what do I provide?" has a single
  answer; a run is reproducible from one saved object. But it blurs *data* (observations)
  and *model choices* (configuration), and it is the largest change to `fit_mod()`'s
  signature.
- **No (keep separate on purpose):** data = observations, `fit_mod()` args = model choices.
  Cleaner conceptual split; smaller blast radius. But the "hidden second surface" problem
  remains — users still discover half the config through `fit_mod()`/`build_*()`.

A middle path: keep `fit_mod()` args authoritative, but let the data object optionally
*carry* a `model_config` slot (defaulted, overridable at fit time) so a saved object is
self-describing without changing `fit_mod()`'s contract.

My weak recommendation: **the middle path** — it captures reproducibility without the API
churn or the conceptual blur. But this is genuinely your call.

---

## 7. Immediate, low-cost wins (do regardless of the larger redesign)

1. Guard the `read_data()` sheet reads (Phase 1 core) — a few `if (%in% sheetnames)`.
2. Ship a minimal single-species Excel template.
3. Regenerate `meta_data_names.xlsx` from the current `switch_check` column set (it predates
   the `Sel_*` columns) — or drop it in favor of the roxygen.
4. Add the missing `fleet_control` presence guards in `rearrange_data()` (the unguarded
   column pulls) so a hand-built table fails in `data_check()` with a clear message, not
   deep in reshaping. (This session already hardened one such path: the switch converters
   now default a missing `Index_loglike` column instead of erroring.)
