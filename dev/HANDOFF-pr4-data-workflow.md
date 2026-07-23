# Handoff — PR 4: requirement table, `data_requirements()`, `build_data()`, `model_config`, spec tree

Start a **fresh session** for this. PR 4 is a large, design-heavy PR (new user-facing
data-construction surface + introspection) with open design questions — it deserves its
own plan-first pass, not a tail-of-marathon start. This handoff is self-contained.

## Where things stand (as of this handoff)

Branch **`dev-data-workflow`**, working tree clean, all suites green, golden **bit-identical**
(pinned in the gitignored `dev/golden-ref.rds`: BS2017SS = 10241.0304275402,
BS2017MS = 10304.2212737156; the plan/`/golden-check` documents the same MS value).

**Completed this session** (the formula-linkage grammar is now productionized + a naming sweep):
- **Grammar densities** shipped earlier (IID/rw/ar1 RE, DM comp priors `build_composition`).
- **Naming sweep** (intuitive-naming priority): `est_M1` silent-loss bug fix + dictionary
  reconciliation; readable string aliases (`estDynamics`, `Estimate_*_sd`, `srr_est_mode`,
  `suitMode`); `Time_varying_{q,sel}_sd_prior → _sd` and `Q_prior/Index_sd_prior/Catch_sd_prior
  → Q_init/Index_sd/Catch_sd` renames (dual-path alias, `.rda` regen); soft-deprecation of
  `Time_varying_q/sel`/`M1_re` toward the grammar. All adversarially reviewed.
- **`M1_re = 6`** separable 2D-AR1 rho-mapping bug fixed (was collapsing to IID).
- **QAR1 (Rogers 2024)** productionized: effect sizes reported (`beta_linkage_obs` + ADREPORT
  SE); env covariates need not span the series (auto-extend + per-slot observed mask +
  fixed-effect NA guard in `materialize_linkage`); `obs_sd` fixed by default, **estimable via
  `linkage_spec(obs_sd_est = TRUE)`** (opt-in — it degenerates on smooth covariates).
- **Selectivity penalty SDs**: `Sel_curve_pen1/2/3` weights expressible as SDs
  (`Sel_shape_sd`+`Sel_shape_dir`, `Sel_curvature_sd`, `Sel_devmag_sd`), `weight = 1/(2·sd²)`,
  for NonParametric (2/9) AND LogisticPM (11); converted in `switch_check`, legacy untouched.

**Deferred / not this PR** (recorded in the plans): the full **GOApollock ADMB bit-identity
bridge** (`dev/PLAN-qar1-productionization.md`, north-star); a **prior on `obs_sd`** to
regularise the identifiability degeneracy; SD treatment for the **2D/3D-AR1** sel forms (they
reuse `Sel_curve_pen` for logit-ρ, not weights); the `Sel_curve_pen` 4-way overloading design
(`dev/PLAN-sel-curve-pen-reformulation.md`); auto-translators for legacy switches (decided
against — soft-deprecation instead).

## PR 4 scope (from `dev/PLAN-data-workflow-and-linkage-grammar.md` §"PR 4", ~line 566)

**Key enabling insight: the required/optional logic already exists internally** — this PR is a
*thin user-facing layer over existing logic*, not a rewrite:
- `clean_data()` already default-fills every optional block (`comp_data`, `caal_data`,
  `emp_sel`, `NByageFixed`, `ration_data`, `diet_data`, `env_data`, `index_cov`, …).
- `data_check()` already encodes the *conditional* requirements (gating on `msmMode`,
  `growth_model`, `estDynamics`, `suitMode`, per-fleet `Catchability`/`Selectivity`/
  `Time_varying_*`, env-index refs).
- `switch_check()` fills switch defaults + derives values (e.g. `Sel_start_year`).

**Sequence (each step enables the next):**
1. **Extract the requirement table** — refactor `data_check()`'s gates into a declarative spec
   (one row per `data_list` element + the condition under which it's required); `data_check()`
   then *consumes* the table. **No behaviour change — must hold golden equivalence.** The
   enabling refactor; do this first and verify bit-identical.
2. **`data_requirements()`** — a reader over the table: Required / Optional-and-defaulted /
   Ignored for a given model spec. Pure introspection, no risk.
3. **`build_data()`** — constructor taking only the blocks a model uses, defaulting the rest via
   `clean_data()`, validating against the same table. `read_data()`/`write_data()` untouched.
   Target ergonomics:
   ```r
   dat <- build_data(nspp = 1, styr = 1977, endyr = 2023, nages = 10,
                     maturity = mat, sex_ratio = sr, fleet_control = fc,
                     catch_data = catch, index_data = survey, comp_data = ages)
   fit_mod(dat)
   ```
4. **Per-process builders** (`build_population() |> add_fishery() |> add_survey() |> …`) — sugar
   over `build_data()`, ONLY if demand justifies the surface (largest new surface to test).

**Also in PR 4:**
- **`model_config`** — a defaulted, validated data-object slot carrying `msmMode`, `initMode`,
  `HCR`, `srr_*`, `M1_*`, `growth_*`, `comp_offset`. **`fit_mod()`'s signature stays unchanged**;
  its args override the slot when explicitly supplied — detect with `missing()`, not a sentinel.
- **`print()`/`summary()` spec tree** — indented tree (PyTorch `nn.Module`-style): dimensions →
  fleets → per-process form → active linkages (with formulas) → estimated/fixed/mapped-out.
  Extends the S3 methods in `R/0-rceattle_class.R` (~L20, L64). Highest-leverage navigability item.
- **`sanity()` + parameter dictionary** (already landed internally on `bugfix-correctness`) —
  `sanity()` must refuse to run on a build-only (`estimateMode = 3`) object, not report on a
  placeholder.

**Invariant any constructor must preserve:** `fleet_control$Fleet_code` must equal the row
number (enforced in `data_check()`; see CLAUDE.md).

## Open design questions (decide first, via AskUserQuestion / EnterPlanMode)
1. Does `build_data()` validate **eagerly** or defer to `data_check()` at fit time? (Deferring =
   one source of truth; eager = better messages.) *Lean: defer for truth, add a thin eager
   pre-check for message quality.*
2. Is the requirement table **internal or exported** (so users can query it)? *Lean: export it —
   `data_requirements()` is more useful if the table is inspectable.*
3. For `msmMode > 0`, does `build_data()` take predation/diet blocks as optional args, or a
   separate `add_predation()`? *Lean: optional args on `build_data()`, defer the builder chain.*
4. Is `write_template()` (emit only the sheets a config needs) worth it for the xlsx user base,
   or does the code-first path supersede it? *Lean: skip initially.*

## How to work it (governing standards for this repo)
- **Plan first, get approval** (EnterPlanMode/ExitPlanMode), then **incremental verified
  sub-commits** — each compiles, `/golden-check` bit-identical (step 1 especially — it's a
  refactor), fast suite green.
- **Federal-rigor standard** (memory `federal-management-rigor`): statistical accuracy verified
  constructively; every change golden/constructively verified; edge cases noted + tested, not
  "it runs". Spawn an **adversarial reviewer** after each substantive step to refute against
  source before calling it done.
- **Released package — preserve the public API.** `read_data()`/`write_data()` untouched; new
  surface is additive; `fit_mod()` signature unchanged.
- **Commits: plain messages, no AI-attribution trailer.** Prefix R with `export PATH=/usr/bin:$PATH`.
  After a schema/jnll change, `R CMD INSTALL .` before local parallel tests (retro/jitter/mse).
  Sweep the sibling repo `../Rceattle-models` after any API change.
- Recommended **starting point: step 1 (requirement-table refactor)** — golden-preserving, and it
  unblocks 2–4. It also forces reading `data_check()` end to end, which is the map you need.

## Key files
- `R/1-data_check.R` (the gates to extract), `R/0-switches.R` (`switch_check`, defaults/derivation),
  `R/0-clean_data.R` (`clean_data` optional-fill), `R/0-read_write_excel_data.R`,
  `R/0-rceattle_class.R` (S3 print/summary), `R/6-fit_mod.R` (arg override + `model_config`).
- Plan: `dev/PLAN-data-workflow-and-linkage-grammar.md` (PR 4 ~L566; PR 5 optional-Excel/schema,
  PR 6 `save_config`/`load_config`, PR 7 contributor docs + C++ legibility follow).
