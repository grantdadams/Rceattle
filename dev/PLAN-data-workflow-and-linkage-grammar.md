# Rceattle: one formula grammar for every process + data-workflow redesign

## Context

As Rceattle has grown (predation, DSEM env linkages, time-varying selectivity, OSA, growth),
the input surface has grown with it — but any single model needs only a subset. Grounded
counts from this session:

- **~100–110 user-facing switches**: ~27 enum/flag columns in `fleet_control` (48 distinct
  names referenced in R; shipped `BS2017SS$fleet_control` has 32), ~13 `data_list` switches
  defaulted in `switch_check()`, plus 19 `fit_mod()` args, 12 `build_srr()`, 9 `build_M1()`,
  11 `build_hcr()`, 16 `fit_control()`.
- **143 `DATA_*` / 48 `PARAMETER_*`** declarations in `ceattle_v01_11.cpp:91-432`.
- Schema drift, an Excel-vs-programmatic asymmetry, and late cryptic validation
  (diagnosed in `Rceattle-models/DATA_WORKFLOW_RFC.md`, 2026-07-05).

The key finding: **the formula-driven linkage system is already generic and already reserves
`q` and `sel`** ([0-linkage_table.R:56](R/0-linkage_table.R#L56),
[linkage.hpp:35-36](src/TMB/linkage.hpp#L35-L36)), and its schema already carries an
**unused `re_group` column** ([0-linkage_table.R:49](R/0-linkage_table.R#L49)). So this is
extension along a designed-for seam, not new architecture.

Intended outcome: a fisheries scientist can express *any* effect on *any* process —
covariate, time block, annual deviate, random walk, AR1 — in **one formula grammar**; see
their whole model in one printout; supply only the data their model needs; and git-diff two
model configurations.

### Decisions locked with Grant this session
1. **Linkages first, data ergonomics follows.**
2. **Config home = middle path**: `fit_mod()` args stay authoritative; the data object
   optionally carries a defaulted `model_config` slot.
3. **Back-compat = deprecate with warning, translate internally** to one TMB code path.
4. **Navigability = all four**: spec tree, `required_inputs()`, diffable text config,
   contributor docs + C++ legibility.
5. **Fixed vs random = lme4 bar notation** (see grammar below).
6. **`~ Year` = linear trend; `~ factor(Year)` / `~ (1|Year)` = annual deviates.** No magic.
7. **Absorb `Time_varying_sel` / `Time_varying_q` / `M1_re` into the grammar**, deprecating
   and translating them.

### Findings from the adversarial verification pass (these corrected the first draft)

Six claims in the draft were wrong or overstated. The corrections shape the plan:

- **`sel` and `q` are not symmetric today.** `LINKAGE_PARAM_CODES$q = c(q = 0L)` exists, so a
  `q` row *encodes fine* and then contributes **nothing** to the model while still being
  estimated and prior-penalized — a silent latent bug. `sel = character(0)` means a `sel` row
  *errors*. ([0-linkage_encode.R:73-74](R/0-linkage_encode.R#L73-L74))
- **`link = "logit"` is a live silent bug.** R validates code 2; no accumulator implements it
  ([linkage.hpp:114,159,228,292,369,416](src/TMB/linkage.hpp#L114)). A logit row is
  estimated, contributes exactly zero, **and** `map_linkage_adjuster` still masks the base
  parameter for a slope-only group — so the process parameter is fixed at init with nothing
  replacing it. Fix this regardless of everything else.
- **`map_linkage_adjuster()` is called FIFTH, not last** —
  [3-build_map.R:44](R/3-build_map.R#L44), *before* `build_map_selectivity` (`:47`) and
  `build_map_catchability` (`:49`). Harmless today (it only writes growth/M/rec params), but
  **adding `q`/`sel` branches at the current position would be silently clobbered.** The
  ordering must change as part of PR 1.
- **`base + dev(yr)` is NOT universal in `selectivity.hpp`.** Additive: cases 1, 2, 3, 4, 5, 8.
  **Case 0 (Fixed/empirical) has no parameter slot at all**; **case 9 NonParametricRPM is a
  random walk in years and renormalizes twice, so a bin-constant offset cancels entirely**;
  **case 11 LogisticPM is multiplicative** (`sel_inf * exp(dev)`, `:516`); cases 6/7 are
  additive but inside a logit. Linkage support must be scoped per form, not assumed.
- **The `Estimate_q` comma-list is parsed for code 5 only**, not 5-or-6
  ([3-build_map.R:1088-1101](R/3-build_map.R#L1088-L1101)). Separately, a real bug:
  `Time_varying_q = "1,3"` becomes `NA_integer_` in the TMB `index_varying_q` vector via
  [5-rearrange_data.R:173-174](R/5-rearrange_data.R#L173-L174).
- **`read_data()` guards are the majority, not the exception**: 14 of 23 `read_xlsx` sites are
  guarded; **9 unguarded sites cover 12 sheets** — `control` (`:193`), the `matrix_data` loop
  over `fleet_control`/`comp_data`/`emp_sel`/`NByageFixed` (`:221`), `age_trans_matrix`
  (`:320`), `age_error` (`:325`), `sex_ratio` (`:361`), `M1_base` (`:365`),
  `bioenergetics_control` (`:375`), `env_data` (`:383`). The last four plus `:221` are the
  real hazards for a minimal workbook.

---

## External validation: the FIMS CIE peer review (April 2026)

The three independent CIE reviewers of FIMS (Holmes, Howell, Ono) plus the summary report
converged on a short list of recommendations. Rceattle is a different tool at a different
maturity, but every one of the reviewers' *cross-cutting* themes is something this plan
already targets — which is a useful independent check that the priorities are right, and in
two places it sharpens them.

**Direct hits — themes all reviewers raised, mapped to where this plan delivers them:**

1. **"Make users aware of default versus other settings, and their implications"** — the
   single most-repeated ask (all four documents, TOR 3). This *is* the PR 4 spec tree +
   `data_requirements()` + parameter dictionary: turn "unknowingly accept defaults" into a
   printout that shows every process's form, what is estimated vs fixed vs mapped-out, and
   which linkages are active. Ono explicitly wants defaults' *consequences* stated, not just
   their values — the spec tree should annotate, e.g., that `~ (1|Year)` recruitment smooths
   episodic recruitment (his worked example of a harmful default).
2. **Estimation-failure warnings must be unavoidable** — every reviewer, emphatically
   ("watermark on figures", "header in tables", "not hidden among routine messages"). This is
   the same principle as PR 0–2's design rule: **a mis-specified linkage errors loudly rather
   than silently contributing zero.** The `.check_sel_linkage_support()` guard and the
   logit/`q`-noop fixes in PR 0 are exactly this instinct applied to model *construction*; the
   reviewers want it applied to model *fitting* too. → adds a diagnostics item below.
3. **"Avoid feature sprawl and overlapping, redundant functionality"** (Howell and Ono both,
   TOR 5) — the entire "one grammar for every process" thesis and the condensing pass. Ono:
   "rename options to improve interpretability, e.g. `constant` → `not estimated`" — a direct
   nudge for the parameter dictionary's plain-language column.
4. **sdmTMB-style tiered documentation** (Ono and the summary name sdmTMB explicitly) — model
   description with equations, getting-started, one topic-based section per component
   (selectivity, growth, likelihoods…), estimation+diagnostics, projections. This is the
   shape PR 7's vignettes should take; the plan already cites sdmTMB's `sanity()`.
5. **Version-controlled, diffable, archivable model configuration** (Ono, TOR 5 "streamlined
   outputs and configuration tracking") — precisely PR 6 `save_config()`/`load_config()`.
6. **A concise reference manual complementing vignettes** (Holmes, Howell, summary) — the
   "diffable text config with doc strings" of PR 6 plus the canonical schema of PR 5 generate
   most of this for free.

**Two sharpenings this plan should absorb:**

- **Add a diagnostics/warnings item to PR 7** (or a small PR 7.5): a single `sanity()`-style
  reporter that runs after `fit_mod()` and flags, unmissably, non-convergence, a
  non-invertible Hessian, and parameters at bounds — the reviewers' universal ask. Rceattle
  already computes a `hessian_conditioning` diagnostic (improved on `dev-ebs-pk`); the gap is
  making its verdict prominent and attaching it to plotted/tabulated output rather than
  leaving it in the console. Cheap, high-leverage, and it is the thing three independent
  experts each singled out.
- **OSA / one-step-ahead residuals for composition data** (Howell, Ono): Rceattle already has
  `comp_osa.hpp` and `0-osa_data.R`, so this is largely a plotting-surface question, not new
  model code — worth confirming it is exposed and documented in the PR 7 diagnostics vignette.

**Deliberately out of scope** (FIMS structural asks that are not Rceattle gaps or not this
work): seasonality, sex/sub-stock partitions, and reference-point machinery are FIMS
catch-up items; Rceattle already has sex, SPR proxies and HCR projection. Noted only so the
reports are not mis-read as a task list for this branch.

## The grammar

One grammar, five (soon six) processes. `env_data` already carries a `Year` column, so the
time terms need no new data.

```r
~ temp                              # fixed covariate effect
~ temp + PDO                        # several covariates
~ Year                              # linear time TREND (standard model.matrix semantics)
~ block(Year, breaks = c(1990, 2005))   # time blocks -> fixed effects
~ factor(Year)                      # annual deviates, FIXED, no sigma parameter
~ factor(Year)  + priors = list(...)    # annual deviates, INPUT sigma (penalized)
~ (1 | Year)                        # annual deviates, ESTIMATED sigma (random effect)
~ rw(Year)                          # random walk, estimated sigma
~ ar1(Year)                         # AR1, estimated sigma + rho
```

**The fixed/random rule — one rule, no new switch.** Anything *outside* a bar is a **fixed
effect**: it has no sigma parameter, and if you want it shrunk you attach a **prior** — that
is the "input sigma" case. Anything *inside* a bar `(1|Year)` or a structured wrapper
`rw()`/`ar1()` is a **random effect**: sigma becomes an estimated parameter and the
coefficients enter `MakeADFun(random=)`. The existing `re_group` column carries the bar's
grouping name; it is currently written but never read, so this is the consumer it was
designed for.

**Identifiability.** This lands exactly on the ADMB `dev_vector` rule in CLAUDE.md:
`~ factor(Year)` uses treatment contrasts, so `model.matrix()` drops the first level and the
vector is identifiable. `~ (1|Year)` / `rw()` / `ar1()` keep all levels because the density
pins them. `~ factor(Year)` with a prior may also keep all levels (`~ 0 + factor(Year)`),
since the penalty pins it.

**Parser — REVISED after the package survey.** Do **not** hand-roll `findbars()`. As of
lme4 2.0.x those helpers were split into the standalone CRAN package **`reformulas`**, now
shared by lme4 *and* glmmTMB, and it already handles `/`, `:`, `||`, `0+`, and the
special-term syntax `ar1(times + 0 | group)`. Add `reformulas` to `Imports:` and call it at
the single existing seam [0-build_linkage.R:483](R/0-build_linkage.R#L483). It is small,
focused, and maintained by Bolker — reimplementing it would be strictly worse.

**Term naming — REVISED.** Adopt **glmmTMB's covariance-structure vocabulary verbatim**
rather than inventing one: `ar1()`, `ou()`, `exp()`, `mat()`, `toep()`, `cs()`, `us()`,
`diag()`, plus the `hom*` prefix for homogeneous-variance variants. Syntax is
`struct(terms + 0 | group)`, so the grouping is explicit in the term and the separate
`by = ~fleet` argument becomes unnecessary for structured terms. Back these with a small
integer registry shared with C++ — exactly the pattern `LINKAGE_PROCESS_CODES` /
`RCEATTLE_PROC_*` already uses.

**THE TRAP TO NOT INHERIT.** glmmTMB's `ar1(times + 0 | group)` requires `times` to be a
**factor**, so lag distance is *column-index difference*, not real elapsed time. Gappy or
mis-ordered levels silently produce a plausible-looking ρ that means something else — the
covstruct vignette calls it "a common mistake." `numFactor()` exists purely to smuggle
numeric coordinates through that channel. A fisheries `Year` is genuinely numeric and
genuinely gappy, so Rceattle must take **numeric `Year`, compute true lags**, and error
loudly on gaps it cannot represent. This is the single most valuable correction the survey
produced.

**Legacy translation table** (deprecate-with-warning, one TMB code path):

| Legacy | Formula equivalent |
|---|---|
| `Time_varying_sel = "Off"` / `Time_varying_q = "Off"` | `~ 1` |
| `= "IID"` | `~ (1 \| Year)` |
| `= "RandomWalk"` | `~ rw(Year)` |
| `= "AR1"` | `~ ar1(Year)` |
| `= "Block"` | `~ block(Year, ...)` (breaks from `Selectivity_block`) |
| `Estimate_q = 5` + `Time_varying_q = "1,3"` | `linkage_spec(~ env1 + env3, by = ~fleet)` on `q` |
| `M1_re = "iid_year"` / `"ar1_year"` / `"ar1_age"` … | `~ (1\|Year)` / `~ ar1(Year)` / `~ ar1(age_bin)` on `M1` |

## Branch

Cut **`dev-data-workflow`** from `dev`. The existing `accessibility-and-code-review` branch
is empty relative to `dev`; its content lands here as PR 7.

## Standing constraint: take runtime wins where they cost no legibility

Grant's rule for this work: reduce run time opportunistically as we pass through code, but
**never at the expense of a fisheries scientist being able to read it.** A faster line that
needs a comment explaining the trick is usually the wrong trade; a faster line that is also
the more obvious line is free. Anything that moves numbers still needs `/golden-check`.

Concrete opportunities the reviews already identified, ordered by (win / risk):

**Free — no numeric change, code gets clearer or stays the same**
- `convert_length_selectivity()` takes `array<Type> sel_length` and `array<Type>
  growth_matrix` **by value** ([selectivity.hpp:156-157](src/TMB/selectivity.hpp#L156)) and is
  called inside a `for(yr) for(sex)` loop. `growth_matrix` is
  `(nspp*2+n_flt, max_sex, max_age, max_nlengths, nyrs)` — a multi-megabyte deep copy per
  call, per objective evaluation. Change to `const array<Type>&`. **Biggest single win in the
  template.** Same pattern, lower frequency: `calculate_selectivity()` (9 arrays by value) and
  `estimate_growth_within_yr()` (3 by value, while its sibling `estimate_growth()` already
  takes the same objects by reference).
- Hoist loop-invariant work: the M1 env product recomputed inside the `age` loop though it
  depends only on `(sp, sex, yr)`; the entire HCR `switch` inside `for(age) for(sex) for(yr)`
  though `proj_F` depends only on `(sp, yr)`; `suit_other` assigned in the innermost of a
  six-deep nest; the plus-group geometric series recomputed per `(sex, yr)`.
- Hoist the `vector<Type> inv_pred_denom(nyrs)` allocation out of the `r_age` loop.
- R side: replace `x <- c(x, ...)` / `rbind`-in-loop growth (empirical-selectivity expansion,
  `mse_summary` accumulating six vectors in a doubly-nested loop — ~18,000 full copies for
  100 sims x 30 years). `R/0-osa_data.R:86-125` already does this correctly; use it as the
  reference.
- `clusterExport(varlist = ls(envir = environment()))` ships the entire function frame —
  `om`, `em`, every intermediate — to every worker. Export an explicit varlist.

**Cheap, needs a golden-check**
- Cache-hostile loop ordering: TMB arrays are column-major, but `predation.hpp` and
  `bioenergetics.hpp` put `yr` (slowest-varying) innermost. `calculate_msvpa_predation()`
  already gets this right with `yr` outermost.
- `process_residuals()` densifies the inverse of the full joint precision matrix to read one
  block — solve against a sparse selector instead.
- `.osa_sdnr_tails()` rebuilds an `nsim x n` matrix per group with a fixed seed, so identical
  work repeats for every group of the same size. Memoize on `n`.

**Structural, only if it pays for itself**
- A `refit()` fast path (lme4's idea): reuse the compiled AD object and start from the fitted
  parameters instead of re-running the whole `fit_mod()` pipeline. Retro peels, jitter and MSE
  loops currently rebuild everything each time. Real speedup for the `9-*`/`10-*` code, but it
  is new public surface — only worth it if those loops are actually a bottleneck in practice.
- The `fit_mod()` re-invocation block is copy-pasted 11 times across ~900 lines and has
  already drifted. Collapsing it to one `.refit_like(model, overrides)` is primarily a
  correctness and legibility win; any speedup is incidental.

**Explicitly not worth it:** micro-optimising the likelihood arithmetic, or anything that
replaces a readable expression with an index-arithmetic trick. `construct_Q()`'s O(n^3) dense
inverse is genuinely slow but only reachable for `Sel_type = 7`; leave it until someone uses
that path.

## The architectural constraint that rules the design (from the RFC, §2)

**A compiled TMB template's `DATA_*` members must all exist.** "Off" cannot be *absent* — it
must be a **present, flagged or zeroed structure**. This is what rules out FIMS's
true-absence module model for Rceattle and points at WHAM's dummy-off, and it is already
Rceattle's own `map -> NA` convention.

Everything in PR 4 has to respect this. `build_data()` may let a *user* omit a block, but
what reaches `MakeADFun` is always a full set of slots — the omission is realised as
`clean_data()` filling a defaulted, inert structure, not as a missing element. The
requirement table's three states map onto this directly: Required (user must supply),
Optional (defaulted to a present-but-inert structure), Ignored (present, never read).

### The worked exemplar already in the tree (RFC §5)

The MVN covariance survey likelihood is this pattern in miniature, already merged, and is the
template to copy rather than reinvent:

- **Optional and keyed by name** — `data_list$index_cov <- list(BTS = Sigma)`; a fleet that
  is not MVN never mentions it.
- **Flag-driven with a back-compatible default** — `fleet_control$Index_loglike` defaults to
  `"Lognormal"`, so existing models are numerically unchanged.
- **Present-but-inert for the compiled model** — TMB always receives a length-`n_flt` list of
  precision matrices; non-MVN fleets get a 1x1 dummy that is never read.
- **Validated early, at the flag** — `data_check()` errors naming the fleet and the exact
  dimension mismatch, rather than failing deep inside `rearrange_data()`.

Under PR 4 this becomes `set_index(d, index_df, cov = list(BTS = Sigma))` with
`Index_loglike = "MVN"` inferred from the presence of `cov` for that fleet.

## Deferred: condense the codebase (TODO, after the data workflow)

**Sequencing (Grant's call): do the data workflow first, then condense as a separate pass.**
Not woven into each PR. The `.refit_like()` investigation showed why — collapsing 11 call
sites that differ in six invisible ways is its own piece of work with its own failure mode,
and bundling it into a feature PR means a behaviour change disguised as a refactor. Keep the
two apart so each can be golden-checked on its own.

The risk of PRs 1-7 is that they *add* surface — `build_data()`, `data_requirements()`,
`model_config`, `save_config()`, the spec tree, `sanity()`. This pass is what makes the net
negative. The requirement table in PR 4 is the one exception that condenses *while* adding,
because it **replaces** the scattered conditionals in `data_check()` rather than sitting
alongside them; prefer that shape where it is available.

Rough sizes are from the review sweep and want confirming before each is attempted.

**Collapse candidates**

| Target | What | Approx |
|---|---|---|
| `9-retro_and_jitter.R`, `10-run_mse.R`, `10-project-no-F.R` | The `fit_mod()` re-invocation block, copy-pasted **11 times** — `build_hcr(11 args)` + `build_srr(13)` + `build_M1(10)` + `build_growth()` + 8 pass-throughs. **Needs the tabulation step first** (see below). | ~900 → ~120 |
| `5-rearrange_data.R` | 26 `pull()` blocks in three shapes; two helpers also put the 1-based→0-based conversion in **one auditable place** instead of four — the same class of bug as the year-index defect | ~110 → ~20 |
| `0-build_srr_and_M.R` | Three `.coerce_*` functions with the same four-branch skeleton (they even disagree on return type), plus four `.warn_*_deprecation()` from one template | ~130 → ~40 |
| `5-rearrange_data.R` | `check_composition_data()` / `check_caal_data()` ~90% identical — and the CAAL one's message says "Composition data have NAs" | ~90 → ~50 |
| `0-read_write_excel_data.R` | The sheet-reading idiom repeated 12×; the CRAN core-cap block 6×; PSOCK dispatch 5× | meaningful |
| `linkage.hpp` | Six accumulators sharing ~90% of their bodies — one templated helper parameterised on an index-mapping functor. **Do this before PR 1 adds two more**, or it becomes eight near-copies. | ~250 → ~90 |
| `growth.hpp` | `estimate_growth()` and `estimate_growth_within_yr()` are ~200-line near-duplicates; two sections are **byte-identical** | ~200 → ~120 |
| `7-plot_ceattle.R` | Six `plot_*` wrappers byte-identical but for two strings, each declaring and forwarding 21 formals; a `.ts_wrapper(output, ylab)` factory | ~300 → ~15 |

**Prerequisite for the `.refit_like()` collapse — do not skip.** The 11 sites are not
homogeneous. Each is ~40 arguments read from the *working* model with a few deliberate
exceptions read from the *pristine* original:

- OM site: 38 x `om_use$`, **5 x `om$`** (`suit_styr`, `suit_endyr`, `srr_mse_switchyr`,
  `srr_hat_styr`, `srr_hat_endyr`) — the OM pinning its suitability and stock-recruit
  reference periods to the original so its biology stays fixed through the projection.
- EM site: 42 x `em_use$`, **1 x `em$`** (`HCRorder`) — origin unclear, worth a decision.
- And one argument is **computed, not forwarded**: `srr_mse_switchyr = em_use$data_list$endyr`
  reads a *different field* and advances every iteration. A pure field-forwarder cannot
  express it; `.refit_like()` needs computed overrides.

So: tabulate all 11 sites x every argument x source object first; classify each exception as
intentional pin or accidental drift; then write the helper with the pins as **named
overrides** so a reader sees "this deliberately uses the original OM" at the call site
instead of inferring it from a variable name. Golden-check each site before and after.

(An earlier draft of this plan repeated a review claim that the OM and EM SRR arguments had
"drifted". They have not — those are two different models each reading its own family. The
real subtlety is the computed `endyr` substitution above.)

**Straight deletion — no replacement needed**

- `R/11-model_average.R`: **181 of 188 lines are commented out** (three abandoned uncertainty
  implementations, 52% of the file). Also a stray `plot_ssb()` call inside a non-plotting
  function that opens a device unbidden — and does it *before* `sdrep$sd` is populated, so it
  plots stale CIs.
- ~350 lines of commented-out Kinzey code in `predation.hpp`, plus the dead triplicate blocks
  in `2-build_params.R`, `4-build_parameter_bounds.R` and `OPT-phaser.R`, and two
  `growth_model == 3` branches nested inside `growth_model < 3` (unreachable).
- The 27-name MSE "quantities to keep" whitelist appears **3×** — and lists
  `"ssb_depletion"` twice in all three, which is the tell that it was hand-copied.
- `t_col()` in `7-plot_stock_recruit.R`: defined, roxygen-documented, called nowhere.
- Six no-op self-assignments (`proj_F(sp, yr) = proj_F(sp, yr);`) in the projection HCR switch.
- `R/dev/`: 12 files, 106 KB, including what look like abandoned forks of shipped code
  (`6-fit_mod_MSE.R`, `11a-mse_run_parallel_fast.R`). Correctly Rbuildignored, but a
  navigation hazard inside a numbered-collation directory. Move to top-level `dev/` or delete.
- Stale empty `tests/testthat/tests-*/` directories left from the flat-layout migration,
  each containing only a `_snaps` folder that `test_check` never reaches.

**The honest arithmetic.** The new surface across PRs 1–7 is roughly 1,500–2,000 lines
(grammar parser, requirement table, `build_data()`, config round-trip, spec tree, dictionary,
tests). The collapses and deletions above are ~2,500–3,000. Net should be **down**, with far
less duplication — but only if the collapses are actually done rather than deferred, which is
why each is tied to a PR that already has that file open.

---

## PR 0 — Fix the three latent bugs (small, independent, ships first)

Standalone because each is a real defect today and each is a prerequisite for the rest.

1. **Logit link.** Add the `linkfn == 2` branch to all six accumulators in `linkage.hpp`.
   Until then, make `link = "logit"` **error** in `validate_linkage_table()` rather than
   silently contribute zero.
2. **`q`/`sel` asymmetry.** Make an unconsumed `q` row error at `pool_linkages()` with a
   clear "not yet supported" message instead of estimating a no-op parameter.
3. **`index_varying_q` NA.** Fix the `as.integer("1,3")` → `NA` coercion at
   [5-rearrange_data.R:173-174](R/5-rearrange_data.R#L173-L174).

Regression tests for each. No numeric change to any valid existing model — `/golden-check`
must be bit-identical.

## PR 0.5 — Fifteen confirmed correctness bugs (approved; blocks all structural work)

Two independent code reviews found 15 verified defects. Grant's call: fix them **before** PR 1,
because that PR adds a fleet dimension and a second parameter vector on top of code where the
selectivity penalty, the M-linkage dot product and array indexing are already wrong — and
golden-reference equivalence only means something against a correct baseline.

**R — silently wrong results**
- [9-retro_and_jitter.R:444](R/9-retro_and_jitter.R#L444) — Mohn's rho writes `mohns[ind, ]`
  but reads `mohns[j, ]`. Retrospective bias wrong for every forecast year > 0. The N counter
  one line up is correctly `[ind, 3]`, which is what hides it.
- [10-run_mse.R:397](R/10-run_mse.R#L397) — projection CAAL built at `:390-396`, then written
  to the wrong field with the wrong payload. **CAAL absent from every MSE.**
- [8-sim_mod.R:90-95](R/8-sim_mod.R#L90-L95) — CAAL branch tests `Comp_loglike` with integer
  codes while the comp block at `:50-55` uses `CAAL_loglike` with string codes. Under string
  control both are FALSE and `:111` writes the previous row's draw.
- [8-sim_mod.R:218](R/8-sim_mod.R#L218), [9-retro_and_jitter.R:292](R/9-retro_and_jitter.R#L292)
  — `log(mean(log(R) - log(R_hat)))`; the inner term is already log-scale and centred near
  zero, so the mean is routinely negative → NaN into `rec_dev`. Sibling branches are correct.
- [0-build_hcr.R:137](R/0-build_hcr.R#L137) — `any()` where the comment and message say `all()`.

**R — errors / corruption**
- [0-combine_data_sets.R:71](R/0-combine_data_sets.R#L71) — merged `env_data` assigned to the
  wrong object and discarded; `:67` loops `mat_names[1:17]` over a 15-element vector, appending
  two `NA`-named elements.
- [0-read_write_excel_data.R:271](R/0-read_write_excel_data.R#L271) — inverted guard **zeroes a
  user-supplied `Month` column**.
- [10-mse_summary.R:58-64](R/10-mse_summary.R#L58-L64) — references `hcr_map`, which exists only
  as a local inside `fit_mod()`. Use the package constant at `R/0-switches.R:170`.
- [OPT-phaser.R:92,97](R/OPT-phaser.R#L92) — `log_M1 = 4` declared twice; `TMBphase()` pairs
  phases to parameters positionally across two differently-ordered lists.
- [6-fit_mod.R:862](R/6-fit_mod.R#L862) — `opt` undefined for `estimateMode = 2` + `NoFishing`.

**C++**
- `ceattle_v01_11.cpp:2762-2771` — type-2 selectivity curvature penalty writes
  `first_difference(first_difference(sel_tmp))` back into `sel_tmp`; only the `age = 0` term is
  the intended second difference. Copy the correct type-9 pattern at `:2843`. **[NUMERIC]**
- `ceattle_v01_11.cpp:1016-1023` — `beta_M1_tmp = M1_beta(sp,sex,i)` broadcasts a scalar across
  the vector, so the env-M1 dot product collapses to `beta_last * sum(env)`. **[NUMERIC]**
- `ceattle_v01_11.cpp:1788-1793` and 4 more sites — two consecutive `if`s where the second must
  be `else if`; a pre-`styr` observation year is transformed twice into a large negative index.
- `selectivity.hpp:101-103, 451-454` — max-normalisation branches on an AD `Type`, freezing the
  argmax at tape time. Use `max2()`, as `:74` already does. **[NUMERIC]**
- `predation.hpp:266-273` — guard tests `diet_prop_sum` but the division is by `denom`.

**Sequencing.** Land the non-numeric fixes first (golden-check must stay bit-identical), then
the numeric ones one at a time, each with its own `/golden-check`. Write regression tests for
`sim_mod()`, `mse_summary()`, the `read_data`/`write_data` round-trip, and `build_map()`
fleet-code sharing **before** touching `adjust_map_shared_params()` — three of those files have
no direct test coverage today, which is exactly why these bugs survived.

**Deferred to a later cleanup PR** (found, not fixed here): the `fit_mod()` re-invocation block
copy-pasted 11 times across ~900 lines and already drifted (`10-run_mse.R:667` pins SRR args to
`om` where `:864` uses `em_use`); the `1:n`-when-n-is-0 sweep (~50 sites); the index-vs-Fleet_code
conflation in `adjust_map_shared_params()`; `call. = FALSE` on 150 of 160 `stop()` calls.

## PR 1 — Fleet dimension + random-effect machinery + catchability linkages

The two structural additions the whole design rests on.

**Fleet stratum.** Linkage strata are `species`/`sex`/`age_bin`; q and selectivity are indexed
by **fleet**. Add a real `fleet` column, do not overload `age_bin`.
- `R/0-linkage_table.R` — `fleet = "integer"` in `LINKAGE_COLS:33`; arg on `linkage_row():183`;
  resolve in `.linkage_row_indices():237` (`NA` = broadcast); add to the `print` columns.
- `R/0-build_linkage.R` — allow `fleet` in `by`: edit the `unknown_by` check at `:488-493` and
  `expand_linkage_strata():521`; add a `fleet =` filter beside `species =`/`sex =` at `:578-585`.
- `R/0-linkage_encode.R` — emit `linkage_fleet` (0-based, sentinel `0` = all).
- `ceattle_v01_11.cpp` — `DATA_IVECTOR(linkage_fleet)` in the linkage block (`:217-228`);
  inject from `R/6-fit_mod.R:485-507`.

**Variances become linkage parameters (Grant's ask, 2026-07-21).** The random-effect forms
introduce `log_sigma_linkage` as the deviate variance, and — crucially — it flows through the
*same* `linkage_spec(init=, bounds=, priors=)` contract as any coefficient. So `init` inputs a
fixed variance, `priors` puts a prior on the variance, and the map decides fixed-vs-estimated —
one mechanism for both of Grant's asks. This is the whole reason to route variances through
the grammar rather than add more `fleet_control` columns.

Two facts this pass must respect / fix:
- **No prior on the deviate SD exists today.** `sel_dev_log_sd` / `index_q_dev_log_sd` are
  estimated-or-fixed but nothing shrinks them; the `jnll_comp(5)` terms are the deviate
  density, not a prior on sigma. The RE machinery adds the prior via the grammar.
- **`Time_varying_sel_sd_prior` / `Time_varying_q_sd_prior` are misnamed** — they are the
  *input value* of the SD ([2-build_params.R:292,370](R/2-build_params.R#L292)), not a prior.
  Rename with a deprecation alias as part of PR 5's schema cleanup; the grammar's `init=`/
  `priors=` split is the correct replacement.

**Random effects.** TMB's `random=` takes whole parameter names
([6-fit_mod.R:395-420](R/6-fit_mod.R#L395-L420)), so a subset of `beta_linkage` cannot be
made random. Split the vector:
- `PARAMETER_VECTOR(beta_linkage)` (fixed) + **`PARAMETER_VECTOR(beta_linkage_re)`** (added to
  `random_vars`), plus **`log_sigma_linkage`** (one per `re_group`) and
  **`trans_rho_linkage`** (one per AR1 group).
- The linkage table gains a `re_index` column saying which vector each row lives in and where.
- The RE density (IID / RW first-differences / `AR1()`) goes in a **new `jnll_comp` row 20** —
  which per CLAUDE.md means adding `"Linkage random effects"` to
  [R/6-rename_output.R:130-151](R/6-rename_output.R#L130-L151) in the same commit.
- Per-spec random/fixed replaces the global `random_rec`/`random_q`/`random_sel` logicals.
  This also fixes the gap the review found: `random_q = TRUE` with IID/RW currently integrates
  the deviates while leaving sigma **fixed** — a Laplace-integrated fixed-sigma model.

**Ordering fix.** Move `map_linkage_adjuster()` from
[3-build_map.R:44](R/3-build_map.R#L44) to run **after** `build_map_selectivity` (`:47`) and
`build_map_catchability` (`:49`), so a linkage genuinely overrides `fleet_control`. Verify the
existing growth/M/rec behavior is unchanged by the move (`/golden-check`).

**Catchability.**
- `rceattle_apply_q_linkages` + `..._natural` in `linkage.hpp`, shape `[n_flt, nyrs]`, copying
  the body shape of `rceattle_apply_M_linkages:198` exactly (process gate, sentinel expansion,
  `min(nyrs, linkage_X.rows())` clamp).
- Consume at [ceattle_v01_11.cpp:614](src/TMB/ceattle_v01_11.cpp#L614):
  `index_q(flt,yr) = exp(index_log_q + index_q_dev + q_offset(flt,yr)) + q_offset_nat(flt,yr)`.
- New exported `build_catchability(linkages = ...)` mirroring `build_M1()`; `fit_mod()` stashes
  `data_list$q_linkages` and adds a 4th `spec_groups` entry at
  [6-fit_mod.R:321-323](R/6-fit_mod.R#L321-L323), with `fleet` in `strata`.
- `q =` branch in `map_linkage_adjuster():1495` masking `index_log_q`.
- Translators for `Estimate_q = 5` (comma-list → covariate spec) and `Time_varying_q` modes
  (→ `(1|Year)`/`rw()`/`ar1()`), each emitting one `deprecate_warn()` naming the equivalent.

**Hard gate:** legacy `Estimate_q = 5` and `Time_varying_q = "AR1"` models must fit
**bit-identically** through the translated path.

**Legibility note (judgment call).** These accumulators already take 14 positional arguments
and this PR adds more. Bundle the encoded table into a plain POD struct
`rceattle_linkage_data<Type>` passed by const reference. This is a parameter object, not the
OOP functor hierarchy the accessibility handoff's guardrails prohibit — flagging it explicitly
because it sits near that line.

## PR 2 — Selectivity linkages

- Populate `LINKAGE_PARAM_CODES$sel` (currently `character(0)`,
  [0-linkage_encode.R:74](R/0-linkage_encode.R#L74)): `slp_asc`, `inf_asc`, `slp_desc`,
  `inf_desc`, `peak`, `sigma_asc`, `sigma_desc`, `right_floor`, `coff`.
- `rceattle_apply_sel_linkages` + `..._natural`, shape `[n_sel, nsex, nyrs, n_sel_params]`.
- **Scope by form, per the review:**
  - *Supported (additive slot):* case 1 Logistic (`:342-343`), 2 NonParametric-Ianelli
    (`:385-386`), 3 DoubleLogistic (`:412-415`), 4 DescendingLogistic (`:426-427`), 5 Hake
    (`:440`), 8 DoubleNormal (`:470-473`).
  - *Supported with care:* case 11 LogisticPM — devs are **multiplicative**
    ([selectivity.hpp:516](src/TMB/selectivity.hpp#L516)), so the offset goes **inside** the
    same `exp()`. Cases 6/7 (2D/3DAR1) — additive but inside a logit, so `link = "logit"`
    (PR 0) is required and the effect is on the logit scale; document that.
  - *Not supported — error in `data_check()` naming the fleet and form:* case 0
    (Fixed/empirical — no parameter slot) and case 9 (NonParametricRPM — a year random walk
    that renormalizes twice, so a bin-constant offset cancels to nothing).
- `build_selectivity(linkages = ...)`, 5th `spec_groups` entry, `sel =` branch in
  `map_linkage_adjuster()` masking `log_sel_slp`/`sel_inf`/`sel_coff`.
- `Time_varying_sel` translator (IID/RW/AR1/Block/Off → formula), deprecate-with-warning.
- `M1_re` translator (7 codes → `(1|Year)`/`ar1(Year)`/`ar1(age_bin)`/separable), same pattern.

## PR 3 — Dirichlet-multinomial priors

There is currently **no prior of any kind** on `comp_weights` / `caal_weights` /
`diet_comp_weights` ([cpp:409,410,414](src/TMB/ceattle_v01_11.cpp#L409)), and no entries for
them in `R/4-build_parameter_bounds.R`. Rather than a bespoke `DM_prior` column, add a
**`"comp"` process to the linkage grammar** with params `theta_comp`, `theta_caal`,
`theta_diet` and `by = ~fleet`:

```r
build_composition(linkages = list(
  theta_comp = linkage_spec(~ 1, by = ~fleet,
                            priors = list(`(Intercept)` = normal(0, 1)))))
```

This costs ~15 lines of C++ because the prior loop already re-targets intercept rows onto the
base parameter ([cpp:3365-3395](src/TMB/ceattle_v01_11.cpp#L3365-L3395)) — add a
`RCEATTLE_PROC_COMP` branch pointing at `comp_weights(flt)`. The parameters are already on the
log scale (`DM_pars_comp = exp(comp_weights)`, `:411`), so `prior_lognormal()` maps directly.
Bonus: time-varying / covariate-driven DM weights come free from the same grammar.

Caveat to handle: `comp_weights` is **dual-purpose** — a plain likelihood multiplier under the
multinomial (`:2610`, `:2617`) and log-DM-theta under DM. Only attach priors when
`Comp_loglike == "DirichletMultinomial"`; error clearly otherwise. Also note it is stored
**un-logged** at [2-build_params.R:387](R/2-build_params.R#L387), unlike its siblings.

## PR 4 — Requirement table, `data_requirements()`, `model_config`, spec tree

Absorbs `HANDOFF-data-workflow-modernization.md` (deleted from the repo root once folded in
here). Its framing supersedes the earlier sketch of this PR.

**The key finding from that handoff: the required/optional logic already exists internally.**
It is only the user-facing surface that is missing, which is what makes this tractable on a
released package — a thin layer over existing logic, not a rewrite.

- `clean_data()` already default-fills every optional block (`comp_data`, `caal_data`,
  `emp_sel`, `NByageFixed`, `ration_data`, `diet_data`, `env_data`, `index_cov`,
  `MSSB0`/`MSB0`) — the engine already tolerates a minimal input.
- `data_check()` already encodes what is *conditionally* required, gating on `msmMode`,
  `growth_model`, `estDynamics`, `suitMode`, per-fleet `Catchability` / `Selectivity` /
  `Time_varying_*`, and environmental-index references.
- `switch_check()` fills switch defaults and derives values from the data (e.g.
  `Sel_start_year` from the first year of observations, grouped by `Selectivity_index`).

**Sequence** (each step enables the next):

1. **Extract the requirement table.** Refactor the `data_check()` gates into a declarative
   spec — one row per `data_list` element with the condition under which it is required.
   `data_check()` then *consumes* the table instead of hard-coding conditions. No behavior
   change; must hold golden equivalence. This is the enabling refactor.
2. **`data_requirements()`** — a reader over that table reporting Required /
   Optional-and-defaulted / Ignored for a given model spec. Pure introspection, no risk.
   (Earlier drafts of this plan called it `required_inputs()`; use the handoff's name.)
3. **`build_data()`** — constructor taking only the blocks a model uses, defaulting the rest
   through `clean_data()`, validating against the same table. `read_data()` / `write_data()`
   untouched.
   ```r
   dat <- build_data(nspp = 1, styr = 1977, endyr = 2023, nages = 10,
                     maturity = mat, sex_ratio = sr,
                     fleet_control = fc, catch_data = catch,
                     index_data = survey, comp_data = ages)
   fit_mod(dat)
   ```
4. **Per-process builders** (`build_population() |> add_fishery() |> add_survey() |> ...`) —
   sugar over `build_data()`, only if demand justifies the surface. Most self-documenting,
   largest new surface to test and maintain.

**Also in this PR:**

- **`model_config`** — a defaulted, validated slot carrying `msmMode`, `initMode`, `HCR`,
  `srr_*`, `M1_*`, `growth_*`, `comp_offset`. `fit_mod()`'s signature is **unchanged**; its
  args override the slot when explicitly supplied (detect via `missing()`, not a sentinel).
- **`print()`/`summary()` spec tree** — indented tree (PyTorch `nn.Module` repr as the model):
  dimensions → fleets → per-process form → active linkages (with their formulas) → estimated
  vs fixed vs mapped-out. Extends the S3 methods at
  [R/0-rceattle_class.R:20,64](R/0-rceattle_class.R#L20). Highest-leverage navigability item:
  it turns "read 600 lines of switch tables" into "read the printout".
- **`sanity()` + the parameter dictionary** (already landed internally on
  `bugfix-correctness`) — see the sdmTMB notes below; `sanity()` must refuse to run on a
  build-only object rather than report on a placeholder.

**Invariant a constructor must preserve:** `Fleet_code` must equal the `fleet_control` row
number (now enforced in `data_check()`).

**Open questions carried over from the handoff:**
1. Does `build_data()` validate eagerly, or defer to `data_check()` at fit time? Deferring
   keeps one source of truth; erroring early gives better messages.
2. Is the requirement table internal, or exported so users can query it?
3. For `msmMode > 0`, does `build_data()` take the predation/diet blocks as optional args, or
   is there a separate `add_predation()`?
4. Is `write_template()` (emit only the sheets a configuration needs) worth doing for the
   existing xlsx user base, or does the code-first path supersede it?

## PR 5 — Optional Excel data + one canonical schema

- Guard the **9 unguarded `read_xlsx` sites** listed in the findings above; route misses
  through the existing `clean_data()`/`switch_check()` defaulting. Ship a **minimal
  single-species template** generated from `required_inputs()`.
- Fix `read_data():207` swallowing non-numeric control entries via
  `suppressWarnings(as.numeric(...))` — a typo currently becomes `NA` and fails much later.
- **Kill the schema drift permanently:** define the `fleet_control` schema **once** (name,
  type, default, allowed values, TMB target, doc string) and generate the `R/data.R` roxygen,
  the `switch_check()` defaults, `meta_data_names.xlsx`, and the Excel template from it. This
  is the fix for the **22 used-but-undocumented** columns, the stale entries
  (`Acuumulation_age_lower/upper`, `Catch_units`, `Sex`, and the `Estimate_survey_sd` /
  `Survey_sd_prior` legacy aliases documented *instead of* `Estimate_index_sd` /
  `Index_sd_prior`), the `weight1_Numbers2` case mismatch, and the
  "must be kept in sync by eye" comment at [0-switches.R:27-30](R/0-switches.R#L27-L30).
- Add the missing `fleet_control` presence guards in `rearrange_data()`.

## PR 6 — `save_config()` / `load_config()`

SAM's `saveConf`/`loadConf` pattern: round-trip the config to documented plain text with each
field's doc string emitted as a comment above it. Makes two assessment runs **git-diffable**
and a configuration archivable alongside an assessment. Reuses PR 5's canonical schema for the
doc strings and PR 4's `model_config` for the content.

## PR 7 — Contributor docs + C++ legibility

The paused accessibility work, unblocked now that PR #98's testing suite and
`.claude/commands/golden-check.md` are merged:

- `vignettes/environmental-linkages.Rmd` — the unified grammar across all six processes, with
  the legacy translation table above as a migration guide.
- `vignettes/extending-rceattle.Rmd` — add an SRR form / a selectivity form / a likelihood
  component. `.github/CONTRIBUTING.md` + pkgdown nav.
- Name the `jnll_comp` magic integer rows via an enum (now 0–20, bare integers across
  `ceattle_v01_11.cpp:2368-3315`, legend only in `R/6-rename_output.R:130-151`).
- Repair the non-monotonic cpp section index (no 5.2, no 6.10, `10.1.1` twice, `12.2.1`
  before `12.1`).
- Doxygen for the three bare headers: `linkage.hpp` (434), `helper_functions.hpp` (215),
  `comp_osa.hpp` (133) — emulating `recruitment.hpp`. Also delete the now-false note at
  `linkage.hpp:23-26` ("currently only exposes the accumulator against the growth parameters").

---

## Verification

Per PR, in order:

1. `pkgload::load_all(".")` with `export PATH=/usr/bin:$PATH`.
2. **`/golden-check`** — `BS2017SS` + a multispecies fit within tolerance. Hard gate for
   PRs 0–3. PR 0 and the PR 1 map-ordering move must be **bit-identical**.
3. `NOT_CRAN=true Rscript -e 'devtools::test()'` green.
4. New tests, flat `test-<area>-<topic>.R` convention, fits behind `skip_on_cran()`:
   - `test-linkage-logit-link.R`, `test-linkage-unsupported-process.R`,
     `test-switches-varying-q-parse.R` — the PR 0 regressions.
   - `test-linkage-formula-grammar.R` — the parser: bars, `rw()`, `ar1()`, `block()`,
     `factor(Year)` vs bare `Year`, and that a bare `Year` gives a *trend* not deviates.
   - `test-linkage-catchability.R` — a `linkage_spec(~temp, by=~fleet)` on `q` reproduces a
     hand-built `Estimate_q = 5` fit bit-for-bit; the deprecation path warns once.
   - `test-linkage-random-effects.R` — `~(1|Year)` estimates sigma and enters `random=`;
     `~factor(Year)` + prior does not; the two give different `jnll_comp` rows.
   - `test-linkage-selectivity.R` — supported forms move as expected; LogisticPM's
     multiplicative path; cases 0 and 9 error in `data_check()` naming the fleet.
   - `test-linkage-tv-translation.R` — every row of the legacy translation table fits
     bit-identically to its pre-change equivalent.
   - `test-likelihood-dirmult-prior.R`, `test-data-required-inputs.R`,
     `test-data-excel-optional-sheets.R`, `test-config-roundtrip.R`.
5. `rcmdcheck::rcmdcheck()` backgrounded before each PR lands.

Per CLAUDE.md, each PR updates **`NEWS.md`**, bumps **`DESCRIPTION` Version:** (minor for
1–4, patch for 0 and 5–7), and touches affected vignettes. After `devtools::document()`,
`git checkout` the spurious roxygen 8.0.0 churn in unrelated `man/` files and `DESCRIPTION`.
Commit messages plain, no AI-attribution trailer; Grant confirms before each commit.
