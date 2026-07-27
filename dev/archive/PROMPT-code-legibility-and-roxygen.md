# Prompt — code legibility (condense/clarify) + roxygen clarity pass

Copy the block below into a fresh session. This is the deferred **"condense the codebase"** pass
plus **PR 7** (contributor docs + C++ legibility) from `dev/PLAN-data-workflow-and-linkage-grammar.md`.
Run it as its own effort, on its own branch — **not** woven into a feature PR.

> **STATUS (essentially COMPLETE, 2026-07-27; branch `pr7-legibility`, DESCRIPTION 5.0.0).**
> Tiers A–C (R-side condense + deletions + the `.ts_wrapper()`/`.coerce_switch_arg()`/`.pull_int`
> factories + two bug fixes), the C++ legibility pass (section-index fixes, Doxygen headers, the
> `JnllRow` enum), **D1** (reframed to the `build_growth(sd_plus_group=)` WHAM/SS3 feature), and
> **E** (the 11-site `.refit_like()` collapse, `R/6-refit_like.R`, verified bit-identical incl. a
> multispecies MSE) are all **shipped, not yet merged**. **D2** (`linkage.hpp` accumulators) was
> **DECLINED** (net-negative). The **roxygen-clarity sweep was skipped** (docs milder than feared).
> **Only remaining: merge `pr7-legibility` → `dev-data-workflow`.** Full detail:
> `dev/HANDOFF-pr7-condense-tiers-DE.md`. Kinzey code preserved per Grant.

---

We're doing a **legibility + documentation pass** on Rceattle (the TMB/C++ + R climate-linked
multispecies stock-assessment model). Data-workflow **PRs 0–6 are shipped** (v4.13.0) on
`dev-data-workflow`; the tree is clean and golden-green at tip `e2e3a110`. **Branch PR 7 off that
tip** (its own branch — do NOT weave it into `dev-data-workflow`). Two goals, no new features:

1. **Condense/clarify the code** — collapse duplication and delete dead code so the package gets
   *smaller and clearer*, with zero behavior change.
2. **Improve the roxygen docs** — make them clearer and more accurate, and **less AI-verbose**.

**Governing standard — this software sets US federal fishing quotas** (Magnuson-Stevens management
advice), so a silently-wrong number can misset a catch limit. A refactor is only correct if it is
**behavior-preserving**: keep golden-reference equivalence (`/golden-check`, all 4 models — BS2017
SS/MS + GOA2018 SS/MS — **bit-identical**) for anything that can move a fit, and **adversarially
verify** each collapse against the source before asserting it is safe (the plan explicitly warns
that "collapsing 11 call sites that differ in six invisible ways is its own failure mode" — a
behavior change disguised as a refactor). Automated tests passing is necessary, not sufficient.

## Where the work is already scoped

**Read `dev/PLAN-data-workflow-and-linkage-grammar.md`, section "Deferred: condense the codebase"**
first — it is a pre-scouted target list from a review sweep, with approximate line savings. Do not
re-derive it from scratch; **confirm each target against the current tree** (sizes/line numbers
have drifted since the sweep — the data-workflow + rename PRs 4–6 landed), then execute the safe
ones. Highlights it already names:

- **Collapse candidates:** the `fit_mod()` re-invocation block copy-pasted **11×** across
  `9-retro_and_jitter.R` / `10-run_mse.R` / `10-project-no-F.R` (→ a `.refit_like()` helper — but
  see the prerequisite below); 26 near-identical `pull()` blocks in `5-rearrange_data.R`; the three
  `.coerce_*` + four `.warn_*_deprecation()` skeletons in `0-build_srr_and_M.R`; the six
  `plot_*` timeseries wrappers in `7-plot_ceattle.R` (→ a `.ts_wrapper()` factory); the six C++
  accumulators in `linkage.hpp` (~90% shared → one templated helper) and the near-duplicate
  `estimate_growth()`/`estimate_growth_within_yr()` in `growth.hpp`.
- **Straight deletion (no replacement):** `11-model_average.R` (181/188 lines commented out);
  ~350 lines of commented Kinzey code in `predation.hpp`; the dead triplicate blocks in
  `2-build_params.R` / `4-build_parameter_bounds.R` / `OPT-phaser.R`; six no-op self-assignments in
  the projection HCR switch; the triplicated 27-name MSE "quantities to keep" whitelist (it lists
  `ssb_depletion` twice — the copy-paste tell); `t_col()` in `7-plot_stock_recruit.R` (defined,
  documented, called nowhere); stale empty `tests/testthat/tests-*/` dirs; and the `R/dev/` fork
  files (Rbuildignored but a navigation hazard).

**Prerequisite for the `.refit_like()` collapse — do not skip (the plan is emphatic):** the 11
sites are NOT homogeneous — each reads ~40 args from the *working* model with a few deliberate
exceptions from the *pristine* original (the OM pins `suit_styr`/`suit_endyr`/`srr_*` to the
original; one arg is *computed*, `srr_mse_switchyr = em_use$data_list$endyr`, not forwarded).
Tabulate all 11 sites × every arg × source object, classify each exception as intentional-pin vs
accidental-drift, then write the helper with the pins as **named overrides** visible at the call
site. Golden-check each site before and after. Don't trust the earlier "the SRR args drifted"
claim — they didn't; two different models each read their own family. **Note:** a v4.13.0 bug fix
(`74f029d6`) added a `getsd` argument to `self_test()` so all four refit dispatchers
(`retrospective`/`jitter`/`self_test`/`profile`) now handle `getsd` identically — fold that into the
tabulation; the `getsd`-inheriting default (`!is.null(<model>$sdrep)`) is a shared idiom worth
hoisting.

## Roxygen clarity — the bar

Rewrite docs to be **precise and terse**, matching the codebase's existing voice. Anchors of the
good voice already in-tree: `src/TMB/recruitment.hpp` (the Doxygen exemplar CLAUDE.md points to),
`R/0-fit_control.R`, `R/0-column_schema.R`. **Avoid AI verbosity:**

- Cut filler openers ("This function is designed to…", "It is important to note that…", "In order
  to…", "Simply…"). Lead with the verb: "@description Fit the CEATTLE model…".
- Don't restate the signature in prose or pad every trivial `@param`. DO document what a reader
  can't infer: **units, allowed values / switch codes, what the default *means*, and side effects**.
- Cut hedging ("should", "may", "generally") where the behavior is definite; state it.
- `@details` only for non-obvious mechanics (an equation, an ordering constraint, a gotcha) — not
  a restatement of `@description`.
- Use the domain terms exactly (per CLAUDE.md): SSB = female spawning-stock biomass; F40%/F35%/B40%
  SPR proxies not "MSY"; name the selectivity form; Francis/McAllister-Ianelli/Dirichlet-multinomial
  weighting; Mohn's rho; OSA residuals. A stock assessor should read a doc string and be right.
- **Keep the "why" comments** the codebase favors (the explanatory section headers, the Doxygen on
  the C++). Condense the "what"; never strip the rationale — do not sacrifice a real explanation to
  brevity.

## Traps (from CLAUDE.md — read it in full)

- **Toolchain:** prefix R compile/check with `export PATH=/usr/bin:$PATH` (system toolchain first).
  `pkgload::load_all(".")` recompiles the TMB DLL after any `.cpp`/`.hpp` edit; add `compile = FALSE`
  for R-only changes. Never commit the `*.o`/`*.so` artifacts.
- **Numbered file prefixes are meaningful — don't renumber/rename wholesale.**
- **Released package — preserve the public API.** Keep deprecation aliases; sweep `../Rceattle-models`
  (+ the sibling assessment repos) after any user-visible change. Ignore the vendored
  `EBS_CEATTLE_TMB` fork.
- **roxygen churn trap (current form):** the repo is now on roxygen2 8.0.0
  (`DESCRIPTION: Config/roxygen2/version: 8.0.0` — the 7.3.3→8.0.0 swap already happened; CLAUDE.md's
  note about `RoxygenNote: 7.3.3` is itself stale). The live nuisance is that `document()` re-touches
  ~10 **stale `man/*.Rd` source-path lines** ("% Please edit documentation in R/OPT-TMBAIC.R" →
  `R/6-tmb_aic.R`, plus the osa/process `.Rd`) from R files renamed without re-documenting. After
  every `document()`, `git checkout` those unrelated `.Rd` and keep only the ones you meant to change
  — OR, since this pass touches docs wholesale, consider a deliberate one-off "resync stale man/
  source paths" commit up front so the noise stops. Give internal helpers `@noRd` (not just
  `@keywords internal`). Never insert a helper between a function's roxygen block and its definition
  (the block binds to whatever follows, silently stealing `@export`).
- `jnll_comp` likelihood-row names are magic integers in `ceattle_v01_11.cpp`, mirrored by hand in
  `R/6-rename_output.R` — keep in sync if you touch either.

## How to run it

- **Separate branch, incremental golden-gated sub-commits**, each: `load_all` → targeted tests →
  `/golden-check` (bit-identical for anything touching a fit; pure deletions of dead code and
  doc-only edits are golden-inert but still run the fast suite) → **spawn an adversarial reviewer**
  (an Explore agent) to try to break the "behavior-preserving" claim → **Grant confirms before each
  commit**. Plain commit messages, no AI-attribution trailer.
- **Deliver a prioritized findings list first** (target → what/why → est. savings → risk →
  golden-relevant y/n), then execute lowest-risk-highest-value first (straight deletions and the
  plot/`pull()`/`.coerce_*` factories before the `.refit_like()` tabulation).
- Bump `DESCRIPTION` patch per change batch; add NEWS bullets only for anything user-visible
  (most of this is internal and gets none). `rcmdcheck` before landing.

Note: PRs 0–6 are all shipped on `dev-data-workflow` (tip `e2e3a110`, golden PASS) — branch PR 7
from there. First action in the session: `export PATH=/usr/bin:$PATH && NOT_CRAN=true` fast suite +
`/golden-check` to confirm the green baseline before touching anything. Then deliver the prioritized
findings list (below) and wait for Grant to pick targets before editing.

Two adjacent deferred items have their own briefs — do NOT fold them into this pass, but know they
exist: `dev/PROMPT-mse-runtime-optimization.md` (decrease MSE run time — a perf pass, overlaps the
`10-run_mse.R` files this pass reads) and `dev/PROMPT-osa-with-composition-accumulation.md`
(OSA-residuals-with-accumulation feature follow-up).

---
