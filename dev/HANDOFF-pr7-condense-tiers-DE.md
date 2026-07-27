# HANDOFF — PR 7 condense pass, deferred Tiers D + E (+ roxygen sweep)

> **STATUS UPDATE (session 2, branch `pr7-legibility`, ~19 commits, DESCRIPTION 4.14.0).**
> Since the Tier A–C header below, a second session took this much further. **DONE:**
> - **D1 (growth.hpp)** — reframed by Grant from a line-extraction into a **feature**: per-species
>   `build_growth(sd_plus_group=)` WHAM/SS3 plus-group SD-at-age switch (`467fc4a8`), plus an
>   NA-inherit fix so refits keep the choice (`4c8de80e`). Golden bit-identical (default WHAM),
>   constructive + refit tests, adversarial SAFE. The originally-planned "extract the ~100-line
>   byte-identical growth tail" was **not** done — it's a 21-arg helper, net-neutral, so skipped.
> - **E (`.refit_like`) — DONE (2026-07-27).** The ~60-line `fit_mod()` refit block copy-pasted at
>   **11 sites** (retrospective ×2, jitter, self_test, profile, run_mse ×5, remove_F) is now one
>   internal helper `R/6-refit_like.R`. It reconstructs the HCR / recFun / M1Fun / growthFun + shared
>   switches from a source `data_list`; every per-site divergence is a named override defaulting to the
>   `data_list` value, so each site shows only what it changes (e.g. site-9's pins to pristine `om$`
>   are now 5 explicit, commented lines). Console suppression left at each call site (differs, and is
>   console-only). **Verified golden-unverifiably** via a bespoke before/after harness
>   (`dev/verify-refit-like.R`) that runs all six entry functions on a tiny seeded model incl. a
>   **multispecies MSE** (exercises the suit-window pin) — **bit-identical across all 9 sections**
>   (obj+ssb+R). Plus `dev/verify-mse-hindcast-invariant.R` confirms the MSE does not perturb the
>   hindcast (OM SSB over styr:endyr fixed under simulate_data F/F and T/T, SS and MS). The six
>   originally-flagged divergences are preserved as the named overrides (none was a bug).
> - **C++ legibility** — cpp section-index dup/misnumber fixes (`85a31543`); Doxygen on all three
>   bare headers linkage/helper_functions/comp_osa (`622a8d71`). The false `linkage.hpp` note was
>   already gone.
> - **Multi-session audit follow-ups** — `plot_logindex` restored as a deprecated shim (`9f637e01`);
>   `mse_summary()` + `sim_mod()` regression tests added (previously zero coverage).
> - **Contributor docs** — already complete (`environmental-linkages-and-priors.Rmd`,
>   `articles/developer-guide.Rmd`); the **roxygen clarity sweep was SKIPPED** (docs milder than
>   feared, low value).
>
> **`jnll_comp` enum — DONE (`56733c6f`):** the 21 magic likelihood-row integers across 112
> `jnll_comp(N,…)` / `unweighted_jnll_comp(N,…)` sites are now named constants (`JNLL_INDEX`..
> `JNLL_LINKAGE_RE`) via a `JnllRow` enum by the declaration, cross-referenced to
> `R/6-rename_output.R`; each enum value equals the integer it replaced, so golden bit-identical
> (incl. per-row `jnll_dev`). Also relabeled `14.3`→`13.2` (Diet likelihood is a section-13 term).
>
> **Also DONE this session:** the `mse_summary()` reshape to a per-entity `list(species, fleet,
> total, meta)` (breaking → DESCRIPTION **5.0.0**); the developer-guide expansion for the PRs 1–7
> data-workflow surface; and a **CLAUDE.md refresh** (column schema / config / linkage pointers +
> `JnllRow` bullet + Active context).
>
> **STILL OPEN:** (1) **Merge** `pr7-legibility` → `dev-data-workflow` — the only remaining item.
> (2) Minor: the orphaned `12.2.1 "Calculate CAAL"` cpp section label (needs a structural call on
> where section 12 begins) was left as-is.
>
> **CLOSED — Tier D2 (`linkage.hpp` accumulators): DECLINED (wontfix), 2026-07-26.** Grant's
> call after review. The five accumulators span **two stratum universes** (population:
> species×sex×age for growth/M/rec; fleet: fleet[×sex] for q/sel) with four tensor ranks and
> `sel`'s three-way param branch — they do **not** factor onto one writer-parameterized helper
> without functors over both universes. The genuinely-shared slice is only the ~9-line
> row-matching harness (~5–15 net lines). Against that: AD-taped code that sets federal quotas,
> already Doxygen'd and individually simple. Verdict: a merge is **net-negative for legibility**.
> (Coverage worry from the original plan is moot: the suite already exercises all five compiled
> accumulators — `test-dynamics-recruitment-linkage.R` *constructively pins the offset tensor*,
> q/sel via real fits — so a future edit is well-netted regardless.)

**Context.** PR 7 (the "condense the codebase" + roxygen pass, brief:
`dev/PROMPT-code-legibility-and-roxygen.md`) was run on branch **`pr7-legibility`** off
`e2e3a110`. **Tiers A–C are done and committed** (9 commits, every fit-touching one golden
bit-identical across all 4 reference models; DESCRIPTION at 4.13.1). The branch is **not merged**
to `dev-data-workflow`. This handoff covers what remains: **Tiers D and E**, and the **roxygen
clarity sweep**, which the first pass barely touched.

Read `dev/PLAN-data-workflow-and-linkage-grammar.md` → "Deferred: condense the codebase" (its
STATUS block records which sweep claims were over-scoped) and CLAUDE.md before starting.

## What shipped in Tiers A–C (so you don't redo it)

`80da7dc5` dead-code deletion + `R/dev/`→`dev/` · `bfc64228` plot_ssb bug fix ·
`e3c6f095` B1 HCR self-assigns (C++) · `b49d4ba3` B2 MSE whitelist constant ·
`5d342b74` B3 man/ source-path resync · `b0a0a7c4` C1 `.ts_wrapper()` factory ·
`64192361` C2 `.coerce_switch_arg()` · `05343873` C3 CAAL message bug · `b22fa021` C3 `.pull_int`.

**Kinzey code (predation.hpp block + R `logH` commented blocks) is PRESERVED per Grant — do not
delete it.**

## How to verify (mandatory for D + E — both touch or re-run fits)

Golden reference is pinned at `dev/golden-ref.rds`; recipe in the `/golden-check` skill. Fit the
4 models (BS + GOA, SS + MS; MS warm-started from SS MLEs) and diff — **must be bit-identical**
(0.000e+00) for any behavior-preserving change. Reference objectives:
ss=10241.0304275402, ms=10267.2478327352, goa_ss=12807.4375258732, goa_ms=12866.2957391829.
Toolchain: `export PATH=/usr/bin:$PATH` before any compile/check (rcmdcheck's build subprocess
drops it → TMB compile fails; use `R CMD build` + `R CMD check <tarball>` instead).

---

## Tier D — C++ numeric-path collapses (highest risk; each individually golden-gated)

**D1 — `src/TMB/growth.hpp` shared tail.** `estimate_growth()` (L59–270) and
`estimate_growth_within_yr()` (L321–505) are ~200-line near-duplicates. The **divergent** part is
the length-projection core (VBGF/Richards `switch`; est_growth uses year-lag cohort recursion +
has a plus-group-correction block within_yr lacks; within_yr uses `fracyr` within-year advance).
The **byte-identical tail** (~100 lines) is the SD-interpolation section + the "Build Growth
Matrix & Weight-at-Age" loop — differing only in one self-referential comment. Extract that tail
into one helper both call. Also drop the **unreachable `growth_model == 3` branches** (L226 and
L461; both already commented, nested inside `growth_model < 3`, and the switch only yields 1/2).
Risk: AD taping — golden bit-identical is mandatory. Est. ~100 lines.

**D2 — `src/TMB/linkage.hpp` accumulators: DECLINED (wontfix, 2026-07-26).** See the STATUS
block at the top for the rationale. The analysis below stands as the record of *why* it doesn't
factor; it is no longer a backlog item.

There are **FIVE**
`rceattle_apply_{growth,M,recruitment,q,sel}_linkages` (not six), sharing ~90% skeleton
(`beta.size()` guard, `yr_hi=min(...)`, the `linkage_process`/`link_code` row filter,
`rceattle_stratum_range()` expansion, `+= b * linkage_X(yr, xc)`). They differ **only** in the
offset-tensor writer + its stratification loops. **BUT `sel` diverges hard** — it writes 3 tensors
(`slp/inf/coff_offset`) with a `param`-code branch — so a single templated helper must special-case
it. **There is NO `linkfn`/logit branch anywhere** (the PR-0 logit item is not in the tree; `link`
here is only a row filter). Biggest concern: **the 4 golden models may not exercise all 5
accumulators** — write dedicated linkage-path regression tests (fit a small model with a
q/sel/M/rec/growth linkage each) BEFORE collapsing, or golden bit-identity proves nothing. Est.
~160 lines. Recommend: do D1 first; attempt D2 only with linkage tests in hand, else leave it.

---

## Tier E — the `.refit_like()` capstone (~900 → ~120 lines; tabulation below is CONFIRMED)

The `fit_mod()` re-invocation block is copy-pasted **11×** across `9-retro_and_jitter.R` (5),
`10-run_mse.R` (5), `10-project-no-F.R` (1). Each calls `build_hcr` + `build_srr` (in
`suppressWarnings`) + `build_M1` + `build_growth` + ~8 pass-throughs, ending in `getsd`/`verbose=0`.
**They are NOT homogeneous** — the whole point of the prerequisite tabulation is that the
divergences are deliberate and must survive as *named, visible overrides at the call site*, not be
silently unified. Write `.refit_like(source_object, ..., srr_mse_switchyr = <override>,
getsd = <override>)` and golden-check **each site** before/after.

| # | File / dispatcher | Lines | `srr_mse_switchyr` | `getsd` | Source obj |
|---|---|---|---|---|---|
| 1–2 | `retrospective`/run_one_peel A,B | 242–300, 350–408 | `min(dl$srr_mse_switchyr, endyr_peel)` | `getsd` var | local `data_list` |
| 3 | `jitter`/run_one_jitter | 611–669 | `min(dl$…, dl$endyr)` | `getsd` var | local |
| 4 | `self_test`/run_one_sim | 785–843 | `min(dl$…, dl$endyr)` | `getsd` var (L842 = v4.13.0 fix `74f029d6`) | local |
| 5 | `profile`/run_one_point | 1197–1254 | `min(dl$…, dl$endyr)` | `getsd` var | local |
| 6 | `om <- fit_mod` | 83–140 | `om$…srr_mse_switchyr` (no `min`) | `FALSE` | pristine `om$` |
| 7 | `em <- fit_mod` | 228–285 | `em$…` (no `min`) | `FALSE` | `em$` |
| 8 | `em` input-F refit | 301–346 | `em$…` | `FALSE` | `em$`, HCR hardcoded `2` |
| 9 | `om_use <- fit_mod` (sim loop) | 635–692 | **`om$…` (pristine)** while HCR/growth read `om_use$` | `FALSE` | mixed — **intentional pin** |
| 10 | `em_use <- fit_mod` (sim loop) | 843–900 | **`em_use$data_list$endyr`** (computed, not `srr_mse_switchyr`) | `FALSE` | `em_use$` |
| 11 | `remove_F` | 26–83 | `min(Rceattle$…, Rceattle$…endyr)` | `TRUE` | `Rceattle$` |

**Do NOT "fix" sites 9 and 10.** Site 9 pins the OM's stock-recruit reference period to the
*pristine* `om$` while everything else reads the working `om_use$` — deliberate (the OM's biology
stays fixed through the projection). Site 10's `srr_mse_switchyr = em_use$data_list$endyr` is a
*computed* override that advances every iteration — a pure field-forwarder can't express it, so
`.refit_like()` must accept computed overrides. (Line numbers are pre-Tier-E; re-confirm against
the tree — this file was not touched by A–C except no MSE-whitelist lines near 929/605.) The
`getsd`-inheriting default `!is.null(<model>$sdrep)` is a shared idiom worth hoisting.
**Note the overlap** with `dev/PROMPT-mse-runtime-optimization.md` (same `10-run_mse.R` files) —
coordinate so the two passes don't collide.

---

## Roxygen clarity sweep (largely UNSTARTED)

Tiers A–C only touched docs opportunistically (the new helpers' roxygen + the B3 path resync).
The brief's "make docs precise/terse, less AI-verbose" sweep of the heavy-doc files
(`0-build_srr_and_M.R` ~480 roxygen lines, `0-build_linkage.R` ~358, `7-plot_ceattle.R` ~323,
`6-fit_mod.R`) is **not done**. **Reality check from the first pass: the roxygen is much milder
than the brief feared** — ~0 classic filler openers ("This function is designed to…") across 4,463
roxygen lines. So this is low-value polish, not a rewrite: prioritize documenting what a reader
can't infer (units, switch codes, what a default *means*, side effects) on the heavy files, and
don't manufacture churn. Voice anchors: `R/0-fit_control.R`, `R/0-column_schema.R`,
`src/TMB/recruitment.hpp`. Give internal helpers `@noRd`; never wedge a helper between a roxygen
block and its function.

## Contributor docs + C++ legibility (the "PR 7" section of the PLAN — UNSTARTED)

These are the *other half* of PR 7 (PLAN "PR 7 — Contributor docs + C++ legibility"), untouched by
the A–C condense pass. All doc/comment-level (golden-inert), but re-confirm line numbers — they
drifted through PRs 4–6 + PR 7 A–C:
- `vignettes/environmental-linkages.Rmd` — the unified grammar across all six processes, with the
  legacy→grammar translation table as a migration guide. `vignettes/extending-rceattle.Rmd` — how to
  add an SRR form / selectivity form / likelihood component; `.github/CONTRIBUTING.md` + pkgdown nav.
- **Name the `jnll_comp` magic-integer rows via an enum** (now 0–20, bare integers in
  `ceattle_v01_11.cpp`; legend only in `R/6-rename_output.R:~130-151`). If you touch the row
  order, update BOTH in lockstep (CLAUDE.md trap).
- **Repair the non-monotonic cpp section index** (missing 5.2 / 6.10, `10.1.1` twice, `12.2.1`
  before `12.1`).
- **Doxygen the three bare headers** — `linkage.hpp`, `helper_functions.hpp`, `comp_osa.hpp` —
  emulating `recruitment.hpp`. Delete the now-false `linkage.hpp` note ("currently only exposes the
  accumulator against the growth parameters") — five accumulators now exist.

## Also noted, out of scope for this pass

- **Two pre-existing suite failures — FIXED (`c0995d12`), and they were NOT order-dependent.** They
  were deterministic stale-test assertions (both reproduce in isolation): `test-dynamics-initial.R`
  asserted `initMode = 5` errors, but 5 is now valid (FishedEquilibrium); and the index-covariance
  retro-peel sim changed `endyr` without clamping `suit_endyr` the way `retrospective()` does. Both
  fixed in the tests (golden-inert). The full suite is now 0 failures. (An earlier read of a
  "shifting pair" of `selectivity-nonparametric` failures was a mistake — those were testthat
  WARNING blocks for the `Time_varying_q` soft-deprecation, not failures.)
- **Per-species `srr_fun`** (Grant's observation). `M1_model` and growth `fun` are per-species
  vectors; `srr_fun` is model-wide (length-1). Making it per-species is a *feature* (touches
  `build_srr()`, `rearrange_data()`, the `srr_fun` DATA scalar in the cpp), not a legibility item.
