# dev/ — working notes (gitignored)

Scratch plans, handoffs, prompts, and verification scripts. **Not shipped** (`dev/` is
gitignored). Durable knowledge — how the subsystems work, the traps — lives in `CLAUDE.md`
and `vignettes/articles/developer-guide.Rmd`; these files are process/status only.

## Active

The data-workflow + linkage-grammar effort (PRs 0–7) is **complete and merged onto
`dev-data-workflow`** (not yet released to `main`). What remains is forward backlog:

| File | What it is |
|---|---|
| `PLAN-data-workflow-and-linkage-grammar.md` | Master roadmap + historical record of the whole effort. Top banner has the current status; the backlog is the prompts below. |
| `PROMPT-doc-quality-review-pr1-7.md` | Documentation-only review of the whole PR 1–7 doc surface (roxygen, vignettes, NEWS, field dictionary, console messages) against three criteria: not AI-verbose, current-capabilities not change-log, scientist-legible. Recommended before the 5.0.0 release. |
| `PLAN-sel-curve-pen-reformulation.md` | Future work: rename the form-dependent `Sel_curve_pen1/2/3` columns to purpose-named SD columns. Not started (`Sel_curve_pen` still in code). |
| `PROMPT-mse-runtime-optimization.md` | Future work: speed up the MSE loop. Not started; its own branch. |
| `PROMPT-osa-with-composition-accumulation.md` | Future work: OSA residuals when composition tail-accumulation is active (deferred from the `Comp_accum_*` feature). |

## Verification scripts

| File | Purpose |
|---|---|
| `golden-ref.rds` | Pinned 4-model golden reference (see the `/golden-check` skill). |
| `verify-refit-like.R` | Before/after equivalence harness for the `.refit_like` collapse (runs every refit entry function incl. SS + MS MSE; compares bit-identical). |
| `verify-mse-hindcast-invariant.R` | Checks the MSE does not perturb the hindcast (OM SSB over styr:endyr fixed under simulate_data/sample_rec, SS + MS). |

## archive/

Completed / superseded / resolved docs, kept for history only. Safe to delete.
The PR 7 handoff + legibility prompt (now merged), PR 4, PR 5 (excel-schema + phase-3),
the linkage random-effect handoff, the QAR1 productionization status, the resolved GOA-MS
M-estimation prompt, and the superseded next-session prompt. Their durable outcomes are
already in `CLAUDE.md` / the developer guide / auto-memory.
