# dev/ — working notes (gitignored)

Scratch plans, handoffs, prompts, and verification scripts. **Not shipped** (`dev/` is
gitignored). Durable knowledge — how the subsystems work, the traps — lives in `CLAUDE.md`
and `vignettes/articles/developer-guide.Rmd`; these files are process/status only.

## Active

| File | What it is |
|---|---|
| `PLAN-data-workflow-and-linkage-grammar.md` | Master roadmap for the data-workflow + linkage-grammar effort (PRs 0–7) and the prioritized backlog. Start here. |
| `HANDOFF-pr7-condense-tiers-DE.md` | PR 7 (legibility/condense) status. Tiers A–C, D1 feature, C++ legibility, `jnll_comp` enum, and the `.refit_like` collapse all done; **D2 declined**. Only open item: **merge `pr7-legibility` → `dev-data-workflow`**. |
| `PROMPT-code-legibility-and-roxygen.md` | Original PR 7 brief (status banner points at the handoff above). |
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
PR 4, PR 5 (excel-schema + phase-3), the linkage random-effect handoff, the QAR1
productionization status, the resolved GOA-MS M-estimation prompt, and the superseded
next-session prompt. Their durable outcomes are already in `CLAUDE.md` / the developer
guide / auto-memory.
