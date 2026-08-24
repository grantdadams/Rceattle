# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**Four stacked PRs are pushed and ready to open**, each a strict prefix of the next, all four on
`dev`. Open and merge them in order; each is independently reviewable and its own CI runs.

| Branch | Commits | Contents | Can it move a number? |
|---|---|---|---|
| `pr1-api-tooling-and-guards` | 8 | Stage 6 `object=` API, Stage 0 Claude tooling + CLAUDE.md, Stage 1 drift guards | **No** — and it carries the guards the rest rely on, so merge it first |
| `pr2-schema-authority` | +5 | Stage 2 docs, Stage 3 schema authority, Stage 4/5, both adversarial-review rounds | **Yes — this is the one that can refuse a workbook** |
| `pr3-schema-order-and-qar1` | +2 | fleet_control column order, `flt_sel_ind` removal, QAR1 deprecation | changes workbook layout only |
| `pr4-analytical-catch-sigma` | +1 | `Estimate_catch_sd = "Analytical"` implemented in the template | new capability; golden bit-identical |

`api-object-argument` is the old integration branch, at the same commit as `pr3`. It can be
deleted once the four are open.

**The single thing to read closely in pr2** is the three-column hard error
(`Selectivity_dimension`, `Sel_shape_dir`, `Sel_shape_mode`). It is the only change in the whole
series that can refuse a workbook that previously loaded. It was pre-flighted over every workbook
in the ecosystem -- 196 with a `fleet_control` sheet, **zero** rejected -- and each column is
validated against exactly the spellings its consumer implements, so no working input is refused.
It could not be split out of pr2 without hunk surgery on two schema files.

## Done & verified

- **Final state**: suite **7248 pass / 0 fail / 3 skip**; golden `identical()` to the `dev`
  baseline; `verify-refit-like.R` bit-identical across all 32 fits; all seven new `@examples`
  execute cleanly through `tools::Rd2ex`.
- **Stage 6** (at the time): suite 6925 pass / 0 fail. Golden 4-model capture `identical()` TRUE, every field
  `0.000e+00`. `verify-refit-like.R` bit-identical across all 32 fits. Ecosystem sweep found 15
  old-name call sites, all shimmed, and zero calls passing both spellings.
- **The golden reference on `dev` is `ss = 10241.0304272585`** (`ms = 10267.2478324443`,
  `goa_ss = 12868.0052289274`, `goa_ms = 12932.7931701136`), measured against a `dev` worktree.
  A different value had been recorded in the local agent file and was wrong.
- **The Stage 3 hard error was pre-flighted, not assumed.** A report-only pass over every
  workbook in the ecosystem found 196 carrying a `fleet_control` sheet and **not one** value the
  new check rejects. That is why it ships as an error rather than a warning, and why this is
  5.12.0 rather than 6.0.0.
- Every drift guard was mutation-tested: renumbering a map value, adding an undispatched code,
  dropping a documented one, and removing the deprecated SRR codes each turn the suite red.

## For Grant's review

1. **`write_data()` now writes `fleet_control` in schema order**, matching what the control and
   bioenergetics sheets already did. Values round-trip identically on all three bundled datasets.
   But **`Month` moves from column 5 to column 48**, because the schema declares it late. That is
   a worse layout for a human reading the workbook. Options: reorder the schema rows so the
   workbook reads sensibly, or drop the fleet_control reordering. **Your call.**
2. **`.PAR_INFO` omits 7 declared parameters** -- `index_q_pow` and the six Kinzey-Punt
   `logH_*`/`H_4`. All belong to stubbed features, so they are pinned as exempt rather than
   documented. If any of those features is being picked up, the exemption should go.
3. **`flt_sel_ind` is dead.** `rearrange_data()` computes it from `Fleet_code` on every fit and
   nothing -- R or C++ -- reads it. Left alone because removing it is a behaviour change, not a
   cleanup.
4. **`inst/dev/CLEANUP_BACKLOG.md` Tier 0 lists 9 known defects** taken from FIXME text in the
   source, including that **QAR1 is inert** -- `Catchability = "AR1"` (the Rogers et al. 2024
   form) can be set and does nothing, because the deviate map is gated on `Time_varying_q`,
   which under QAR1 holds an `env_data` column index rather than a mode
   (`R/3-build_map.R:1181`). Note this is **not** `Time_varying_q = "AR1"`, which is a
   different switch sharing the same string and works correctly. Also that `run_mse()`'s catch
   fill-in does not work for `assessment_period > 1`. Each needs an issue if it is real.

## Known flags

- **5.10.0 moved three predation figures' numbers** -- `plot_m2_at_age_prop()` (a share now, not
  a contribution), `plot_ration()` (x average numbers-at-age) and `plot_b_eaten()` (million mt).
  Any figure regenerated from them differs from earlier runs. `plot_selectivity()` also renames
  `p$data$Age` to `Bin`.
- **Result-changing changes on this line that are not labelled breaking**: the mode-5 selectivity
  penalty fix (GOA Pacific cod SSB 2050 -14.1%), parameter bounds previously applied to the wrong
  parameters, composition weights warm-starting from `inits`, failed `run_mse()` simulations
  returning only a marker, the `mse_summary()` reshape, the recruitment fixes, and `sim_mod()`
  drawing the index under the fleet's own `Index_distribution`. **A model carrying GOA numbers
  forward needs a refit.**
- The `nages` age-vs-index defect in the three predation plotters **is fixed** (`d4cdfb2f`), and
  the class is now closed by `minage != 1` fixtures. CLAUDE.md's "still wrong" note was stale and
  has been corrected. The same class survives in the template at `src/TMB/ceattle.cpp:1943`.

## Blocked

Nothing.

## Resume here

Open the four PRs above, in order. Then work `inst/dev/CLEANUP_BACKLOG.md` using
`inst/dev/BACKLOG-PLAN.md`, which sequences it by who is exposed rather than by tier label and
says what covers each item -- the short version being that `/golden-check` will be green through
almost all of it, because none of the four reference models reaches these inputs.

Queued and not started: `inst/dev/CONTRIBUTOR-EXPERIENCE.md` — making the codebase navigable to
a fisheries scientist, drawn from a review of FIMS. **Item 0 comes first and is a conversation,
not code**: ask the sibling-repo authors where they actually stopped, because the ordering of
A-H rests on an assumption nobody has tested. A-E and H are doc and tooling; G is additive; F is
the only item that can break a caller and needs `/golden-check` + `/ecosystem-sweep`. Each recipe
ships with a drift guard or it is not done. It carries its own start prompt.

## Older paused work

**The `accessibility-and-code-review` refactor is superseded, not outstanding.** The branch is
gone from local and remote, its plan
(`~/Downloads/HANDOFF-accessibility-refactor-implementation.md`) no longer exists, and no commit
in any branch mentions it. Verified 2026-08-22 against `dev`: three of its four locked "chosen
extras" have landed by other routes -- the `JnllRow` enum (`ceattle.cpp:2963`, 148 uses), the
repaired cpp section index (47 numbered sections), and real Doxygen `@file`/`@brief` blocks on
the previously bare headers.

**The one concrete leftover is splitting `R/0-build_srr_and_M.R`** -- still 1,497 lines and 29
functions in one file, and still the grab-bag the plan named. It is in
`inst/dev/CLEANUP_BACKLOG.md`.

Anything else that plan contained is unrecoverable from the repo. Do not resume from the branch
name or the handoff path; both are dead references.
