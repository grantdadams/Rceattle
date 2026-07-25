# Next session — starting prompt

> **CURRENT STATE (2026-07):** the data-workflow PRs have progressed well past this prompt.
> **PR 5 is DONE (v4.12.0)** — the canonical column schema + the full intuitive-naming rename
> pass + read robustness + `write_template()`. **PR 6 (`save_config`/`load_config`) is IN
> PROGRESS** — sub-commit A (YAML serialization + `save_config`/`load_config`/`config()`) is
> built + reviewed; remaining is the additive `fit_mod(config=)` overlay + finalize. A parallel
> session added **composition tail-accumulation** (`Comp_accum_young`/`Comp_accum_old` + the C++
> fold), being finished to schema standard + committed standalone. The authoritative status +
> next steps live in `dev/PLAN-data-workflow-and-linkage-grammar.md` (the "Status" section).
> The linkage random-effect-density brief below is retained for the record; check it against the
> current tree before assuming it is still the open thread.

Copy the block below into the next session.

---

We're continuing the Rceattle formula-linkage grammar work on branch **`dev-data-workflow`**
(shipped v4.9.0: catchability + selectivity environmental linkages, all parametric forms).
The branch is green (full suite 4059/0) and every commit is golden-verified.

**Governing standard — this software sets US federal fishing quotas.** Rceattle output drives
management advice under Magnuson-Stevens, so a silently-wrong number can misset a catch limit.
Hold every change to that bar: verify *statistical accuracy*, not just that a fit converges;
keep golden-reference equivalence for anything that can move an existing fit; and for anything
without a golden reference (the RE density), verify constructively — simulation
self-consistency (recover a known sigma), limiting cases (sigma→0 = the fixed fit),
cross-checks against the M1_re analogue and the real reference models. **Note and test edge
cases explicitly** (gappy/mis-ordered years, single- vs two-sex, series ending in different
years, mirrored/shared params, mapped-out "Off" structures, non-convergence); an unsupported
path must error loudly naming why, never silently do nothing. Automated tests passing is
necessary but not sufficient — reason through the math, and adversarially verify correctness
claims against the source before asserting them (several review "findings" this project did not
survive that check).

**Primary task — the random-effect linkage density.** This is the deferred, most delicate
piece: Laplace-approximation AD taping with NO golden reference to catch a subtle bug, so it
was banked for a fresh session. **Before writing any C++, read
`dev/HANDOFF-linkage-random-effects.md` in full** — it has the step order, the M1_re verbatim
template (`ceattle_v01_11.cpp:~3388`), the constructive-verification recipe, the
`load_all`-vs-installed parallel-worker gotcha (`R CMD INSTALL .` before local parallel
tests), and — importantly — **four real validation models in `../Rceattle-models`** (EBS/GOA
pollock) that use the existing `Time_varying_sel/q` + `random_sel/q` + `Estimate_q=6`
mechanisms. Reproduce each bit-identically through the new grammar as you go.

Build it in the same incremental, verified sub-commits used for q and selectivity:
encode the RE group registry → split the accumulators (read `beta_linkage_re` vs
`beta_linkage` by `re_index`) → the row-20 IID density (`dnorm(beta_re, 0, exp(log_sigma))`) →
populate the vectors → remove the fit-time guard → **route sigma through
`linkage_spec(init=, priors=)`** (Grant's explicit ask: input a variance value AND put a prior
on it — one mechanism) → then rw()/ar1() → then the legacy translators
(`Time_varying_*`/`M1_re` → grammar, bit-identical) and the `Time_varying_*_sd_prior` → `_sd`
rename. Row 20 must stay exactly 0 when no RE is present (bit-identity still holds).

**Second task, when the density lands — Dirichlet-multinomial priors (plan PR 3).** There is
currently NO prior of any kind on `comp_weights` / `caal_weights` / `diet_comp_weights`. Add a
`"comp"` process to the same linkage grammar (params `theta_comp`/`theta_caal`/`theta_diet`,
`by = ~fleet`) so a prior falls out of `linkage_spec(priors = ...)` — ~15 lines of C++ because
the prior loop already re-targets intercept rows onto the base parameter. Only attach when
`Comp_loglike == "DirichletMultinomial"`; error otherwise. Note `comp_weights` is stored
un-logged (`2-build_params.R:387`), unlike its siblings. Full spec in the plan.

**Standing design goal across everything — intuitive names.** Grant and the FIMS CIE reviews
both stress it: pick self-explanatory names for objects/switches (a stock assessor should read
them correctly on sight), prefer readable strings over positional integer codes, and when
touching a misleading legacy name (e.g. `Time_varying_*_sd_prior`, which is an input value not
a prior) rename with a deprecation alias — it's a released package, so keep back-compat.

Plan: `dev/PLAN-data-workflow-and-linkage-grammar.md` (canonical
`~/.claude/plans/shiny-hugging-yeti.md`). Commit messages plain, no AI trailer; confirm before
committing. Prefix R compile/test commands with `export PATH=/usr/bin:$PATH`.

Start by confirming the branch is still `dev-data-workflow` and green, then read the handoff.

---
