# Session prompt — GOA multispecies M-estimation convergence

> **RESOLVED (2026-07-24) — no code change; fixed-M golden kept.**
> The premise below ("`obj$fn` evaluates to NaN / `NA` objective") did **not** reproduce on
> HEAD. The M-estimating GOA-MS fit (`build_M1(M1_model = c(1,2,1))`) **converges**: finite
> objective **13440.28**, max gradient **~7e-4**, no NaN in `obj$fn`/`report()`/`jnll_comp`.
> The historical "`NA` objective" was a **reporting artifact**: with `getsd = TRUE`, `sdreport`
> fails on a non-invertible Hessian, and `TMBhelper::fit_tmb` then returns an `opt` **without**
> `$objective` (reads back as `NULL`/NA). With `getsd = FALSE` the objective is a clean 13440.28.
> Root cause of the non-PD Hessian: **M1 is confounded with fleet-8 (GOA_pollock_fishery)
> selectivity** — the classic M-vs-selectivity ridge. Verified levers:
> - `getsd = FALSE` → populated `opt$objective = 13440.28` (finite likelihood confirmed).
> - Disabling fleet-8 time-varying selectivity → **still** non-invertible (2 non-identifiable
>   *static* `sel_inf`), so the confound is with selectivity itself, not just its deviations.
> - Ruled out: `log_M1_dev` free (mapped `NA` at `M1_re=0`), M1 bounds, predation-M2, priors,
>   and any recent-merge regression.
> A PD-Hessian/identified M-estimating fit is unreachable without pinning M (e.g. an M prior),
> which would deviate from the published no-prior model. **Decision:** keep the fixed-M `goa_ms`
> golden (12866.2957391829, PD Hessian, reproducible) and treat M-estimation on the bundled
> `GOA2018SS` as an inherent identifiability limit, not a bug. Nothing further to do here.
> (Diagnostic scripts: `oracle.R`/`diag2.R`/`compare.R`/`confirm.R` in the session scratchpad.)

---

Paste this to start a fresh session. Self-contained.

---

We're on branch `dev-data-workflow` in the Rceattle package. The `/golden-check`
golden reference now covers four models (BS2017SS/MS, GOA2018 SS, GOA2018 MS). The GOA
multispecies reference currently uses **fixed M**, because the **published** GOA
multispecies recipe (Adams et al. 2022, `examples/Fit_2018_GOA_multi-species_model.R`),
which **estimates M** via `M1_model = c(1, 2, 1)`, returns an **`NA` objective** on the
bundled `GOA2018SS` data. Your task: **figure out why the M-estimating GOA multispecies fit
fails to converge, and make it converge to a stable finite objective** — then it can replace
the fixed-M reference as the `goa_ms` golden.

## What we already know (don't re-derive)
- **The failure:** the M-estimating fit's optimizer lands at a parameter set where
  `obj$fn(par)` evaluates to **NaN**, so `.fit_tmb()` records `opt$objective = NA`
  (`R/0-tmb_helpers.R`). It is **not**: a `getsd`/`sdreport` artifact (that's `tryCatch`'d →
  only NAs the SEs, `R/0-tmb_helpers.R:95-102`); a Newton-step overshoot (`newtonsteps`
  defaults to 0, so that block is skipped); or a merge regression (it fails on `main` too —
  it predates all the recent work).
- **Fixed-M GOA-MS converges fine** (obj = 12866.2957391829, grad ~0.01). So predation-M2 on
  fixed M1 is stable; it's the **joint M1 + M2 estimation** that's ill-conditioned.
- `GOA2018SS` structure is consistent with the recipe: 3 species (Pollock 1-sex, Arrowtooth
  **2-sex**, Cod 1-sex); `M1_model = c(1,2,1)` (ATF sex-specific, matches its 2 sexes);
  `other_food` all positive (no suitability divide-by-zero); `Ceq = c(2,2,2)`
  (temperature-dependent consumption → leans on `env_data` adequacy).

## The exact recipe that fails
```r
suppressMessages(pkgload::load_all(".", quiet = TRUE))   # export PATH=/usr/bin:$PATH
ss <- Rceattle::fit_mod(GOA2018SS, file = NULL, inits = NULL, estimateMode = 0,
        random_rec = FALSE, msmMode = 0, fit_control = fit_control(phase = TRUE, verbose = 0))
ms <- Rceattle::fit_mod(GOA2018SS, inits = ss$estimated_params, file = NULL,
        estimateMode = 0, M1Fun = build_M1(M1_model = c(1, 2, 1)), niter = 3,
        random_rec = FALSE, msmMode = 1, suitMode = 0,
        fit_control = fit_control(phase = TRUE, verbose = 0))
# ms$opt$objective is NA
```

## Suggested investigation angles
1. **Localize the NaN.** Instrument: build with `estimateMode = 3` (real objective, no
   optimize) at the SS-MLE inits, then `obj$fn(par)` and `obj$report()` — is the objective
   already NaN at the *starting* point, or only after the optimizer moves? Walk the reported
   quantities (M2/predation mortality, suitability, `NByage`, `Z`, log-terms) to find which
   goes NaN first. A `log`/`sqrt` of a non-positive intermediate is the usual culprit.
2. **Identifiability.** Estimating M1 *and* predation M2 simultaneously is often
   under-identified. Check whether a **prior on M** pins it: `build_M1(M1_model = c(1,2,1),
   M1_use_prior = ..., M_prior = ..., M_prior_sd = ...)`. Also check the M1 **bounds** in
   `R/4-build_parameter_bounds.R` — an unbounded/degenerate M can wander to a NaN region.
3. **Phasing / niter.** Try `phase = FALSE` (does the base optimization diverge, or is it a
   phase transition?) and different `niter` (the predation fixed-point loop). Does a slower/
   more-phased schedule stabilize it?
4. **Ground truth — the sibling repos.** `../Rceattle-models/GOA CEATTLE/` and `../GOAceattle/`
   hold the *real* Adams-et-al-2022 GOA multispecies runs. Read their exact `fit_mod()` call:
   the published fit converges, so they likely pass **bounds, an M prior, a different
   `initMode`, or a specific phasing** the bundled example omits. Mirror it.
5. **Data adequacy.** `Ceq = 2` needs an `env_data` temperature index with enough columns
   (`data_check()` gates `Cindex`); confirm `GOA2018SS$env_data` is adequate. Confirm the diet
   data and bioenergetics scalars are sane for all 3 predators.

## When it converges
- Confirm a **finite objective + reasonable max gradient**, reproducible across re-runs.
- Re-pin `goa_ms` in `dev/golden-ref.rds` (regenerate via `/golden-check`) with the
  M-estimating recipe, update the `goa_ms` recipe in `.claude/commands/golden-check.md`, and
  update the `golden-reference-state` memory. Keep the **other three** models bit-reproducible.
- Federal-rigor standard: verify constructively (it converges *and* is identified — e.g. the
  Hessian is PD / a profile is well-behaved), not just "it ran".
```
