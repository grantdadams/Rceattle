# TODO: stock-recruit curves in multispecies models

Status: **open**. Opened 2026-09-11 with 5.30.0 (`8c5fe1ba`) and 5.31.0, after
two adversarial reviews of each.

Spawning biomass per recruit (SPR) is undefined under predation: total mortality
carries M2, which scales with predator abundance, so the template leaves `SPR0`
and `SPRFinit` at 0 when `msmMode > 0`. Everything below follows from that.

## Where it stands

- **5.30.0**
  - The Ianelli penalty is allowed under predation (`srr_fun = "mean"`, curve as
    `srr_pred_fun`).
  - A Beverton-Holt steepness prior is refused. An alpha prior goes through a
    linkage, which acts on `rec_pars` itself: `normal()` is on natural-scale
    alpha; `prior_lognormal()` is `dnorm(log alpha)`, centred so `exp(meanlog)`
    is alpha's mean under `bias_adjust_proc` (5.33.0) and its median without.
  - Species with input numbers-at-age (`estDynamics > 0`) carry no curve terms.
  - Dynamic B0 keeps those species at their input numbers.
  - `fit_mod()` warns when a supplied map fixes the curve.
- **5.31.0**
  - A curve fitted in the hindcast is allowed under predation, in every
    initMode (see item 1).
  - `R_init` is the free level `exp(rec_pars[, "R0"])`; steepness is 0; the
    Ricker `posfun(alpha * SPR0 - 1)` penalty is skipped.
  - `sample_rec(sample_rec = FALSE)` and `retrospective()` project
    `log(mean(R / R_hat))` = `log(mean(exp(rec_dev)))`, which is
    mean-unbiased.

## Measured on hake (4 species, `MSE_yr2024.R` setup, MSVPA, 2026-09-11)

| Configuration | Objective | Hake alpha, beta | Terminal SSB vs baseline |
|---|---|---|---|
| NonEquilibrium, mean recruitment (baseline) | 2663.81 | — | — |
| Ianelli BH penalty, alpha/beta free | 2732.20 | 18.8, 2.7e-6 | +0.04% |
| FreeParams, mean recruitment (control) | 2520.14 | — | −34.7% |
| FreeParams, BH in the hindcast | 2520.08 | 4.5e8, 57 (flat) | −34.7% |
| NonEquilibrium, BH in the hindcast, free `R_init` | 2654.13 | 5.2e8, 58.5 (flat) | −24.9% |

Measured on 5.31.0. 5.33.0 mean-centred the lognormal priors and the Ianelli
penalty, so these objectives cannot be compared with 5.33.0 fits.

- **The penalty row's objective** excludes the fixed predators' curve penalty
  (5,137.9 nats), which 5.30.0 removed. Objectives are not comparable across
  initModes: FreeParams drops the initial-abundance density (−28.8 nats on
  hake).
- **The hake data do not inform a curve fitted in the hindcast.** Both hindcast
  fits run to the flat ridge, recruitment ≈ alpha/beta ≈ 8e6 at any SSB (the
  predicted/asymptote ratio is 1.000 at median SSB). Their Hessians are
  positive definite but near-singular (condition 1.5e14 and 6.5e13).
- **FreeParams drives the shift, not the curve.** It drops M1 from 0.265 to
  0.190 and mean recruitment by 43%, and moves SSB −56% in 1980 and −35% at the
  end. The 1980 ages 16–20 fall to 259–450 fish, with log-scale standard errors
  of 2.7–4.0. The composition fit improves only 6.7 nats.

## Open items, most important first

1. **Done in 5.31.0: every initMode is allowed.** Under a hindcast curve under
   predation, `R_init = R0` for every initMode (`ceattle.cpp` 6.3), and none
   reads SPR. Still open:
   - Under 3/4, `R_init` is confounded with Finit.
   - NonEquilibrium converged on the fixture when started from the
     mean-recruitment fit (gradient 0.05), but collapsed from default starts.
     Add a fitted test and document warm-starting.
2. **NonEquilibrium decays the initial ages with M1 only**, so under predation
   the initial deviations absorb the missing M2. Hake's mean `init_dev` is −1.6,
   against age-1 M2 of 0.81 plus σ²/2 of 0.77.
   - Decaying with `M_at_age(..., 0)` (M1 plus year-1 M2) under `msmMode > 0`
     would fix it. It is a starting age structure, not a reference point, so it
     stays SPR-free.
   - It moves every multispecies NonEquilibrium fit, including the hake
     baseline (2663.8053181057). Re-record baselines.
3. **Anchor the curve.** Under predation nothing ties alpha and beta to an
   equilibrium, so they run to the flat ridge.
   - On the 2-species fixture log alpha reached 702, next to double-precision
     overflow at 709.8, and a 1-peel retrospective failed with NaN gradients.
     Holding alpha at its single-species value cost 2.5 nats.
   - `rec_pars` has no bounds.
   - Options:
     - bounds on `rec_pars[, 2:3]`;
     - an alpha prior through a linkage as the documented configuration;
     - item 4's diagnostic.
4. **Add an SPR-free degenerate-curve check to `.check_stock_recruit()`.** Under
   predation it returns a NOTE before it looks at alpha or beta, so the flat
   hake curve passes. Check:
   - the predicted/asymptote ratio across the observed SSB range;
   - alpha or beta at a boundary, or with NaN standard errors;
   - optionally, a clearly labelled M1-only steepness as a diagnostic.
5. **Projected-deviation conventions: done in 5.33.0.** The Ianelli branch now
   projects `log(mean(R / R_hat))` over the penalty years, and its
   penalty is mean-centred under `bias_adjust_proc`. The median form ran hake
   projections at about 45% of the curve's mean. Single-species keeps
   `log(mean R) - log R0`. Still open: a test that expected projected R agrees
   with and without resampling.
6. **The Ianelli dynamic B0 mixes bases.** It applies `curve × exp(rec_dev)`,
   with `rec_dev` about the mean (`ceattle.cpp` ~2052–2055), rather than
   `log R - log R_hat`. The fitted `log R_hat - log R0` reaches 0.61, so
   no-fishing recruitment can be off by up to 1.84× at equal SSB.
7. **Multispecies SB0 (`MSSB0`).**
   - It is derived only when `HCR != "NoFishing"` (`fit_mod.R`); otherwise it is
     the 999 mt placeholder.
   - It comes from a mean-recruitment projection to `projyr`, which is not an
     equilibrium: at M ≈ 0.27, cohorts alive at `endyr` still carry about 33%
     of SSB seven years on.
   - Verify whether HCRs 3/6 can reach the objective with the placeholder
     (reviewer: `ceattle.cpp` ~2136–2183, 4543–4548). The other session's
     `test-switches-hcr-multispecies.R` covers the same area.
   - Options:
     - a converged no-fishing projection;
     - reporting depletion against dynamic B0;
     - an in-model multispecies equilibrium (fixed-point M2), which item 8
       needs.
8. **Steepness under predation.**
   - `h = 0.2(1 + beta S0)/(1 + 0.2 beta S0)` is steepness only when `S0` is
     the curve's own equilibrium.
   - If it is needed, define it from an in-model equilibrium S0 and report it
     as conditional steepness.
   - Alpha (recruits per unit female SSB, the slope at the origin) is what
     carries over from a single-species fit, not h. It must be restated in the
     multispecies model's SSB units and definition.
9. **Simulation recovery.**
   - No test shows the hindcast curve or the Ianelli penalty recovering a known
     alpha/beta. The fixture's SSB spans only 1.7×, too little contrast; build
     one with real contrast.
   - Also measure the σR inflation the Ianelli double density implies under
     `random_rec` (up to √2 by the algebra; `review4/re.R` did not finish).
10. **Per-species `srr_fun` / `srr_pred_fun`**, so hake can carry a curve while
    the predators stay on mean recruitment. Both are global today.
11. **MSE resampling.** `sample_rec()` resamples independent hindcast
    deviations, including year 1 and years outside the `srr_hat` window. A
    block bootstrap within the window would keep autocorrelated recruitment.
12. **Smaller.**
    - `retrospective()`'s predation guard is untested, and the peel refit
      overwrites the column it writes.
    - A year-varying R0 linkage is inert on the hindcast-curve path.
    - `R_init` starts at exp(9) = 8103, and the alpha/beta defaults are
      off-scale for tonnes; document seeding with `srr_alpha_init` /
      `srr_beta_init`.
    - The 5.31.0 tests are build-only (`estimateMode = 3`); add one fitted
      case.
13. **What to report with a multispecies stock-recruit fit:**
    - the SSB-recruitment pairs with the curve;
    - the predicted/asymptote ratio;
    - alpha and beta with their standard errors;
    - sensitivity to initMode;
    - how σR was configured, and the projection convention;
    - the dynamic B0 series;
    - a note that objectives are not comparable across these configurations.

    Document what the Ianelli penalty fits: recruitment centred halfway between
    the mean and the curve, at half the variance.
