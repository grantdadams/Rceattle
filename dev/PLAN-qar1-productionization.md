# QAR1 (Rogers et al. 2024) productionization — status

> **North-star:** the end goal is to **match the afsc-assessments/GOApollock**
> assessment (https://github.com/afsc-assessments/GOApollock), the official ADMB
> GOA pollock model, not merely the Rceattle legacy `Estimate_q = 6` path. That
> model **estimates** the QAR1 observation SD, so the grammar must too
> (decision below), and its env-covariate handling is the parity target.

## Locked design decisions (this session)

- **Estimable `obs_sd`** — implement it (GOApollock estimates the obs SD). Contract:
  the `obs_sd` value is the *start*; estimated by default; fixable (for a known
  measurement error) via the same `est_phase`/fix mechanism as `sigma`/`rho`.
  Needs a **simulation self-consistency** check (β and `obs_sd` are only jointly
  identified when `obs_sd` is informative).
- **env_data mismatch = skip un-observed years.** Auto-extend `env_data` to the
  model year range with `NA` for missing years; apply the QAR1 observation ONLY
  where the covariate is present (no fabricated observations — do NOT replicate the
  legacy mean-fill). The latent AR1 still spans all years. Verify density
  equivalence (bit-identity vs legacy) on a *full-coverage* env_data case instead.


The grammar expresses the Rogers-2024 environmental-index-driven QAR1 catchability
model as `ar1(1 | Year)` with `observe = "<env column>", obs_sd = <value>`: an AR1
latent q-deviation observed against an env_data column with a fixed measurement SD,
entering q through an estimated effect size (`beta_linkage_obs`). The controlled
density-math was verified exact (~1e-13) in the RE-density session.

## Done

- **(b) Effect sizes are reported** (commit 7ea7e1f9). `beta_linkage` (coefficients),
  `beta_linkage_re` (deviations), and `beta_linkage_obs` (QAR1 effect size) are
  `REPORT`'d into `fit$quantities`; coefficients + effect size are `ADREPORT`'d so
  they carry a standard error. Golden bit-identical.

- **(c) Estimable `obs_sd`** (implemented). The QAR1 observation SD is now an
  estimated parameter `log_obs_sd_linkage` (one per observed group), started from
  the spec `obs_sd`, replacing the fixed DATA value — matching the reference
  `Estimate_q = 6` / GOApollock. Golden bit-identical (length 0 with no QAR1).
  **Validated on the real pollock QcovPol data**: `obs_sd` estimates to a sensible
  0.60 (β = 0.65, σ = 0.58) — non-degenerate. On a *smooth synthetic* covariate it
  collapses toward 0 (the documented β/obs_sd joint-identifiability caveat: only
  identified when the covariate is informative). **Still TODO:** a prior / explicit
  fixed-`obs_sd` option (the "fixable" half of the decision) for weakly-informative
  cases; a full simulate-and-recover check against GOApollock.

- **env_data skip-un-observed years** — NOT yet implemented (next). The materialize
  already yields `re_obs_value = NA` for years absent from env_data; the remaining
  work is: auto-extend env_data to `styr:endyr`, encode a per-slot observed mask,
  and have the cpp observation skip masked slots (see the machinery map below).

## GOA-pollock cross-check (part d) — findings

Reference: `../Rceattle-models/GOA pollock/2024 pollock model.R` + its data
(`Data/GOA_24_pollock_single_species_1970-2024.xlsx`), which carries the `QcovPol`
catchability covariate in `env_data`. (The base model itself fits **RW** q on the
acoustic survey (fleet 1, `Time_varying_q = 4`), not the QAR1 — the QAR1 is the
`Estimate_q = 6` capability the covariate enables.)

- **The grammar QAR1 fits the real pollock data end-to-end** — after extending
  `env_data` back to `styr` (see gap 1), a grammar `ar1(observe="QcovPol")` q-linkage
  on fleet 1 builds and returns a finite objective, the effect size, and 55 yearly
  deviations (1970–2024).
- **It does NOT yet bit-identically reproduce the legacy `Estimate_q = 6`.** At
  matched initial process σ (`Time_varying_q_sd[1] = 0.038`) and obs SD
  (`Q_sd_prior[1] = 1`), the estimateMode-3 objectives differ by ~129
  (legacy 430409.7 vs grammar 430539.0). The divergence is **not** in the core
  density (that was verified exact on a controlled fixture) but in the real-config
  env handling below.

## Open gaps (issues to file / future work)

1. **env_data must span `styr`; the legacy path auto-fills missing years with the
   mean.** `QcovPol` starts at 1983 but the model starts at 1970. The grammar
   *errors* ("env_data must start at styr"); the legacy `Estimate_q = 5/6` path fills
   pre-data years with the covariate mean internally. To reproduce the legacy fit the
   user currently has to hand-extend `env_data`. **And the mean-fill / observed-year
   set almost certainly differs between the manual extension and the legacy internal
   fill — the prime suspect for the ~129 objective gap.** Productionizing means the
   grammar should either auto-extend (mean-fill, matching legacy) or the observation
   should skip un-observed years; decide and match the legacy exactly, then re-run the
   bit-identity check.

2. **Estimable `obs_sd`** *(deferred at Grant's request — file as an issue)*. The
   grammar fixes `obs_sd` (a DATA value, `linkage_re_obs_sd`); the legacy
   `Estimate_q = 6` **estimates** the observation SD (`index_q_log_sd`) alongside the
   process SD, effect size, and rho. So a fully-estimated legacy QAR1 cannot be
   reproduced by the grammar without adding an estimable-`obs_sd` option: a
   `PARAMETER` (log obs SD) per observed group, mapped estimable when requested, with
   the observation density using the estimated value. Note (from the RE-density work):
   β and `obs_sd` are only jointly identified when `obs_sd` is informative — so this
   needs a simulation self-consistency check on recovery, not just "it converged".

3. **read/write round-trip.** `env_data` (the `observe` column) already round-trips
   through `write_data`/`read_data`. Persisting the `observe`/`obs_sd` *spec* fields
   (a `build_catchability()` argument, not xlsx data) is the `save_config()` /
   `load_config()` feature (unbuilt), not an xlsx concern.

## Next session

Resolve gap 1 (env_data span/mean-fill parity with the legacy path) and re-run the
GOA-pollock bit-identity cross-check; then, if wanted, gap 2 (estimable `obs_sd`)
with a simulation-recovery check. Gap 3 folds into `save_config()`.
