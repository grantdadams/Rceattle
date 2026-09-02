# The sibling assessment repos

`../Rceattle-models` and `../GOA-ATF-ESP` are live consumers of this package's API. A breaking
change here breaks scripts that produce federal catch advice.

- **`../Rceattle-models`** — EBS/GOA pollock, sablefish, arrowtooth, plaice, POP, hake.
- **`../GOA-ATF-ESP`** — GOA arrowtooth and its multispecies (cannibalism) run: **the only live
  two-sex, `suitMode = 0` model**, so it is what exercises the sexed and predation paths.
- **`../GOA-multispecies-assessment`** — the 2026 GOA three-species assessment chapter (pollock,
  arrowtooth, cod), submitted to the September Plan Team. **The widest single exercise of the API
  in the ecosystem**: `fit_mod()` in both modes, `retrospective()`, `jitter()`, `self_test()`,
  `osa_residuals()` / `osa_diagnostics()`, `run_mse()` / `check_mse()` / `mse_summary()`,
  `build_hcr()` (rules 3 and 5), `save_config()` / `load_config()`, `residuals()`, and most of
  `plot_*()`. `run_all.R` is the entry point and gates each stage on its inputs.
  Sixteen fleets, three species (nages 10/21/10), one two-sex species.
- **Ignore `EBS_CEATTLE_TMB`** — a vendored fork, not a consumer.

Fitted `*.rds` are ~50 MB each. Keep them out of git.

## Sweeping

```
grep -rn "<symbol>" --include=*.R \
  "../Rceattle-models" "../GOA-ATF-ESP" "../GOA-multispecies-assessment"
```

Two limits, and both need saying out loud when you report a sweep:

- **It catches removed or renamed API. It does not catch behavioural drift.** A script that
  still parses can still produce different numbers.
- **Some models there are partially implemented and do not run regardless**, so not every hit
  needs chasing.

The v4.7.0 ggplot migration is the cautionary case: `plot_*()` began returning ggplot objects,
so every `plot_x(...); mtext(...)` chain in those scripts failed with "plot.new has not been
called yet". A grep for the function names would have found them; nobody ran one.

## Running them, when a sweep is not enough

Neither repo caches a fitted object, so verifying against a real assessment means refitting.
That is cheaper than it sounds: the terminal fit is under a minute for either pollock model,
~13 min for the ATF chain.

| Model | Entry point |
|---|---|
| GOA pollock 2025 | `../Rceattle-models/GOA pollock/2025/04-fit-and-diagnostics.R` |
| EBS pollock 2024 | `../Rceattle-models/EBS pollock/2024/04-fit-and-diagnostics.R` |
| GOA arrowtooth | `../GOA-ATF-ESP/R/Run_2025_ceattle.R` |
| Pacific hake MSE | `../Rceattle-models/Pacific hake/04-mse.R` |
| GOA 3-species assessment | `../GOA-multispecies-assessment/run_all.R` |

**The hake MSE is the one script that runs `run_mse()` end to end**, and the only routine
exercise of three-species predation with estimated suitability, of `suitMode` differing per
predator, and of Dirichlet-multinomial comps carrying a prior on their own weight. Golden
covers none of that: its four models are single- and multi-species Bering Sea and Gulf of
Alaska hindcasts, with no MSE and no estimated suitability. Run it after touching predation,
suitability, the DM likelihood, `sim_mod()`, or `run_mse()`.

Its four fits take ~3.5 min together on an M-series Mac, plus ~2 min for `nsim = 2, cores = 2`.
Reference objectives (verified equal on `dev` and `pr3-schema-order-and-qar1`, 2026-08-22):

| Stage | -log L |
|---|---|
| single-species | 2133.8207228717 |
| single-species + category-1 HCR | 2134.4713926593 |
| MSVPA, estimated M | 2137.4433306648 |
| estimated suitability | 2260.7063099135 |

All four were bit-identical across the two branches (delta 0.000e+00), as was the vulnerability
matrix: 0.8172 (arrowtooth to hake) and 0.7686 (sablefish to hake). The script's own
inline comments give the first three ~5 higher and the fourth as 2262.318: those are Rceattle
5.6.1 numbers that still carried the `theta_diet` prior constants. `README.md` in that folder
records the "clean" values, which are what the current package reproduces.

Traps:

- **`run_mse(cores > 1)` works under `pkgload::load_all()`** because `.parallel_lapply()` forks;
  a PSOCK fallback would not see the loaded tree. Do not assume a parallel MSE failure is the
  model.
- `mse_summary()` returns a ragged list (`species`, `fleet`, `total`, `meta`), not a data frame.
  `as.data.frame()` on it errors -- that is the caller's bug, not a broken MSE.

- The pollock scripts' `Data/` paths are relative to the **project** root, not the year folder.
- **`../GOA-multispecies-assessment` runs entirely off saved fits.** `R/03`–`R/07` load
  `models/GOA_26_mod_list.RData` rather than refitting, so a breaking change surfaces there as a
  script error, not as a moved number, and the saved objects can lag the package. Fits are
  ~1 min each; `run_all.R` rebuilds them from `R/02`. Its `R/07_figures_tables.R` renders every
  figure the chapter uses, so it is a cheap end-to-end check of the `plot_*()` surface: source it
  and confirm the manifest it prints reports no new skips.
- **The ATF script cannot be sourced straight through on any version** — it references three
  objects it never assigns (`:364`, `:480`, `:570`, the last gating the whole final figure
  block). This is a property of the script, not of your change.
- Its `file =` arguments write **into the assessment repo**. Run it from a sandbox that
  symlinks `Data/`.
- Force plots through `ggplot2::ggplot_build()`. A figure that assembles but cannot render is
  not a pass.
