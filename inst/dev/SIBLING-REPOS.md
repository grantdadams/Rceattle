# The sibling assessment repos

`../Rceattle-models` and `../GOA-ATF-ESP` are live consumers of this package's API. A breaking
change here breaks scripts that produce federal catch advice.

- **`../Rceattle-models`** — EBS/GOA pollock, sablefish, arrowtooth, plaice, POP, hake.
- **`../GOA-ATF-ESP`** — GOA arrowtooth and its multispecies (cannibalism) run: **the only live
  two-sex, `suitMode = 0` model**, so it is what exercises the sexed and predation paths.
- **Ignore `EBS_CEATTLE_TMB`** — a vendored fork, not a consumer.

Fitted `*.rds` are ~50 MB each. Keep them out of git.

## Sweeping

```
grep -rn "<symbol>" --include=*.R "../Rceattle-models" "../GOA-ATF-ESP"
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

Traps:

- The pollock scripts' `Data/` paths are relative to the **project** root, not the year folder.
- **The ATF script cannot be sourced straight through on any version** — it references three
  objects it never assigns (`:364`, `:480`, `:570`, the last gating the whole final figure
  block). This is a property of the script, not of your change.
- Its `file =` arguments write **into the assessment repo**. Run it from a sandbox that
  symlinks `Data/`.
- Force plots through `ggplot2::ggplot_build()`. A figure that assembles but cannot render is
  not a pass.
