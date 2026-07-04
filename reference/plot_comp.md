# Plot composition fits and residuals

Diagnostic plots for age / length composition data from a fitted
[Rceattle](https://grantdadams.github.io/Rceattle/reference/Rceattle-package.md)
model, drawn with ggplot2 for a consistent look with
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md).
With the default `residual_type = "pearson"` three figures are produced
per fleet x composition type:

- **Pearson residual bubbles** by year and bin, faceted by fleet (and,
  for joint-sex data, by sex); red = positive, blue = negative, sized by
  magnitude. The Pearson residual is \\(p - \hat p)/\sqrt{\hat p (1 -
  \hat p)/N}\\, the same form used by
  [`residuals.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/residuals.Rceattle.md).

- **Annual composition** – observed (shaded area) vs fitted (line)
  proportion at age / length, one panel per year. Joint-sex data are
  mirrored (females up, males down).

- **Aggregated composition** – the same, summed over (hindcast) years.

The shaded area and fitted line span only the observed bins (they do not
extend past the first/last bin), and bins with zero observed proportion
are retained (only `NA` bins are dropped), so the curves are not
interpolated across gaps.

## Usage

``` r
plot_comp(
  Rceattle,
  file = NULL,
  model_names = NULL,
  species = NULL,
  cex = 3,
  lwd = 3,
  right_adj = 0,
  residual_type = c("pearson", "osa")
)
```

## Arguments

- Rceattle:

  A single fitted `Rceattle` model.

- file:

  Optional filename stem; when supplied, each figure is also saved as a
  PNG.

- model_names:

  Unused (kept for back-compatibility).

- species:

  Optional species code(s) to plot (matched against the `Species`
  column); `NULL` (default) plots all. Mirrors the `species` argument of
  [`residuals.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/residuals.Rceattle.md)
  and
  [`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md).

- cex:

  Unused (kept for back-compatibility).

- lwd:

  Width of the fitted-composition line. Default `3`.

- right_adj:

  Unused (kept for back-compatibility).

- residual_type:

  `"pearson"` (default) for the ggplot2 Pearson-residual and
  composition-fit figures drawn here, or `"osa"` to instead draw the
  one-step-ahead residual diagnostics via
  [`osa_residuals()`](https://grantdadams.github.io/Rceattle/reference/osa_residuals.md)
  and
  [`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md)
  – a Q-Q plot (with SDNR / tail annotation) alongside signed OSA- and
  Pearson-residual bubbles. The `"osa"` path builds its observation data
  on demand, so it works with any fit.

## Value

Invisibly, a named list of the `ggplot` objects. Called for its side
effect of drawing (and optionally saving) the figures.
