# CPUE residual

Plot of residuals CPUE indices on log-scale (r4ss-style)

## Usage

``` r
plot_indexresidual(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  species = NULL,
  right_adj = 0,
  top_adj = 0.05,
  incl_proj = FALSE,
  single.plots = FALSE,
  width = NULL,
  height = NULL,
  residual_type = c("pearson", "osa")
)
```

## Arguments

- Rceattle:

  Single or list of Rceattle model objects exported from `Rceattle`

- file:

  name of a file to identified the files exported by the function.

- model_names:

  Names of models to be used in legend

- line_col:

  Colors of models to be used for line color

- species:

  Species names for legend

- right_adj:

  How much right side of the x-axis for fitting the legend. As
  percentage.

- top_adj:

  How much top side of the y-axis for fitting the legend. As percentage.

- incl_proj:

  TRUE/FALSE include projections years

- single.plots:

  if TRUE plot invidual fits else make multiplot

- width:

  plot width

- height:

  plot hight

- residual_type:

  `"pearson"` (default) for the legacy index/catch residual plots, or
  `"osa"` to draw one-step-ahead residual diagnostics (Q-Q plot with
  SDNR / tail annotation and residual-by-year) for a single fitted model
  via
  [`osa_residuals()`](https://grantdadams.github.io/Rceattle/reference/osa_residuals.md)
  and
  [`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md).
