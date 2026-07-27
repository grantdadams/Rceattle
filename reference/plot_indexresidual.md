# Survey index residuals

Plots log residuals `log(predicted) - log(observed)` of the survey index
by year, faceted by survey fleet (`residual_type = "pearson"`, the
default), or one-step-ahead (OSA) residual diagnostics for a single fit
(`residual_type = "osa"`, via
[`osa_residuals()`](https://grantdadams.github.io/Rceattle/reference/osa_residuals.md)
/
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md)).

## Usage

``` r
plot_indexresidual(
  Rceattle,
  file = NULL,
  model_names = NULL,
  species = NULL,
  width = 7,
  height = 6.5,
  residual_type = c("pearson", "osa"),
  line_col = NULL,
  right_adj = 0,
  top_adj = 0.05,
  incl_proj = FALSE,
  single.plots = FALSE
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; if given the figure is written to
  `<file>_survey_indices.png`.

- model_names:

  Legend labels for the models.

- species:

  Species (indices) to include. Default all.

- width, height:

  Saved figure size (inches).

- residual_type:

  `"pearson"` (log index residuals) or `"osa"`.

- line_col, right_adj, top_adj, single.plots:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

- incl_proj:

  Include projection years.

## Value

A `ggplot` object.
