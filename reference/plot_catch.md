# Landings fits

Plot of fitted landings data on natural-scale (r4ss-style)

## Usage

``` r
plot_catch(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  species = NULL,
  right_adj = 0,
  top_adj = 1.2,
  incl_proj = FALSE,
  single.plots = FALSE,
  width = NULL,
  height = NULL,
  alpha = 0.4,
  lwd = 2,
  ymax = NULL,
  maxyr = NULL,
  mse = FALSE,
  error = TRUE,
  fleets = NULL
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

  How much top side of the y-axis for fitting the legend. As percentage
  (default = 1.2).

- incl_proj:

  TRUE/FALSE include projections years

- single.plots:

  if TRUE plot invidual fits else make multiplot

- width:

  plot width

- height:

  plot hight

- alpha:

  Transparency of lines for MSE plots

- lwd:

  Line width

- ymax:

  Fleet specific upper ylim

- maxyr:

  Max year to plot

- mse:

  Is if an MSE object from
  [`load_mse`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  or
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md).
  Will plot data from OMs.

- error:

  Include observed data with error bars?

- fleets:

  Which fishing fleets to include (defaults to all = NULL)
