# Plot biomass depletion

Plots the mean biomass depletion and 95\\

## Usage

``` r
plot_depletion(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  species = NULL,
  spnames = NULL,
  add_ci = FALSE,
  lwd = 3,
  save = FALSE,
  right_adj = 0,
  legend.pos = "topright",
  width = 7,
  height = 6.5,
  minyr = NULL,
  maxyr = NULL,
  incl_proj = FALSE,
  mod_cex = 1,
  lty = rep(1, length(Rceattle)),
  alpha = 0.4,
  mod_avg = rep(FALSE, length(Rceattle)),
  mse = FALSE,
  OM = TRUE,
  reference = NULL,
  ylab = NULL
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

  What species to include 1:nspp

- spnames:

  Species names for legend

- add_ci:

  If the confidence interval is to be added (see Details for how it is
  constructed)

- lwd:

  Line width as specified by user

- save:

  Save derived quantity?

- right_adj:

  Multiplier for to add to the right side of the figure for fitting the
  legend.

- legend.pos:

  Position of the legend as used by
  [`legend`](https://rdrr.io/r/graphics/legend.html) (default =
  "topright").

- width:

  Figure width in inches

- height:

  Figure height in inches

- minyr:

  First year to plot

- maxyr:

  max year to plot

- incl_proj:

  TRUE/FALSE, include projection years

- mod_cex:

  Cex of text for model name legend

- lty:

  line type

- alpha:

  shading for confidence intervals

- mod_avg:

  TRUE/FALSE

- mse:

  Is an MSE object from
  [`load_mse`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  or
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)

- OM:

  if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?

- reference:

  Reference model

- ylab:

  Y-axis label. `NULL` (default) derives one from `output` and the
  model's `minage`.

## Value

Returns and saves a figure with the biomass depletion trajectory.
