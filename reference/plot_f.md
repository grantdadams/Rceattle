# plot F

Function that plots the F time series per species from Rceattle

## Usage

``` r
plot_f(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  species = NULL,
  spnames = NULL,
  add_ci = FALSE,
  lwd = 3,
  right_adj = 0,
  width = 7,
  height = 6.5,
  minyr = NULL,
  maxyr = NULL,
  incl_proj = FALSE,
  mod_cex = 1,
  alpha = 0.4,
  mod_avg = rep(FALSE, length(Rceattle)),
  mse = FALSE,
  OM = TRUE
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

  Which species to plot e.g. c(1,4). Default = NULL plots them all

- spnames:

  Species names for legend

- add_ci:

  NOT WORKING If the confidence interval is to be added

- lwd:

  Line width as specified by user

- right_adj:

  Multiplier for to add to the right side of the figure for fitting the
  legend.

- width:

  plot width

- height:

  plot height

- minyr:

  First year to plot

- maxyr:

  max year to plot

- incl_proj:

  TRUE/FALSE, include projection years

- mod_cex:

  Cex of text for model name legend

- alpha:

  shading for confidence intervals

- mod_avg:

  is the list a model average? (DEPRECATED)

- mse:

  Is an MSE object from
  [`load_mse`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  or
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)

- OM:

  if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?

## Value

Returns and saves a figure with the population trajectory.
