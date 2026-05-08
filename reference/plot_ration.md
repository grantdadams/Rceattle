# Plot ration

Function that plots the ration across ages (minage:nages) as estimated
from Rceattle. Returns and saves a figure with the ration trajectory.
Ration is multiplied by biomass-at-age/sex to get population level
estimates

## Usage

``` r
plot_ration(
  Rceattle,
  file = NULL,
  minage = 1,
  model_names = NULL,
  line_col = NULL,
  spnames = NULL,
  species = NULL,
  lwd = 3,
  lty = 1,
  right_adj = 0,
  minyr = NULL,
  width = 7,
  height = 6.5,
  incl_proj = FALSE,
  incl_mean = FALSE,
  add_ci = FALSE
)
```

## Arguments

- Rceattle:

  Single or list of Rceattle model objects exported from `Rceattle`

- file:

  name of a file to identified the files exported by the function.

- minage:

  minage to plot ration (i.e. age "minage"+)

- model_names:

  Names of models to be used in legend

- line_col:

  Colors of models to be used for line color

- spnames:

  Species names for legend

- species:

  Which species to plot e.g. c(1,4). Default = NULL plots them all

- lwd:

  Line width as specified by user

- lty:

  Line type

- right_adj:

  Multiplier for to add to the right side of the figure for fitting the
  legend.

- minyr:

  first year to plot

- width:

  Figure width in inches

- height:

  Figure height in inches

- incl_proj:

  TRUE/FALSE include projections years

- incl_mean:

  TRUE/FALSE include time series mean as horizontal line

- add_ci:

  TRUE/FALSE, includes 95 percent confidence interval
