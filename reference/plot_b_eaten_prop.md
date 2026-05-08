# Plot biomass consumed of each prey species by predator

Function that plots the biomass consumed trends as estimated from
Rceattle. Returns and saves a figure with the biomass eaten trajectory.

## Usage

``` r
plot_b_eaten_prop(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  spnames = NULL,
  species = NULL,
  lwd = 3,
  right_adj = 0,
  top_adj = 0.15,
  minyr = NULL,
  mohns = NULL,
  width = 7,
  height = 6.5,
  incl_proj = FALSE,
  incl_mean = FALSE,
  add_ci = FALSE,
  mod_cex = 1
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

- spnames:

  Species names for legend

- species:

  Which species to plot e.g. c(1,4). Default = NULL plots them all

- lwd:

  Line width as specified by user

- right_adj:

  Multiplier for to add to the right side of the figure for fitting the
  legend.

- top_adj:

  Adjustment for top margin

- minyr:

  first year to plot

- mohns:

  data.frame of mohn's rows extracted from
  [`retrospective`](https://grantdadams.github.io/Rceattle/reference/retrospective.md)

- width:

  Figure width in inches

- height:

  Figure height in inches

- incl_proj:

  TRUE/FALSE include projections years

- incl_mean:

  TRUE/FALSE include horizontal long term mean

- add_ci:

  TRUE/FALSE, includes 95 percent confidence interval

- mod_cex:

  Cex of text for model name legend
