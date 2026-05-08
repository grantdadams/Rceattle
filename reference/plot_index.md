# CPUE fits

Plot of fitted CPUE indices on natural-scale (r4ss-style)

## Usage

``` r
plot_index(
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
  error = TRUE
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

- error:

  include observed data and error bars?
