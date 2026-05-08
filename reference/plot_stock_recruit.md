# Plot stock recruit function

Function the plots the stock recruit function as estimated from Rceattle

## Usage

``` r
plot_stock_recruit(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  width = 7,
  height = 6.5,
  species = NULL,
  spnames = NULL,
  lwd = 3,
  lty = 1,
  incl_proj = FALSE,
  plot_env = FALSE,
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

- width:

  Figure width in inches

- height:

  Figure height in inches

- species:

  Which species to plot e.g. c(1,4). Default = NULL plots them all

- spnames:

  Species names for legend

- lwd:

  Line width as specified by user

- lty:

  Line type

- incl_proj:

  TRUE/FALSE, include projection years for environmental relationship

- plot_env:

  TRUE/FALSE, plot environmental covariate relationship

- mod_cex:

  Cex of text for model name legend
