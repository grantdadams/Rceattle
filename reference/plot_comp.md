# Plot time series of comp data

Function the plots the comp data as estimated from Rceattle

## Usage

``` r
plot_comp(
  Rceattle,
  file = NULL,
  model_names = NULL,
  species = NULL,
  cex = 3,
  lwd = 3,
  right_adj = 0
)
```

## Arguments

- Rceattle:

  Single or list of Rceattle model objects exported from `Rceattle`

- file:

  name of a file to identified the files exported by the function.

- model_names:

  Names of models to be used in legend

- species:

  Species names for legend

- cex:

  Line width as specified by user

- lwd:

  Line width for observed data lines

- right_adj:

  How many units of the x-axis to add to the right side of the figure
  for fitting the legend.

## Value

Returns and saves a figure
