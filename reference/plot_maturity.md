# Plot maturity

Function that plots the maturity of each species

## Usage

``` r
plot_maturity(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  species = NULL,
  width = 4,
  height = 5.5,
  lwd = 3
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

- width:

  Figure width in inches

- height:

  Figure height in inches

- lwd:

  Line width as specified by user
