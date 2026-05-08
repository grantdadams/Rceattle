# Plot selectivity

Function that plots the fishery and survey selectivity as estimated from
Rceattle

## Usage

``` r
plot_selectivity(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  width = 7,
  height = 6.5,
  species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
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

- width:

  Figure width in inches

- height:

  Figure height in inches

- species:

  Species names for legend

- lwd:

  Line width as specified by user
