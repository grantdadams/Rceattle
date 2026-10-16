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

  Path/prefix for the saved figure. `NULL` (default) draws without
  saving.

- model_names:

  Names of models to be used in legend

- line_col, lwd:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored; colour comes from the shared model palette and the line
  width is fixed.

- species:

  Species labels for the facet strips, length `nspp`. Default is the
  model's own names. It labels panels and does not subset.

- width:

  Figure width in inches

- height:

  Figure height in inches
