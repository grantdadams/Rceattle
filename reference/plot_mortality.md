# Plot M1 + M2

Function that plots the M1 and M2 as estimated from Rceattle

## Usage

``` r
plot_mortality(
  Rceattle,
  file = NULL,
  incl_proj = FALSE,
  zlim = NULL,
  type = 0,
  width = 8,
  height = 5.5,
  title = NULL,
  log = FALSE,
  minyr = NULL,
  theta = 155,
  species = NULL,
  maxage = NULL,
  title_cex = 10,
  M2 = TRUE
)
```

## Arguments

- Rceattle:

  Single or list of Rceattle model objects exported from `Rceattle`

- file:

  name of a file to identified the files exported by the function.

- incl_proj:

  Include the projection years (TRUE/FALSE)

- zlim:

  zlim for M1 + M2 plots. Character - use max range across species in
  model. NULL - use species specific ranges. Vector of two.

- type:

  0 = Tiles, 1 = contour, 2 = facet lines, 3 = persp

- width:

  Plot width when saved "inches"

- height:

  Plot height when saved "inches"

- title:

  Additional title to add. Will also add species names if not NULL

- log:

  TRUE/FALSE use log M1 + M2

- minyr:

  First year to plot

- theta:

  theta for persp plot

- species:

  Species to plot. Plots all if null.

- maxage:

  Plot up to this age. Plots all ages if NULL

- title_cex:

  Font size for title

- M2:

  TRUE/FALSE Use M2 only (True) or total M (False)
