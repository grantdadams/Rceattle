# Plot selectivity

Selectivity over the hindcast, one panel per fleet. Each fleet is drawn
on the dimension it was estimated on: age for an age-based fleet, length
bin for a length-based one (`Selectivity_dimension`). A model mixing the
two gets one figure per dimension.

## Usage

``` r
plot_selectivity(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  width = 7,
  height = 6.5,
  species = NULL,
  lwd = 3,
  spnames = NULL,
  lty = 1,
  minyr = NULL,
  maxyr = NULL,
  alpha = 0.25,
  colour_by = c("auto", "year", "model")
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; the figure is written to `<file>_<suffix>.png` if
  given.

- model_names:

  Legend labels for the models.

- line_col:

  Line colours; names, hex, or base-graphics integers. `NULL` uses the
  colorblind-safe Okabe-Ito palette. Applied to whichever variable the
  figure separates by colour, in legend order. Too few colours are
  recycled, with a warning.

- width, height:

  Saved figure size in inches.

- species:

  Species to include, as indices (`c(1, 3)`), names, a logical mask, or
  `"all"`. Default `NULL` plots every species. Species **labels** belong
  in `spnames`; a character `species` that matches no species name is
  read as labels for back-compatibility.

- lwd:

  Line width on the base-graphics scale: the default `3` renders as a
  standard-weight ggplot line. A vector varies it across series.

- spnames:

  Species labels, length `nspp`. Default: the model's own.

- lty:

  Line type. A vector varies it across the levels of whatever the figure
  separates by line type.

- minyr, maxyr:

  First / last year to plot.

- alpha:

  Transparency of the faintest year when the fan is drawn by
  transparency, i.e. under `colour_by = "model"`.

- colour_by:

  What colour separates: `"year"` (a fan), `"model"`, or `"auto"` (the
  default) for year with a single model and model with several.

## Value

A `ggplot`, or for a model mixing age- and length-based fleets a named
list of them, one per dimension.

## What the colours mean

With one model, colour is the year, so time-varying selectivity reads as
a fan and a time-invariant fleet collapses to a single line. With
several models, colour separates the models and the year fan moves to
transparency, faintest at the earliest year drawn. The transparency
scale spans the years shown across all models, so the same year is the
same shade everywhere and a short retrospective peel stops short of
solid. `colour_by` overrides the choice either way, and `alpha` sets the
faintest end.

Line type separates the sexes, and `lty` supplies its values. Panels are
fleets, so `spnames` does not label anything here – it only lets
`species` select by name.
