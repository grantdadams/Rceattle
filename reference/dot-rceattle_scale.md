# Add the standard colorblind-safe scales to a plot

Discrete aesthetics (series identity, e.g. model) use the Okabe-Ito
qualitative palette; continuous aesthetics (ordered magnitude, e.g.
year) use viridis. Both are colorblind-safe.

## Usage

``` r
.rceattle_scale(
  p,
  discrete = TRUE,
  aesthetics = c("colour", "fill"),
  line_col = NULL
)
```

## Arguments

- p:

  A ggplot.

- discrete:

  Use the discrete (Okabe-Ito) or continuous (viridis) scale.

- aesthetics:

  Which of `"colour"`/`"fill"` to add.

- line_col:

  User-supplied colours, or `NULL` for the package default. Recycled
  over the mapped variable's levels in plotting order.

## Details

`line_col` overrides the default, supplying the palette for whichever
variable the plot maps to colour (see `?rceattle-plot-args`). On a
continuous colour aesthetic the colours are ramp anchors instead.
