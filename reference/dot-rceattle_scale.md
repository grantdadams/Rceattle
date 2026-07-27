# Add the standard colourblind-safe scales to a plot

Discrete aesthetics (series identity, e.g. model) use the Okabe-Ito
qualitative palette; continuous aesthetics (ordered magnitude, e.g.
year) use viridis. Both are colourblind-safe.

## Usage

``` r
.rceattle_scale(p, discrete = TRUE, aesthetics = c("colour", "fill"))
```

## Arguments

- p:

  A ggplot.

- discrete:

  Use the discrete (Okabe-Ito) or continuous (viridis) scale.

- aesthetics:

  Which of `"colour"`/`"fill"` to add.
