# Bubble plot of composition residuals (afscOSA styling)

Bubble plot of composition residuals (afscOSA styling)

## Usage

``` r
.osa_bubble_plot(osa, ylab = "Bin", title = "OSA residuals")
```

## Arguments

- osa:

  A data frame with `source`, `year`, `age_length_bin`, and `residual`
  columns. Bubbles are placed at (year, age/length bin); red = positive,
  blue = negative; size and transparency scale with the absolute
  residual; outliers (`|resid| > 3`) are drawn as triangles.

- ylab:

  Y-axis label (e.g. `"Age bin"` or `"Length bin"`).

- title:

  Panel title.

## Value

A `ggplot` object.
