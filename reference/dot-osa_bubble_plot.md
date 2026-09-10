# Bubble plot of composition residuals (afscOSA styling)

The size scale is pinned to `[0, .RCE_BUBBLE_MAX]` so two figures
compare by eye, with larger residuals truncated onto it by
`.rce_truncate_resid()`.

## Usage

``` r
.osa_bubble_plot(osa, ylab = "Bin", title = "OSA residuals")
```

## Arguments

- osa:

  A data frame with `source`, `year`, `age_length_bin`, and `residual`
  columns. Bubbles are placed at (year, age/length bin); red = positive,
  blue = negative; size scales with the absolute residual; outliers
  (`|resid| > 3`) are drawn as triangles.

- ylab:

  Y-axis label (e.g. `"Age bin"` or `"Length bin"`).

- title:

  Panel title.

## Value

A `ggplot` object.
