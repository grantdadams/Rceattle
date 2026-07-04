# One bin-side (age or length) composition column: Q-Q + OSA + Pearson bubbles

One bin-side (age or length) composition column: Q-Q + OSA + Pearson
bubbles

## Usage

``` r
.osa_comp_side(comp, pear, side)
```

## Arguments

- comp:

  Composition rows with `source` and `.side` columns.

- pear:

  Reshaped Pearson rows with `source` and `.side` columns, or `NULL`.

- side:

  `"age"` or `"length"`.

## Value

A stacked `cowplot`/`ggplot` object, or `NULL` if no rows on that side.
