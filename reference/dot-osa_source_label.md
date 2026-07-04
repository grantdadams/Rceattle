# Build a two-row data-source label for residual panel headers

The first row is the fleet name from `fleet_control` (falling back to
the fleet code); the second row is the data type, tagged with whether
composition bins are ages or lengths. The two rows are separated by a
newline so they render as a two-line facet strip.

## Usage

``` r
.osa_source_label(osa)
```

## Arguments

- osa:

  A data frame with `source`, `fleet`, and (optionally) `fleet_name` and
  `index_label` columns.

## Value

Character vector of labels (one per row).
