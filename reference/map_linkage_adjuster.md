# Helper to turn off base linked parameters

This function automatically sets parameters to NA (not estimated) if
they include linkaged.

## Usage

``` r
map_linkage_adjuster(map_list, data_list)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  an Rceattle data_list (with the pooled `linkage_table` from
  [`pool_linkages()`](https://grantdadams.github.io/Rceattle/reference/pool_linkages.md)).

## Value

Updated `map_list`.
