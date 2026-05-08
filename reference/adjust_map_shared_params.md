# Helper to adjust map for shared catchability/selectivity indices

Enforces parameter sharing by mapping parameters for fleets with a
common `Selectivity_index` or `Q_index` to the same value as the initial
index.

## Usage

``` r
adjust_map_shared_params(map_list, data_list)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

## Value

Updated `map_list`.
