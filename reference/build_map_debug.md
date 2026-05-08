# Helper to set map for debug mode

Sets all parameters in the map list to NA, except the `dummy` parameter,
for use in debug or testing modes.

## Usage

``` r
build_map_debug(map_list, debug)
```

## Arguments

- map_list:

  The current TMB map list.

- debug:

  Logical. If TRUE, debug mode is activated.

## Value

Updated `map_list`.
