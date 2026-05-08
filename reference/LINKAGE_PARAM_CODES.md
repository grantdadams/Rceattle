# Integer codes for the `param` column, per process

Each process has its own independent parameter namespace (`log_K` means
something different to growth than `log_alpha` does to recruitment), so
the param column is encoded against a per-process table. Unknown
parameter names error.

## Usage

``` r
LINKAGE_PARAM_CODES
```

## Format

An object of class `list` of length 5.
