# Rerun with F = 0.

Function to update hindcast and set F to 0. Useful for determining
dynamic reference points for multi-species models under climate-change.

## Usage

``` r
remove_F(object = NULL, Rceattle = NULL)
```

## Arguments

- object:

  A fitted Rceattle model object

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.
