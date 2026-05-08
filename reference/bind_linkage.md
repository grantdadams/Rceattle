# Row-bind one or more linkage tables, preserving the schema

Wraps `rbind` and re-applies validation and class. Accepts either a list
of tables or several tables passed as separate arguments. Empty inputs
are skipped.

## Usage

``` r
bind_linkage(...)
```

## Arguments

- ...:

  linkage tables, or a single list of them.

## Value

An `Rceattle_linkage_table`.
