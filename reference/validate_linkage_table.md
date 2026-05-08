# Validate a linkage table against the schema

Checks that all required columns are present with the correct R type,
that no required values are `NA`, and that enum-like columns (`process`,
`link`) only contain allowed values.

## Usage

``` r
validate_linkage_table(x)
```

## Arguments

- x:

  object to validate.

## Value

`x` invisibly, on success. Throws an error otherwise.
