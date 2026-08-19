# Default-fill and tidy a data list before fitting

Fills optional blocks with correctly-shaped empties, filters
observations to the model's year range, extends catch to the projection
years, and re-keys any survey-index covariance matrices. Whether a
missing block is actually a problem is left to
[`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md),
which knows the model configuration.

## Usage

``` r
clean_data(data_list)
```

## Arguments

- data_list:

  Rceattle data list
