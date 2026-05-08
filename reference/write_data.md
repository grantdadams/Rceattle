# Write data file

Write data file

## Usage

``` r
write_data(data_list, file = "Rceattle_data.xlsx")
```

## Arguments

- data_list:

  Rceattle data_list object

- file:

  Filname to be used. Must end with '.xlsx'

## Examples

``` r

library(Rceattle)
data(BS2017SS)
out_file <- file.path(tempdir(), "BS2017SS.xlsx")
write_data(data_list = BS2017SS, file = out_file)
file.remove(out_file)
#> [1] TRUE
```
