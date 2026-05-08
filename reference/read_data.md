# Read a CEATTLE excel data file

Read a CEATTLE excel data file

## Usage

``` r
read_data(file = "Rceattle_data.xlsx")
```

## Arguments

- file:

  Filname to be used. Must end with '.xlsx'

## Examples

``` r

library(Rceattle)
data(BS2017SS)
out_file <- file.path(tempdir(), "BS2017SS.xlsx")
write_data(data_list = BS2017SS, file = out_file)
data_list <- read_data(file = out_file)
#> Adding 'Month' column to 'fleet_control' with default value of 0
file.remove(out_file)
#> [1] TRUE
```
