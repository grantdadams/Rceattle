# Check and Clean Composition Data

This function checks the composition data for zero sum rows, handles NA
values, and verifies that the composition data spans the required range
of ages/lengths.

## Usage

``` r
check_composition_data(data_list)
```

## Arguments

- data_list:

  A list containing the following components:

  - comp_obs: A matrix or data frame of composition observations.

  - comp_data: A data frame with metadata including species, sex, and
    age/length information.

  - nsex: A vector of sex counts by species.

  - nages: A vector of age counts by species.

  - nlengths: A vector of length counts by species. throws An error if
    any rows of `comp_obs` sum to 0 or if the composition data does not
    span the range of ages/lengths.

## Value

The modified `data_list` with NA values in `comp_obs` converted to 0.

## Examples

``` r
data_list <- list(
  comp_obs = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2),
  comp_data = data.frame(Species = c(1, 2), Sex = c(1, 1), Age0_Length1 = c(0, 0)),
  nsex = c(1, 1),
  nages = c(3, 3),
  nlengths = c(3, 3)
)
cleaned_data_list <- check_composition_data(data_list)
```
