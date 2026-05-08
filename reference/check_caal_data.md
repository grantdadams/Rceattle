# Check and Clean CAAL Data

This function checks the CAAL data for zero sum rows, handles NA values,
and verifies that the CAAL data spans the required range of ages

## Usage

``` r
check_caal_data(data_list)
```

## Arguments

- data_list:

  A list containing the following components:

  - caal_obs: A matrix or data frame of composition observations.

  - caal_data: A data frame with metadata including species, sex, and
    age/length information.

  - nsex: A vector of sex counts by species.

  - nages: A vector of age counts by species.

  - nlengths: A vector of length counts by species. throws An error if
    any rows of `caal_obs` sum to 0 or if the composition data does not
    span the range of ages

## Value

The modified `data_list` with NA values in `caal_obs` converted to 0.
cleaned_data_list \<- check_caal_data(data_list)
