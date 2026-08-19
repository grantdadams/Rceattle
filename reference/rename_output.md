# Label the derived quantities reported by a CEATTLE fit

Attaches dimension names (species, sex, age, length bin, year, fleet) to
the report objects returned by the TMB model, and names the rows/columns
of the likelihood-component matrices, so the fitted quantities read as
labeled arrays rather than bare numbers.

## Usage

``` r
rename_output(data_list = NULL, quantities = NULL)
```

## Arguments

- data_list:

  an Rceattle data_list

- quantities:

  list of "report" objects from Rceattle.
