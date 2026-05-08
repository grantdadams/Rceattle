# Function to load .RDs files from MSE runs

Function to load .RDs files from MSE runs

## Usage

``` r
load_mse(dir = NULL, file = NULL, exclude = NULL, include_em = TRUE)
```

## Arguments

- dir:

  Directory used to save files from
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)

- file:

  file name used to save files from
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)

- exclude:

  index of MSE simulations not to load

- include_em:

  whether the EMs should be loaded or not (default = TRUE)

## Value

list of MSE simulations/run
