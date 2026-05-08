# Model average of derived quantities

Model average of derived quantities

## Usage

``` r
model_average(Rceattle, weights = NULL, uncertainty = FALSE, nboot = 10000)
```

## Arguments

- Rceattle:

  list of Rceattle model objects

- weights:

  vector of weights to be used for weighting models

- uncertainty:

  TRUE/FALSE Sample uncertainty across derived quantities using weighted
  bootstrap from the asymptotic distribution of MLEs

- nboot:

  Number of bootstraps taken from asymptotic distribution of MLEs.
  Default = 10000

## Value

an Rceattle object with derived quantities weighted by the specified
weights. The length of the derived quantities spans the years which
overlap across all models.
