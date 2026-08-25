# Model average of derived quantities

Model average of derived quantities

## Usage

``` r
model_average(
  object = NULL,
  weights = NULL,
  uncertainty = FALSE,
  nboot = 10000,
  Rceattle = NULL
)
```

## Arguments

- object:

  list of Rceattle model objects

- weights:

  vector of weights to be used for weighting models

- uncertainty:

  TRUE/FALSE Sample uncertainty across derived quantities using weighted
  bootstrap from the asymptotic distribution of MLEs

- nboot:

  Number of bootstraps taken from asymptotic distribution of MLEs.
  Default = 10000

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

an Rceattle object with derived quantities weighted by the specified
weights. The length of the derived quantities spans the years which
overlap across all models.

## Examples

``` r
if (FALSE) { # \dontrun{
# Equal weights across two candidate models.
model_average(list(fit1, fit2), weights = c(1, 1))
# Weighted by AIC, with bootstrapped uncertainty.
model_average(list(fit1, fit2), weights = c(0.7, 0.3), uncertainty = TRUE)
} # }
```
