# Sample historical recruitment deviations into the projection

Sample historical recruitment deviations into the projection

## Usage

``` r
sample_rec(
  object = NULL,
  sample_rec = TRUE,
  update_model = TRUE,
  rec_trend = 0,
  Rceattle = NULL
)
```

## Arguments

- object:

  CEATTLE model object exported from `Rceattle`

- sample_rec:

  Include resampled recruitment deviations from the hindcast in the OM
  projection. Resampled deviations are used rather than drawing from
  N(0, sigmaR) because the initial deviations bias R0 low. If FALSE,
  uses the mean recruitment deviation.

- update_model:

  Update model dynamics. Default = TRUE

- rec_trend:

  Linear increase or decrease in mean recruitment from `endyr` to
  `projyr`. This is the terminal multiplier
  `mean rec * (1 + (rec_trend/projection years) * 1:projection years)`.
  Can be of length 1 or of length nspp. If length 1, all species get the
  same trend.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

Rceattle model
