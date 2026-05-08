# Sample historical recruitment deviates and place in the projection

Sample historical recruitment deviates and place in the projection

## Usage

``` r
sample_rec(Rceattle, sample_rec = TRUE, update_model = TRUE, rec_trend = 0)
```

## Arguments

- Rceattle:

  CEATTLE model object exported from `Rceattle`

- sample_rec:

  Include resampled recruitment deviates from the"hindcast" in the
  projection of the OM. Resampled deviates are used rather than sampling
  from N(0, sigmaR) because initial deviates bias R0 low. If false, uses
  mean of recruitment deviates.

- update_model:

  Update model dynamics. Default = TRUE

- rec_trend:

  Linear increase or decrease in mean recruitment from `endyr` to
  `projyr`. This is the terminal multiplier
  `mean rec * (1 + (rec_trend/projection years) * 1:projection years)`.
  Can be of length 1 or of length nspp. If length 1, all species get the
  same trend.

## Value

Rceattle model
