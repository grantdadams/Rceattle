# Run a management strategy evaluation

Runs a forward-projecting management strategy evaluation (MSE).
Projected selectivity, catchability, foraging days, and weight-at-age
are held at the operating model's terminal hindcast year. Survey SD is
set to the average over the historical time series, and composition
sample size is held at the last year. There is no implementation error
and no observation error on catch.

## Usage

``` r
run_mse(
  om,
  em,
  nsim = 10,
  start_sim = 1,
  assessment_period = 1,
  sampling_period = 1,
  simulate_data = TRUE,
  regenerate_past = FALSE,
  sample_rec = TRUE,
  rec_trend = 0,
  fut_sample = 1,
  cap = NULL,
  catch_mult = NULL,
  seed = 666,
  regenerate_seed = seed,
  loopnum = 1,
  file = NULL,
  dir = NULL,
  timeout = 999,
  endyr = NA,
  cores = NULL
)
```

## Arguments

- om:

  CEATTLE model object exported from `Rceattle`

- em:

  CEATTLE model object exported from `Rceattle`

- nsim:

  Number of simulations to run (default 10)

- start_sim:

  First simulation number to start at. Useful if the code stops at
  specific seed/sim (default = 1).

- assessment_period:

  Period of years that each assessment is taken

- sampling_period:

  Period of years data sampling is conducted. Single value or vector the
  same length as the number of fleets.

- simulate_data:

  Include simulated random error proportional to that estimated/provided
  for the data from the OM.

- regenerate_past:

  Refits the EM to historical/conditioning data prior to the MSE, where
  the data are generated from the OM with `simulate_data = TRUE` or
  without `simulate_data = FALSE` sampling error.

- sample_rec:

  Include resampled recruitment deviations from the hindcast in the OM
  projection. Resampled deviations are used rather than drawing from
  N(0, sigmaR) because the initial deviations bias R0 low. If FALSE,
  uses the mean recruitment deviation.

- rec_trend:

  Linear increase or decrease in mean recruitment from `endyr` to
  `projyr`. This is the terminal multiplier
  `mean rec * (1 + (rec_trend/projection years) * 1:projection years)`.
  Can be of length 1 or of length nspp. If length 1, all species get the
  same trend.

- fut_sample:

  future sampling effort relative to last year.
  ` Log_sd * 1 / fut_sample` for index and ` Sample_size * fut_sample`
  for comps

- cap:

  A cap on the catch in the projection. Can be a single number applied
  to all species (proportional to recommended catch) or vector of length
  `nspp` applied to each species. Default = NULL

- catch_mult:

  A multiplier for the catch in the projection. Can be a single number
  or vector of length nspp. Default = NULL

- seed:

  seed for the simulation

- regenerate_seed:

  seed for regenerating data

- loopnum:

  number of times to re-start optimization (where `loopnum=3` sometimes
  achieves a lower final gradient than `loopnum=1`)

- file:

  (Optional) Filename where each OM simulation with EMs will be saved.
  If NULL, no files are saved.

- dir:

  (Optional) Directory where each OM simulation is saved

- timeout:

  length of time (minutes) estimation will run before stopping a sim
  (default 999 minutes)

- endyr:

  Terminal year of the MSE projection. Default = NA uses `projyr` from
  the operating model.

- cores:

  Number of cores to use for parallel simulations. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

## Value

A list of operating models (differ by simulated recruitment determined
by `nsim`) and estimation models fit to each operating model (differ by
terminal year).
