# Management strategy evaluation performance metric summary

Management strategy evaluation performance metric summary

## Usage

``` r
mse_summary(mse, om_only = FALSE)
```

## Arguments

- mse:

  MSE runs from
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  or
  [`load_mse`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)

- om_only:

  only include performance metrics from OMs

## Value

Alist of two data.frames with MSE summary statistics of performance
metrics including: data.frame 1

1.  Average annual catch across projection years and simulations per
    fleet and across fleets

2.  Average interannual variation in catch (IAV) across projection
    years (n) per fleet and across fleets

3.  % of years in which the fishery is closed across simulations (s)

4.  Average relative mean squared error in estimate of spawning biomass
    in the terminal year across simulations

5.  % of years in which the population is perceived as undergoing
    overfishing as determined from F_Limit across simulations via
    [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
    in the EM

6.  % of years in which the population is perceived to be overfished as
    determined from B_Limit across simulations via
    [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
    in the EM

7.  % of years in which the population is undergoing overfishing as
    determined from the “true” F_Limit across simulations via
    [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
    in the OM

8.  % of years in which the population is overfished as determined from
    the “true” B_Limit across simulations via
    [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
    in the OM

9.  Average ratio of spawning biomass over B_target in the terminal year
    across simulations in the OM 10-14. Terminal biomass, SSB, SSB
    depletion (relative to equilibrium), SSB depletion (relative to
    dynamic SB0)
