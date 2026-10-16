# Helper to set map for Selectivity parameters

Maps base selectivity parameters (`log_sel_slp`, `sel_inf`, `sel_coff`)
and time-varying deviations, based on `Selectivity` and
`Time_varying_sel` settings in `fleet_control`.

## Usage

``` r
build_map_selectivity(map_list, data_list, nyrs_hind, random_sel)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

- nyrs_hind:

  Number of historical years.

- random_sel:

  Logical indicating if selectivity deviations are random effects.

## Value

Updated `map_list`.

## Details

`Selectivity` sets the curve's shape and `Time_varying_sel` whether it
varies over time; for the allowed values of both, see "3. Selectivity
(per fleet)" in
[`vignette("model-options-and-functionality")`](https://afsc-assessments.github.io/Rceattle/articles/model-options-and-functionality.md),
or
[rceattle_switch_maps](https://afsc-assessments.github.io/Rceattle/reference/rceattle_switch_maps.md)
for the code-to-name maps themselves. `N_sel_bins` is how many age or
length bins the non-parametric and AR1 forms estimate (`Selectivity` 2,
5, 6, 7, 9 or 13); it is unused otherwise. `random_sel` in
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
treats the deviations as random effects and estimates their variance.

The double-normal forms blend Gaussian ascending and descending limbs at
the peak. `DoubleNormal` reads its four parameters from the shared
arrays: `sel_inf[1]` = peak, `sel_inf[2]` = logit(right_floor) (the
right-tail floor, SS3 P6/`end_logit`), `log_sel_slp[1]` = log(sigma_asc)
and `log_sel_slp[2]` = log(sigma_desc). right_floor -\> 0 is
dome-shaped; right_floor -\> 1 is logistic ascending only.
`DoubleNormalSS3` instead carries its own six-parameter array
(`sel_dn6`).
