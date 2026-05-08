# Plot diet composition fits

If year == 0, diet data are averaged from suit_styr to suit_endyr If
prey_age \>= 0 diet data are diet proportion of prey-at-age in
predator-at-age If prey_age \< 0 diet data are diet proportion of
prey-spp in predator-at-age (sum across prey ages) If prey_age \< 0 and
pred_age \< 0, diet data are mean diet proportion of prey-spp in
predator-spp (sum across prey ages and take mean across predator ages)
If prey_age \< 0 and pred_age \< -500, diet data are weighted mean diet
proportion of prey-spp in predator-spp (sum across prey ages and take
weighted mean across predator ages)

## Usage

``` r
plot_diet_comp(Rceattle, file = NULL, species = NULL)
```

## Arguments

- Rceattle:

  Single or list of Rceattle model objects exported from `Rceattle`

- file:

  name of a file to identified the files exported by the function.

- species:

  Species names for legend

## Value

Returns and saves a figure
