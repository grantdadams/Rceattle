# Selectivity specification

Carries environmental linkages on selectivity parameters. The effect on
a parameter is written as a formula and composes additively with any
`Time_varying_sel` process error on the same fleet (the two are separate
mechanisms: a covariate effect versus a deviation).

The parameter names are the shape parameters of the parametric
selectivity forms:

- `slp_asc`, `slp_desc`:

  ascending / descending logistic slope (log scale); for a double-normal
  the ascending / descending width, aliased `sigma_asc` / `sigma_desc`.

- `inf_asc`, `inf_desc`:

  ascending / descending inflection age/length (natural scale); for a
  double-normal the peak and the logit right-floor, aliased `peak` /
  `right_floor`.

- `coff`:

  non-parametric selectivity-at-bin coefficients.

Every parameter accepts `link = "log"` (multiplicative on the natural
parameter) or `link = "identity"` (additive), like the other processes.

**Priors on a selectivity parameter.** An intercept-only formula (`~ 1`)
with a `priors` entry places a prior on the selectivity parameter itself
(no year-to-year offset is added). Read the prior on the parameter's own
scale: the slopes (`slp_asc` / `slp_desc`) are on the log scale (use
`lognormal()`), the inflections (`inf_asc` / `inf_desc`) on the natural
scale (use `normal()`). See Examples for a normal prior on the ascending
inflection. This mirrors the prior-only
[`build_composition()`](https://grantdadams.github.io/Rceattle/reference/build_composition.md)
path.

A selectivity prior targets one parameter, so in a two-sex model an
unstratified `~ 1` prior constrains sex 1 only – use `by = ~ sex` for a
per-sex prior. An `init` on a selectivity intercept has no effect (the
starting value comes from the data), and a prior on the double-normal
`right_floor` is not supported. For a fleet that mirrors another fleet's
selectivity (shared `Selectivity_index`), place the prior on the lead
fleet so the shared parameter block is not penalized more than once.

## Usage

``` r
build_selectivity(linkages = NULL)
```

## Arguments

- linkages:

  Optional named list of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by selectivity parameter. Coefficients are per fleet by
  default (`by = ~ fleet`); use the `fleet` argument of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  to restrict a spec to particular fleets.

## Value

A list of selectivity settings for
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).

## Examples

``` r
# \donttest{
# A cold-pool effect on the ascending inflection of a logistic fleet
build_selectivity(linkages = list(
  inf_asc = linkage_spec(~ cold_pool, by = ~ fleet)))
#> $linkages
#> $linkages$inf_asc
#> <Rceattle linkage spec>
#>   param:   inf_asc
#>   formula: ~cold_pool
#>   link:    log
#> 
#> 

# A normal prior on the ascending inflection (intercept-only formula)
build_selectivity(linkages = list(
  inf_asc = linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 3)))))
#> $linkages
#> $linkages$inf_asc
#> <Rceattle linkage spec>
#>   param:   inf_asc
#>   formula: ~1
#>   prior:    (Intercept) ~ normal(0, 3)
#>   link:    log
#> 
#> 
# }
```
