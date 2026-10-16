# Specify the harvest control rule (HCR) used for Rceattle

Defines the harvest control rule and its associated reference points.
The rules are tabulated, with which ones work under predation, in
[`vignette("projections-and-reference-points")`](https://afsc-assessments.github.io/Rceattle/articles/projections-and-reference-points.md);
the formulations of the three ramped rules are below.

## Usage

``` r
build_hcr(
  HCR = 0,
  DynamicHCR = FALSE,
  Ftarget = 0.4,
  Flimit = 0.35,
  Ptarget = 0.4,
  Plimit = 0,
  Alpha = 0.05,
  Pstar = 0.45,
  Sigma = 0.5,
  Fmult = 1,
  HCRorder = 1
)
```

## Arguments

- HCR:

  Harvest control rule, as an integer or string alias: `0`/`"NoFishing"`
  (no catch, hindcast only), `1`/`"CMSY"` (maximize joint catch,
  optionally holding depletion at or above `Plimit`), `2`/`"ConstantF"`
  (constant input F at `Ftarget`), `3`/`"ConstantFSSB"` (the F reaching
  `Ftarget` percent of SSB0 by the end of the projection),
  `4`/`"ConstantFSPR"` (constant Fspr, scalable by `Fmult` following
  NEFSC), `5`/`"NPFMC"`, `6`/`"PFMC"` or `7`/`"SESSF"`. Default `0`.
  **Only 0, 1, 2, 3 and 6 work in multispecies mode.**

- DynamicHCR:

  TRUE/FALSE. Use dynamic rather than static reference points (default
  FALSE).

- Ftarget:

  Target fishing mortality rate (yr^-1), SPR- or depletion-based, or the
  input F for projections. For SPR F40%, enter 0.40.

- Flimit:

  Limit fishing mortality rate (yr^-1), SPR- or depletion-based. For SPR
  F35%, enter 0.35.

- Ptarget:

  Target spawning-stock biomass as a proportion of static or dynamic
  spawning-stock biomass at F = 0. For B40%, enter 0.40.

- Plimit:

  Limit spawning-stock biomass on the same basis. For B10%, enter 0.10.

- Alpha:

  Parameter used in the NPFMC Tier 3 rule.

- Pstar:

  Quantile for the uncertainty buffer, given `Flimit` and `Sigma`.

- Sigma:

  Standard deviation of the normal uncertainty buffer, given `Flimit`
  and `Pstar`.

- Fmult:

  Multiplier on the target F (default 1), scaling it following NEFSC
  convention.

- HCRorder:

  For multispecies models, the order in which to project fishing
  (predators first, then prey, for instance).

## Value

A `list` containing the harvest control rule and associated biological
reference points.

## Details

`hcr = 5` or `"NPFMC"`, the North Pacific Fishery Management Council
Tier 3 spawner-per-recruit rule: Stock status: \\SB \> SB at Ftarget\\
\$\$Fofl = Flimit\$\$ \$\$Fuse = Ftarget\$\$ Stock status: \\Alpha \< SB
/ SB at Ftarget \le 1\\ \$\$Fofl = Flimit \* (SB / SB\_{Ftarget} -
Alpha) / (1 - Alpha)\$\$ \$\$Fuse = Ftarget \* (SB / SB\_{Ftarget} -
Alpha) / (1 - Alpha)\$\$ Stock status: \\SB / SB at Ftarget \le Alpha\\
or \\SB \< Plimit \* SB0\\ \$\$Fofl = 0\$\$ \$\$Fuse = 0\$\$

`hcr = 6` or `"PFMC"`, the Pacific Fishery Management Council category 1
40-10 annual catch limit rule (PFMC 2020), using Fspr in single species
or the F reaching X% of SSB0 under predation. The uncertainty buffer is
`qnorm(Pstar, Flimit, Sigma)`, which equals
`Flimit + qnorm(Pstar, 0, Sigma)` and so sits *below* `Flimit` when
\\Pstar \< 0.5\\. The taper runs between \\SB0 \* Plimit\\ and \\SB0 \*
Ptarget\\, so the 40-10 shape needs `Ptarget = 0.40` and
`Plimit = 0.10`; the default `Plimit = 0` gives a 40-0 rule. Stock
status: \\SB \> SB0 \* Ptarget\\ \$\$Fofl = Flimit\$\$ \$\$Fuse =
qnorm(Pstar, Flimit, Sigma)\$\$ Stock status: \\SB0 \* Plimit \< SB \le
SB0 \* Ptarget\\ \$\$Fofl = Flimit\$\$ \$\$Fuse = qnorm(Pstar, Flimit,
Sigma) \* \frac{SB0 \* Ptarget \* (SB - SB0 \* Plimit)}{SB \* SB0 \*
(Ptarget - Plimit)}\$\$ Stock status: \\SB \< SB0 \* Plimit\\ \$\$Fofl =
0\$\$ \$\$Fuse = 0\$\$

`hcr = 7` or `"SESSF"`, the Southern and Eastern Scalefish and Shark
Fishery Tier 1 spawner-per-recruit rule, with F_Limit = F\_(20%) and
B_Limit = SB_20 (AFMA 2017): Stock status: \\SB \> SB0 \* Ptarget\\
\$\$Fofl = Flimit\$\$ \$\$Fuse = Ftarget\$\$ Stock status: \\SB0 \*
Ptarget \> SB \> SB0 \* Plimit\\ \$\$Fofl = Flimit \* (SB / (SB0 \*
Plimit) - 1)\$\$ \$\$Fuse = Ftarget \* (SB / (SB0 \* Plimit) - 1)\$\$
Stock status: \\SB \< SB0 \* Plimit\\ \$\$Fofl = 0\$\$ \$\$Fuse = 0\$\$

## Examples

``` r
# Tier 3 NPFMC control rule: F40% as FABC, F35% as FOFL.
build_hcr(HCR = "NPFMC", DynamicHCR = FALSE, Ftarget = 0.4, Flimit = 0.35)
#> $HCR
#> [1] "NPFMC"
#> 
#> $DynamicHCR
#> [1] FALSE
#> 
#> $Ftarget
#> [1] 0.4
#> 
#> $Flimit
#> [1] 0.35
#> 
#> $Ptarget
#> [1] 0.4
#> 
#> $Plimit
#> [1] 0
#> 
#> $Alpha
#> [1] 0.05
#> 
#> $Pstar
#> [1] 0.45
#> 
#> $Sigma
#> [1] 0.5
#> 
#> $Fmult
#> [1] 1
#> 
#> $HCRorder
#> [1] 1
#> 

# No fishing: the projection a B0 reference point is read from.
build_hcr(HCR = "NoFishing")
#> $HCR
#> [1] "NoFishing"
#> 
#> $DynamicHCR
#> [1] FALSE
#> 
#> $Ftarget
#> [1] 0.4
#> 
#> $Flimit
#> [1] 0.35
#> 
#> $Ptarget
#> [1] 0.4
#> 
#> $Plimit
#> [1] 0
#> 
#> $Alpha
#> [1] 0.05
#> 
#> $Pstar
#> [1] 0.45
#> 
#> $Sigma
#> [1] 0.5
#> 
#> $Fmult
#> [1] 1
#> 
#> $HCRorder
#> [1] 1
#> 
```
