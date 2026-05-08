# Specify the harvest control rule (HCR) used for Rceattle

Defines the harvest control rule and associated reference points used in
the Rceattle model.

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

  Harvest control rule to use. Accepts an integer or equivalent string
  alias. Default = 0. See Details for the full list of available
  options.

- DynamicHCR:

  TRUE/FALSE. Whether to use static or dynamic reference points (default
  = FALSE).

- Ftarget:

  Target fishing mortality rate (yr^-1) (SPR or depletion based) or
  input F for projections. For example, if Ftarget is spr F40%, enter
  0.40.

- Flimit:

  Limit fishing mortality rate (yr^-1) (SPR or depletion based). For
  example, if Flimit is spr F35%, enter 0.35.

- Ptarget:

  Target spawning-stock biomass as a percentage of static or dynamic
  spawning-stock-biomass at F = 0 (accounts for recruitment).

- Plimit:

  Limit spawning-stock biomass as a percentage of static or dynamic
  spawning-stock-biomass at F = 0 (accounts for recruitment).

- Alpha:

  Parameter used in NPFMC Tier 3 HCR.

- Pstar:

  Quantile used for uncertainty buffer given `Flimit` and `Sigma`.

- Sigma:

  Standard deviation used for normally distributed uncertainty buffer
  given `Flimit` and `Pstar`.

- Fmult:

  Multiplier on the target F (default = 1). Used to scale the target
  fishing mortality following NEFSC convention.

- HCRorder:

  For multi-species models, the order in which to project fishing (e.g.,
  predators first, then prey).

## Value

A `list` containing the harvest control rule and associated biological
reference points.

## Details

**Harvest control rule formulations currently implemented in Rceattle:**

`hcr = 0` or `"NoFishing"`: No catch. Estimate the hindcast.

`hcr = 1` or `"CMSY"`: CMSY. Maximize catch across all species
simultaneously. CMSY can be constrained such that depletion does not
fall below `Plimit`.

`hcr = 2` or `"ConstantF"`: Constant input F set at `Ftarget` for each
species (vector or single F). SPR (single-species only) based Flimit is
specified via `Flimit`.

`hcr = 3` or `"ConstantFSSB"`: F that achieves `Ftarget`% of SSB0 in the
end of the projection.

`hcr = 4` or `"ConstantFSPR"`: Constant Fspr set at `Ftarget` for each
species. Can be multiplied by `Fmult` following NEFSC.

`hcr = 5` or `"NPFMC"`: The North Pacific Fishery Management Council
(NPFMC) Tier 3 spawner-per-recruit-based harvest control rule: Stock
status: \\SB \> SB at Ftarget\\ \$\$Fofl = Flimit\$\$ \$\$Fuse =
Ftarget\$\$ Stock status: \\Alpha \< SB / SB at Ftarget \le 1\\ \$\$Fofl
= Flimit \* (SB / Ptarget - Alpha) / (1 - Alpha)\$\$ \$\$Fuse = Ftarget
\* (SB / Ptarget - Alpha) / (1 - Alpha)\$\$ Stock status: \\SB / SB at
Ftarget \le Alpha\\ or \\SB \< Plimit \* SB0\\ \$\$Fofl = 0\$\$ \$\$Fuse
= 0\$\$

`hcr = 6` or `"PFMC"`: An HCR based on the The Pacific Fishery
Management Council (PFMC) category 1 40-10 annual catch limit (ABC)
harvest control rule assuming Fofl is normally distributed with a
standard deviation (sigma) = 0.5 and an uncertainty quantile buffer
(P\*) of 0.45 (PFMC 2020). The model uses Fspr if single-species or F
that achieves X% of SSB0 for multi-species. Target biological reference
points are calculated based on the normal cumulative distribution
function (Phi(quantile,mean,standard deviation)), P\*, and sigma as
follows: Stock status: \\SB \> SB0 \* Ptarget\\ \$\$Fofl = Flimit\$\$
\$\$Fuse = Phi(Pstar, Flimit, Sigma)\$\$ Stock status: \\SB0 \* Plimit
\< SB \le SB0 \* Ptarget\\ \$\$Fofl = Flimit\$\$ \$\$Fuse = Phi(Pstar,
Flimit, Sigma) \* \frac{SB0 \* Ptarget \* (SB - SB0 \* Plimit)}{SB \*
SB0 \* (Ptarget - Plimit)}\$\$ Stock status: \\SB \< SB0 \* Plimit\\
\$\$Fofl = 0\$\$ \$\$Fuse = 0\$\$

`hcr = 7` or `"SESSF"`: An HCR based on the The Southern and Eastern
Scalefish and Shark Fishery (SESSF) spawner-per-recruit-based Tier 1
harvest control rule where F_Limit=F\_(20%), B_Limit=SB_20, F_Target
(AFMA 2017) calculated as follows: Stock status: \\SB \> SB0 \*
Ptarget\\ \$\$Fofl = Flimit\$\$ \$\$Fuse = Ftarget\$\$ Stock status:
\\Ptarget \> SB \> SB0 \* Plimit\\ \$\$Fofl = Flimit \* (SB / (SB0 \*
Ptarget) - 1)\$\$ \$\$Fuse = Ftarget \* (SB / (SB0 \* Ptarget) - 1)\$\$
Stock status: \\SB \< SB0 \* Plimit\\ \$\$Fofl = 0\$\$ \$\$Fuse = 0\$\$

NOTE: only HCRs 1, 2, 3, and 6 will work in multi-species mode.
