# Discrete palette function: Okabe-Ito for up to 8 series, else interpolate

Mirrors the arity contract of `scale_*_viridis_d` (accepts any `n`) so
the scale never errors on a many-model overlay; beyond 8 series it
interpolates the Okabe-Ito anchors (rare enough that CVD-optimality is
not guaranteed).

## Usage

``` r
.okabe_ito_pal(n)
```
