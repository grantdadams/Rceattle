# String\<-\>integer mapping for `srr_fun` / `srr_pred_fun` in [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)

Either form is accepted; the canonical integer code is what the TMB
template ultimately consumes. Only the structural codes (0, 2, 4) get
string aliases. The historical env-driven codes (1, 3, 5) still work
with a soft-deprecation warning – their structural part is identical to
0 / 2 / 4 respectively, and the env effect is now expressed via the
`linkages` argument to
[`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md).

## Usage

``` r
.SRR_FUNS
```

## Format

An object of class `integer` of length 3.
