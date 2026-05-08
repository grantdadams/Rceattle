# Internal: data mask exposing short prior names to NSE-aware callers

Returned as a named list rather than an environment so it can be passed
straight to
[`rlang::eval_tidy()`](https://rlang.r-lib.org/reference/eval_tidy.html)
as a `data` mask. Names are bound to the corresponding `prior_*`
constructors so user code such as `priors = list(temp = normal(0, 1))`
resolves correctly without masking
[`base::gamma()`](https://rdrr.io/r/base/Special.html) or
[`base::beta()`](https://rdrr.io/r/base/Special.html) at the package
level.

## Usage

``` r
.prior_dispatch_mask()
```
