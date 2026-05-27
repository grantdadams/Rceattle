# Set or override the target parameter name on a linkage spec

Used by `build_*()` helpers that infer the parameter name from the list
key under which a spec is registered (e.g.
`linkages = list(K = linkage_spec(~temp))` -\> set `param = "K"`). If
the spec already names a different parameter the function errors to
surface user mistakes.

## Usage

``` r
.set_linkage_param(spec, param)
```

## Arguments

- spec:

  an `Rceattle_linkage_spec`.

- param:

  target parameter name (character scalar).

## Value

The spec, with `param` set.
