# cran-comments.md

This file is the message that accompanies a CRAN submission via
`devtools::release()` or `devtools::submit_cran()`. It is excluded from
the package tarball via `.Rbuildignore`.

## Test environments

- macOS Sonoma (local), R 4.5.x release
- Ubuntu 22.04 (GitHub Actions, r-lib/actions), R release / devel / oldrel-1
- Windows Server (GitHub Actions, r-lib/actions), R release
- win-builder (devel and release): TBD before submission
- R-hub: TBD before submission (`rhub::rhub_check()`)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the `installed size is 11.3 Mb` flag, driven by the bundled
example datasets in `data/` (Bering Sea, Gulf of Alaska, Georges Bank,
and Atka mackerel applications shipped to make the vignettes runnable
without an external download). We are tracking moving the largest of
these to `inst/extdata/` for the next release.

## Downstream dependencies

There are currently no reverse dependencies on CRAN.

## Submission notes

This is the first CRAN submission of `Rceattle`. The package implements
the multispecies CEATTLE assessment model
(Holsman et al. 2016, *Deep Sea Res. II*; Adams et al. 2022,
*Fisheries Research*) on top of TMB.

`TMBhelper` (a GitHub-only package) is **not** required: it is in
`Suggests:` and `fit_mod()` falls back to plain `nlminb` + `sdreport`
when it is not installed, so CRAN's check will exercise the
no-`TMBhelper` path automatically.
