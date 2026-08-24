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

Measured at 5.15.0 with `rcmdcheck(args = c("--no-manual", "--as-cran"))`
on macOS, R 4.5.x. The note is `checking CRAN incoming feasibility`, and
it flags three things, all expected for a first submission:

- **New submission.**
- **`Remotes` in `DESCRIPTION`.** A devtools field for installing
  `TMBhelper` from GitHub during development. It must be removed before
  submission; CRAN does not accept it.
- **`TMBhelper` not in a mainstream repository.** See the submission note
  below: it is in `Suggests:` and the package works without it.

Installed size is reported as INFO rather than a note on this toolchain
(14.6 Mb: `data` 5.1 Mb, `libs` 4.3 Mb, `R` 2.0 Mb, `extdata` 1.3 Mb),
driven by the bundled example datasets in `data/` (Bering Sea, Gulf of
Alaska, Georges Bank, and Atka mackerel applications shipped to make the
vignettes runnable without an external download). It may still be raised
as a note on other platforms. We are tracking moving the largest of these
to `inst/extdata/`.

**The source tarball is 33.9 Mb, which is well over what CRAN accepts.**
Shrinking `data/` is a prerequisite for submission, not just a tidy-up.

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
