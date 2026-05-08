# Pre-release checklist

Steps to run *before* tagging a new version of Rceattle. This is a
research / management-facing package, so the priority is reproducibility
and clear release notes — not speed.

## 1. Source state

- [ ] All in-progress changes committed; working tree clean
      (`git status`).
- [ ] `DESCRIPTION` `Version:` field bumped (semver: bump patch for
      bug fixes / docs, minor for new features, major for breaking
      API changes).
- [ ] `NEWS.md` top section heading matches the new version.
      Convert any "Unreleased" placeholder to the version number and
      add the release date.

## 2. Verify

```r
# Regenerate man/ and NAMESPACE from roxygen comments
roxygen2::roxygenise()

# Locally
devtools::document()
devtools::test()

# Vignettes (rendering, not just chunk parsing)
devtools::build_vignettes()

# CRAN-equivalent check
devtools::check(args = c("--no-manual", "--as-cran"))

# URL rot
urlchecker::url_check()

# Spell check (DESCRIPTION + man/)
devtools::spell_check()
```

All four should be 0 errors / 0 warnings. The `installed size`
NOTE is acceptable until we move large `data/` objects to
`inst/extdata/`.

## 3. Tag

```bash
# Replace X.Y.Z with the DESCRIPTION version
git tag -a vX.Y.Z -m "Rceattle vX.Y.Z"
git push origin vX.Y.Z
```

The `pkgdown` GitHub Actions workflow rebuilds the website on the
`release` event, so a GitHub Release should be drafted from the tag
(Releases → "Draft a new release" → select tag → paste the relevant
`NEWS.md` section as the body).

## 4. Confirm reproducible install

From a clean R session on a different machine (or in a `renv` sandbox):

```r
remotes::install_github("grantdadams/Rceattle@vX.Y.Z")
library(Rceattle)
packageVersion("Rceattle")  # should match X.Y.Z
citation("Rceattle")        # should list all references
example("fit_mod", package = "Rceattle", give.lines = TRUE)
```

## 5. Operational consumers

For users running Rceattle inside an assessment / MSE pipeline, communicate
that they should pin to the new tag in their lockfile (`renv.lock`) or
`DESCRIPTION` `Remotes:` field rather than tracking `master`. The
maintainer email in `DESCRIPTION` is the formal contact.

## 6. (Optional) CRAN submission

When ready:

```r
devtools::release()
```

This reads `cran-comments.md` and uploads to CRAN via the official
submission form.
