---
description: Run the testthat suite (optionally one file/filter) with NOT_CRAN=true
argument-hint: "[file path or name filter, e.g. test-likelihood-osa-residuals]"
---

Run the Rceattle tests with `NOT_CRAN=true` (so `skip_on_cran()` blocks actually run) and
`options(testthat.max_fails = Inf)`. Always prefix with `export PATH=/usr/bin:$PATH`.

**No argument** → run the whole suite:

```
export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript -e 'options(testthat.max_fails = Inf); devtools::test()'
```

**With `$ARGUMENTS`** → resolve it to a repo-relative path under `tests/testthat/` (accept a
partial name; the snippet takes the resolved path as-is, so do not prefix it again) and run just that file with the shared helpers loaded — the project idiom:

```
export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript -e '
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  e <- new.env(parent = asNamespace("Rceattle"))   # a plain new.env() fails: "could not find function \"data_check\""
  for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, envir = e)
  options(testthat.max_fails = Inf)
  testthat::test_file("<resolved-path>", env = e, reporter = "progress")'
```

**If `src/TMB/*.cpp` or `*.hpp` changed since the last load**, recompile first and run the
suite serially -- `DESCRIPTION` sets `Config/testthat/parallel: true`, and the workers cannot
load a freshly rebuilt DLL. The suite then aborts with `getDLLRegisteredRoutines.DLLInfo`,
which is a toolchain failure, not a test failure:

```
export PATH=/usr/bin:$PATH && Rscript -e 'pkgload::load_all(".")' && \
  NOT_CRAN=true TESTTHAT_PARALLEL=false Rscript -e 'options(testthat.max_fails = Inf); devtools::test()'
```

Report pass/fail counts and show the first failures with file:line. Don't fix anything
unless I ask.
