
testthat::skip_on_cran()

# Tests are organized into subdirectories under tests/testthat/ for easier
# navigation (tests-Dynamics, tests-Mortality, ...). testthat does not
# discover subdirectory tests automatically, so this master file pulls
# them in.
#
# We deliberately use source() instead of testthat::test_dir() because
# nested test_dir() invocations under test_check() can trigger an
# 'evaluation nested too deeply: infinite recursion' error inside
# rlang's trace deparser whenever any test fails (the deparser hits a
# recursion limit when formatting deeply-nested test_that backtraces
# from a nested reporter). Direct sourcing keeps test_that() blocks in
# the outer reporter where they aggregate cleanly.

source(testthat::test_path("helpers-make-msm-data.R"))
source(testthat::test_path("helpers.R"))

subdirs <- c(
  "tests-Dynamics",
  "tests-Mortality",
  "tests-Likelihoods",
  "tests-Data-processing",
  "tests-Selectivity",
  "tests-Functions",
  "tests-Growth"
)

for (subdir in subdirs) {
  files <- list.files(
    testthat::test_path(subdir),
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
  for (f in files) {
    source(f, local = TRUE)
  }
}
