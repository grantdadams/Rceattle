testthat::skip_on_cran()
testthat::test_that("TMB model compiles and library exists", {
  testthat::skip_if_not_installed("TMB")

  # Use test helpers to compile
  source(file.path("tests", "testthat", "helpers.R"))
  compile_tmb_if_needed()

  # Expect the compiled shared library to be present
  lib <- paste0("ceattle_v01_11", .Platform$dynlib.ext)
  testthat::expect_true(file.exists(lib), info = lib)
})
