testthat::skip_on_cran()
testthat::test_that("TMB dynamic library loads", {
  testthat::skip_if_not_installed("TMB")
  source(file.path("tests", "testthat", "helpers.R"))

  lib <- paste0("ceattle_v01_11", .Platform$dynlib.ext)

  # Compile if missing
  compile_tmb_if_needed()

  testthat::expect_true(file.exists(lib), info = lib)

  # Load/unload safely using helper
  with_loaded_dll(lib, {
    loaded <- getLoadedDLLs()
    testthat::expect_true(any(grepl("ceattle_v01_11", names(loaded))))
  })
})
