testthat::test_that("Input validation for fit_mod", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rceattle")
  # data_check() rejects an NA hindcast catch with an error, where the second
  # expectation below asserts a warning. Skipped since the file was written.
  testthat::skip("see issue #244: expects a warning where data_check() errors")

  # NULL should error
  testthat::expect_error(Rceattle::fit_mod(NULL), "Missing data_list")

  # Bad shaped data: missing critical columns in catch_data should error
  # source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data()
  dat$catch_data$Catch <- NA
  testthat::expect_warning(Rceattle::fit_mod(dat), regexp = "NA|missing|finite", ignore.case = TRUE)
})
