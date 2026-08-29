# The parameter dictionary is the user-facing explanation of the TMB parameter
# vector, so it has to stay in step with the parameters build_params() actually
# creates. A name that drifts out of the dictionary shows up to the user as a
# warning from parameter_dictionary(), not as an error.

test_that("parameter_dictionary() returns the whole table by default", {
  dict <- parameter_dictionary()
  expect_s3_class(dict, "data.frame")
  expect_named(dict, c("internal", "natural", "process", "meaning", "dims"))
  expect_gt(nrow(dict), 0)
  expect_false(any(duplicated(dict$internal)))
  expect_false(any(is.na(dict$meaning)))
})

test_that("parameter_dictionary() looks up by name and by process", {
  expect_identical(parameter_dictionary("R_log_sd")$natural, "sigma_R")
  expect_identical(nrow(parameter_dictionary("R_log_sd")), 1L)

  sel <- parameter_dictionary(process = "selectivity")
  expect_true(all(sel$process == "selectivity"))
  expect_true("sel_coff" %in% sel$internal)

  # Both filters together intersect.
  expect_identical(nrow(parameter_dictionary("R_log_sd", process = "selectivity")), 0L)
})

test_that("parameter_dictionary() rejects an unknown process but only warns on an unknown name", {
  expect_error(parameter_dictionary(process = "wishful"), "Unknown process")
  expect_warning(parameter_dictionary("not_a_parameter"), "Not in the parameter dictionary")
})

test_that("every parameter build_params() creates is in the dictionary", {
  data_list <- make_test_data()
  pars <- suppressWarnings(build_params(data_list))
  missing <- setdiff(names(pars), parameter_dictionary()$internal)
  expect_identical(missing, character(0),
                   info = paste("Undocumented parameters:", paste(missing, collapse = ", ")))
})
