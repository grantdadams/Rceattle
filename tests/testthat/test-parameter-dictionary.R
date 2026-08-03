# The parameter dictionary is only useful if it stays in sync with the actual
# parameter list. These tests fail loudly when a parameter is added, removed or
# renamed without updating R/0-parameter_dictionary.R -- which is the drift that
# makes this kind of lookup table rot.

.dict_params <- function() {
  dat <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))
  suppressWarnings(Rceattle::build_params(data_list = dat))
}

testthat::test_that("every parameter has a dictionary entry", {
  pars <- .dict_params()
  missing <- setdiff(names(pars), Rceattle:::.PAR_INFO$internal)
  testthat::expect_equal(
    missing, character(0),
    info = paste("parameters absent from .PAR_INFO:",
                 paste(missing, collapse = ", "))
  )
})


testthat::test_that("the dictionary has no stale entries", {
  pars <- .dict_params()
  stale <- setdiff(Rceattle:::.PAR_INFO$internal, names(pars))
  testthat::expect_equal(
    stale, character(0),
    info = paste(".PAR_INFO entries with no matching parameter:",
                 paste(stale, collapse = ", "))
  )
})


testthat::test_that("dictionary is well formed", {
  d <- Rceattle:::.PAR_INFO
  testthat::expect_equal(anyDuplicated(d$internal), 0L)
  testthat::expect_setequal(names(d),
                            c("internal", "natural", "process", "meaning", "dims"))
  # No empty cells -- a blank meaning is worse than no entry, because it
  # silently produces an uninformative message.
  for (col in names(d)) {
    testthat::expect_true(all(nzchar(d[[col]])), info = col)
    testthat::expect_false(anyNA(d[[col]]), info = col)
  }
  # Process labels are a closed set so downstream grouping can rely on them.
  testthat::expect_true(all(d$process %in% c(
    "recruitment", "mortality", "growth", "fishing", "catchability",
    "selectivity", "observation", "predation", "linkage", "internal")))
})


testthat::test_that(".par_label() is message-safe", {
  # Known parameter: name, natural-scale name and meaning.
  lbl <- Rceattle:::.par_label("log_M1")
  testthat::expect_match(lbl, "^log_M1 \\(M1\\): ")

  # Unknown parameter must fall back rather than error -- a diagnostic that
  # dies because of a stale name is worse than one that is terse.
  testthat::expect_identical(Rceattle:::.par_label("no_such_param"),
                             "no_such_param")

  # Vectorised lookup returns one row per requested name, in order.
  info <- Rceattle:::.par_info(c("rec_dev", "no_such_param", "sel_coff"))
  testthat::expect_equal(nrow(info), 3L)
  testthat::expect_equal(info$internal,
                         c("rec_dev", "no_such_param", "sel_coff"))
  testthat::expect_true(is.na(info$meaning[2]))
})
