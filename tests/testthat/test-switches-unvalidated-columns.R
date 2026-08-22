# Three fleet_control columns reached the fit unvalidated until 5.12.0:
# Selectivity_dimension, Sel_shape_dir and Sel_shape_mode. A typo in any of them
# resolved to NA rather than erroring -- Selectivity_dimension became a missing
# selectivity dimension, the Sel_shape_* columns a missing penalty mode -- and
# nothing downstream re-checked them. The model then fitted, and reported a
# number, on an input nobody had accepted.
#
# A report-only pass over the 196 workbooks in the ecosystem that carry a
# fleet_control sheet found no value the new check rejects, so no real
# assessment is affected. That is why this ships as an error rather than a
# warning-then-error.

fc_base <- function() {
  d <- make_test_data(nyrs = 4, nages = 3)
  d <- suppressMessages(suppressWarnings(switch_check(d)))
  # switch_check() does not create Sel_shape_dir (it is only defaulted where a
  # non-parametric fleet needs it), so without this one of the three columns
  # this file exists for was silently skipped by the `next` below.
  if (is.null(d$fleet_control$Sel_shape_dir)) {
    d$fleet_control$Sel_shape_dir <- "Decreasing"
  }
  d
}

test_that("an invalid value in any of the three is an error, naming the fleet", {
  base <- fc_base()
  cases <- list(
    Selectivity_dimension = c("Age", "Length"),
    Sel_shape_dir         = c("Decreasing", "Increasing"),
    Sel_shape_mode        = c("Directional", "Smooth")
  )
  checked <- 0L
  for (col in names(cases)) {
    testthat::expect_false(is.null(base$fleet_control[[col]]),
                           info = paste(col, "missing from the fixture -- it would be skipped"))
    if (is.null(base$fleet_control[[col]])) next
    checked <- checked + 1L
    d <- base
    d$fleet_control[[col]][1] <- "Agee"     # a plausible typo, not gibberish
    err <- validate_switches(d)

    testthat::expect_gt(length(err), 0)
    testthat::expect_true(any(grepl(col, err, fixed = TRUE)),
                          info = paste(col, "error does not name the column"))
    # The message must say what IS allowed, or the user cannot act on it.
    for (v in cases[[col]]) {
      testthat::expect_true(any(grepl(v, err, fixed = TRUE)),
                            info = paste(col, "error omits the valid value", v))
    }
  }
  # All three, or the loop passed by skipping.
  testthat::expect_equal(checked, 3L)
})

test_that("every valid value passes, in the spellings the consumers accept", {
  base <- fc_base()
  # rearrange_data() matches Selectivity_dimension on the exact strings only, so
  # the integer spelling is deliberately NOT valid -- it would reach the
  # template as NA. Pin that, since the map carries an integer side.
  d0 <- base; d0$fleet_control$Selectivity_dimension <- 1
  testthat::expect_gt(length(validate_switches(d0)), 0)

  for (v in c("Smooth", "smooth", 1)) {
    d <- base; d$fleet_control$Sel_shape_mode <- v
    testthat::expect_length(validate_switches(d), 0)
  }
  for (v in c("Increasing", "increasing", "-1", "Decreasing")) {
    d <- base; d$fleet_control$Sel_shape_dir <- v
    testthat::expect_length(validate_switches(d), 0)
  }
  for (v in c("Age", "Length")) {
    d <- base
    d$fleet_control$Selectivity_dimension <- v
    testthat::expect_length(validate_switches(d), 0)
  }
  # A clean fixture raises nothing at all.
  testthat::expect_length(validate_switches(base), 0)
})

test_that("an Off fleet is not validated, and NA is allowed where the column is optional", {
  base <- fc_base()
  d <- base
  d$fleet_control$Fleet_type[1] <- "Off"
  d$fleet_control$Selectivity_dimension[1] <- "Agee"
  testthat::expect_length(validate_switches(d), 0)

  # Sel_shape_dir / Sel_shape_mode only bite where a value was supplied.
  for (col in c("Sel_shape_dir", "Sel_shape_mode")) {
    if (is.null(base$fleet_control[[col]])) next
    d2 <- base
    d2$fleet_control[[col]] <- NA
    testthat::expect_length(validate_switches(d2), 0)
  }
})

test_that("validate_switches reads its maps from the schema", {
  # The point of the change: which values a column may take is a schema fact.
  # Adding a switch column means adding a row, not another hardcoded map
  # reference inside validate_switches().
  for (col in c("Selectivity_dimension", "Sel_shape_dir", "Sel_shape_mode",
                "Selectivity", "Catchability", "Index_distribution")) {
    m <- Rceattle:::.rce_allowed_map(col)
    testthat::expect_false(is.null(m), info = paste(col, "has no map in the schema"))
    testthat::expect_gt(length(m), 1)
  }
  testthat::expect_setequal(names(Rceattle:::.rce_allowed_map("Selectivity_dimension")),
                            c("Age", "Length"))
})
