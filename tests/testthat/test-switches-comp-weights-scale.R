# A composition weight is read on a scale set by its distribution: a natural-
# scale multiplier under a multinomial, but the LOG of the starting weight under
# a Dirichlet-multinomial. write_template() seeds these columns with 1, so a
# template-built model switched to a Dirichlet-multinomial silently starts at e.
# fit_mod() now says so, once per fit. These tests pin that it reports the right
# fleets and, above all, that it changes nothing.

# Two fleets, both switched on and both carrying fitted composition rows, so the
# only thing under test is the distribution/weight pair. A fitted row is one with
# Year > 0; the helper skips fleets whose weight is never read.
fc_with <- function(dist, wt, dist_col = "Comp_distribution", wt_col = "Comp_weights",
                    fleet_type = c("Survey", "Fishery"), data_years = c(2000, 2000)) {
  fc <- data.frame(Fleet_name = c("Survey", "Fishery"),
                   Fleet_code = 1:2, Fleet_type = fleet_type,
                   stringsAsFactors = FALSE)
  fc[[dist_col]] <- dist
  fc[[wt_col]]   <- wt
  obs <- data.frame(Fleet_code = 1:2, Year = data_years)
  d <- list(fleet_control = fc)
  d[[if (dist_col == "CAAL_distribution") "caal_data" else "comp_data"]] <- obs
  d
}

flag <- function(d) Rceattle:::.rce_flag_dm_weight_scale(d)

test_that("a Dirichlet-multinomial fleet left at the template weight is reported", {
  expect_message(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1))),
                 "starting weight of e")
  # Named, so the user knows which fleet to edit.
  expect_message(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1))),
                 "Survey")
  # The multinomial fleet is also at 1, but for it 1 means 1. Not reported.
  expect_false(grepl("Fishery",
    paste(capture_messages(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1)))),
          collapse = " ")))
})

test_that("the integer spelling of the distribution is recognised too", {
  # Comp_distribution is a canonical string until convert_switches(), but
  # Diet_distribution is resolved to an integer early, so both must work.
  expect_message(flag(fc_with(c(1, 0), c(1, 1))), "starting weight of e")
})

test_that("a deliberate weight is left alone", {
  # Only an exact 1 -- the template value -- is worth a comment.
  expect_silent(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(0, 1))))
  expect_silent(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(2.5, 1))))
  # A multinomial fleet at 1 is correct and must never be reported.
  expect_silent(flag(fc_with(c("Multinomial", "Multinomial"), c(1, 1))))
  # Missing values are not a scale problem.
  expect_silent(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(NA, NA))))
})

test_that("CAAL and diet weights carry the same overload, and are reported the same way", {
  expect_message(
    flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1),
                 "CAAL_distribution", "CAAL_weights")),
    "'CAAL_weights' is 1")
  expect_message(
    flag(list(Diet_distribution = c(1, 0), Diet_comp_weights = c(1, 1),
                   spnames = c("Pollock", "Cod"))),
    "'Diet_comp_weights' is 1")
})

test_that("the report is inert -- it changes no value", {
  # This is the property that matters. A message that also edited the weights
  # would move every Dirichlet-multinomial fit in the package.
  for (d in list(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1)),
                 fc_with(c("Multinomial", "Multinomial"), c(1, 1)),
                 fc_with(c("DirichletMultinomial", "Multinomial"), c(0, 1)),
                 list(Diet_distribution = c(1, 0), Diet_comp_weights = c(1, 1)))) {
    expect_identical(suppressMessages(flag(d)), d)
  }
})

test_that("a data list without the columns is handled", {
  expect_silent(flag(list()))
  expect_identical(suppressMessages(flag(list())), list())
  expect_silent(flag(list(fleet_control = data.frame(Fleet_code = 1:2))))
})

test_that("fleets whose weight is never read are not reported", {
  # An Off fleet's composition likelihood is not evaluated...
  expect_silent(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1),
                             fleet_type = c("Off", "Fishery"))))
  # ...and neither is one with no fitted composition rows. A negative year marks
  # a row that is carried but not fitted.
  expect_silent(flag(fc_with(c("DirichletMultinomial", "Multinomial"), c(1, 1),
                             data_years = c(-2000, 2000))))
})

test_that("a ragged diet pair is left for data_check() to report", {
  # This runs before data_check(), so recycling here would raise a warning ahead
  # of the check that names the real problem.
  expect_silent(flag(list(Diet_distribution = c(1, 1, 0), Diet_comp_weights = c(1, 1))))
})
