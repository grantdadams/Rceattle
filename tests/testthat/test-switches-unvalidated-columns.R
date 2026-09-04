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


# build_map() is reached directly by run_mse() and retrospective(), so it can be
# handed a fleet_control that data_check() never validated. Selectivity and
# Time_varying_sel are then NA, and the first `==` against them fails with
# "missing value where TRUE/FALSE needed", naming neither column nor fleet.
# switch_check() is not the validator -- build_map() already calls it, and it
# passes NA through -- so the message names data_check().
#
# NA is legal on a fleet that is Off, so the guard keeps that exemption, and the
# random-effects sigma loop that runs before it is not gated on Fleet_type: it
# has to tolerate the NA the guard just allowed.

.unvalidated <- function(tv = "Off", sel = "Logistic", ftype = NULL, q = NULL) {
  d <- make_test_data(nyrs = 6, nprojyrs = 2, nages = 5)
  # fit_mod() normally supplies growth_model from growthFun(); build_map() is
  # called directly here, so set the empirical-growth default explicitly.
  d$growth_model <- rep(0, d$nspp)
  d$fleet_control$Selectivity      <- sel
  d$fleet_control$Time_varying_sel <- tv
  if (!is.null(ftype)) d$fleet_control$Fleet_type  <- ftype
  if (!is.null(q))     d$fleet_control$Catchability <- q
  d
}

.build_map_on <- function(d, random_sel = FALSE) {
  p <- suppressMessages(suppressWarnings(Rceattle::build_params(d)))
  suppressMessages(suppressWarnings(
    Rceattle::build_map(d, p, random_sel = random_sel)))
}

test_that("build_map() names the column and fleet when a switch reached it as NA", {
  for (rs in c(FALSE, TRUE)) {
    testthat::expect_error(.build_map_on(.unvalidated(tv = NA), rs),
                           "Time_varying_sel is NA for fleet")
    testthat::expect_error(.build_map_on(.unvalidated(sel = NA), rs),
                           "Selectivity is NA for fleet")
    # Names the fleet the way data_check() does, and points at the validator
    # that can actually explain it.
    testthat::expect_error(.build_map_on(.unvalidated(tv = NA), rs), "Fishery")
    testthat::expect_error(.build_map_on(.unvalidated(tv = NA), rs), "data_check")
    # Not the message a bare comparison would have given.
    testthat::expect_error(
      .build_map_on(.unvalidated(tv = NA), rs),
      "^(?!.*missing value where TRUE/FALSE needed).*$", perl = TRUE)
  }
})

test_that("an Off fleet may still carry NA, under random_sel too", {
  # switch_check() gates its Selectivity / Time_varying_sel checks on
  # `.rce_is_on`, so a fleet that is off may leave them unset. The sigma loop
  # runs before the guard and is not gated on Fleet_type, so it is the one that
  # has to tolerate the exemption -- every mode, since it keys on Time_varying_sel.
  for (rs in c(FALSE, TRUE)) {
    testthat::expect_no_error(.build_map_on(.unvalidated(tv = NA, ftype = "Off"), rs))
    testthat::expect_no_error(.build_map_on(.unvalidated(sel = NA, ftype = "Off"), rs))
    for (tv in c("IID", "AR1", "RandomWalk", "RandomWalkAscending")) {
      testthat::expect_no_error(
        .build_map_on(.unvalidated(tv = tv, sel = NA, ftype = "Off"), rs))
    }
  }
})

test_that("an all-Off fleet_control with NA switches fits end to end", {
  # The exemption is only real if data_check() also tolerates it: its
  # selectivity-form checks compare Selectivity to strings with no Fleet_type
  # gate, so an unset column there failed before build_map() was ever reached.
  testthat::skip_on_cran()
  d <- .unvalidated(tv = NA, sel = NA, ftype = "Off")
  testthat::expect_no_error(suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, msmMode = 0, estimateMode = "DebugBuild",
                      fit_control = Rceattle::fit_control(verbose = 0)))))
})

test_that("Catchability tolerates the NA switch_check() explicitly permits", {
  # Unlike the selectivity columns, NA Catchability is legal on ANY fleet
  # (R/0-switches.R allows `%in% c(NA, q_map, ...)`), so build_map() must not
  # error on it -- only avoid comparing it as though it were a string.
  testthat::expect_no_error(.build_map_on(.unvalidated(q = NA)))
})

test_that("a validated data list is unaffected", {
  for (tv in c("Off", "IID")) {
    for (rs in c(FALSE, TRUE)) {
      testthat::expect_no_error(.build_map_on(.unvalidated(tv = tv), rs))
    }
  }
})
