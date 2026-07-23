# build_data() -- the code-first data-list constructor. Fast unit behaviour is
# unguarded; the golden-equivalence fit is skipped on CRAN.

testthat::test_that("base + no overrides equals clean_data(base)", {
  testthat::skip_if_not_installed("Rceattle")
  a <- suppressWarnings(suppressMessages(build_data(base = BS2017SS, .check = FALSE)))
  b <- suppressWarnings(suppressMessages(Rceattle::clean_data(BS2017SS)))
  testthat::expect_identical(a, b)
})

testthat::test_that("named overrides are applied", {
  testthat::skip_if_not_installed("Rceattle")
  d <- suppressWarnings(suppressMessages(
    build_data(base = BS2017SS, projyr = 2060, .check = FALSE)))
  testthat::expect_identical(d$projyr, 2060)
})

testthat::test_that("legacy top-level aliases map to canonical names", {
  testthat::skip_if_not_installed("Rceattle")
  d <- suppressWarnings(suppressMessages(
    build_data(base = BS2017SS, wt = BS2017SS$weight, .check = FALSE)))
  testthat::expect_false(is.null(d$weight))
  testthat::expect_null(d$wt)
})

testthat::test_that("a typo (near-miss of a known name) errors with a suggestion", {
  testthat::skip_if_not_installed("Rceattle")
  testthat::expect_error(
    build_data(base = BS2017SS, maturty = 1),
    "Did you mean.*maturity")
  testthat::expect_error(
    build_data(base = BS2017SS, fleet_controll = data.frame()),
    "Did you mean.*fleet_control")
})

testthat::test_that("case-variant names are caught, not passed as junk", {
  testthat::skip_if_not_installed("Rceattle")
  # Case-variant of a canonical name -> typo, with the correctly-cased suggestion
  # (not silently accepted as a second junk element).
  testthat::expect_error(
    build_data(base = BS2017SS, Maturity = 1),
    "Did you mean 'maturity'")
  # Case-variant of a legacy alias -> suggest the alias, not a distant scalar.
  testthat::expect_error(
    build_data(base = BS2017SS, WT = BS2017SS$weight),
    "Did you mean 'wt'")
})

testthat::test_that("alias collisions are rejected (post-mapping duplicate)", {
  testthat::skip_if_not_installed("Rceattle")
  # wt maps to weight; supplying both must error, not silently last-wins.
  testthat::expect_error(
    build_data(base = BS2017SS, wt = BS2017SS$weight, weight = BS2017SS$weight),
    "resolve to the same data element")
  # Two aliases for diet_data likewise collide.
  testthat::expect_error(
    build_data(base = BS2017SS,
               UobsWtAge = data.frame(), stom_prop_data = data.frame()),
    "resolve to the same data element")
})

testthat::test_that("empty-scalar required element is flagged, not passed", {
  testthat::skip_if_not_installed("Rceattle")
  # present(numeric(0)) must be FALSE so the pre-check fires here rather than
  # letting clean_data() crash with a cryptic dplyr error.
  d <- BS2017SS
  d$styr <- numeric(0)
  testthat::expect_error(build_data(base = d), "required data element")
})

testthat::test_that("a lone positional data.frame is rejected as base", {
  testthat::skip_if_not_installed("Rceattle")
  testthat::expect_error(
    build_data(data.frame(x = 1)), "must be an Rceattle data list")
})

testthat::test_that("a genuinely novel element name passes through", {
  testthat::skip_if_not_installed("Rceattle")
  expect_ok <- tryCatch({
    suppressWarnings(suppressMessages(
      build_data(base = BS2017SS, my_custom_note = 1, .check = FALSE)))
    TRUE
  }, error = function(e) FALSE)
  testthat::expect_true(expect_ok)
})

testthat::test_that("base and file are mutually exclusive", {
  testthat::expect_error(
    build_data(base = list(), file = "x.xlsx"), "at most one")
})

testthat::test_that("unnamed / mis-typed arguments are rejected", {
  testthat::skip_if_not_installed("Rceattle")
  # An unnamed positional arg is captured by `file`; alongside `base` that trips
  # the mutual-exclusion guard (a valid rejection either way).
  testthat::expect_error(
    build_data(base = BS2017SS, data.frame(x = 1)),
    "at most one|must be a single path")
  # A non-path `file` on its own is rejected with a steering message.
  testthat::expect_error(build_data(file = 42), "must be a single path")
})

testthat::test_that("pre-check flags missing required core, passes when complete", {
  testthat::skip_if_not_installed("Rceattle")
  # Minimal partial object: missing weight/fleet_control/... -> friendly error.
  testthat::expect_error(
    build_data(nspp = 1, styr = 1977, endyr = 2000),
    "required data element")
  # A complete object passes.
  testthat::expect_no_error(
    suppressWarnings(suppressMessages(build_data(base = BS2017SS))))
  # .check = FALSE lets a partial object through.
  testthat::expect_no_error(
    suppressWarnings(suppressMessages(
      build_data(nspp = 1, styr = 1977, .check = FALSE))))
})

testthat::test_that("clean_data is idempotent (build_data relies on it)", {
  testthat::skip_if_not_installed("Rceattle")
  # build_data() runs clean_data(); fit_mod() runs it again. The second pass
  # must be a no-op or a build_data() object would fit differently than the
  # same blocks fed straight to fit_mod().
  c1 <- suppressWarnings(suppressMessages(Rceattle::clean_data(BS2017SS)))
  c2 <- suppressWarnings(suppressMessages(Rceattle::clean_data(c1)))
  testthat::expect_identical(c1, c2)
})

testthat::test_that("build_data(base=) fits bit-identically to the source", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rceattle")
  bd <- suppressWarnings(suppressMessages(build_data(base = BS2017SS)))
  f_bd <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    bd, file = NULL, inits = NULL, estimateMode = 0, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  f_src <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    BS2017SS, file = NULL, inits = NULL, estimateMode = 0, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  testthat::expect_equal(f_bd$opt$objective, f_src$opt$objective, tolerance = 1e-10)
  testthat::expect_equal(f_bd$opt$par, f_src$opt$par, tolerance = 1e-10)
})
