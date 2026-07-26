# Direct tests for the shared .coerce_switch_arg() helper that backs
# .coerce_srr_fun() and .coerce_M1_arg(). These validate user-supplied switch
# codes that feed the fit, so behaviour must be exact.

test_that(".coerce_switch_arg maps strings and integers to canonical codes", {
  map <- c(a = 0L, b = 1L, c = 2L)
  expect_identical(Rceattle:::.coerce_switch_arg("b", map, "x"), 1L)
  expect_identical(Rceattle:::.coerce_switch_arg(c("a", "c"), map, "x"), c(0L, 2L))
  expect_identical(Rceattle:::.coerce_switch_arg(2L, map, "x"), 2L)
})

test_that(".coerce_switch_arg errors on unknown or out-of-range values", {
  map <- c(a = 0L, b = 1L, c = 2L)
  expect_error(Rceattle:::.coerce_switch_arg("z", map, "x"), "unknown `x`")
  expect_error(Rceattle:::.coerce_switch_arg(9L, map, "x"), "must be one of")
  expect_error(Rceattle:::.coerce_switch_arg(list(1), map, "x"),
               "must be a string or integer")
})

test_that(".coerce_switch_arg enforces the length rule per switch", {
  map <- c(a = 0L, b = 1L)
  # length_exact_one = TRUE (srr_fun style)
  expect_error(Rceattle:::.coerce_switch_arg(c("a", "b"), map, "x",
                                             length_exact_one = TRUE),
               "must be length 1")
  # default (M1_model style) rejects only length 0
  expect_error(Rceattle:::.coerce_switch_arg(character(0), map, "x"),
               "length >= 1")
  expect_identical(Rceattle:::.coerce_switch_arg(c("a", "b"), map, "x"),
                   c(0L, 1L))
})

test_that(".coerce_switch_arg warns via warn_fn only for deprecated codes", {
  map <- c(a = 0L, b = 1L)
  expect_warning(
    Rceattle:::.coerce_switch_arg(7L, map, "x", deprecated = 7L,
                                  warn_fn = function(int) warning("dep!")),
    "dep!")
  # no deprecation set -> no warning
  expect_silent(Rceattle:::.coerce_switch_arg(1L, map, "x"))
  # legacy_note is appended inside the range error
  expect_error(
    Rceattle:::.coerce_switch_arg(9L, map, "x", legacy_note = ", plus legacy"),
    "plus legacy")
})

test_that("srr_fun / M1_model wrappers preserve their documented behaviour", {
  expect_identical(Rceattle:::.coerce_srr_fun("BevertonHolt", "srr_fun"), 2L)
  expect_warning(Rceattle:::.coerce_srr_fun(3, "srr_fun"), "soft-deprecated")
  expect_error(Rceattle:::.coerce_srr_fun(c(2, 2), "srr_fun"), "must be length 1")
})
