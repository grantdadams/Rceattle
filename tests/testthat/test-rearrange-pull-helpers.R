# .pull_int() / .pull_int0() replace the repeated
# `fleet_control %>% pull(col) %>% as.integer()` idiom in rearrange_data(), with
# .pull_int0() centralising the 1-based -> 0-based (R -> C++) index shift.

test_that(".pull_int coerces a fleet_control column to integer", {
  fc <- data.frame(A = c(3L, 1L, 2L), B = c("2", "0", "1"),
                   stringsAsFactors = FALSE)
  expect_identical(Rceattle:::.pull_int(fc, "A"), c(3L, 1L, 2L))
  expect_identical(Rceattle:::.pull_int(fc, "B"), c(2L, 0L, 1L)) # character -> integer
})

test_that(".pull_int0 shifts to 0-based and keeps the original result type", {
  fc <- data.frame(A = c(3L, 1L, 2L))
  expect_identical(Rceattle:::.pull_int0(fc, "A"), c(2, 0, 1))
  # `as.integer(x) - 1` is double; the helper must preserve that (not integer),
  # matching the exact type the old `pull() %>% as.integer() - 1` produced.
  expect_identical(typeof(Rceattle:::.pull_int0(fc, "A")), "double")
})
