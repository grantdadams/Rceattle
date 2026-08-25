# The spec tree must show a parameter that carries more than one linkage.
#
# Provenance: the grammar lets one parameter hold a LIST of specs -- a shared
# prior plus a fleet-specific random walk, which is how GOA pollock 2025 is
# configured. The renderer read `$formula` off that list, got NULL, and printed
# "NULL". Three of that model's six linkage rows displayed that way, and they
# were the three carrying the random-effect structure; the rows that rendered
# were the trivial intercept-only priors.

testthat::skip_on_cran()

stacked <- function() {
  data("BS2017SS", envir = environment())
  d <- suppressWarnings(suppressMessages(Rceattle::switch_check(BS2017SS)))
  d$sel_linkages <- Rceattle::build_selectivity(linkages = list(
    slp_asc = list(
      Rceattle::linkage_spec(~ 1, fleet = 1),
      Rceattle::linkage_spec(~ rw(1 | Year), fleet = 1)),
    slp_desc = Rceattle::linkage_spec(~ 1, fleet = 1)))$linkages
  d$q_linkages <- Rceattle::build_catchability(linkages = list(
    q = list(
      Rceattle::linkage_spec(~ ar1(1 | Year), fleet = 1),
      Rceattle::linkage_spec(~ rw(1 | Year), fleet = 2))))$linkages
  d
}

# Only the linkage rows. Anchored on the tree glyph followed immediately by the
# process slot, so the fleets block -- whose rows also carry "sel:" and "q:",
# as `[1] Name | Survey | sel: Logistic | q: Estimated` -- is not swept up.
tree_lines <- function(d) {
  tr <- Rceattle:::.rce_spec_tree(d)
  tr[grepl("^\\s*[├└]─ (q|sel|growth|M1|srr|comp): ", tr)]
}


testthat::test_that("a stacked parameter renders its formulas, not NULL", {
  lines <- tree_lines(stacked())
  testthat::expect_false(any(grepl("NULL", lines, fixed = TRUE)))
  testthat::expect_true(any(grepl("rw(1 | Year)", lines, fixed = TRUE)))
  testthat::expect_true(any(grepl("ar1(1 | Year)", lines, fixed = TRUE)))
})


testthat::test_that("stacked specs are tagged [i/n] and unstacked ones are not", {
  lines <- tree_lines(stacked())
  # slp_asc holds two specs, so both its rows say so.
  testthat::expect_equal(sum(grepl("slp_asc [1/2]", lines, fixed = TRUE)), 1L)
  testthat::expect_equal(sum(grepl("slp_asc [2/2]", lines, fixed = TRUE)), 1L)
  # slp_desc holds one, and carries no tag.
  desc <- grep("slp_desc", lines, fixed = TRUE, value = TRUE)
  testthat::expect_length(desc, 1L)
  testthat::expect_false(grepl("[1/1]", desc, fixed = TRUE))
})


testthat::test_that("the count is the number of linkages, not of parameters", {
  tr <- Rceattle:::.rce_spec_tree(stacked())
  # 2 (q) + 2 (slp_asc) + 1 (slp_desc) = 5 specs across 3 parameters.
  testthat::expect_true(any(grepl("linkages (5)", tr, fixed = TRUE)))
  testthat::expect_length(tree_lines(stacked()), 5L)
})


testthat::test_that("a malformed slot costs one row, not the whole tree", {
  # `$` on an atomic vector is an error, not NULL, so a stray non-spec in a
  # linkage slot must not propagate: both print methods catch a tree failure by
  # replacing the ENTIRE tree with "(spec tree unavailable)", which is a worse
  # outcome than dropping the one row that cannot be described.
  data("BS2017SS", envir = environment())
  d <- suppressWarnings(suppressMessages(Rceattle::switch_check(BS2017SS)))
  d$sel_linkages <- list(bad = c(1, 2),
                         good = Rceattle::linkage_spec(~ 1, fleet = 1))

  tr <- testthat::expect_no_error(Rceattle:::.rce_spec_tree(d))
  testthat::expect_false(any(grepl("spec tree unavailable", tr, fixed = TRUE)))
  testthat::expect_true(any(grepl("sel: good ~1", tr, fixed = TRUE)))
  testthat::expect_true(any(grepl("linkages (1)", tr, fixed = TRUE)))
})


testthat::test_that("the flattener is total on hostile input", {
  f <- Rceattle:::.rce_linkage_specs
  for (x in list(c(1, 2), "x", NULL, NA, list(), list(1, 2), list(NULL), ~ 1)) {
    testthat::expect_length(f(x), 0L)
  }
})


testthat::test_that("a model with no linkages prints no linkage block", {
  data("BS2017SS", envir = environment())
  d <- suppressWarnings(suppressMessages(Rceattle::switch_check(BS2017SS)))
  tr <- Rceattle:::.rce_spec_tree(d)
  testthat::expect_false(any(grepl("linkages (", tr, fixed = TRUE)))
})
