# Provenance: adversarial review of PR #117 (Tier 0 cleanup), defect 1.
#
# Arrays carrying a sex dimension are dimensioned max_sex = imax(nsex), so when
# EVERY species is one-sex that dimension has length 1 and sex index 1 does not
# exist. Ten lines wrote it unconditionally -- the "males get 1 - sex_ratio"
# half of each recruitment and sex-ratio assignment.
#
# The value written is 0: section 6.4 sets sex_ratio(sp, 0) = 1 for a one-sex
# species before any of them run. But the INDEX is out of range, and Eigen does
# not bounds-check in a release build (src/TMB/compile.R passes
# safebounds = FALSE), so it is a silent write into whatever the flat offset
# lands on. With max_sex = 1 that offset is exactly (sp, 0, age + 1, yr) -- the
# next age of the same species and year -- which the age loop overwrites
# immediately afterwards. Benign by execution order, not by construction: any
# reordering of those loops turns it into live corruption of numbers-at-age.
#
# Measured on BS2017SS at estimateMode = 3, compiling the template with
# TMB::compile(..., safebounds = TRUE):
#   before: "TMB has received an error from Eigen. The following condition was
#            not met: ( (dim*0 <= tup) && (tup < dim) ).all()"
#   after:  clean, objective 1537036.287629372 (unchanged)
#
# There is no reproduction from R against the shipped DLL -- an unchecked write
# is not observable -- so the invariant is pinned statically instead.
#
# The census is pinned by NAME AND COUNT, not floored at "more than zero". A
# floor is satisfied by whatever survives a rename or a reformat: measured
# 2026-10-12, an `if(true) N_at_age(sp, 1, 0, 0) = Type(0.0);` and a respelled
# `N_at_age(sp,1,0,yr)` both left the old anchored pattern at 10 hits and passed,
# and respelling five sites took it to 5 and still passed. The arrays are now
# discovered from their declarations rather than listed by hand, which found
# `N_eq` -- 4 writes that were outside the census entirely -- taking the real
# count from 10 to 15.

test_that("every write to sex index 1 is guarded by a two-sex condition", {
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)

  # Every write to sex index 1, per array. Measured on dev 2026-10-12.
  EXPECTED_SEX1_WRITES <- c(N_at_age = 7, N_eq = 4, NByage0 = 1, NByageF = 1,
                            N_at_age_dB0 = 1, N_at_age_dBF = 1)

  # The hazard class is "declared with a max_sex dimension", so a new array of
  # that shape joins the census instead of silently sitting outside it.
  decl <- grep("^\\s*array<Type>\\s+\\w+\\(\\s*nspp\\s*,\\s*max_sex\\s*,", src,
               value = TRUE)
  arrays <- unique(sub("^\\s*array<Type>\\s+(\\w+)\\(.*$", "\\1", decl))
  expect_gt(length(arrays), 0)

  # Counted with a SEPARATE grepl per array: one alternation capture is
  # leftmost-longest, so a line holding two of these arrays -- `N_eq(sp, 1, age)
  # = N_at_age(sp, 1, age, 0);` is one -- is attributed to only the last.
  wpat <- function(a) paste0("\\b", a, "\\(\\s*sp\\s*,\\s*1\\s*,")
  found <- vapply(arrays, function(a) sum(grepl(wpat(a), src)), integer(1))
  found <- found[found > 0L]

  expect_setequal(names(found), names(EXPECTED_SEX1_WRITES))
  expect_equal(found[names(EXPECTED_SEX1_WRITES)], EXPECTED_SEX1_WRITES,
               info = paste("sex-index-1 census moved:",
                            paste(sprintf("%s %d", names(found), found),
                                  collapse = ", ")))

  pat  <- paste0("\\b(", paste(arrays, collapse = "|"), ")\\(\\s*sp\\s*,\\s*1\\s*,")
  hits <- grep(pat, src)

  # 15 writes sit on 14 LINES, because one line carries two of these arrays:
  # `N_eq(sp, 1, age) = N_at_age(sp, 1, age, 0);`. Both numbers are pinned, since
  # the lookback below is per line while the census above is per write, and a
  # merge or split of that line moves one without the other.
  expect_equal(length(hits), 14L)

  # Either an explicit nsex test, or inside the `if(sex == 1)` arm of a loop
  # bounded by nsex(sp) -- both mean the slot exists. Comment lines are NOT
  # stripped before this lookback: the window is positional, and removing lines
  # lets a neighbouring block's guard drift into range and satisfy it.
  unguarded <- Filter(function(n) {
    ctx <- src[max(1, n - 6):(n - 1)]
    !any(grepl("if\\(nsex\\(sp\\) > 1\\)", ctx)) && !any(grepl("if\\(sex == 1\\)", ctx))
  }, hits)

  expect_equal(unguarded, integer(0),
               info = paste("unguarded:", paste(src[unguarded], collapse = " | ")))
})


test_that("a one-sex model reports only the sex dimension it has", {
  testthat::skip_on_cran()

  # The consequence the guards protect: BS2017SS is one-sex throughout, so
  # max_sex is 1 and nothing may claim a second sex.
  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = BS2017SS, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = fit_control(verbose = 0))))

  expect_true(all(m$obj$env$data$nsex == 1))
  expect_equal(dim(m$obj$report(m$obj$env$last.par)$N_at_age)[2], 1L)

  # And the objective is the value measured under bounds checking above, so the
  # guards did not move it.
  expect_equal(as.numeric(m$obj$fn(m$obj$par)), 1537036.287629372, tolerance = 1e-8)
})
