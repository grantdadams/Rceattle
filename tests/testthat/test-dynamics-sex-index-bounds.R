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

test_that("every write to sex index 1 is guarded by a two-sex condition", {
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)

  # The population arrays that carry a sex dimension.
  arrays <- c("N_at_age", "NByage0", "NByageF", "N_at_age_dB0", "N_at_age_dBF")
  pat    <- paste0("^\\s*(", paste(arrays, collapse = "|"), ")\\(sp, 1,")
  hits   <- grep(pat, src)
  expect_gt(length(hits), 0)

  # Either an explicit nsex test, or inside the `if(sex == 1)` arm of a loop
  # bounded by nsex(sp) -- both mean the slot exists.
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
