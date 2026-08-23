# Provenance: inst/dev/CLEANUP_BACKLOG.md Tier 0, the ceattle.cpp marker
# "will blow up if nlengths is less than nages".
#
# age_hat and age_obs_hat are written at AGE indices -- up to nages*2 for
# joint-sex composition -- whatever dimension the fleet's data are on, because a
# length composition is built at age and then converted through the age-length
# transition. Both used to be constructed as `= comp_obs`, taking their width
# from the workbook's Comp_ columns. For a model whose composition data are all
# length-based with nlengths < nages that width is nlengths, so the writes ran
# past the matrix.
#
# There is no end-to-end reproduction here on purpose. Eigen does not
# bounds-check in a release build, so the pre-fix failure was an out-of-bounds
# WRITE into adjacent memory -- it does not raise, and what it corrupts is not
# reproducible from R. The invariant is what can be checked, so the invariant is
# what is pinned: the two matrices are wide enough for the age indices they are
# written at, and are not sized from comp_obs.

test_that("age_hat and age_obs_hat are wide enough for the ages written into them", {
  testthat::skip_on_cran()

  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)
  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = fit_control(verbose = 0))))

  r <- m$obj$report(m$obj$env$last.par)
  need <- 2L * max(d$nages)          # joint-sex writes age + nages*sex

  expect_gte(ncol(r$age_hat), need)
  expect_gte(ncol(r$age_obs_hat), need)

  # And never narrower than the observations they are compared against, so the
  # REPORTed objects cannot shrink for anyone reading them.
  expect_gte(ncol(r$age_hat), ncol(r$comp_hat))
  expect_gte(ncol(r$age_obs_hat), ncol(r$comp_hat))

  expect_true(all(is.finite(r$age_hat)))
  expect_true(all(is.finite(r$age_obs_hat)))
})


test_that("the template does not size the age matrices from comp_obs", {
  # The static half. `matrix<Type> age_hat = comp_obs;` is the defect itself, so
  # pin its absence: a later edit that reinstates it puts the out-of-bounds write
  # back, and nothing else in the suite would notice.
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)

  expect_false(any(grepl("age_hat = comp_obs", src, fixed = TRUE)))
  expect_false(any(grepl("age_obs_hat = comp_obs", src, fixed = TRUE)))
  # Sized by age instead.
  expect_true(any(grepl("max_age * 2", src, fixed = TRUE)))
})
