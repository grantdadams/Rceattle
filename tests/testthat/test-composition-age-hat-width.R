# Provenance: inst/dev/CLEANUP_BACKLOG.md Tier 0, the ceattle.cpp marker
# "will blow up if nlengths is less than nages".
#
# age_hat and age_obs_hat are written at AGE indices -- nages for a
# sexes-combined or sex-specific composition row, nages*2 for a joint-sex one
# (comp_data Sex = 3) -- whatever dimension the fleet's data are on, because a
# length composition is built at age and then converted through the age-length
# transition. Both used to be constructed as `= comp_obs`, taking their width
# from the workbook's Comp_ columns. For a model whose composition data are all
# length-based with nlengths < nages that width is nlengths, so the writes ran
# past the matrix.
#
# There is no end-to-end reproduction of the ORIGINAL overrun here, and none is
# possible from the suite: Eigen does not bounds-check in a release build, so
# the failure was a silent out-of-bounds WRITE whose effect is not observable
# from R, and no fixture can even set it up -- helpers.R fixes nlengths = nages,
# and comp_obs is one matrix across all rows, so on the bundled data its width
# is set by the widest row and no species is left short. That half is held by
# the source pin in the last test.
#
# What IS reproducible, and what the fits below pin, is the opposite error.
# Reserving nages*2 unconditionally also closes the overrun, but it hands every
# model trailing all-zero columns it never had: on BS2017SS that is age_hat at
# 42 columns against comp_obs's 25. age_hat is REPORTed and assessment scripts
# cbind it into their output (GOA_22.1.pollock.fixed.check.R), so its shape is
# not free.

# Widest age index the model's composition rows will be written at. Mirrors the
# max_age_cols loop in section 4.8 of ceattle.cpp; pinned against it below.
# comp_ctl is read from the TMB data, not from fit_mod()'s returned data_list --
# that one is the workbook as supplied and carries no rearranged comp_ctl, so
# reading it there makes every check here vacuously true.
# Columns are Fleet_code, Species, Sex, Age0_Length1, Year (R/5-rearrange_data.R).
comp_ctl_of <- function(m) {
  ctl <- m$obj$env$data$comp_ctl
  testthat::expect_gt(nrow(ctl), 0)
  ctl
}

comp_age_width <- function(m) {
  ctl   <- comp_ctl_of(m)
  nages <- as.integer(m$obj$env$data$nages)
  max(ifelse(ctl[, 3] == 3, 2L, 1L) * nages[ctl[, 2]])
}

build_only <- function(data_list) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = data_list, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = fit_control(verbose = 0))))
}


test_that("age_hat and age_obs_hat are exactly as wide as the ages written into them", {
  testthat::skip_on_cran()

  m <- build_only(make_test_data(nyrs = 8, nages = 5, seed = 1))
  r <- m$obj$report(m$obj$env$last.par)

  # Never narrower than the observations either, so the REPORTed objects cannot
  # shrink for anyone reading them.
  need <- max(comp_age_width(m), ncol(r$comp_hat))

  expect_equal(ncol(r$age_hat), need)
  expect_equal(ncol(r$age_obs_hat), need)
  expect_true(all(is.finite(r$age_hat)))
  expect_true(all(is.finite(r$age_obs_hat)))
})


test_that("a bundled single-sex model gains no trailing columns", {
  testthat::skip_on_cran()

  # BS2017SS is the discriminating case: every comp row is Sex = 0, comp_obs is
  # 25 columns wide, and max_age is 21 -- so an unconditional nages*2 reservation
  # gives age_hat 42 columns, 17 of them permanently zero, and changes the shape
  # of a REPORTed object for a model whose numbers have not moved since 2017.
  m <- build_only(BS2017SS)
  r <- m$obj$report(m$obj$env$last.par)

  expect_true(all(comp_ctl_of(m)[, 3] != 3))
  expect_equal(ncol(r$age_hat), max(comp_age_width(m), ncol(r$comp_hat)))
  expect_equal(ncol(r$age_hat), ncol(r$comp_hat))
  expect_lt(ncol(r$age_hat), 2L * max(as.integer(m$obj$env$data$nages)))
})


test_that("a joint-sex model reserves the second sex's ages", {
  testthat::skip_on_cran()

  # GOA2018SS carries Sex = 3 rows, so it is the only bundled exercise of the
  # nages*2 branch. Without it the joint-sex reservation is covered by the
  # source pin alone.
  m <- build_only(GOA2018SS)
  r <- m$obj$report(m$obj$env$last.par)

  ctl <- comp_ctl_of(m)
  expect_true(any(ctl[, 3] == 3))

  expect_equal(ncol(r$age_hat), max(comp_age_width(m), ncol(r$comp_hat)))
  expect_equal(ncol(r$age_obs_hat), ncol(r$age_hat))

  # Wide enough for every joint-sex row's own second-sex ages.
  nages <- as.integer(m$obj$env$data$nages)
  joint <- ctl[ctl[, 3] == 3, , drop = FALSE]
  expect_gte(ncol(r$age_hat), max(2L * nages[joint[, 2]]))

  expect_true(all(is.finite(r$age_hat)))
  expect_true(all(is.finite(r$age_obs_hat)))
})


test_that("the template sizes the age matrices per composition row, not from comp_obs", {
  # The static half, and the only cover for the original out-of-bounds write.
  # `matrix<Type> age_hat = comp_obs;` is the defect itself, so pin its absence:
  # a later edit that reinstates it puts the write back, and nothing else in the
  # suite would notice.
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)

  expect_false(any(grepl("age_hat = comp_obs", src, fixed = TRUE)))
  expect_false(any(grepl("age_obs_hat = comp_obs", src, fixed = TRUE)))

  # And that the width is not reserved from max_age, which would widen every
  # model. Scoped to the sizing statement rather than banning the identifier.
  sizing <- grep("max_age_cols", src, value = TRUE)
  expect_gt(length(sizing), 0)
  expect_false(any(grepl("max_age", gsub("max_age_cols", "", sizing), fixed = TRUE)))

  # Sized from the per-row age width instead, joint-sex rows included. Keep in
  # lockstep with comp_age_width() above.
  expect_true(any(grepl("int comp_row_cols = (comp_ctl(comp_row, 2) == 3) ? nages(comp_row_sp) * 2 : nages(comp_row_sp);",
                        src, fixed = TRUE)))
})
