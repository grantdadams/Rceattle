# Provenance: inst/dev/CLEANUP_BACKLOG.md Tier 0, two markers.
#
# 1. R/10-run_mse.R "does not work for assessments that don't occur annually".
#    `dat_fill_ind` indexes every filled catch row in the assessment interval,
#    which spans assessment_period years. The scalar-cap branch tested
#    sum(Catch[dat_fill_ind]) > cap, so an annual ceiling was applied to a
#    multi-year total: at assessment_period = 2 a two-year sum was scaled to one
#    year's cap, roughly halving projected catch. The species-specific branch is
#    per row and was always right.
#
# 2. R/10-mse_summary.R "EM uses fixed-depletion proxy for HCR 2 (0.5*0.35)".
#    Every other arm reads the configured depletion; HCR 2 read a literal.
#
# Both are tested against the arithmetic rather than a fitted model: these are
# bookkeeping rules, and a full run_mse() would take minutes to exercise one
# `ifelse`.

# The scalar-cap rule as run_mse() now applies it, lifted to a function so the
# per-year grouping can be checked directly. Keep in lockstep with the "Apply
# cap" block in R/10-run_mse.R.
apply_scalar_cap <- function(catch, year, cap) {
  ind <- seq_along(catch)
  for (cap_yr in unique(year)) {
    yr_ind <- ind[year == cap_yr]
    yr_tot <- sum(catch[yr_ind])
    if (yr_tot > cap) catch[yr_ind] <- cap * catch[yr_ind] / yr_tot
  }
  catch
}


test_that("a scalar catch cap is an annual ceiling, not an interval one", {
  # Two species, two projection years, 30 t each -- 60 t in a year, 120 t over
  # the interval. A 100 t annual cap must not bind.
  catch <- c(30, 30, 30, 30)
  year  <- c(2020, 2020, 2021, 2021)

  expect_equal(apply_scalar_cap(catch, year, cap = 100), catch)

  # The old rule summed across both years (120 > 100) and scaled everything by
  # 100/120, cutting each year's total to 50 t against a 100 t cap.
  old <- if (sum(catch) > 100) 100 * catch / sum(catch) else catch
  expect_equal(sum(old[year == 2020]), 50)
  expect_lt(sum(old[year == 2020]), 100)
})


test_that("the cap still binds within a year, and is unchanged at period 1", {
  catch <- c(80, 80, 10, 10)
  year  <- c(2020, 2020, 2021, 2021)
  got   <- apply_scalar_cap(catch, year, cap = 100)

  # 2020 is over (160 > 100) and is scaled to exactly the cap; 2021 is under
  # (20) and is untouched.
  expect_equal(sum(got[year == 2020]), 100)
  expect_equal(got[year == 2021], c(10, 10))
  # Scaling is proportional, so the split between species is preserved.
  expect_equal(got[1] / got[2], catch[1] / catch[2])

  # assessment_period = 1: one year per interval, so the grouped rule and the
  # old ungrouped one agree exactly. This is why no existing MSE result moves.
  one_yr <- c(70, 70)
  expect_equal(apply_scalar_cap(one_yr, c(2020, 2020), cap = 100),
               100 * one_yr / sum(one_yr))
})


test_that("mse_summary()'s HCR 2 overfished threshold reads the configured limit", {
  # The rule as mse_summary() now applies it. HCR 2 (constant catch / constant
  # F) carries no target of its own, so it takes the configured limit, like the
  # final arm.
  em_thresh <- function(HCR, Plimit, Ptarget) {
    if (HCR == 2)                            Plimit
    else if (HCR == 4)                       0.5 * Ptarget
    else if (HCR == 5)                       0.5 * Plimit
    else if (HCR == 6 && Ptarget == 0.40)    0.25
    else if (HCR == 6)                       0.5 * Ptarget
    else                                     Plimit
  }

  # A run configured with Plimit = 0.20 is now scored against 0.20, not the
  # literal 0.175 that 0.5 * 0.35 produced regardless of configuration.
  expect_equal(em_thresh(2, Plimit = 0.20, Ptarget = 0.35), 0.20)
  expect_false(isTRUE(all.equal(em_thresh(2, 0.20, 0.35), 0.5 * 0.35)))

  # HCR 2 and the fallback arm now agree, which is the point.
  expect_equal(em_thresh(2, 0.20, 0.35), em_thresh(99, 0.20, 0.35))

  # The arms that already read their configuration are untouched.
  expect_equal(em_thresh(4, 0.20, 0.35), 0.175)
  expect_equal(em_thresh(5, 0.20, 0.35), 0.10)
  expect_equal(em_thresh(6, 0.20, 0.40), 0.25)
})


test_that("run_mse() and mse_summary() still carry the rules these tests copy", {
  # The two tests above check arithmetic against a local copy, so pin the
  # originals: a copy that silently stops matching its source tests nothing.
  mse   <- readLines(testthat::test_path("..", "..", "R", "10-run_mse.R"), warn = FALSE)
  summ  <- readLines(testthat::test_path("..", "..", "R", "10-mse_summary.R"), warn = FALSE)

  expect_true(any(grepl("cap * new_catch_data$Catch[yr_ind] / yr_tot",
                        mse, fixed = TRUE)))
  expect_true(any(grepl("for(cap_yr in unique(new_catch_data$Year[dat_fill_ind]))",
                        mse, fixed = TRUE)))
  expect_true(any(grepl("if(HCR == 2)                              Plimit[sp]",
                        summ, fixed = TRUE)))
  # And that the literal proxy is gone.
  expect_false(any(grepl("0.5 * 0.35", summ, fixed = TRUE)))
})
