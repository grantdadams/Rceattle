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
#    The old line had a SECOND defect, which is why the fix moves results at
#    assessment_period = 1 as well: it was ifelse(sum(x) > cap, cap*x/sum(x), x)
#    with a length-1 test, and ifelse() returns the shape of its test. So the
#    scaled vector was truncated to its first element and R recycled that one
#    value over every capped row -- every species got species 1's share, and the
#    total came out over the cap unless the split happened to be equal.
#
# 2. R/10-mse_summary.R "EM uses fixed-depletion proxy for HCR 2 (0.5*0.35);
#    OM uses 0.5*SBF". HCR 2 read a literal on the depletion scale while the OM
#    scored the same rule as an absolute biomass, so the cross-tab compared two
#    criteria. Plimit is not the answer: it defaults to 0, which would make the
#    EM's P(SSB < SSBlimit) identically zero for a default ConstantF run. Under
#    msmMode > 0 both sides do read Plimit and so do report 0 at the default --
#    that is the OM's own multispecies rule, left alone here because changing it
#    would desynchronise the two sides again.
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


# The line run_mse() used to run, transcribed exactly -- ifelse() and all, since
# the recycling is the point. Assignment into a subset is what recycled the
# length-1 result, so it is reproduced here too.
apply_old_scalar_cap <- function(catch, year, cap) {
  out <- catch
  out[seq_along(catch)] <- ifelse(sum(catch) > cap,
                                  cap * catch / sum(catch),
                                  catch)
  out
}


test_that("a scalar catch cap is an annual ceiling, not an interval one", {
  # Two species, two projection years, 30 t each -- 60 t in a year, 120 t over
  # the interval. A 100 t annual cap must not bind.
  catch <- c(30, 30, 30, 30)
  year  <- c(2020, 2020, 2021, 2021)

  expect_equal(apply_scalar_cap(catch, year, cap = 100), catch)

  # The old rule summed across both years (120 > 100) and scaled by 100/120,
  # cutting each year's total to 50 t against a 100 t cap.
  old <- apply_old_scalar_cap(catch, year, cap = 100)
  expect_equal(sum(old[year == 2020]), 50)
})


test_that("the cap still binds within a year and is applied per row", {
  catch <- c(80, 80, 10, 10)
  year  <- c(2020, 2020, 2021, 2021)
  got   <- apply_scalar_cap(catch, year, cap = 100)

  # 2020 is over (160 > 100) and is scaled to exactly the cap; 2021 is under
  # (20) and is untouched.
  expect_equal(sum(got[year == 2020]), 100)
  expect_equal(got[year == 2021], c(10, 10))
  # Scaling is proportional, so the split between species is preserved.
  expect_equal(got[1] / got[2], catch[1] / catch[2])
})


test_that("the ifelse recycling moved catch at assessment_period = 1 too", {
  # One year, two species, unequal catch. This is assessment_period = 1, where
  # the interval grouping is a no-op -- so anything that differs here is the
  # recycling, not the regrouping.
  catch <- c(80, 20)
  year  <- c(2020, 2020)

  got <- apply_scalar_cap(catch, year, cap = 50)
  expect_equal(got, c(40, 10))
  expect_equal(sum(got), 50)

  # ifelse() returns the shape of its length-1 test, so the old line kept only
  # species 1's scaled catch and recycled it: both species got 40 t and the
  # total came to 80 t against a 50 t cap.
  old <- apply_old_scalar_cap(catch, year, cap = 50)
  expect_equal(old, c(40, 40))
  expect_gt(sum(old), 50)

  # An equal split is the one case where recycling coincides with the right
  # answer, which is why it cannot stand in for this test.
  even <- c(70, 70)
  expect_equal(apply_scalar_cap(even, year, cap = 100),
               apply_old_scalar_cap(even, year, cap = 100))
})


test_that("a cap that never binds left catch alone", {
  # The worse half of the same defect, and the easier one to miss: ifelse()
  # truncates whichever branch it takes, so the untouched-catch branch was
  # recycled too. A cap set generously enough to be inert still rewrote every
  # species to species 1's catch -- inventing removals rather than misallocating
  # them. A projection is not free to gain catch because a ceiling was declared.
  catch <- c(80, 20)
  year  <- c(2020, 2020)

  expect_equal(apply_scalar_cap(catch, year, cap = 500), catch)

  old <- apply_old_scalar_cap(catch, year, cap = 500)
  expect_equal(old, c(80, 80))
  expect_equal(sum(old), 160)   # true total is 100
})


test_that("HCR 2's EM overfished threshold is the OM's own absolute rule", {
  # mse_summary() routes HCR 2 through ssb_limit_thresh(), the same helper the
  # OM arm uses, so the two sides of the cross-tab cannot drift apart. The
  # helper is a local inside mse_summary(); its HCR 2 arm is copied here and
  # pinned against the source in the last test of this file.
  ssb_limit_thresh_hcr2 <- function(SBF_val) list(val = 0.5 * SBF_val, is_dep = FALSE)

  # ConstantF: 0.5 * SBF on the ABSOLUTE scale, not a depletion.
  thr <- ssb_limit_thresh_hcr2(SBF_val = 300)
  expect_equal(thr$val, 150)
  expect_false(thr$is_dep)

  # Plimit is NOT the threshold. It defaults to 0, so reading it would make
  # P(SSB < SSBlimit) identically FALSE for a default ConstantF run -- an
  # unfishable-looking stock reported as never overfished.
  expect_equal(formals(Rceattle::build_hcr)$Plimit, 0.0)
  expect_false(isTRUE(all.equal(thr$val, formals(Rceattle::build_hcr)$Plimit)))

  # And not the old literal depletion proxy either.
  expect_false(isTRUE(all.equal(thr$val, 0.5 * 0.35)))

  # A stock at 140 t against SBF = 300 is below the limit; 160 t is not. Under
  # the Plimit reading both would score FALSE.
  expect_true(140 < thr$val)
  expect_false(160 < thr$val)
})


test_that("run_mse() and mse_summary() still carry the rules these tests copy", {
  # The two tests above check arithmetic against a local copy, so pin the
  # originals: a copy that silently stops matching its source tests nothing.
  # Skipped where the sources are absent -- under R CMD check the tests run
  # against the INSTALLED package, whose R/ holds Rceattle.rdb rather than the
  # .R files, so reading them unguarded is an error rather than a check.
  find_src <- function(f) {
    p <- c(file.path("R", f), testthat::test_path("..", "..", "R", f))
    p[file.exists(p)]
  }
  mse_p  <- find_src("10-run_mse.R")
  summ_p <- find_src("10-mse_summary.R")
  testthat::skip_if(length(mse_p) == 0 || length(summ_p) == 0,
                    "R/ sources not available")

  mse  <- readLines(mse_p[1],  warn = FALSE)
  summ <- readLines(summ_p[1], warn = FALSE)

  expect_true(any(grepl("cap * new_catch_data$Catch[yr_ind] / yr_tot",
                        mse, fixed = TRUE)))
  expect_true(any(grepl("for(cap_yr in unique(new_catch_data$Year[dat_fill_ind]))",
                        mse, fixed = TRUE)))
  # HCR 2's arm in ssb_limit_thresh(), copied by the test above.
  expect_true(any(grepl("if(HCR == 2)                                  list(val = 0.5 * SBF_val,        is_dep = FALSE)",
                        summ, fixed = TRUE)))
  # The EM side reads that helper rather than a rule of its own.
  expect_true(any(grepl("if(HCR == 2 & msmMode == 0){", summ, fixed = TRUE)))
  expect_true(any(grepl("em_below <- if(em_thr$is_dep) em_q$ssb_depletion[sp, end_yr_col] < em_thr$val",
                        summ, fixed = TRUE)))
  # And that neither the literal proxy nor the degenerate Plimit reading is back.
  expect_false(any(grepl("0.5 * 0.35", summ, fixed = TRUE)))
  expect_false(any(grepl("if(HCR == 2)                              Plimit[sp]",
                        summ, fixed = TRUE)))
})
