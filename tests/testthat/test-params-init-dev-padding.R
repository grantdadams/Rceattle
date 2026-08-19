# =============================================================================
# init_dev padding in build_params().
#
# The cpp indexes init_dev as init_dev(sp, age - 1) over age = 1:(nages - 1) in
# every initMode -- the free-parameter branch, the equilibrium / non-equilibrium
# branches, and the init_dev penalty alike -- so columns nages:max_age are
# unused padding for ALL modes and carry the -999 sentinel.
#
# This was previously split in two, with the nages column gated on
# `data_list$initMode > 0`. By that point switch_check() has resolved initMode
# to its canonical *string*, so the test compared "FreeParams" > "0" -- TRUE by
# lexicographic collation -- and the gate was always open. These tests pin the
# invariant the gate was silently providing, so a future alias beginning with a
# digit (or any reintroduced string comparison) cannot quietly change it.
# =============================================================================

.idev <- function(mode, dat) {
  d <- dat
  d$initMode <- mode
  suppressMessages(suppressWarnings(Rceattle::build_params(d)))$init_dev
}

testthat::test_that("init_dev padding is identical for every initMode", {
  data(BS2017SS, package = "Rceattle")
  ref <- .idev(names(Rceattle:::initMode_map)[1], BS2017SS)

  # Every string alias agrees with the first.
  for (m in names(Rceattle:::initMode_map)) {
    testthat::expect_identical(.idev(m, BS2017SS), ref, info = m)
  }
  # ... and so does the equivalent integer code, which switch_check() resolves.
  for (i in unname(Rceattle:::initMode_map)) {
    testthat::expect_identical(.idev(i, BS2017SS), ref, info = paste("int", i))
  }
})

testthat::test_that("init_dev is 0 up to nages-1 and -999 from nages on", {
  data(BS2017SS, package = "Rceattle")
  idev <- .idev("FreeParams", BS2017SS)          # mode 0: the case the old gate got wrong
  max_age <- ncol(idev)

  for (sp in seq_len(BS2017SS$nspp)) {
    nages_sp <- BS2017SS$nages[sp]
    # Columns the cpp actually reads.
    testthat::expect_true(all(idev[sp, seq_len(nages_sp - 1)] == 0))
    # Columns it never reads.
    testthat::expect_true(all(idev[sp, nages_sp:max_age] == -999))
  }
})

testthat::test_that("padding respects ragged nages in a multispecies model", {
  data(BS2017MS, package = "Rceattle")
  idev <- .idev("FreeParams", BS2017MS)
  max_age <- ncol(idev)
  # BS2017MS carries species of differing nages, so the -999 boundary must move
  # per species rather than sitting at a single column.
  for (sp in seq_len(BS2017MS$nspp)) {
    nages_sp <- BS2017MS$nages[sp]
    testthat::expect_true(all(idev[sp, seq_len(nages_sp - 1)] == 0))
    testthat::expect_true(all(idev[sp, nages_sp:max_age] == -999))
  }
})
