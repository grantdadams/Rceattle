# Provenance: inst/dev/CLEANUP_BACKLOG.md Tier 0, the ceattle.cpp marker
# "will bomb if minage > 1".
#
# Recruits arriving at age minage in year yr were spawned in yr - minage, so for
# yr < minage the stock-recruit relationship indexed ssb() at a negative column.
# Eigen does not bounds-check in a release build, so this did not bomb: it read
# adjacent memory. Measured on dev before the fix, BevertonHolt, nages = 5:
#
#   minage   objective        R[1, 1:4]
#   1        14407.3797486    1.2124872  1.2124872  1.1094781  1.0204317
#   2        15162.5962515    1.2124872  8.5948466e-314  ...
#   3        17532.154116     1.2124872  0  8.5948466e-314  ...
#
# is.finite() was TRUE throughout -- a denormalized double is finite -- so
# nothing flagged it, and the objective moved by thousands of units purely from
# reading memory that was not SSB. Only reachable with an active stock-recruit
# relationship: srr_fun = 0 never calls the SRR with an SSB.
#
# Fix: for yr < minage the SRR has no modelled SSB, so those years take the
# equilibrium mean R0 with their own deviation (srr_fun 0). This follows Stock
# Synthesis, which covers the pre-start period with an equilibrium age
# composition plus early recruitment deviations rather than extending the
# relationship backwards. WHAM avoids the boundary entirely by fixing the
# recruitment lag at one year (helper_functions.hpp: SSB(y-1), loop from y = 1).

minage_fit <- function(minage) {
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, minage = minage)
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    recFun = build_srr(srr_fun = "BevertonHolt", srr_pred_fun = "BevertonHolt",
                       srr_est_mode = 1, Bmsy_lim = 0),
    fit_control = fit_control(verbose = 0))))
}


test_that("recruitment before the first spawning year is a real number, not memory", {
  testthat::skip_on_cran()

  for (ma in c(2L, 3L)) {
    m <- minage_fit(ma)
    r <- m$obj$report(m$obj$env$last.par)

    # The years whose spawning year predates styr.
    early <- r$R[1, seq_len(ma - 1L) + 1L]

    expect_true(all(is.finite(early)))
    # A denormalized read passes is.finite(), so require the values be of the
    # order of recruitment itself rather than 1e-314.
    expect_true(all(early > 1e-6))

    # And that they are exactly the equilibrium mean with their own deviation,
    # R0 * exp(rec_dev), rather than any stock-recruit prediction. rec_dev is a
    # parameter, not a REPORTed object.
    yrs <- seq_len(ma - 1L) + 1L
    rdev <- m$obj$env$parList()$rec_dev[1, yrs]
    expect_equal(unname(early), unname(r$R0[1, yrs] * exp(rdev)), tolerance = 1e-12)
  }
})


test_that("a minage > 1 objective is reproducible across builds", {
  testthat::skip_on_cran()

  # The sharpest symptom of the old behaviour: reading uninitialized memory made
  # the objective depend on whatever happened to be there.
  objs <- vapply(1:3, function(i) {
    m <- minage_fit(3L)
    as.numeric(m$obj$fn(m$obj$par))
  }, numeric(1))

  expect_length(unique(objs), 1L)
  expect_true(is.finite(objs[1]))
})


test_that("minage = 1 is untouched -- the guard cannot fire there", {
  testthat::skip_on_cran()

  # Every bundled dataset and all three live assessments are minage = 1, so this
  # is what makes the fix safe to ship: at minage = 1 the spawning year is
  # yr - 1 >= 0 for every yr >= 1 and the guard never runs.
  m <- minage_fit(1L)
  expect_equal(as.numeric(m$obj$fn(m$obj$par)), 14407.3797486, tolerance = 1e-8)
})


test_that("no recruitment site indexes the spawning year without guarding it", {
  # Drift guard. The defect was four separate unguarded `yr - minage(sp)` reads;
  # a fifth added later would be just as silent, so require every one of them to
  # be a named local that is range-checked, not an inline array index.
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)
  src <- src[!grepl("^\\s*//", src)]

  uses <- grep("yr\\s*-\\s*minage\\(sp\\)", src, value = TRUE)
  expect_gt(length(uses), 0)

  # Each must be an assignment to a spawn-year local, never a subscript such as
  # ssb(sp, yr - minage(sp)).
  expect_true(all(grepl("int\\s+\\w*spawn_yr\\s*=", uses)),
              info = paste("unguarded:", paste(setdiff(uses, grep("int\\s+\\w*spawn_yr\\s*=", uses, value = TRUE)), collapse = " | ")))
})
