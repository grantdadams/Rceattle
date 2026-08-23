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
# Fix: for yr < minage the SRR has no modelled SSB, so those years take R_init,
# equilibrium recruitment at F = Finit, with their own deviation -- exactly what
# year 0 already gets at R(sp, 0) = R_init * exp(rec_dev(sp, 0)). This follows
# Stock Synthesis, which covers the pre-start period with an equilibrium age
# composition plus early recruitment deviations rather than extending the
# relationship backwards. WHAM avoids the boundary entirely by fixing the
# recruitment lag at one year (helper_functions.hpp: SSB(y-1), loop from y = 1).
#
# NOT R0: under a Beverton-Holt or Ricker hindcast build_map() maps the
# mean-recruit parameter out and the template overwrites only R0[, 1] with the
# derived (alpha - 1/SPR0)/Beta, so R0[, yr] stays at the build_params() start of
# exp(9) = 8103.08. A guard reading R0[, yr] therefore replaces a denormalized
# memory read with a number ~6700x the neighbouring recruitments, and an
# assertion written as R == R0[, yrs] * exp(rdev) pins that rather than catching
# it. The checks below are against R_init and against the neighbouring years.

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

    # A denormalized read passes is.finite(), so bound both sides against the
    # hindcast years that ARE stock-recruit predictions. 1e-314 fails low; the
    # exp(9) an R0[, yr] guard would give fails high.
    nyrs_hind <- length(m$data_list$styr:m$data_list$endyr)
    fitted    <- r$R[1, seq(ma + 1L, nyrs_hind)]
    expect_true(all(early > 0.1 * min(fitted)))
    expect_true(all(early < 10  * max(fitted)))

    # They are the initial equilibrium recruitment with their own deviation,
    # R_init * exp(rec_dev) -- the same expression year 1 is built from, so the
    # whole pre-styr-spawned block is one equilibrium. rec_dev is a parameter,
    # not a REPORTed object.
    yrs  <- seq_len(ma - 1L) + 1L
    rdev <- m$obj$env$parList()$rec_dev[1, yrs]
    expect_equal(unname(early), unname(r$R_init[1] * exp(rdev)), tolerance = 1e-12)

    # Year 1 is that same equilibrium, so the guarded years are continuous with
    # it rather than a separate rule.
    expect_equal(unname(r$R[1, 1]),
                 unname(r$R_init[1] * exp(m$obj$env$parList()$rec_dev[1, 1])),
                 tolerance = 1e-12)

    # And they are NOT R0[, yr] * exp(rec_dev). Under a stock-recruit hindcast
    # the mean-recruit parameter is mapped out and only R0[, 1] is overwritten
    # with the derived equilibrium, so R0[, yr] is a starting value. Pinning the
    # trap directly: a guard reading it would satisfy every check above.
    expect_false(isTRUE(all.equal(unname(early),
                                  unname(r$R0[1, yrs] * exp(rdev)),
                                  tolerance = 1e-8)))
  }
})


test_that("expected and realised recruitment share one anchor before styr", {
  testthat::skip_on_cran()

  # The stock-recruit penalty scores log(R) against log(R_hat). In the guarded
  # years neither side has a spawning biomass, so the two must rest on the SAME
  # equilibrium and the residual must be rec_dev alone -- otherwise a level
  # difference with no data behind it is read as information about the curve.
  #
  # The Ianelli configuration is where the two anchors can diverge: srr_fun is
  # mean recruitment and srr_pred_fun is a Beverton-Holt curve, so R_init (the
  # mean) and R_hat's own first-year value (the curve at F = Finit) are far
  # apart -- 8103.08 against 1.212487 on this fixture. Anchoring R_hat on the
  # curve cost a fixed 0.83 nats per guarded year.
  ianelli_fit <- function(minage) {
    d <- make_test_data(nyrs = 8, nages = 5, seed = 1, minage = minage)
    suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
      random_rec = FALSE,
      recFun = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                         srr_est_mode = 1, Bmsy_lim = 0, proj_mean_rec = FALSE),
      fit_control = fit_control(verbose = 0))))
  }

  for (ma in c(2L, 3L)) {
    m   <- ianelli_fit(ma)
    r   <- m$obj$report(m$obj$env$last.par)
    yrs <- seq_len(ma - 1L) + 1L

    # Both sides on R_init, so realised and expected agree exactly (rec_dev is 0
    # at estimateMode = 3) and the penalty reads no level difference.
    expect_equal(unname(r$R_hat[1, yrs]), rep(unname(r$R_init[1]), length(yrs)),
                 tolerance = 1e-12)
    expect_equal(unname(r$R[1, yrs]), unname(r$R_hat[1, yrs]), tolerance = 1e-12)

    # And specifically NOT the curve's own first-year equilibrium, which is what
    # the two anchors diverging looked like.
    expect_false(isTRUE(all.equal(unname(r$R_hat[1, yrs[1]]),
                                  unname(r$R_hat[1, 1]), tolerance = 1e-6)))
  }
})


test_that("the guarded years carry their own recruitment deviation", {
  testthat::skip_on_cran()

  # At estimateMode = 3 every rec_dev starts at 0, so the checks above compare
  # against a bare R_init and a guard that dropped exp(rec_dev) entirely would
  # pass them all. Re-report at a parameter vector with the deviations set to
  # known, distinct values: the pre-styr-spawned years must move with their own
  # deviation, not with year 1's and not at all.
  ma <- 3L
  m  <- minage_fit(ma)

  pars <- m$obj$env$last.par
  idx  <- which(names(pars) == "rec_dev")
  expect_gt(length(idx), 0)

  # rec_dev is (nspp x nyrs) in column-major order, so species 1's year yr sits
  # at idx[(yr - 1) * nspp + 1]. This fixture is single-species.
  expect_equal(as.integer(m$obj$env$data$nspp), 1L)
  dev_set <- c(0.30, -0.20, 0.45)          # years 1, 2, 3
  pars[idx[1:3]] <- dev_set

  r <- m$obj$report(pars)

  yrs   <- seq_len(ma - 1L) + 1L           # years 2 and 3: spawned before styr
  early <- r$R[1, yrs]

  expect_equal(unname(early), unname(r$R_init[1] * exp(dev_set[yrs])),
               tolerance = 1e-12)
  expect_equal(unname(r$R[1, 1]), unname(r$R_init[1] * exp(dev_set[1])),
               tolerance = 1e-12)

  # The deviations are distinct, so each year is genuinely its own -- dropping
  # exp(rec_dev), or reusing year 1's, changes these numbers.
  expect_false(isTRUE(all.equal(unname(early[1]), unname(early[2]))))
  expect_false(isTRUE(all.equal(unname(early[1]), unname(r$R_init[1]))))
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

  # Each must bind the lag to a named local, never index with it inline as
  # ssb(sp, yr - minage(sp)). This checks the binding only -- that the local is
  # then clamped or branched on is not something a line-scoped grep can see, so
  # a new site still needs reading. Its value is catching the inline form, which
  # is what all four original defects were.
  named <- "int\\s+\\w+\\s*=\\s*yr\\s*-\\s*minage\\(sp\\)"
  expect_true(all(grepl(named, uses)),
              info = paste("unguarded:", paste(setdiff(uses, grep(named, uses, value = TRUE)), collapse = " | ")))
})
