# The diagnostics share one display contract.
#
# Provenance: convergence_diagnostics() was the only member of the family with a
# class, a print method, a severity vocabulary and an overall status. The rest
# returned bare frames and lists. osa_diagnostics() in particular printed sixteen
# columns at seven significant figures -- three screen-widths in an 80-column
# terminal -- with no line saying how many sources had failed; and Mohn's rho
# came back as a bare number while osa_diagnostics() shipped its null intervals,
# so two diagnostics took opposite conventions on the same question.
#
# These tests pin BOTH halves: that the display exists and carries a verdict, and
# that adding it did not change any return value. The second half is the one that
# matters for back-compatibility -- every object is still the frame or list it
# was, and downstream code indexes it unchanged.

testthat::skip_on_cran()

# A severity-tagged header line, e.g. "<Rceattle jitter>  status: OK".
expect_header <- function(txt, what) {
  hdr <- grep("^<Rceattle ", txt, value = TRUE)
  testthat::expect_length(hdr, 1L)
  testthat::expect_match(hdr, what, fixed = TRUE)
  testthat::expect_match(hdr, "status: (OK|NOTE|WARN|FAIL)")
}

shown <- function(x) utils::capture.output(print(x))


testthat::test_that("osa_diagnostics carries a class, a verdict, and its columns", {
  testthat::skip_if_not_installed("TMB")
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))
  osa <- suppressWarnings(suppressMessages(
    Rceattle::osa_residuals(fit, source = c("index", "catch"))))
  d <- Rceattle::osa_diagnostics(osa)

  testthat::expect_s3_class(d, "rceattle_osa_diagnostics")
  testthat::expect_s3_class(d, "data.frame")

  # Return value unchanged: every documented column still present and reachable.
  for (nm in c("group", "source", "fleet", "n", "sdnr", "sdnr_lo", "sdnr_hi",
               "lower", "lower_lo", "lower_hi", "upper", "upper_lo", "upper_hi",
               "sdnr_ok", "lower_ok", "upper_ok")) {
    testthat::expect_true(nm %in% names(d), info = nm)
  }
  testthat::expect_true(is.numeric(d$sdnr))
  testthat::expect_true(is.data.frame(as.data.frame(d)))

  txt <- shown(d)
  expect_header(txt, "OSA diagnostics")
  # The count that was missing: how many sources failed, not just per-row flags.
  testthat::expect_true(any(grepl("source\\(s\\) outside a null interval", txt)))
  # A per-source line for each fleet, tagged.
  testthat::expect_true(sum(grepl("^  \\[(OK|NOTE|WARN|FAIL)", txt)) >= 1)
  # And it must be narrow enough to read.
  testthat::expect_true(max(nchar(txt)) <= 100)
})


testthat::test_that("retrospective reports Mohn's rho against a band", {
  testthat::skip_if_not_installed("TMB")
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))
  retro <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 2)))

  testthat::expect_s3_class(retro, "Rceattle_retro")
  # Return value unchanged.
  testthat::expect_true(is.data.frame(retro$mohns))
  testthat::expect_true(all(c("Object", "Forecast year", "N") %in%
                             names(retro$mohns)))
  testthat::expect_true(is.list(retro$Rceattle_list))

  txt <- shown(retro)
  expect_header(txt, "retrospective")
  # The band is stated, which the bare number never did.
  testthat::expect_true(any(grepl("outside \\+/-", txt)))
  # And it is an argument, not a constant: a band of 0 flags everything.
  wide <- utils::capture.output(print(retro, band = 1e6))
  testthat::expect_true(any(grepl("none of", wide)))
})


testthat::test_that("jitter reports the fraction that reached the optimum", {
  testthat::skip_if_not_installed("TMB")
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))
  jit <- suppressWarnings(suppressMessages(Rceattle::jitter(fit, njitter = 3)))

  testthat::expect_s3_class(jit, "Rceattle_jitter")
  # Return value unchanged, plus njitter so the denominator is knowable: the
  # converged runs alone cannot say how many starts were attempted.
  testthat::expect_true(is.numeric(jit$nll))
  testthat::expect_true(is.list(jit$Rceattle_list))
  testthat::expect_equal(jit$njitter, 3)

  txt <- shown(jit)
  expect_header(txt, "jitter")
  testthat::expect_true(any(grepl("reached the best optimum", txt)))
  testthat::expect_true(any(grepl("best -log L", txt)))
})


testthat::test_that("self_test reports the fraction that converged", {
  testthat::skip_if_not_installed("TMB")
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))
  sims <- suppressWarnings(suppressMessages(
    Rceattle::self_test(fit, nsim = 2, start = "estimated", cores = 1)))

  testthat::expect_s3_class(sims, "Rceattle_selftest")
  # Return value unchanged: still the list of fits every downstream call reads.
  # `nsim` is carried because dropped runs make the returned length the wrong
  # denominator -- the same reason jitter() carries njitter.
  testthat::expect_true(is.list(sims))
  testthat::expect_equal(attr(sims, "nsim"), 2)
  if (length(sims)) testthat::expect_s3_class(sims[[1]], "Rceattle")

  # The two operations the live assessment scripts perform on this object.
  testthat::expect_type(c(sims, list(fit)), "list")
  testthat::expect_length(c(sims, list(fit)), length(sims) + 1L)
  testthat::expect_false(inherits(sims[seq_along(sims)], "Rceattle_selftest"))

  txt <- shown(sims)
  expect_header(txt, "self-test")
  testthat::expect_true(any(grepl("simulation\\(s\\) converged", txt)))
  # What generated the replicates, since it decides what the spread means.
  testthat::expect_true(any(grepl("processes held at their fitted values", txt)))
})


testthat::test_that("a self-test names the processes redrawn, not their masks", {
  # sim_mod() returns each deviation beside a same-shaped `_drawn` mask. That is
  # bookkeeping, not a second process, so listing it would report "rec_dev,
  # rec_dev_drawn" as two things redrawn -- and the line exists to say what the
  # replicates were generated under.
  mk <- function(s) structure(list(convergence = list(status = s)),
                              class = "Rceattle")
  sims <- structure(
    list(Sim_1 = mk("OK")), nsim = 1L,
    process_sim = list(Sim_1 = list(rec_dev = 1, rec_dev_drawn = TRUE,
                                    log_M1_dev = 2, log_M1_dev_drawn = TRUE)),
    class = "Rceattle_selftest")

  txt <- shown(sims)
  testthat::expect_true(any(grepl("processes redrawn: log_M1_dev, rec_dev", txt)))
  testthat::expect_false(any(grepl("_drawn", txt)))
  # And it points at the truth to compare against, which is not the OM.
  testthat::expect_true(any(grepl("not the operating model", txt)))
})


testthat::test_that("a profile says whether the grid brackets the minimum", {
  # Constructed rather than fitted: the states worth pinning are an edge
  # minimum and a grid where nothing converged, and reaching either through a
  # real profile means choosing `values` that make the optimizer fail on
  # purpose. The shipped function's own class and shape are asserted by the
  # live profile in test-functions-profile-param.R.
  mk <- function(nll, v = seq(0.1, 1.5, by = 0.2)) structure(
    list(Rceattle_list = list(), grid = data.frame(slot_1 = v), nll = nll,
         param = "sigmaR", slots = list(1)), class = "Rceattle_profile")

  interior <- mk(c(20, 12, 6, 4, 4.5, 7, 15, 30))
  txt <- shown(interior)
  expect_header(txt, "profile")
  testthat::expect_true(any(grepl("status: OK", txt)))
  testthat::expect_true(any(grepl("minimum -log L", txt)))
  # The profile-likelihood interval, read off the grid: 0.7 is the minimum and
  # 0.9 is within 1.92 of it, 1.1 is not.
  testthat::expect_true(any(grepl("within 1.92 of the minimum : 0.7 to 0.9", txt)))

  # A minimum at the last grid point has not been bracketed -- the grid ran
  # out. It plots as an ordinary curve, which is why it is said out loud.
  edge <- shown(mk(c(30, 22, 16, 12, 9, 7, 5, 4)))
  testthat::expect_true(any(grepl("status: WARN", edge)))
  testthat::expect_true(any(grepl("does not bracket it", edge)))
  testthat::expect_true(any(grepl("to >=1.5", edge)))

  # A failed point is a NOTE, not a WARN: the minimum is still bracketed.
  gappy <- shown(mk(c(20, NA, 6, 4, 4.5, NA, 15, 30)))
  testthat::expect_true(any(grepl("status: NOTE", gappy)))
  testthat::expect_true(any(grepl("2 did not converge", gappy)))

  # Nothing converged: report it and stop, rather than min() over an empty set.
  none <- shown(mk(rep(NA_real_, 8)))
  testthat::expect_true(any(grepl("status: FAIL", none)))
  testthat::expect_true(any(grepl("no grid point converged", none)))

  # A cross-profile gets no interval -- the cutoff would be chisq_k/2 and the
  # region is not an interval -- but every cell is still edge-tested.
  g <- expand.grid(slot_1 = c(0.1, 0.2, 0.3), slot_2 = c(0.1, 0.2, 0.3))
  cross <- structure(list(Rceattle_list = list(), grid = g,
                          nll = c(9, 8, 9, 8, 4, 8, 9, 8, 9), param = "M1",
                          slots = list(c(1, 1, 1), c(1, 2, 1))),
                     class = "Rceattle_profile")
  txt2 <- shown(cross)
  testthat::expect_true(any(grepl("slot_1 = 0.2, slot_2 = 0.2", txt2)))
  testthat::expect_false(any(grepl("within 1.92", txt2)))
})


testthat::test_that("the severity vocabulary is the convergence one", {
  # One vocabulary across the family, so a reader meets the same four words.
  testthat::expect_equal(Rceattle:::.CONV_SEVERITY,
                         c("OK", "NOTE", "WARN", "FAIL"))
  testthat::expect_equal(Rceattle:::.rce_worst(c("OK", "WARN", "NOTE")), "WARN")
  testthat::expect_equal(Rceattle:::.rce_worst(c("OK", "OK")), "OK")
  testthat::expect_equal(Rceattle:::.rce_worst(character(0)), "OK")
  testthat::expect_equal(Rceattle:::.rce_worst(c(NA, "FAIL")), "FAIL")
})


testthat::test_that("the table helper is total on an empty frame", {
  # A diagnostic with nothing to report must print its header and stop, not error.
  out <- utils::capture.output(
    Rceattle:::.rce_diag_table(data.frame(a = character(0)), c("A" = "a")))
  testthat::expect_length(out, 0L)
})
