# retrospective(), jitter() and self_test() accept fit_control() (5.11.0), but
# they refit through .refit_like(), which forwards only `phase` and `getsd`.
# Two things have to hold: a field this path cannot reach must be an error
# rather than a silent no-op, and a field the caller never touched must not be
# applied just because fit_control() gave it a default.

ctl <- function(...) Rceattle:::.rce_refit_control(fit_control(...), "retrospective")

test_that("no fit_control means no override", {
  expect_null(Rceattle:::.rce_refit_control(NULL, "retrospective"))
})

test_that("a plain fit_control() overrides nothing", {
  # Every field is at its default, so the caller asked for nothing.
  expect_identical(ctl(), fit_control()[character(0)])
  expect_length(ctl(), 0)
})

test_that("only the fields the caller changed are returned", {
  expect_identical(names(ctl(phase = TRUE)), "phase")
  expect_true(ctl(phase = TRUE)$phase)

  expect_identical(names(ctl(getsd = FALSE)), "getsd")
  expect_false(ctl(getsd = FALSE)$getsd)

  expect_setequal(names(ctl(phase = TRUE, getsd = FALSE)), c("phase", "getsd"))
})

test_that("fit_control(getsd = FALSE) cannot silently un-phase the refits", {
  # The defect this guards: fit_control() defaults `phase` to FALSE, but
  # retrospective() phases its peels because otherwise the parameters barely
  # move from the full-model fit and Mohn's rho is biased towards zero. If the
  # whole bundle were applied, a request about standard errors would quietly
  # flatten the retrospective -- a diagnostic reporting clean when it is not.
  got <- ctl(getsd = FALSE)
  expect_false("phase" %in% names(got))
})

test_that("a field the refit path cannot reach is an error, not a no-op", {
  # .refit_like() never forwards these, so accepting them would leave a user
  # believing every peel had used them.
  expect_error(ctl(newtonsteps = 3), "cannot apply")
  expect_error(ctl(newtonsteps = 3), "newtonsteps")
  expect_error(ctl(bias.correct = TRUE), "cannot apply")
  expect_error(ctl(rel_tol = 2), "cannot apply")
  # The message says where to set it instead.
  expect_error(ctl(newtonsteps = 3), "fit_mod")
})

test_that("nested and NULL-valued fields compare without spurious errors", {
  # comp_offset and TMBfilename default to NULL, and nlminb_control is a nested
  # list -- all three go through identical() in the changed-field test.
  expect_length(ctl(), 0)
  expect_length(ctl(nlminb_control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)), 0)
  expect_error(ctl(comp_offset = 1e-4), "cannot apply")
  expect_error(ctl(nlminb_control = list(trace = 1)), "cannot apply")
})

test_that("fit_control must come from fit_control()", {
  expect_error(Rceattle:::.rce_refit_control(list(phase = TRUE), "retrospective"),
               "must come from fit_control")
})

test_that("the three diagnostics take fit_control, and retrospective takes phase", {
  for (fn in c("retrospective", "jitter", "self_test")) {
    expect_true("fit_control" %in% names(formals(get(fn, envir = asNamespace("Rceattle")))),
                info = paste(fn, "should accept fit_control"))
  }
  # retrospective() used to hardcode phase = TRUE with no way to change it.
  # The default must still be TRUE, or every retrospective in the package moves.
  fm <- formals(retrospective)
  expect_true("phase" %in% names(fm))
  expect_true(isTRUE(fm$phase))
})

test_that("profile() is knowingly excluded", {
  # It refits too, but takes only a bare getsd=. Recorded in the developer
  # guide; pinned here so a future fit_control= addition updates both.
  expect_false("fit_control" %in% names(formals(getS3method("profile", "Rceattle"))))
})
