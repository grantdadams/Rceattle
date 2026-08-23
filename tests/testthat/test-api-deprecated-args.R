# The fitted-model argument of the ten diagnostics was renamed to `object` in
# 5.11.0. The old spellings (`Rceattle`, `fit`) are still accepted so that the
# assessment scripts in ../Rceattle-models and ../GOA-ATF-ESP keep running
# unchanged. These tests pin the whole contract: the old name still reaches the
# model, it is silent, both names at once is an error, and positional calls are
# untouched.

# fn -> the name it used to take.
deprecated_model_args <- c(
  retrospective     = "Rceattle",
  jitter            = "Rceattle",
  self_test         = "Rceattle",
  sim_mod           = "Rceattle",
  sample_rec        = "Rceattle",
  remove_F          = "Rceattle",
  model_average     = "Rceattle",
  reweight_comps    = "fit",
  osa_residuals     = "fit",
  process_residuals = "fit"
)

# Deliberately NOT a fitted model, so every function stops at its own input
# guard. Old and new spellings must then produce the identical error, which is
# what proves the value was plumbed to the same place -- without running a real
# retrospective or 50 jitters to find out.
not_a_fit <- function() "not a fitted model"

# Warnings raised by a call, as a character vector.
warnings_from <- function(expr) {
  w <- character()
  withCallingHandlers(
    tryCatch(force(expr), error = function(e) invisible(NULL)),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  w
}

test_that("both formals carry a default, so clusterExport can export the frame", {
  # retrospective(), jitter() and self_test() hand the caller's frame to PSOCK
  # workers with clusterExport(cl, ls(envir = frame), frame). ls() lists formals
  # that were never supplied, and clusterExport get()s each one, so a formal
  # with no default makes the export throw. .parallel_lapply() catches that and
  # quietly falls back to serial -- on Windows only, since FORK never calls
  # clusterExport. Guarding it here because nothing else would show it.
  for (fn in names(deprecated_model_args)) {
    fm <- formals(get(fn, envir = asNamespace("Rceattle")))
    old <- deprecated_model_args[[fn]]

    expect_true(names(fm)[1] == "object",
                info = paste(fn, "must take `object` first"))
    expect_true(old %in% names(fm),
                info = paste(fn, "must still accept", old))

    # A default that is missing shows up as an empty symbol.
    for (nm in c("object", old)) {
      expect_false(identical(fm[[nm]], quote(expr = )),
                   info = paste0(fn, "(", nm, ") must have a default"))
    }
  }
})

test_that("the deprecated name is declared last, so positional calls still bind object", {
  for (fn in names(deprecated_model_args)) {
    fm  <- formals(get(fn, envir = asNamespace("Rceattle")))
    old <- deprecated_model_args[[fn]]
    nms <- setdiff(names(fm), "...")
    expect_identical(nms[length(nms)], old,
                     info = paste(fn, "must declare", old, "last"))
  }
})

test_that("a new argument is appended, so an existing positional call still means what it did", {
  # The leading formals as 5.10.0 shipped them, with only the fitted model
  # renamed. An assessment script calling retrospective(mod, 5, FALSE, 3, 4)
  # means cores = 4; inserting an argument ahead of `cores` would rebind that 4
  # to something else with no error and no warning. New arguments go on the end.
  pre_rename <- list(
    retrospective = c("object", "peels", "rescale", "nyrs_forecast", "cores", "getsd"),
    jitter        = c("object", "njitter", "sd", "phase", "seed", "cores", "getsd", "timeout"),
    self_test     = c("object", "nsim", "simulate", "seed", "cores", "getsd", "phase",
                      "start", "debug", "timeout", "process")
  )
  for (fn in names(pre_rename)) {
    nms <- names(formals(get(fn, envir = asNamespace("Rceattle"))))
    expect_identical(nms[seq_along(pre_rename[[fn]])], pre_rename[[fn]],
                     info = paste(fn, "reordered an argument that predates 5.11.0"))
  }
})

test_that("the deprecated formal is cleared, so a PSOCK worker gets the model once", {
  # clusterExport() serializes each name in the exported frame separately, so
  # while both `object` and the old name are bound to the fit, every Windows
  # worker receives a whole fitted model twice. The old name is the spelling the
  # assessment scripts use, so that is the common path, not the rare one.
  for (fn in c("retrospective", "jitter", "self_test")) {
    body_txt <- paste(deparse(body(get(fn, envir = asNamespace("Rceattle")))),
                      collapse = "\n")
    cleared  <- regexpr("Rceattle <- NULL", body_txt, fixed = TRUE)
    exported <- regexpr(".parallel_lapply(", body_txt, fixed = TRUE)

    expect_gt(cleared, 0)
    expect_gt(exported, 0)
    expect_lt(cleared, exported)
  }
})

test_that("the old name still reaches the model, and does so silently", {
  m <- not_a_fit()
  for (fn in names(deprecated_model_args)) {
    f   <- get(fn, envir = asNamespace("Rceattle"))
    old <- deprecated_model_args[[fn]]

    old_err <- tryCatch(do.call(f, stats::setNames(list(m), old)),
                        error = function(e) conditionMessage(e))
    new_err <- tryCatch(do.call(f, list(object = m)),
                        error = function(e) conditionMessage(e))

    # Never the signal that the shim failed to fire.
    expect_false(grepl("is missing, with no default", old_err, fixed = TRUE),
                 info = paste(fn, "lost its", old, "back-compatibility"))
    # Same failure point => the argument was plumbed to the same place.
    expect_identical(old_err, new_err,
                     info = paste(fn, "treats", old, "differently from object"))

    # Silent in this release; opt in with the option below. Other warnings a
    # function may raise are not this test's business -- only the rename.
    w <- warnings_from(do.call(f, stats::setNames(list(m), old)))
    expect_false(any(grepl("has been renamed to", w)),
                 info = paste(fn, "warned about the rename before 5.12.0"))
  }
})

test_that("supplying both spellings is an error, including an explicit NULL", {
  m <- not_a_fit()
  for (fn in names(deprecated_model_args)) {
    f   <- get(fn, envir = asNamespace("Rceattle"))
    old <- deprecated_model_args[[fn]]

    expect_error(
      do.call(f, stats::setNames(list(m, m), c("object", old))),
      "give the fitted model once",
      info = paste(fn, "must reject both spellings at once"))

    # `object = NULL` is still the caller naming the argument. Resolving it in
    # favour of the old name would be a guess.
    expect_error(
      do.call(f, stats::setNames(list(NULL, m), c("object", old))),
      "give the fitted model once",
      info = paste(fn, "must reject an explicit NULL plus", old))
  }
})

test_that("the deprecation warning can be switched on ahead of time", {
  m <- not_a_fit()
  op <- options(Rceattle.warn_deprecated_args = TRUE)
  on.exit(options(op), add = TRUE)
  for (fn in names(deprecated_model_args)) {
    f   <- get(fn, envir = asNamespace("Rceattle"))
    old <- deprecated_model_args[[fn]]
    w <- warnings_from(do.call(f, stats::setNames(list(m), old)))
    expect_true(any(grepl("has been renamed to", w)),
                info = paste(fn, "should warn under the option"))
  }
})

test_that("a model-less call names the argument at fault", {
  # These used to fail with "argument is of length zero" or
  # "invalid 'type' (list) of argument" several lines in.
  for (fn in names(deprecated_model_args)) {
    f <- get(fn, envir = asNamespace("Rceattle"))
    err <- suppressWarnings(tryCatch(f(), error = function(e) conditionMessage(e)))
    expect_true(
      grepl("object", err, fixed = TRUE) || grepl("not of class", err, fixed = TRUE),
      info = paste0(fn, "() gave an unhelpful error: ", err))
  }
})

test_that("profile() keeps `fitted`, because stats::profile() defines it", {
  # Renaming a method's first formal away from its generic is an R CMD check
  # failure ("S3 generic/method consistency"), not a style choice.
  expect_identical(
    names(formals(getS3method("profile", "Rceattle")))[1], "fitted")
})

test_that("functions that take something other than a fitted model kept their names", {
  # mse_summary() takes an MSE result and compare_sim() a set of simulations;
  # `object` would say less than the names they already have.
  expect_identical(names(formals(mse_summary))[1], "mse")
  expect_identical(names(formals(compare_sim))[1], "operating_mod")
})
