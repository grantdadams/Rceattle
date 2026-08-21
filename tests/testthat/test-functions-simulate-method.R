# simulate() -- the stats generic for an Rceattle fit. sim_mod() does the work
# and is tested elsewhere; these cover the method's own contract: replicate
# count, the always-a-list return, seed locality, and process pass-through.

testthat::test_that("simulate() validates its arguments", {
  fit <- structure(list(), class = "Rceattle")
  testthat::expect_error(simulate.Rceattle(structure(list(), class = "lm")),
                         "must be an Rceattle")
  testthat::expect_error(simulate(fit, nsim = 0), "positive integer")
  testthat::expect_error(simulate(fit, nsim = c(1, 2)), "positive integer")
  testthat::expect_error(simulate(fit, nsim = NA), "positive integer")
})

testthat::test_that("simulate() returns nsim data sets, always as a list", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  one <- simulate(fit, nsim = 1)
  # A list even at nsim = 1, so callers never special-case the length.
  testthat::expect_type(one, "list")
  testthat::expect_length(one, 1L)
  testthat::expect_named(one, "Sim_1")
  testthat::expect_identical(class(one[[1]]), class(fit$data_list))

  three <- simulate(fit, nsim = 3)
  testthat::expect_length(three, 3L)
  testthat::expect_named(three, c("Sim_1", "Sim_2", "Sim_3"))
  # Independent draws, not one repeated.
  testthat::expect_false(identical(three[[1]]$catch_data$Catch,
                                   three[[2]]$catch_data$Catch))
})

testthat::test_that("simulate() seeds reproducibly without disturbing the stream", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  a <- simulate(fit, nsim = 2, seed = 42)
  b <- simulate(fit, nsim = 2, seed = 42)
  testthat::expect_equal(a[[1]]$catch_data$Catch, b[[1]]$catch_data$Catch)
  testthat::expect_equal(a[[2]]$catch_data$Catch, b[[2]]$catch_data$Catch)
  testthat::expect_equal(attr(a, "seed"), 42)

  # stats::simulate() restores the caller's random state, so a stream in
  # progress is not displaced by having simulated from a model.
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(stats::runif(1))
  state <- .Random.seed
  invisible(simulate(fit, nsim = 1, seed = 7))
  testthat::expect_identical(.Random.seed, state)

  # ...and without a seed it advances the stream, as an unseeded draw should.
  set.seed(99)
  invisible(stats::runif(1))
  invisible(simulate(fit, nsim = 1))
  testthat::expect_false(identical(.Random.seed, state))
  testthat::expect_true(is.numeric(before))
})

testthat::test_that("simulate() passes process through to sim_mod()", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  plain <- simulate(fit, nsim = 1, seed = 3)
  testthat::expect_null(attr(plain[[1]], "process_sim"))

  withp <- simulate(fit, nsim = 2, seed = 3, process = "recruitment")
  tr <- attr(withp[[1]], "process_sim")
  testthat::expect_true(all(c("rec_dev", "init_dev") %in% names(tr)))
  # Each replicate is its own draw of the process, not a shared one.
  testthat::expect_false(identical(tr$rec_dev,
                                   attr(withp[[2]], "process_sim")$rec_dev))
})
