# sample_rec(sample_rec = FALSE) under a single-species curve fitted in the hindcast
# (5.33.0): log(mean R) - log(R0) applied the curve's SSB effect twice.

testthat::test_that("a single-species hindcast curve projects the mean ratio to the curve", {
  bld <- function(inits = NULL) suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = inits, estimateMode = 3, msmMode = 0,
    recFun = Rceattle::build_srr(srr_fun = "BevertonHolt"),
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  m  <- bld()
  nh <- m$data_list$endyr - m$data_list$styr + 1

  # Non-zero deviations, so the mean ratio differs from zero and from log(mean R) - log(R0).
  p <- m$estimated_params
  set.seed(1)
  p$rec_dev[, 1:nh] <- rnorm(length(p$rec_dev[, 1:nh]), 0, 0.5)
  m <- bld(p)
  q <- m$quantities

  s   <- suppressMessages(suppressWarnings(
    Rceattle::sample_rec(m, sample_rec = FALSE, update_model = FALSE)))
  dev <- unname(s$estimated_params$rec_dev[, nh + 1])
  testthat::expect_equal(dev, unname(log(rowMeans((q$R / q$R_hat)[, 1:nh]))), tolerance = 1e-10)

  # Under a hindcast curve R / R_hat is exp(rec_dev), year 0 included.
  testthat::expect_equal(dev, unname(log(rowMeans(exp(p$rec_dev[, 1:nh])))), tolerance = 1e-10)
  testthat::expect_false(isTRUE(all.equal(
    dev, unname(log(rowMeans(q$R[, 1:nh])) - log(q$R0[, 1])))))
})
