# Under the model OSA residuals are i.i.d. N(0, 1), so a fixed |resid| > 3 flag
# expects 0.0027 n outliers per panel: 0.3 at n = 100, 13.5 at n = 5,000, so a
# long length panel always looked worse than a short age panel. The OSA panel
# flags above qnorm(1 - 0.05 / (2 n)) per panel (Bonferroni; 0.05 expected flags
# whatever n). The Pearson panel keeps 3: those residuals are not N(0, 1).

bubble_frame <- function(n, source, r) {
  data.frame(source = source, year = seq_len(n), age_length_bin = 1L,
             residual = c(r, rep(0.1, n - 1)))
}

testthat::test_that("the OSA flag is per panel and does not grow with the panel's size", {
  osa <- rbind(bubble_frame(30, "small", 3.5), bubble_frame(3000, "large", 3.5))
  p <- Rceattle:::.osa_bubble_plot(osa)
  shape <- p$data$shape
  testthat::expect_equal(shape[p$data$source == "small"][1], "outlier")   # 3.5 > 3.14
  testthat::expect_equal(shape[p$data$source == "large"][1], "normal")    # 3.5 < 4.32
  # Expected flags under the null are 0.05 per panel at any n: over 400
  # simulated N(0, 1) panels of 100 and of 5,000, the mean count is about 0.05
  # for both (the fixed 3 would give 0.27 and 13.5).
  set.seed(1)
  flags <- function(n) replicate(400, {
    f <- data.frame(source = "p", year = seq_len(n), age_length_bin = 1L,
                    residual = stats::rnorm(n))
    sum(suppressWarnings(Rceattle:::.osa_bubble_plot(f))$data$shape == "outlier")
  })
  testthat::expect_lt(mean(flags(100)), 0.12)
  testthat::expect_lt(mean(flags(5000)), 0.12)
})

testthat::test_that("the Pearson panel keeps the fixed 3, and NA residuals do not count", {
  osa <- bubble_frame(3000, "large", 3.5)
  osa$residual[2:1500] <- NA
  p <- Rceattle:::.osa_bubble_plot(osa, outlier = "fixed")
  testthat::expect_equal(p$data$shape[1], "outlier")
  p2 <- Rceattle:::.osa_bubble_plot(osa)                      # n = 1501 finite
  testthat::expect_equal(p2$data$shape[1], "normal")           # cut 4.16
  # A residual beyond the bubble scale is flagged on its untruncated value.
  osa$residual[1] <- 8
  p3 <- suppressWarnings(Rceattle:::.osa_bubble_plot(osa))   # warns that 8 is truncated
  testthat::expect_equal(p3$data$shape[1], "outlier")
})
